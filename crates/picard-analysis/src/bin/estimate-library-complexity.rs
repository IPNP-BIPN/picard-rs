//! `EstimateLibraryComplexity` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.markduplicates.EstimateLibraryComplexity.doWork` at tag 3.4.0. The quality
//! check, the grouping, the match rule and the metrics live in
//! `picard_analysis::estimate_library_complexity`.
//!
//! The tool never looks at where a read aligned. It takes both ends of a template in READ order,
//! sorts every template by the first `MIN_IDENTICAL_BASES` of both ends, and inside each run of
//! equal prefixes counts how many templates are the same sequence give or take `MAX_DIFF_RATE`.
//! The size of each such set is a bin of a histogram, one histogram per library, and the metrics
//! are read off those bins.
//!
//! Three things decide the numbers and none of them is the alignment.
//!
//! A template is dropped before the sort if either end fails the quality check, which is an
//! INTEGER mean over the read and a no-call anywhere in the seed.
//!
//! `MIN_GROUP_COUNT` drops a bin from the METRICS and not from the histogram, so a file whose
//! duplicates are all in bins of one reports zeros beside a histogram that says otherwise.
//!
//! The optical duplicates are found per SET, with the template that claimed the others as the
//! keeper, and a set of four takes the finder's graph path rather than its pairwise one.

use std::collections::HashMap;
use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricsFile};
use picard_analysis::collect_duplicate_metrics::DuplicationMetrics;
use picard_analysis::estimate_library_complexity::{
    groups, metrics, passes_quality_check, search_duplicates, seed_order, LibraryHistograms,
    PairedRead, DEFAULT_MAX_DIFF_RATE, DEFAULT_MAX_GROUP_RATIO, DEFAULT_MIN_GROUP_COUNT,
    DEFAULT_MIN_IDENTICAL_BASES, DEFAULT_MIN_MEAN_QUALITY,
};
use picard_analysis::java_hash_map::{string_hash_code, JavaHashMap};
use picard_analysis::mark_duplicates::{location, Location, DEFAULT_OPTICAL_DUPLICATE_DISTANCE};

const READ_PAIRED_FLAG: u16 = 0x1;
const READ_STRAND_FLAG: u16 = 0x10;
const FIRST_OF_PAIR_FLAG: u16 = 0x40;
const SECOND_OF_PAIR_FLAG: u16 = 0x80;
const SECONDARY_FLAG: u16 = 0x100;
const SUPPLEMENTARY_FLAG: u16 = 0x800;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn args_of(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

/// The `@RG` lines of one header, as `(id, library)` in the order they were written.
fn read_groups(header: &htsjdk_bam::header::SamHeader) -> Vec<(String, Option<String>)> {
    header
        .read_groups
        .iter()
        .map(|group| {
            (
                group.id.clone(),
                group.attributes.get("LB").map(str::to_string),
            )
        })
        .collect()
}

fn tag_string(record: &BamRecord, name: &str) -> Option<String> {
    let bytes = name.as_bytes();
    if bytes.len() != 2 {
        return None;
    }
    match record.tags.get(Tag::new(&[bytes[0], bytes[1]])) {
        Some(TagValue::Str(value)) => Some(value.clone()),
        _ => None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let inputs = {
        let mut inputs = args_of(&args, "INPUT=");
        inputs.extend(args_of(&args, "I="));
        inputs
    };
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let number = |key: &str, default: i64| -> i64 {
        arg(&args, key)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let seed = number("MIN_IDENTICAL_BASES=", DEFAULT_MIN_IDENTICAL_BASES as i64).max(0) as usize;
    let min_quality = number("MIN_MEAN_QUALITY=", i64::from(DEFAULT_MIN_MEAN_QUALITY)) as i32;
    let min_group_count = number("MIN_GROUP_COUNT=", DEFAULT_MIN_GROUP_COUNT);
    let max_group_ratio = number("MAX_GROUP_RATIO=", DEFAULT_MAX_GROUP_RATIO);
    let max_read_length = number("MAX_READ_LENGTH=", 0).max(0) as usize;
    let optical_distance = number(
        "OPTICAL_DUPLICATE_PIXEL_DISTANCE=",
        i64::from(DEFAULT_OPTICAL_DUPLICATE_DISTANCE),
    ) as i32;
    let max_diff_rate = arg(&args, "MAX_DIFF_RATE=")
        .and_then(|value| value.parse().ok())
        .unwrap_or(DEFAULT_MAX_DIFF_RATE);
    let barcode_tag = arg(&args, "BARCODE_TAG=");
    let read_one_barcode_tag = arg(&args, "READ_ONE_BARCODE_TAG=");
    let read_two_barcode_tag = arg(&args, "READ_TWO_BARCODE_TAG=");
    let use_barcodes =
        barcode_tag.is_some() || read_one_barcode_tag.is_some() || read_two_barcode_tag.is_some();

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    // The read groups accumulate ACROSS the input files, and a template's read group is an index
    // into that growing list, so the same group in two files is two entries unless the headers
    // agree.
    let mut all_read_groups: Vec<(String, Option<String>)> = Vec::new();
    let mut pairs: Vec<PairedRead> = Vec::new();
    let mut records_read: i64 = 0;

    for input in &inputs {
        let mut raw = Vec::new();
        std::fs::File::open(input)?.read_to_end(&mut raw)?;
        let (header, records) = if raw.starts_with(&[0x1f, 0x8b]) {
            let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
            let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
            let header = reader.header.text.clone();
            let records: Vec<BamRecord> = reader
                .map(|r| r.map_err(|e| format!("{e:?}")))
                .collect::<Result<_, _>>()?;
            (header, records)
        } else {
            let text = String::from_utf8(raw)?;
            read_sam(&text).map_err(|e| format!("{e:?}"))?
        };
        all_read_groups.extend(read_groups(&header));

        // `pendingByName`, which pairs the two ends however far apart they lie.
        struct Pending {
            read1: Option<Vec<u8>>,
            read2: Option<Vec<u8>>,
            quality_ok: bool,
            read_group: i32,
            location: Location,
            barcodes: (i32, i32, i32),
            added: bool,
        }
        let mut pending: HashMap<String, Pending> = HashMap::new();

        for record in &records {
            if record.flags & READ_PAIRED_FLAG == 0 {
                continue;
            }
            if record.flags & (FIRST_OF_PAIR_FLAG | SECOND_OF_PAIR_FLAG) == 0 {
                continue;
            }
            if record.flags & (SECONDARY_FLAG | SUPPLEMENTARY_FLAG) != 0 {
                continue;
            }
            records_read += 1;

            let entry = pending.entry(record.read_name.clone()).or_insert_with(|| {
                let parsed = location(&record.read_name);
                // `addLocationInformation` returning false leaves the read group at -1, and the
                // library then falls back to "Unknown" however the record was tagged.
                let read_group = if parsed.known {
                    let id = tag_string(record, "RG");
                    id.and_then(|id| {
                        all_read_groups
                            .iter()
                            .position(|(group, _)| *group == id)
                            .map(|index| index as i32)
                    })
                    .unwrap_or(-1)
                } else {
                    -1
                };
                Pending {
                    read1: None,
                    read2: None,
                    quality_ok: true,
                    read_group,
                    location: parsed,
                    barcodes: (0, 0, 0),
                    added: false,
                }
            });

            entry.quality_ok = entry.quality_ok
                && passes_quality_check(
                    &record.read_bases,
                    &record.base_qualities,
                    seed,
                    min_quality,
                    max_read_length,
                );

            let mut bases = record.read_bases.clone();
            if record.flags & READ_STRAND_FLAG != 0 {
                htsjdk_bam::sequence::reverse_complement(&mut bases);
            }
            let hash_of = |tag: &Option<String>| -> i32 {
                tag.as_ref()
                    .and_then(|tag| tag_string(record, tag))
                    .map_or(0, |value| string_hash_code(&value))
            };
            if record.flags & FIRST_OF_PAIR_FLAG != 0 {
                entry.read1 = Some(bases);
                if use_barcodes {
                    entry.barcodes.0 = hash_of(&barcode_tag);
                    entry.barcodes.1 = hash_of(&read_one_barcode_tag);
                }
            } else {
                entry.read2 = Some(bases);
                if use_barcodes {
                    entry.barcodes.2 = hash_of(&read_two_barcode_tag);
                }
            }

            if entry.read1.is_some() && entry.read2.is_some() && entry.quality_ok && !entry.added {
                entry.added = true;
                let library = if entry.read_group == -1 {
                    "Unknown".to_string()
                } else {
                    all_read_groups[entry.read_group as usize]
                        .1
                        .clone()
                        .unwrap_or_else(|| "Unknown".to_string())
                };
                pairs.push(PairedRead {
                    read1: entry.read1.clone().expect("read one"),
                    read2: entry.read2.clone().expect("read two"),
                    library,
                    read_group: entry.read_group,
                    location: entry.location,
                    barcodes: entry.barcodes,
                });
            }
        }
    }

    // `SortingCollection` with `PairedReadComparator`. Both sorts are stable, so templates whose
    // seeds are equal keep the order the file gave them.
    pairs.sort_by(|left, right| seed_order(left, right, seed));

    // `meanGroupSize`, whose divisor is four to the power of twice the seed: a long seed makes the
    // expected group one and every real group suspicious.
    let mean_group_size = std::cmp::max(
        1,
        (records_read / 2) / 4i64.saturating_pow((seed as u32).saturating_mul(2)).max(1),
    );

    let mut histograms: JavaHashMap<LibraryHistograms> = JavaHashMap::new();
    for range in groups(&pairs, seed) {
        let group = &pairs[range.clone()];
        if group.len() as i64 > mean_group_size * max_group_ratio {
            // The reference warns and drops the group whole.
            continue;
        }
        // `splitByLibrary`, whose map is a `HashMap`: the libraries of one group are searched in
        // that map's order, which is the order their names hash into.
        let mut by_library: JavaHashMap<Vec<usize>> = JavaHashMap::new();
        for (index, pair) in group.iter().enumerate() {
            match by_library.remove(&pair.library) {
                Some(mut existing) => {
                    existing.push(index);
                    by_library.put(&pair.library, existing);
                }
                None => by_library.put(&pair.library, vec![index]),
            }
        }
        let libraries: Vec<(String, Vec<usize>)> = by_library
            .iter()
            .map(|(name, members)| (name.to_string(), members.clone()))
            .collect();
        for (library, members) in libraries {
            if !histograms.contains_key(&library) {
                histograms.put(&library, LibraryHistograms::default());
            }
            let mut entry = histograms.remove(&library).expect("the library");
            let sequences: Vec<&PairedRead> = members.iter().map(|index| &group[*index]).collect();
            search_duplicates(
                &sequences,
                seed,
                max_diff_rate,
                max_read_length,
                optical_distance,
                use_barcodes,
                &mut entry,
            );
            histograms.put(&library, entry);
        }
    }

    let mut file = MetricsFile::new();
    file.add_header("EstimateLibraryComplexity <command line>");
    file.add_header("Started on: <timestamp>");
    // `for (final String library : duplicationHistosByLibrary.keySet())`: the rows and the
    // histogram columns come out in the map's own order, not in the order the libraries were met.
    let libraries: Vec<(String, LibraryHistograms)> = histograms
        .iter()
        .map(|(name, value)| (name.to_string(), value.clone()))
        .collect();
    let mut rows: Vec<DuplicationMetrics> = Vec::new();
    let mut out_histograms: Vec<OutHistogram> = Vec::new();
    for (library, entry) in &libraries {
        let bins: Vec<(i64, i64, i64)> = entry
            .duplication
            .iter()
            .map(|(size, count)| (*size, *count, entry.optical.get(size).copied().unwrap_or(0)))
            .collect();
        let computed = metrics(&bins, min_group_count);
        rows.push(DuplicationMetrics {
            library: library.clone(),
            read_pairs_examined: computed.read_pairs_examined,
            read_pair_duplicates: computed.read_pair_duplicates,
            read_pair_optical_duplicates: computed.read_pair_optical_duplicates,
            ..DuplicationMetrics::default()
        });
        out_histograms.push(OutHistogram {
            bin_label: "duplication_group_count".to_string(),
            value_label: library.clone(),
            key_class: "java.lang.Integer".to_string(),
            bins: entry
                .duplication
                .iter()
                .map(|(size, count)| (size.to_string(), *count as f64))
                .collect(),
        });
    }
    for row in &rows {
        file.add_metric(row);
    }
    file.histograms = out_histograms;
    std::fs::write(&output, file.write())?;
    Ok(())
}
