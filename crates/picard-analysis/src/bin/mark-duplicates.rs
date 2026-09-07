//! `MarkDuplicates` as a runnable binary: the covering array's port side.
//!
//! The library decides which records are duplicates; this is the file around that decision, and
//! the file is where the arguments live.
//!
//! `ASSUME_SORT_ORDER` is the one that shapes everything. It is not a hint: the reference reads it
//! instead of the header, refuses anything that is neither coordinate nor queryname, and, for
//! queryname, rewrites the OUTPUT header to `SO:unknown GO:query` because a queryname-grouped file
//! is not a sorted one. The writer is always presorted, so assuming coordinate over a
//! queryname-sorted file does not re-sort anything -- it reaches htsjdk's own check and fails
//! there, with htsjdk's wording.
//!
//! What is not ported here is stated rather than implied: `TAG_DUPLICATE_SET_MEMBERS` needs the
//! representative-read index of every duplicate set, which is a second sorting collection in the
//! reference and a surface of its own, so the repository declares no value for it and the array
//! does not vary it.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::mark_duplicates::{mark, Options, Record, ScoringStrategy, TaggingPolicy};

const DUPLICATE_READ: u16 = 0x400;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let metrics_file = arg(&args, "METRICS_FILE=").or_else(|| arg(&args, "M="));
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);

    let options = Options {
        scoring: match arg(&args, "DUPLICATE_SCORING_STRATEGY=").as_deref() {
            None | Some("SUM_OF_BASE_QUALITIES") => ScoringStrategy::SumOfBaseQualities,
            Some("TOTAL_MAPPED_REFERENCE_LENGTH") => ScoringStrategy::TotalMappedReferenceLength,
            Some("RANDOM") => ScoringStrategy::Random,
            Some(other) => {
                return Err(format!("unknown DUPLICATE_SCORING_STRATEGY: {other}").into())
            }
        },
        remove_duplicates: flag("REMOVE_DUPLICATES=", false),
        remove_sequencing_duplicates: flag("REMOVE_SEQUENCING_DUPLICATES=", false),
        tagging_policy: match arg(&args, "TAGGING_POLICY=").as_deref() {
            None | Some("DontTag") => TaggingPolicy::DontTag,
            Some("OpticalOnly") => TaggingPolicy::OpticalOnly,
            Some("All") => TaggingPolicy::All,
            Some(other) => return Err(format!("unknown TAGGING_POLICY: {other}").into()),
        },
        clear_dt: flag("CLEAR_DT=", true),
        optical_duplicate_pixel_distance: arg(&args, "OPTICAL_DUPLICATE_PIXEL_DISTANCE=")
            .map(|v| v.parse::<i32>())
            .transpose()?
            .unwrap_or(100),
        parse_read_names: true,
        barcode_tag: arg(&args, "BARCODE_TAG="),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

    // `openInputs` sets the header's own sort order to ASSUME_SORT_ORDER, so by the time the check
    // below reads "the header", the header is already saying what it was told to say -- which is
    // why the refusal prints the assumed order twice.
    let assumed = arg(&args, "ASSUME_SORT_ORDER=").or_else(|| arg(&args, "ASO="));
    let mut header = header;
    if let Some(order) = &assumed {
        header.set_sort_order(order);
    }
    let header_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    if !matches!(header_order.as_str(), "coordinate" | "queryname") {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: This program requires input that \
             are either coordinate or query sorted (according to the header, or at least \
             ASSUME_SORT_ORDER and the content.) Found ASSUME_SORT_ORDER={} and header \
             sortorder={header_order}",
            assumed.as_deref().unwrap_or("null")
        );
        std::process::exit(1);
    }

    // The read group's library, which is what a duplicate set is cut by, and the group's index in
    // the header, which is what `closeEnough` compares.
    let library_of = |record: &BamRecord| -> (String, i32) {
        let id = string_tag(record, b"RG");
        match id.and_then(|id| {
            header
                .read_groups
                .iter()
                .position(|group| group.id == id)
                .map(|position| (position, &header.read_groups[position]))
        }) {
            Some((position, group)) => (
                group
                    .attributes
                    .get("LB")
                    .unwrap_or("Unknown Library")
                    .to_string(),
                position as i32,
            ),
            None => ("Unknown Library".to_string(), -1),
        }
    };

    let marked_records: Vec<Record> = records
        .iter()
        .map(|record| {
            let (library, read_group) = library_of(record);
            Record {
                name: record.read_name.clone(),
                flags: record.flags,
                reference_index: record.reference_index,
                alignment_start: record.alignment_start,
                cigar: record.cigar.clone(),
                qualities: record.base_qualities.clone(),
                mate_reference_index: record.mate_reference_index,
                library,
                read_group,
                barcode: options.barcode_tag.as_ref().and_then(|tag| {
                    match record
                        .tags
                        .get(Tag::new(tag.as_bytes().try_into().unwrap_or(b"RX")))
                    {
                        Some(TagValue::Str(value)) => Some(value.clone()),
                        _ => None,
                    }
                }),
                existing_dt: string_tag(record, b"DT"),
                mate_cigar: None,
                mate_alignment_start: record.mate_alignment_start,
            }
        })
        .collect();

    let marking = mark(&marked_records, &options);

    // `createOutHeader`: assuming queryname says the file is queryname GROUPED, which is not a
    // sort order, so the output says so.
    let mut out_header = header.clone();
    if assumed.as_deref() == Some("queryname") {
        out_header.set_sort_order("unknown");
        out_header.set_group_order("query");
    }
    let add_pg_tag = flag("ADD_PG_TAG_TO_READS=", true);

    // htsjdk's own check in `SAMFileWriterImpl.addAlignment`, which is where a file assumed to be
    // coordinate-sorted and not is caught. The writer is presorted, so nothing re-sorts and the
    // record that goes backwards is the one that throws.
    let coordinate_out = out_header.attributes.get("SO") == Some("coordinate");
    let mut previous: Option<(i32, i32)> = None;

    let mut writer = BamWriter::new(Vec::new(), &out_header).map_err(|e| format!("{e:?}"))?;
    for (index, record) in records.iter().enumerate() {
        if !marking.written[index] {
            continue;
        }
        if coordinate_out {
            // `SAMSortOrderChecker.isSorted` under `coordinate`, whose key is the contig NAME and
            // the start, and whose message prints that key for both records.
            let here = (record.reference_index, record.alignment_start);
            if let Some(before) = previous {
                // `SAMRecordCoordinateComparator.fileOrderCompare`, in which an UNMAPPED record
                // (reference index -1) sorts LAST rather than first: comparing the pair as plain
                // tuples finds a violation one record early, and the message names the record it
                // found rather than the record htsjdk finds.
                let file_order = |left: (i32, i32), right: (i32, i32)| -> i32 {
                    if left.0 == -1 {
                        return if right.0 == -1 { 0 } else { 1 };
                    }
                    if right.0 == -1 {
                        return -1;
                    }
                    if left.0 != right.0 {
                        return left.0 - right.0;
                    }
                    left.1 - right.1
                };
                if file_order(before, here) > 0 {
                    let key = |reference_index: i32, start: i32| -> String {
                        let name = usize::try_from(reference_index)
                            .ok()
                            .and_then(|i| out_header.sequences.get(i))
                            .map(|s| s.name.as_str())
                            .unwrap_or("*");
                        format!("{name}:{start}")
                    };
                    // Both keys are the OFFENDING record, and that is htsjdk's, not a slip here:
                    // `isSorted` sets `prev = rec` before returning false, so the message's
                    // `getPreviousRecord()` returns the record that has just been rejected.
                    eprintln!(
                        "Exception in thread \"main\" java.lang.IllegalArgumentException: \
                         Alignments added out of order in SAMFileWriterImpl.addAlignment for \
                         file://{output}. Sort order is coordinate. Offending records are at [{}] \
                         and [{}]",
                        key(here.0, here.1),
                        key(here.0, here.1)
                    );
                    std::process::exit(1);
                }
            }
            previous = Some(here);
        }

        let mut written = record.clone();
        written.flags &= !DUPLICATE_READ;
        if marking.duplicate[index] {
            written.flags |= DUPLICATE_READ;
        }
        written.tags.remove(Tag::new(b"DT"));
        if let Some(code) = &marking.duplicate_type[index] {
            written
                .tags
                .insert(Tag::new(b"DT"), TagValue::Str(code.clone()));
        }
        if add_pg_tag {
            written
                .tags
                .insert(Tag::new(b"PG"), TagValue::Str("MarkDuplicates".to_string()));
        }
        writer.write(&written).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;

    // The metrics file is written because the tool requires the argument, and its content is not
    // what this array compares: the OUTPUT is.
    if let Some(path) = metrics_file {
        let mut text = String::from("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n");
        text.push_str(
            "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\t\
             UNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t\
             READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n",
        );
        for m in &marking.metrics {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                m.library,
                m.unpaired_reads_examined,
                m.read_pairs_examined,
                m.secondary_or_supplementary,
                m.unmapped_reads,
                m.unpaired_read_duplicates,
                m.read_pair_duplicates,
                m.read_pair_optical_duplicates,
                m.percent_duplication,
                m.estimated_library_size.unwrap_or(0),
            ));
        }
        std::fs::write(path, text)?;
    }
    Ok(())
}
