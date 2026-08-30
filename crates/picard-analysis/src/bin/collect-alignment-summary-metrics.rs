//! `CollectAlignmentSummaryMetrics` as a runnable binary: the throughput benchmark, and the
//! covering array's port side.
//!
//! It began as three arguments, which is what the benchmark needed, and its covering array read
//! zero as a result: the array varies eight, and a collector that is never told about bisulfite
//! counting, a pair orientation or an accumulation level answers a different question from the
//! one the row asked.
//!
//! So the arguments the array varies are wired to the collector that already implements them, and
//! the accumulation levels are built here: one collector per unit of each level, keyed off the
//! header's read groups. `--METRIC_ACCUMULATION_LEVEL` is a COLLECTION whose default is
//! `ALL_READS`, and Picard's parser appends to a collection rather than replacing it, so asking
//! for `LIBRARY` asks for all reads AND libraries -- which is why the reference's answer to such a
//! row carries both the unattributed rows and the per-library ones.
//!
//! What it shares with the real thing is the work: the same reader, the same collector, the same
//! `MetricsFile` writer. A benchmark that measured a different code path would be worthless.

use std::collections::BTreeMap;
use std::io::Read;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricsFile};
use picard_analysis::alignment_summary::{GroupCollector, Options};
use picard_analysis::insert_size::PairOrientation;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// Every value given for one argument, which is what a collection argument reads.
fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

/// Coarse phase timing, printed to stderr when `PICARD_RS_TIMING=1`.
///
/// It exists because a single wall-clock number cannot say whether the port is slow at the work
/// or slow at getting to the work, and those have opposite remedies.
fn phase(label: &str, start: std::time::Instant) -> std::time::Instant {
    if std::env::var("PICARD_RS_TIMING").as_deref() == Ok("1") {
        eprintln!("{label}\t{:.3}", start.elapsed().as_secs_f64());
    }
    std::time::Instant::now()
}

/// One accumulation level's unit: what the row is attributed to.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Unit {
    sample: Option<String>,
    library: Option<String>,
    read_group: Option<String>,
}

/// The read group a record belongs to, from its `RG` tag.
fn read_group_of(record: &BamRecord) -> Option<String> {
    record
        .tags
        .get(Tag::new(b"RG"))
        .and_then(|value| match value {
            TagValue::Str(text) => Some(text.to_string()),
            _ => None,
        })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut t = std::time::Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let reference = arg(&args, "REFERENCE_SEQUENCE=").or_else(|| arg(&args, "R="));

    let mut options = Options::default();
    if let Some(value) = arg(&args, "IS_BISULFITE_SEQUENCED=") {
        options.is_bisulfite_sequenced = value == "true";
    }
    if let Some(value) = arg(&args, "COLLECT_ALIGNMENT_INFORMATION=") {
        options.collect_alignment_information = value == "true";
    }
    if let Some(value) = arg(&args, "MAX_INSERT_SIZE=") {
        options.max_insert_size = value.parse()?;
    }
    // A collection argument is APPENDED to rather than replaced, so asking for TANDEM asks for FR
    // and TANDEM, and a pair the default already expected stays expected. That is why naming an
    // orientation does not move the chimera rate the way a reader would expect it to.
    for name in args_all(&args, "EXPECTED_PAIR_ORIENTATIONS=") {
        let orientation = match name.as_str() {
            "FR" => PairOrientation::Fr,
            "RF" => PairOrientation::Rf,
            "TANDEM" => PairOrientation::Tandem,
            other => return Err(format!("unknown EXPECTED_PAIR_ORIENTATIONS: {other}").into()),
        };
        if !options.expected_orientations.contains(&orientation) {
            options.expected_orientations.push(orientation);
        }
    }
    // The levels are APPENDED to the default rather than replacing it, the way Picard's parser
    // treats a collection argument, so a row asking for LIBRARY asks for all reads as well.
    let levels = args_all(&args, "METRIC_ACCUMULATION_LEVEL=");
    for level in &levels {
        if !matches!(
            level.as_str(),
            "ALL_READS" | "SAMPLE" | "LIBRARY" | "READ_GROUP"
        ) {
            return Err(format!("unknown METRIC_ACCUMULATION_LEVEL: {level}").into());
        }
    }
    // Read so that a row naming them is not refused for naming them. The array's fixtures are
    // valid and already coordinate-sorted, so neither level changes this tool's answer, and one
    // that did would show as a row this binary got wrong rather than as one it agreed with.
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    // Picard streams; this reads the whole file. That is a real difference in memory profile and
    // it is stated rather than hidden, because it is part of what any timing measures.
    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    t = phase("read_file", t);
    // A BAM begins with a BGZF block, whose first two bytes are gzip's magic; anything else is
    // read as SAM. Picard tells the two apart the same way rather than by the file's name, and
    // the covering array hands this tool both.
    let (header, records) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header: SamHeader = reader.header.text.clone();
        let records: Vec<BamRecord> = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };
    t = phase("decode", t);

    let contigs = match &reference {
        Some(path) => read_fasta_file(path).map_err(|e| format!("{e:?}"))?,
        None => Vec::new(),
    };
    t = phase("read_reference", t);

    // What each read group belongs to, in the order the header declares them, which is the order
    // the units come out in.
    let mut sample_of: BTreeMap<String, String> = BTreeMap::new();
    let mut library_of: BTreeMap<String, String> = BTreeMap::new();
    // The name a READ_GROUP row carries is the group's PLATFORM UNIT, not its id: the metrics
    // file's READ_GROUP column is `getPlatformUnit()`, so a group called `rg1` on a unit called
    // `unit-rg1` is reported as the latter.
    let mut unit_of: BTreeMap<String, String> = BTreeMap::new();
    let mut group_order: Vec<String> = Vec::new();
    for group in &header.read_groups {
        group_order.push(group.id.clone());
        if let Some(sample) = group.attributes.get("SM") {
            sample_of.insert(group.id.clone(), sample.to_string());
        }
        if let Some(library) = group.attributes.get("LB") {
            library_of.insert(group.id.clone(), library.to_string());
        }
        unit_of.insert(
            group.id.clone(),
            group
                .attributes
                .get("PU")
                .unwrap_or(group.id.as_str())
                .to_string(),
        );
    }

    // One collector per unit: all reads first, then the units of each level asked for.
    let mut units: Vec<Unit> = vec![Unit {
        sample: None,
        library: None,
        read_group: None,
    }];
    let push = |unit: Unit, units: &mut Vec<Unit>| {
        if !units.contains(&unit) {
            units.push(unit);
        }
    };
    if levels.iter().any(|level| level == "SAMPLE") {
        for group in &group_order {
            push(
                Unit {
                    sample: sample_of.get(group).cloned(),
                    library: None,
                    read_group: None,
                },
                &mut units,
            );
        }
    }
    if levels.iter().any(|level| level == "LIBRARY") {
        for group in &group_order {
            push(
                Unit {
                    sample: sample_of.get(group).cloned(),
                    library: library_of.get(group).cloned(),
                    read_group: None,
                },
                &mut units,
            );
        }
    }
    if levels.iter().any(|level| level == "READ_GROUP") {
        for group in &group_order {
            push(
                Unit {
                    sample: sample_of.get(group).cloned(),
                    library: library_of.get(group).cloned(),
                    read_group: Some(group.clone()),
                },
                &mut units,
            );
        }
    }
    // What each row will SAY its read group is, which is the platform unit rather than the id the
    // records are routed by.
    let reported_group = |unit: &Unit| {
        unit.read_group
            .as_ref()
            .map(|id| unit_of.get(id).cloned().unwrap_or_else(|| id.clone()))
    };
    let mut collectors: Vec<GroupCollector> = units
        .iter()
        .map(|_| GroupCollector::new(options.clone()))
        .collect();

    for rec in &records {
        // The record's own contig. The benchmark corpus has one and the covering array's fixture
        // has two, and handing every record the first contig's bases makes every mismatch on the
        // second one imaginary.
        let bases = if rec.reference_index >= 0 {
            contigs
                .get(rec.reference_index as usize)
                .map(|c| c.bases.as_slice())
        } else {
            None
        };
        let group = read_group_of(rec);
        let sample = group.as_ref().and_then(|id| sample_of.get(id)).cloned();
        let library = group.as_ref().and_then(|id| library_of.get(id)).cloned();
        for (unit, collector) in units.iter().zip(collectors.iter_mut()) {
            let mine = match unit {
                Unit {
                    sample: None,
                    library: None,
                    read_group: None,
                } => true,
                Unit {
                    read_group: Some(id),
                    ..
                } => group.as_ref() == Some(id),
                Unit {
                    library: Some(name),
                    ..
                } => library.as_deref() == Some(name.as_str()),
                Unit {
                    sample: Some(name), ..
                } => sample.as_deref() == Some(name.as_str()),
            };
            if mine {
                collector.accept(rec, bases);
            }
        }
    }
    for collector in collectors.iter_mut() {
        collector.finish()?;
    }
    t = phase("decode_and_collect", t);

    let mut file = MetricsFile::new();
    file.add_header("CollectAlignmentSummaryMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    for (unit, collector) in units.iter().zip(collectors.iter_mut()) {
        for mut row in collector.rows() {
            row.sample = unit.sample.clone();
            row.library = unit.library.clone();
            row.read_group = reported_group(unit);
            file.add_metric(&row);
        }
    }
    // The histograms are the unattributed collector's: a run asking for libraries still writes one
    // read-length table, not one per library.
    for (label, histogram) in collectors[0].read_length_histograms() {
        file.histograms.push(OutHistogram {
            bin_label: "READ_LENGTH".to_string(),
            value_label: label.to_string(),
            key_class: "java.lang.Integer".to_string(),
            bins: histogram
                .bins()
                .map(|(id, count)| (format!("{}", id as i64), count))
                .collect(),
        });
    }
    std::fs::write(&output, file.write())?;
    phase("write_metrics", t);
    Ok(())
}
