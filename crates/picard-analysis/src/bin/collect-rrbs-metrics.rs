//! `CollectRrbsMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectRrbsMetrics.doWork` and `RrbsMetricsCollector` at tag 3.4.0:
//!
//! * `customCommandLineValidation`, every message at once;
//! * the prefix that gains a `.` unless it ends in one, the sort refusal unless `ASSUME_SORTED`
//!   (whose default is false here), and the reference walker that refuses to rewind once it is
//!   assumed;
//! * every mapped record on a contig `SEQUENCE_NAMES` keeps, unfiltered otherwise: secondary,
//!   supplementary and duplicate records all count;
//! * `MultiLevelCollector`'s units, levels in declaration order and units in header order;
//! * per read: the length and mismatch filters, then each alignment block read on the strand it
//!   was sequenced from. The CpG branch reads the block's own qualities and the non-CpG branch the
//!   whole read's, with the same index, which is Picard's and is kept;
//! * each unit's CpG sites in `CpgLocation` order: contig name as a string, then the 0-based
//!   position.
//!
//! The chart is R's and is not drawn; the reference ignores R's exit status.

use std::collections::BTreeMap;

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::header::ReadGroup;
use htsjdk_bam::sequence::{
    bases_equal, bisulfite_bases_equal, complement, count_mismatches, is_bisulfite_converted,
};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::java_hash_map::JavaHashMap;
use picard_analysis::metrics_cli::{
    absolute, fail, read_group, read_input, refuse_validation, thrown, Args, ReferenceWalker,
};

const TOOL: &str = "CollectRrbsMetrics";
const UNKNOWN: &str = "unknown";
const UNMAPPED: u16 = 0x4;
const REVERSE: u16 = 0x10;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Level {
    AllReads,
    Sample,
    Library,
    ReadGroup,
}

/// `PerUnitRrbsMetricsCollector`.
#[derive(Default)]
struct Unit {
    sample: Option<String>,
    library: Option<String>,
    read_group: Option<String>,
    cyto_converted: i64,
    cyto_total: i64,
    /// `Histogram<CpgLocation>`, a `TreeMap` on (sequence, position).
    cpg_total: BTreeMap<(String, i32), i64>,
    cpg_converted: BTreeMap<(String, i32), i64>,
    mapped: i64,
    small: i64,
    mismatched: i64,
    no_cpg: i64,
}

fn opt(value: &Option<String>) -> Value {
    match value {
        Some(s) => Value::Str(s.clone()),
        None => Value::Null,
    }
}

struct Summary<'a>(&'a Unit);

impl MetricBean for Summary<'_> {
    fn class_name(&self) -> &str {
        "picard.analysis.RrbsSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "READS_ALIGNED",
            "NON_CPG_BASES",
            "NON_CPG_CONVERTED_BASES",
            "PCT_NON_CPG_BASES_CONVERTED",
            "CPG_BASES_SEEN",
            "CPG_BASES_CONVERTED",
            "PCT_CPG_BASES_CONVERTED",
            "MEAN_CPG_COVERAGE",
            "MEDIAN_CPG_COVERAGE",
            "READS_WITH_NO_CPG",
            "READS_IGNORED_SHORT",
            "READS_IGNORED_MISMATCHES",
            "SAMPLE",
            "LIBRARY",
            "READ_GROUP",
        ]
    }
    fn values(&self) -> Vec<Value> {
        let u = self.0;
        // finish().
        let cyto_rate = if u.cyto_total == 0 {
            0.0
        } else {
            u.cyto_converted as f64 / u.cyto_total as f64
        };
        let seen: i64 = u.cpg_total.values().sum();
        let converted: i64 = u.cpg_converted.values().sum();
        let seen = seen as i32;
        let converted = converted as i32;
        let cpg_rate = if seen == 0 {
            0.0
        } else {
            f64::from(converted) / f64::from(seen)
        };
        // Histogram.getMeanBinSize and getMedianBinSize.
        let mean = u.cpg_total.values().sum::<i64>() as f64 / u.cpg_total.len() as f64;
        let median = if u.cpg_total.is_empty() {
            0.0
        } else {
            let mut bins: Vec<f64> = u.cpg_total.values().map(|&v| v as f64).collect();
            bins.sort_by(f64::total_cmp);
            let mid = bins.len() / 2;
            let mut m = bins[mid];
            if bins.len().is_multiple_of(2) {
                m = (m + bins[mid - 1]) / 2.0;
            }
            m
        };
        vec![
            Value::Long(u.mapped),
            Value::Long(u.cyto_total),
            Value::Long(u.cyto_converted),
            Value::Double(cyto_rate),
            Value::Long(i64::from(seen)),
            Value::Long(i64::from(converted)),
            Value::Double(cpg_rate),
            Value::Double(mean),
            Value::Long(median as i32 as i64),
            Value::Long(u.no_cpg),
            Value::Long(u.small),
            Value::Long(u.mismatched),
            opt(&u.sample),
            opt(&u.library),
            opt(&u.read_group),
        ]
    }
}

struct Detail<'a> {
    unit: &'a Unit,
    sequence: &'a str,
    position: i32,
    total: i64,
    converted: i64,
}

impl MetricBean for Detail<'_> {
    fn class_name(&self) -> &str {
        "picard.analysis.RrbsCpgDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "SEQUENCE_NAME",
            "POSITION",
            "TOTAL_SITES",
            "CONVERTED_SITES",
            "PCT_CONVERTED",
            "SAMPLE",
            "LIBRARY",
            "READ_GROUP",
        ]
    }
    fn values(&self) -> Vec<Value> {
        let pct = if self.converted == 0 {
            0.0
        } else {
            self.converted as f64 / self.total as f64
        };
        vec![
            Value::Str(self.sequence.to_string()),
            Value::Long(i64::from(self.position)),
            Value::Long(self.total),
            Value::Long(self.converted),
            Value::Double(pct),
            opt(&self.unit.sample),
            opt(&self.unit.library),
            opt(&self.unit.read_group),
        ]
    }
}

/// One level's units, keyed the way `Distributor` keys them (`None` for a read group without the
/// attribute), in the order they were made.
type Units = Vec<(Option<String>, Unit)>;

struct Thresholds {
    min_read_length: i64,
    max_mismatch_rate: f64,
    c_quality: i64,
    next_quality: i64,
}

/// `isAboveCytoQcThreshold`.
fn above_threshold(qualities: &[u8], index: usize, t: &Thresholds) -> bool {
    index + 1 < qualities.len()
        && i64::from(qualities[index] as i8) >= t.c_quality
        && i64::from(qualities[index + 1] as i8) >= t.next_quality
}

/// `isC`: a C in the reference, read as itself or as its bisulfite conversion.
fn is_c(reference: u8, read: u8) -> bool {
    bases_equal(reference, b'C') && bisulfite_bases_equal(false, read, reference)
}

fn reverse_complement(bases: &mut [u8]) {
    bases.reverse();
    for b in bases.iter_mut() {
        *b = complement(*b);
    }
}

/// `PerUnitRrbsMetricsCollector.acceptRecord`.
fn accept(
    unit: &mut Unit,
    record: &htsjdk_bam::record::BamRecord,
    sequence: &str,
    reference: &[u8],
    t: &Thresholds,
) {
    unit.mapped += 1;
    let length = record.read_bases.len() as i64;
    let negative = record.flags & REVERSE != 0;
    let blocks = alignment_blocks(&record.cigar, record.alignment_start);
    if length < t.min_read_length {
        unit.small += 1;
        return;
    }
    let bound = (length as f64 * t.max_mismatch_rate + 0.5).floor() as i64;
    if i64::from(count_mismatches(
        &record.read_bases,
        &blocks,
        reference,
        0,
        negative,
        true,
    )) > bound
    {
        unit.mismatched += 1;
        return;
    }
    let read_qualities = &record.base_qualities;
    let mut record_cpgs = 0;
    for block in &blocks {
        let block_length = block.length as usize;
        let ref_start = (block.reference_start - 1) as usize;
        let read_start = (block.read_start - 1) as usize;
        let mut ref_fragment = reference[ref_start..ref_start + block_length].to_vec();
        let mut read_fragment = record.read_bases[read_start..read_start + block_length].to_vec();
        let mut quality_fragment = read_qualities[read_start..read_start + block_length].to_vec();
        if negative {
            reverse_complement(&mut ref_fragment);
            reverse_complement(&mut read_fragment);
            quality_fragment.reverse();
        }
        let mut i = 0usize;
        while i + 1 < block_length {
            let current = if negative {
                ref_start + (block_length - 1) - i - 1
            } else {
                ref_start + i
            } as i32;
            if bases_equal(ref_fragment[i], b'C') && bases_equal(ref_fragment[i + 1], b'G') {
                let valid = is_c(ref_fragment[i], read_fragment[i])
                    && bases_equal(ref_fragment[i + 1], read_fragment[i + 1])
                    && above_threshold(&quality_fragment, i, t);
                if valid {
                    record_cpgs += 1;
                    let key = (sequence.to_string(), current);
                    *unit.cpg_total.entry(key.clone()).or_insert(0) += 1;
                    if is_bisulfite_converted(read_fragment[i], ref_fragment[i], false) {
                        *unit.cpg_converted.entry(key).or_insert(0) += 1;
                    }
                }
                i += 1;
            } else if is_c(ref_fragment[i], read_fragment[i])
                && above_threshold(read_qualities, i, t)
                && bisulfite_bases_equal(false, read_fragment[i + 1], ref_fragment[i + 1])
            {
                unit.cyto_total += 1;
                if is_bisulfite_converted(read_fragment[i], ref_fragment[i], false) {
                    unit.cyto_converted += 1;
                }
            }
            i += 1;
        }
    }
    if record_cpgs == 0 {
        unit.no_cpg += 1;
    }
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("M", "METRICS_FILE_PREFIX"),
        ("R", "REFERENCE"),
        ("AS", "ASSUME_SORTED"),
        ("LEVEL", "METRIC_ACCUMULATION_LEVEL"),
    ]);
    let input = args.required("INPUT");
    let mut prefix = args.required("METRICS_FILE_PREFIX");
    let reference = args.required("REFERENCE");
    let t = Thresholds {
        min_read_length: args.int("MINIMUM_READ_LENGTH", 5),
        max_mismatch_rate: args.double("MAX_MISMATCH_RATE", 0.1),
        c_quality: args.int("C_QUALITY_THRESHOLD", 20),
        next_quality: args.int("NEXT_BASE_QUALITY_THRESHOLD", 10),
    };
    let mut sequence_names: JavaHashMap<()> = JavaHashMap::new();
    for name in args.all("SEQUENCE_NAMES") {
        sequence_names.put(&name, ());
    }
    let assume_sorted = args.bool("ASSUME_SORTED", false);
    let mut levels = Vec::new();
    for value in args.collection("METRIC_ACCUMULATION_LEVEL", &["ALL_READS"]) {
        let level = match value.as_str() {
            "ALL_READS" => Level::AllReads,
            "SAMPLE" => Level::Sample,
            "LIBRARY" => Level::Library,
            "READ_GROUP" => Level::ReadGroup,
            other => fail(&format!(
                "Argument 'METRIC_ACCUMULATION_LEVEL' cannot be set to '{other}'"
            )),
        };
        if !levels.contains(&level) {
            levels.push(level);
        }
    }

    // customCommandLineValidation.
    let mut messages = Vec::new();
    if t.max_mismatch_rate < 0.0 || t.max_mismatch_rate > 1.0 {
        messages.push("MAX_MISMATCH_RATE must be in the range of 0-1".to_string());
    }
    if t.c_quality < 0 {
        messages.push("C_QUALITY_THRESHOLD must be >= 0".to_string());
    }
    if t.next_quality < 0 {
        messages.push("NEXT_BASE_QUALITY_THRESHOLD must be >= 0".to_string());
    }
    if t.min_read_length <= 0 {
        messages.push("MINIMUM_READ_LENGTH must be > 0".to_string());
    }
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }

    // doWork.
    if !prefix.ends_with('.') {
        prefix.push('.');
    }
    let summary_out = format!("{prefix}rrbs_summary_metrics");
    let details_out = format!("{prefix}rrbs_detail_metrics");
    let (header, records) = read_input(&input);
    if !assume_sorted && header.attributes.get("SO") != Some("coordinate") {
        thrown(&format!(
            "picard.PicardException: The input file {} does not appear to be coordinate sorted",
            absolute(&input)
        ));
    }
    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| fail(&format!("{e:?}")));

    // MultiLevelCollector.setup: units keyed per level, in header order.
    let key_of = |level: Level, group: &ReadGroup| -> Option<String> {
        let attribute = match level {
            Level::AllReads => return None,
            Level::Sample => "SM",
            Level::Library => "LB",
            Level::ReadGroup => "PU",
        };
        group.attributes.get(attribute).map(str::to_string)
    };
    let make = |level: Level, group: Option<&ReadGroup>| -> Unit {
        let get = |a: &str| match group {
            Some(g) => g.attributes.get(a).map(str::to_string),
            None => Some(UNKNOWN.to_string()),
        };
        match level {
            Level::AllReads => Unit::default(),
            Level::Sample => Unit {
                sample: get("SM"),
                ..Unit::default()
            },
            Level::Library => Unit {
                sample: get("SM"),
                library: get("LB"),
                ..Unit::default()
            },
            Level::ReadGroup => Unit {
                sample: get("SM"),
                library: get("LB"),
                read_group: get("PU"),
                ..Unit::default()
            },
        }
    };
    let mut distributors: Vec<(Level, Units)> = Vec::new();
    for level in [
        Level::AllReads,
        Level::Sample,
        Level::Library,
        Level::ReadGroup,
    ] {
        if !levels.contains(&level) {
            continue;
        }
        let mut units: Vec<(Option<String>, Unit)> = Vec::new();
        if level == Level::AllReads {
            units.push((None, Unit::default()));
        } else {
            for group in &header.read_groups {
                let key = key_of(level, group);
                if !units.iter().any(|(k, _)| *k == key) {
                    units.push((key, make(level, Some(group))));
                }
            }
        }
        distributors.push((level, units));
    }

    let mut walker = ReferenceWalker::default();
    for record in &records {
        if record.flags & UNMAPPED != 0 {
            continue;
        }
        let sequence = header
            .sequences
            .get(record.reference_index.max(0) as usize)
            .map(|s| s.name.as_str())
            .unwrap_or("*");
        if !sequence_names.is_empty() && !sequence_names.contains_key(sequence) {
            continue;
        }
        walker.get(record.reference_index);
        let bases = &contigs[record.reference_index as usize].bases;
        let group = read_group(&header, record);
        for (level, units) in distributors.iter_mut() {
            let index = if *level == Level::AllReads {
                0
            } else {
                let key = group
                    .and_then(|g| key_of(*level, g))
                    .unwrap_or_else(|| UNKNOWN.to_string());
                match units.iter().position(|(k, _)| k.as_deref() == Some(&key)) {
                    Some(i) => i,
                    None => {
                        if key != UNKNOWN {
                            thrown(&format!(
                                "picard.PicardException: Could not find collector for {key}"
                            ));
                        }
                        units.push((Some(key), make(*level, None)));
                        units.len() - 1
                    }
                }
            };
            accept(&mut units[index].1, record, sequence, bases, &t);
        }
    }

    let mut summary_file = MetricsFile::new();
    let mut details_file = MetricsFile::new();
    for file in [&mut summary_file, &mut details_file] {
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
    }
    for (_, units) in &distributors {
        for (_, unit) in units {
            summary_file.add_metric(&Summary(unit));
            for ((sequence, position), &total) in &unit.cpg_total {
                details_file.add_metric(&Detail {
                    unit,
                    sequence,
                    position: *position,
                    total,
                    converted: unit
                        .cpg_converted
                        .get(&(sequence.clone(), *position))
                        .copied()
                        .unwrap_or(0),
                });
            }
        }
    }
    for (path, file) in [(&summary_out, &summary_file), (&details_out, &details_file)] {
        if let Err(e) = std::fs::write(path, file.write()) {
            fail(&format!("{e}"));
        }
    }
}
