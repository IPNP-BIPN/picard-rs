//! `CollectGcBiasMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectGcBiasMetrics` and `GcBiasMetricsCollector` at tag 3.4.0 around
//! the window arithmetic in `picard_analysis::gc`:
//!
//! * `SinglePassSamProgram.makeItSo`: the sort check, the reference walker, `STOP_AFTER`, and no
//!   stop at the unmapped reads (`usesNoRefReads` is true here), which still count as clusters;
//! * `MultiLevelCollector`: the levels in the fixed order all reads, sample, library, read group,
//!   whatever order they were asked in, `METRIC_ACCUMULATION_LEVEL` appending to its default, and
//!   each level's units in header order, keyed by sample, library and platform unit;
//! * one `GcObject` per unit, and a second one without duplicates under `ALSO_IGNORE_DUPLICATES`;
//!   a unit no aligned read reached writes nothing. Every record is counted: the collector has no
//!   filter, so secondary, supplementary and duplicate records all start a read;
//! * the per-contig GC array the collectors share, recomputed whenever the contig changes.
//!
//! The chart is R's and is not drawn. The reference ignores R's exit status, so a chart that fails
//! (an empty summary) does not fail the run, and the port's run does not either.

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{count_deleted_bases, count_inserted_bases, count_mismatches};
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::gc::{
    calculate_all_gcs, calculate_ref_windows_by_gc, dropout_metrics,
    phred_score_from_obs_and_errors, GcBiasDetailMetrics, GcBiasSummaryMetrics, GcObject, BINS,
};
use picard_analysis::metrics_cli::{
    check_coordinate_sorted, fail, read_group, read_input, thrown, Args, ReferenceWalker,
};

const TOOL: &str = "CollectGcBiasMetrics";
const UNKNOWN: &str = "unknown";

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const REVERSE: u16 = 0x10;
const FIRST: u16 = 0x40;
const DUPLICATE: u16 = 0x400;

/// `MetricAccumulationLevel`, in declaration order, which is the order `setup` builds them in.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Level {
    AllReads,
    Sample,
    Library,
    ReadGroup,
}

impl Level {
    /// `ACCUMULATION_LEVEL_*`.
    fn group(self) -> &'static str {
        match self {
            Level::AllReads => "All Reads",
            Level::Sample => "Sample",
            Level::Library => "Library",
            Level::ReadGroup => "Read Group",
        }
    }
}

/// `PerUnitGcBiasMetricsCollector`: one unit's counters, with and without duplicates.
struct Unit {
    level: Level,
    /// The `gcData` key: the read group's platform unit, library, sample, or `All_Reads`.
    name: String,
    all: GcObject,
    unique: Option<GcObject>,
    logged: usize,
}

struct Distributor {
    level: Level,
    /// `LinkedHashMap<String, collector>`: insertion order is header order.
    units: Vec<(String, Unit)>,
}

fn add_read(
    object: &mut GcObject,
    rec: &BamRecord,
    gc: &[i8],
    reference: &[u8],
    window: i32,
    bisulfite: bool,
) {
    if rec.flags & PAIRED == 0 || rec.flags & FIRST != 0 {
        object.total_clusters += 1;
    }
    let pos = if rec.flags & REVERSE != 0 {
        rec.alignment_end() - window
    } else {
        rec.alignment_start
    };
    object.total_aligned_reads += 1;
    if pos > 0 {
        let window_gc = gc[pos as usize];
        if window_gc >= 0 {
            let bin = window_gc as usize;
            object.reads_by_gc[bin] += 1;
            object.bases_by_gc[bin] += rec.read_bases.len() as i64;
            let blocks = alignment_blocks(&rec.cigar, rec.alignment_start);
            object.errors_by_gc[bin] += i64::from(count_mismatches(
                &rec.read_bases,
                &blocks,
                reference,
                0,
                rec.flags & REVERSE != 0,
                bisulfite,
            )) + i64::from(count_inserted_bases(&rec.cigar))
                + i64::from(count_deleted_bases(&rec.cigar));
        }
    }
}

fn clusters(object: &mut GcObject, rec: &BamRecord) {
    if rec.flags & PAIRED == 0 || rec.flags & FIRST != 0 {
        object.total_clusters += 1;
    }
}

/// `addGcDataToFile`: 101 detail rows and a summary, or nothing when no read aligned.
fn rows(
    object: &GcObject,
    level: Level,
    name: &str,
    reads_used: &str,
    windows_by_gc: &[i32],
    window: i32,
) -> Option<(Vec<GcBiasDetailMetrics>, GcBiasSummaryMetrics)> {
    if object.total_aligned_reads <= 0 {
        return None;
    }
    let total_windows: f64 = windows_by_gc.iter().map(|&w| w as f64).sum();
    let total_reads: f64 = object.reads_by_gc.iter().map(|&r| r as f64).sum();
    let mean = total_reads / total_windows;
    let who = |l: Level| (level == l).then(|| name.to_string());
    let mut details = Vec::with_capacity(BINS);
    for (i, &windows) in windows_by_gc.iter().enumerate() {
        let read_starts = i64::from(object.reads_by_gc[i]);
        let mut detail = GcBiasDetailMetrics {
            accumulation_level: level.group().to_string(),
            reads_used: reads_used.to_string(),
            gc: i as i32,
            windows,
            read_starts,
            mean_base_quality: 0,
            normalized_coverage: 0.0,
            error_bar_width: 0.0,
            sample: who(Level::Sample),
            library: who(Level::Library),
            read_group: who(Level::ReadGroup),
        };
        if object.errors_by_gc[i] > 0 {
            detail.mean_base_quality = phred_score_from_obs_and_errors(
                object.bases_by_gc[i] as f64,
                object.errors_by_gc[i] as f64,
            );
        }
        if windows != 0 {
            detail.normalized_coverage = (read_starts as f64 / windows as f64) / mean;
            detail.error_bar_width = ((read_starts as f64).sqrt() / windows as f64) / mean;
        }
        details.push(detail);
    }
    // calculateGcNormCoverage, whose window total is an int.
    let norm = |start: usize, end: usize| {
        let mut windows_total = 0i32;
        let mut sum = 0.0;
        for i in start..=end {
            if windows_by_gc[i] != 0 {
                sum += f64::from(object.reads_by_gc[i]);
                windows_total = windows_total.wrapping_add(windows_by_gc[i]);
            }
        }
        if windows_total == 0 {
            0.0
        } else {
            sum / (f64::from(windows_total) * mean)
        }
    };
    let (at_dropout, gc_dropout) = dropout_metrics(&details);
    let summary = GcBiasSummaryMetrics {
        accumulation_level: level.group().to_string(),
        reads_used: reads_used.to_string(),
        window_size: window,
        total_clusters: object.total_clusters,
        aligned_reads: object.total_aligned_reads,
        at_dropout,
        gc_dropout,
        gc_nc_0_19: norm(0, 19),
        gc_nc_20_39: norm(20, 39),
        gc_nc_40_59: norm(40, 59),
        gc_nc_60_79: norm(60, 79),
        gc_nc_80_100: norm(80, 100),
        sample: who(Level::Sample),
        library: who(Level::Library),
        read_group: who(Level::ReadGroup),
    };
    Some((details, summary))
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("R", "REFERENCE_SEQUENCE"),
        ("CHART", "CHART_OUTPUT"),
        ("S", "SUMMARY_OUTPUT"),
        ("WINDOW_SIZE", "SCAN_WINDOW_SIZE"),
        ("MGF", "MINIMUM_GENOME_FRACTION"),
        ("BS", "IS_BISULFITE_SEQUENCED"),
        ("LEVEL", "METRIC_ACCUMULATION_LEVEL"),
        ("AS", "ASSUME_SORTED"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let reference = args.required("REFERENCE_SEQUENCE");
    let _chart = args.required("CHART_OUTPUT");
    let summary_output = args.required("SUMMARY_OUTPUT");
    let window = args.int("SCAN_WINDOW_SIZE", 100) as i32;
    let _ = args.double("MINIMUM_GENOME_FRACTION", 0.00001);
    let bisulfite = args.bool("IS_BISULFITE_SEQUENCED", false);
    let ignore_duplicates = args.bool("ALSO_IGNORE_DUPLICATES", false);
    let assume_sorted = args.bool("ASSUME_SORTED", true);
    let stop_after = args.int("STOP_AFTER", 0);

    // A collection argument appends to its default; `null` empties it first.
    let asked = args.collection("METRIC_ACCUMULATION_LEVEL", &["ALL_READS"]);
    let mut levels = Vec::new();
    for value in &asked {
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

    let (header, records) = read_input(&input);
    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| fail(&format!("{e:?}")));
    check_coordinate_sorted(&input, &header, assume_sorted);

    // setup.
    let bases: Vec<Vec<u8>> = contigs.iter().map(|c| c.bases.clone()).collect();
    let windows_by_gc = calculate_ref_windows_by_gc(BINS, &bases, window.max(0) as usize);
    let new_unit = |level: Level, name: &str| Unit {
        level,
        name: name.to_string(),
        all: GcObject::default(),
        unique: ignore_duplicates.then(GcObject::default),
        logged: 0,
    };
    let key_of = |level: Level, group: &htsjdk_bam::header::ReadGroup| -> Option<String> {
        let attribute = match level {
            Level::AllReads => return None,
            Level::Sample => "SM",
            Level::Library => "LB",
            Level::ReadGroup => "PU",
        };
        group.attributes.get(attribute).map(str::to_string)
    };
    let mut distributors: Vec<Distributor> = Vec::new();
    for level in [
        Level::AllReads,
        Level::Sample,
        Level::Library,
        Level::ReadGroup,
    ] {
        if !levels.contains(&level) {
            continue;
        }
        let mut units: Vec<(String, Unit)> = Vec::new();
        if level == Level::AllReads {
            units.push((String::new(), new_unit(level, "All_Reads")));
        } else {
            for group in &header.read_groups {
                let key = key_of(level, group).unwrap_or_else(|| "null".to_string());
                if !units.iter().any(|(k, _)| *k == key) {
                    units.push((key.clone(), new_unit(level, &key)));
                }
            }
        }
        distributors.push(Distributor { level, units });
    }

    // The walk. The GC array and the upper-cased contig are the collector's, shared by its units.
    let mut walker = ReferenceWalker::default();
    let mut gc: Option<(i32, Vec<i8>, Vec<u8>)> = None;
    let mut count = 0i64;
    for record in &records {
        let mut reference_bases: Option<&[u8]> = None;
        if record.reference_index != -1 {
            walker.get(record.reference_index);
            reference_bases = Some(&contigs[record.reference_index as usize].bases);
        }
        let group = read_group(&header, record);
        for distributor in distributors.iter_mut() {
            let index = if distributor.level == Level::AllReads {
                0
            } else {
                let key = group
                    .and_then(|g| key_of(distributor.level, g))
                    .unwrap_or_else(|| UNKNOWN.to_string());
                match distributor.units.iter().position(|(k, _)| *k == key) {
                    Some(i) => i,
                    None => {
                        if key != UNKNOWN {
                            thrown(&format!(
                                "picard.PicardException: Could not find collector for {key}"
                            ));
                        }
                        let unit = new_unit(distributor.level, UNKNOWN);
                        distributor.units.push((key, unit));
                        distributor.units.len() - 1
                    }
                }
            };
            let unit = &mut distributor.units[index].1;
            if record.read_bases.is_empty() && unit.logged < 100 {
                unit.logged += 1;
                continue;
            }
            if record.flags & UNMAPPED == 0 {
                if gc
                    .as_ref()
                    .is_none_or(|(i, _, _)| *i != record.reference_index)
                {
                    let upper = match reference_bases {
                        Some(b) => b.to_ascii_uppercase(),
                        None => thrown("java.lang.NullPointerException"),
                    };
                    let last_window_start = upper.len() as i64 - i64::from(window);
                    let gcs = calculate_all_gcs(
                        &upper,
                        last_window_start.max(0) as usize,
                        window.max(0) as usize,
                    );
                    gc = Some((record.reference_index, gcs, upper));
                }
                let (_, gcs, upper) = gc.as_ref().expect("gc");
                add_read(&mut unit.all, record, gcs, upper, window, bisulfite);
                if let Some(unique) = unit.unique.as_mut() {
                    if record.flags & DUPLICATE == 0 {
                        add_read(unique, record, gcs, upper, window, bisulfite);
                    }
                }
            } else {
                clusters(&mut unit.all, record);
                if let Some(unique) = unit.unique.as_mut() {
                    if record.flags & DUPLICATE == 0 {
                        clusters(unique, record);
                    }
                }
            }
        }
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
    }

    // finish: every level's units, each with its ALL rows and then its UNIQUE ones.
    let mut details_file = MetricsFile::new();
    let mut summary_file = MetricsFile::new();
    for file in [&mut details_file, &mut summary_file] {
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
    }
    for distributor in &distributors {
        for (_, unit) in &distributor.units {
            let mut objects = vec![(&unit.all, "ALL")];
            if let Some(unique) = &unit.unique {
                objects.push((unique, "UNIQUE"));
            }
            for (object, used) in objects {
                if let Some((details, summary)) =
                    rows(object, unit.level, &unit.name, used, &windows_by_gc, window)
                {
                    for d in &details {
                        details_file.add_metric(d);
                    }
                    summary_file.add_metric(&summary);
                }
            }
        }
    }
    for (path, file) in [(&output, &details_file), (&summary_output, &summary_file)] {
        if let Err(e) = std::fs::write(path, file.write()) {
            fail(&format!("{e}"));
        }
    }
}
