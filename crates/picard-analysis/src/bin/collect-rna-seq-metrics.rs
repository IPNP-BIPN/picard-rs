//! `CollectRnaSeqMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectRnaSeqMetrics` at tag 3.4.0 around the per-unit collector in
//! `picard_analysis::rnaseq_metrics`:
//!
//! * `customCommandLineValidation`: no ribosomal intervals and a fragment percentage of nought is
//!   refused with a `PicardException` thrown out of the parser;
//! * `SinglePassSamProgram.makeItSo`: the sort check, the optional reference walker that refuses
//!   to rewind, `STOP_AFTER`, and no stop at the unmapped reads (`usesNoRefReads` is true);
//! * `setup`: the refFlat read against the input's dictionary, the ribosomal intervals checked
//!   against it (`Sequence dictionaries differ`), uniqued, and `IGNORE_SEQUENCE` resolved to
//!   indices (`Unrecognized sequence`), in that order;
//! * `MultiLevelCollector`: the levels in the fixed order all reads, sample, library, read group,
//!   each level's units in header order, keyed by sample, library and platform unit, and an
//!   `unknown` unit made the first time a record names none;
//! * `finish`: every unit's row and its normalized-coverage histogram, in unit order.
//!
//! The chart is R's and is not drawn.

use std::collections::HashSet;

use htsjdk_bam::header::ReadGroup;
use htsjdk_bam::interval::{Interval, IntervalList};
use htsjdk_bam::overlap::OverlapDetector;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::metrics_cli::{
    absolute, check_coordinate_sorted, fail, read_group, read_input, thrown, Args, ReferenceWalker,
};
use picard_analysis::refflat;
use picard_analysis::rnaseq_metrics::{RnaSeqMetricsCollector, StrandSpecificity};

const TOOL: &str = "CollectRnaSeqMetrics";
const UNKNOWN: &str = "unknown";

/// One level's units, a `LinkedHashMap` from the unit's key to its collector.
type Units<'a> = Vec<(Option<String>, RnaSeqMetricsCollector<'a>)>;

/// `MetricAccumulationLevel`, in declaration order, which is the order `setup` builds them in.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Level {
    AllReads,
    Sample,
    Library,
    ReadGroup,
}

/// The `@SQ` names and lengths of a SAM-style header, in order.
fn dictionary_of(text: &str) -> Vec<(String, i64)> {
    text.lines()
        .filter(|l| l.starts_with("@SQ"))
        .map(|l| {
            let mut name = String::new();
            let mut length = 0;
            for field in l.split('\t') {
                if let Some(v) = field.strip_prefix("SN:") {
                    name = v.to_string();
                } else if let Some(v) = field.strip_prefix("LN:") {
                    length = v.parse().unwrap_or(0);
                }
            }
            (name, length)
        })
        .collect()
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("R", "REFERENCE_SEQUENCE"),
        ("STRAND", "STRAND_SPECIFICITY"),
        ("CHART", "CHART_OUTPUT"),
        ("LEVEL", "METRIC_ACCUMULATION_LEVEL"),
        ("AS", "ASSUME_SORTED"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let ref_flat = args.required("REF_FLAT");
    let reference = args.get("REFERENCE_SEQUENCE").map(str::to_string);
    let ribosomal_path = args.get("RIBOSOMAL_INTERVALS").map(str::to_string);
    let strand = match args.required("STRAND_SPECIFICITY").as_str() {
        "NONE" => StrandSpecificity::None,
        "FIRST_READ_TRANSCRIPTION_STRAND" => StrandSpecificity::FirstReadTranscriptionStrand,
        "SECOND_READ_TRANSCRIPTION_STRAND" => StrandSpecificity::SecondReadTranscriptionStrand,
        other => fail(&format!(
            "Argument 'STRAND_SPECIFICITY' cannot be set to '{other}'"
        )),
    };
    let minimum_length = args.int("MINIMUM_LENGTH", 500) as i32;
    let ignore_sequence = args.all("IGNORE_SEQUENCE");
    let rrna_fragment_percentage = args.double("RRNA_FRAGMENT_PERCENTAGE", 0.8);
    let end_bias_bases = args.int("END_BIAS_BASES", 100) as i32;
    let assume_sorted = args.bool("ASSUME_SORTED", true);
    let stop_after = args.int("STOP_AFTER", 0);

    // A collection argument appends to its default; `null` empties it first.
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
    if ribosomal_path.is_none() && rrna_fragment_percentage == 0.0 {
        thrown(
            "picard.PicardException: Must use a RIBOSOMAL_INTERVALS file if \
             RRNA_FRAGMENT_PERCENTAGE = 0.0",
        );
    }

    // makeItSo.
    let (header, records) = read_input(&input);
    let mut walker = reference.as_ref().map(|_| ReferenceWalker::default());
    check_coordinate_sorted(&input, &header, assume_sorted);

    // setup.
    let seq_names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let recognized: HashSet<&str> = seq_names.iter().map(String::as_str).collect();
    let ref_flat_text =
        std::fs::read_to_string(&ref_flat).unwrap_or_else(|e| fail(&format!("{e}")));
    let genes = refflat::load(&ref_flat_text, |c| recognized.contains(c))
        .unwrap_or_else(|e| fail(&format!("{e:?}")));

    let mut ribosomal: OverlapDetector<Interval> = OverlapDetector::new(0, 0);
    if let Some(path) = &ribosomal_path {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
        let theirs = dictionary_of(&text);
        let ours: Vec<(String, i64)> = header
            .sequences
            .iter()
            .map(|s| (s.name.clone(), i64::from(s.length)))
            .collect();
        if ours != theirs {
            thrown(&format!(
                "picard.PicardException: Sequence dictionaries differ in {} and {}",
                absolute(&input),
                absolute(path)
            ));
        }
        let names: Vec<String> = theirs.into_iter().map(|(n, _)| n).collect();
        let list = IntervalList::parse_body(names, &text)
            .unwrap_or_else(|e| fail(&format!("{e:?}")))
            .uniqued(true);
        for interval in list.intervals {
            let (contig, start, end) = (interval.contig.clone(), interval.start, interval.end);
            ribosomal.add(&contig, start, end, interval);
        }
    }
    let ribosomal_initial = ribosomal_path.as_ref().map(|_| 0i64);

    let mut ignored: Vec<i32> = Vec::new();
    for name in &ignore_sequence {
        match seq_names.iter().position(|s| s == name) {
            Some(i) => ignored.push(i as i32),
            None => thrown(&format!(
                "picard.PicardException: Unrecognized sequence {name} passed as argument to \
                 IGNORE_SEQUENCE"
            )),
        }
    }

    let new_collector = |sample: Option<&str>, library: Option<&str>, read_group: Option<&str>| {
        RnaSeqMetricsCollector::new(
            &seq_names,
            &genes,
            &ribosomal,
            ribosomal_initial,
            &ignored,
            minimum_length,
            strand,
            rrna_fragment_percentage,
            end_bias_bases,
        )
        .for_unit(
            sample.map(str::to_string),
            library.map(str::to_string),
            read_group.map(str::to_string),
        )
    };
    let key_of = |level: Level, group: &ReadGroup| -> Option<String> {
        let attribute = match level {
            Level::AllReads => return None,
            Level::Sample => "SM",
            Level::Library => "LB",
            Level::ReadGroup => "PU",
        };
        group.attributes.get(attribute).map(str::to_string)
    };
    let unit_for = |level: Level, group: &ReadGroup| {
        let sm = group.attributes.get("SM");
        let lb = group.attributes.get("LB");
        let pu = group.attributes.get("PU");
        match level {
            Level::AllReads => new_collector(None, None, None),
            Level::Sample => new_collector(sm, None, None),
            Level::Library => new_collector(sm, lb, None),
            Level::ReadGroup => new_collector(sm, lb, pu),
        }
    };
    let unknown_for = |level: Level| match level {
        Level::AllReads => unreachable!("the all-reads level has no unknown unit"),
        Level::Sample => new_collector(Some(UNKNOWN), None, None),
        Level::Library => new_collector(Some(UNKNOWN), Some(UNKNOWN), None),
        Level::ReadGroup => new_collector(Some(UNKNOWN), Some(UNKNOWN), Some(UNKNOWN)),
    };

    // One distributor per level, its units in a `LinkedHashMap` keyed as the level keys them; a
    // read group without the attribute is the `null` key, which no record ever reaches.
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
        let mut units: Units = Vec::new();
        if level == Level::AllReads {
            units.push((None, new_collector(None, None, None)));
        } else {
            for group in &header.read_groups {
                let key = key_of(level, group);
                if !units.iter().any(|(k, _)| *k == key) {
                    units.push((key, unit_for(level, group)));
                }
            }
        }
        distributors.push((level, units));
    }

    let mut count = 0i64;
    for record in &records {
        if let Some(walker) = walker.as_mut() {
            if record.reference_index != -1 {
                walker.get(record.reference_index);
            }
        }
        let group = read_group(&header, record);
        for (level, units) in distributors.iter_mut() {
            let index = if *level == Level::AllReads {
                0
            } else {
                let key = group
                    .and_then(|g| key_of(*level, g))
                    .unwrap_or_else(|| UNKNOWN.to_string());
                match units
                    .iter()
                    .position(|(k, _)| k.as_deref() == Some(key.as_str()))
                {
                    Some(i) => i,
                    None => {
                        if key != UNKNOWN {
                            thrown(&format!(
                                "picard.PicardException: Could not find collector for {key}"
                            ));
                        }
                        units.push((Some(key), unknown_for(*level)));
                        units.len() - 1
                    }
                }
            };
            units[index].1.accept(record);
        }
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
    }

    // finish, then addAllLevelsToFile.
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for (_, units) in distributors {
        for (_, collector) in units {
            let (metrics, histogram) = collector.finish();
            file.add_metric(&metrics);
            file.histograms.push(histogram);
        }
    }
    if let Err(e) = std::fs::write(&output, file.write()) {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Could not write to file {output}: {e}"
        ));
    }
}
