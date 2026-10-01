//! `UmiAwareMarkDuplicatesWithMateCigar` as a runnable binary: the covering array's port side.
//!
//! `SimpleMarkDuplicatesWithMateCigar` with its duplicate sets split by UMI; the decisions are
//! [`picard_analysis::umi_aware_duplicates`]'s. What is here is the file around them: the
//! arguments, the sort-order refusal (made after `ASSUME_SORT_ORDER` has been written into the
//! header), the writer, which re-sorts what it is given into coordinate order because the tool
//! opens it with `presorted = false`, and the `UMI_METRICS_FILE`.
//!
//! `METRICS_FILE` is written because the tool requires it; its duplication metrics are not what
//! the covering array compares (the output BAM and the UMI metrics are), and they are written as
//! an empty table rather than as a guess.

use htsjdk_bam::writer::BamWriter;
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::metrics_cli::{sort_order, thrown, Args};
use picard_analysis::umi_aware_duplicates::{run, Scoring, UmiArgs, UmiMetricsRow};

const COLUMNS: [&str; 11] = [
    "LIBRARY",
    "MEAN_UMI_LENGTH",
    "OBSERVED_UNIQUE_UMIS",
    "INFERRED_UNIQUE_UMIS",
    "OBSERVED_BASE_ERRORS",
    "DUPLICATE_SETS_IGNORING_UMI",
    "DUPLICATE_SETS_WITH_UMI",
    "OBSERVED_UMI_ENTROPY",
    "INFERRED_UMI_ENTROPY",
    "UMI_BASE_QUALITIES",
    "PCT_UMI_WITH_N",
];

struct Row(UmiMetricsRow);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.sam.markduplicates.UmiMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let m = &self.0;
        vec![
            m.library.clone().map_or(Value::Null, Value::Str),
            Value::Double(m.mean_umi_length),
            Value::Long(m.observed_unique_umis),
            Value::Long(m.inferred_unique_umis),
            Value::Long(m.observed_base_errors),
            Value::Long(m.duplicate_sets_ignoring_umi),
            Value::Long(m.duplicate_sets_with_umi),
            Value::Double(m.observed_umi_entropy),
            Value::Double(m.inferred_umi_entropy),
            Value::Double(m.umi_base_qualities),
            Value::Double(m.pct_umi_with_n),
        ]
    }
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("M", "METRICS_FILE"),
        ("UMI_METRICS", "UMI_METRICS_FILE"),
        ("DS", "DUPLICATE_SCORING_STRATEGY"),
        ("ASO", "ASSUME_SORT_ORDER"),
        ("PG", "PROGRAM_RECORD_ID"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let metrics_file = args.required("METRICS_FILE");
    let umi_metrics_file = args.required("UMI_METRICS_FILE");

    let scoring = match args.get("DUPLICATE_SCORING_STRATEGY") {
        None | Some("SUM_OF_BASE_QUALITIES") => Scoring::SumOfBaseQualities,
        Some("TOTAL_MAPPED_REFERENCE_LENGTH") => Scoring::TotalMappedReferenceLength,
        Some("RANDOM") => Scoring::Random,
        Some(other) => picard_analysis::metrics_cli::fail(&format!(
            "Argument 'DUPLICATE_SCORING_STRATEGY' cannot be set to '{other}'"
        )),
    };
    let umi_args = UmiArgs {
        scoring,
        max_edit_distance_to_join: args.int("MAX_EDIT_DISTANCE_TO_JOIN", 1) as i32,
        umi_tag: args.get("UMI_TAG_NAME").unwrap_or("RX").to_string(),
        molecular_identifier_tag: args.get("MOLECULAR_IDENTIFIER_TAG").map(str::to_string),
        allow_missing_umis: args.bool("ALLOW_MISSING_UMIS", false),
        duplex_umi: args.bool("DUPLEX_UMI", false),
        remove_duplicates: args.bool("REMOVE_DUPLICATES", false),
        program_record_id: match args.get("PROGRAM_RECORD_ID") {
            Some(id) => Some(id.to_string()),
            None if args_say_null(&args, "PROGRAM_RECORD_ID") => None,
            None => Some("MarkDuplicates".to_string()),
        },
    };

    let (mut header, records) = picard_analysis::metrics_cli::read_input(&input);
    // `openInputs` writes the assumed order into the header before anything looks at it.
    if let Some(order) = args.get("ASSUME_SORT_ORDER") {
        header.set_sort_order(order);
    } else if args.bool("ASSUME_SORTED", false) {
        header.set_sort_order("coordinate");
    }
    if sort_order(&header) != "coordinate" {
        thrown("picard.PicardException: This program requires inputs in coordinate SortOrder");
    }

    let result = match run(&header, records, &umi_args) {
        Ok(result) => result,
        Err(failure) => thrown(&format!("{}: {}", failure.class, failure.message)),
    };

    // `SAMFileWriterFactory.makeWriter(header, false, ...)`: not presorted, so the writer sorts
    // what it is given by the header's order, stably.
    let mut written = result.written;
    written.sort_by(htsjdk_bam::coordinate::compare);
    let mut writer = BamWriter::new(Vec::new(), &header)
        .unwrap_or_else(|e| picard_analysis::metrics_cli::fail(&format!("{e}")));
    for record in &written {
        writer
            .write(record)
            .unwrap_or_else(|e| picard_analysis::metrics_cli::fail(&format!("{e:?}")));
    }
    let bytes = writer
        .finish()
        .unwrap_or_else(|e| picard_analysis::metrics_cli::fail(&format!("{e}")));
    let bytes = if output.ends_with(".sam") {
        let plain = htsjdk_bgzf::decompress_all(&bytes)
            .unwrap_or_else(|e| picard_analysis::metrics_cli::fail(&format!("{e:?}")));
        picard_analysis::sam_format_converter::bam_to_sam(&plain)
            .unwrap_or_else(|e| picard_analysis::metrics_cli::fail(&e))
            .into_bytes()
    } else {
        bytes
    };
    write_or_fail(&output, &bytes);

    let mut duplication = MetricsFile::new();
    duplication.add_header("UmiAwareMarkDuplicatesWithMateCigar <command line>");
    duplication.add_header("Started on: <timestamp>");
    write_or_fail(&metrics_file, duplication.write().as_bytes());

    let mut umi = MetricsFile::new();
    umi.add_header("UmiAwareMarkDuplicatesWithMateCigar <command line>");
    umi.add_header("Started on: <timestamp>");
    for row in result.metrics {
        umi.add_metric(&Row(row));
    }
    write_or_fail(&umi_metrics_file, umi.write().as_bytes());
}

/// Barclay reads `null` as "set to null", which for `PROGRAM_RECORD_ID` turns the PG rewrite off.
fn args_say_null(_args: &Args, name: &str) -> bool {
    std::env::args().any(|raw| {
        let raw = raw.trim_start_matches('-');
        raw == format!("{name}=null") || raw == "PG=null"
    })
}

fn write_or_fail(path: &str, bytes: &[u8]) {
    if let Err(e) = std::fs::write(path, bytes) {
        picard_analysis::metrics_cli::fail(&format!("{e}"));
    }
}
