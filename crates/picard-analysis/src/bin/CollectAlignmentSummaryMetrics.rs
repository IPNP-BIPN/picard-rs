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

use htsjdk_bam::fasta::read_fasta_file;
use picard_analysis::alignment_summary::{LevelCollectors, Options};
use picard_analysis::args::arg;
use picard_analysis::insert_size::PairOrientation;

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CollectAlignmentSummaryMetrics");
    picard_analysis::output::install();
    let mut t = std::time::Instant::now();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let reference = arg(&args, "REFERENCE_SEQUENCE=").or_else(|| arg(&args, "R="));
    // Picard 3.5.0 refuses a chart where there is no R.
    picard_analysis::metrics_cli::refuse_chart_without_r(
        "CollectAlignmentSummaryMetrics",
        arg(&args, "HISTOGRAM_FILE=").is_some() || arg(&args, "H=").is_some(),
    );

    let mut options = Options::default();
    if let Some(value) = arg(&args, "IS_BISULFITE_SEQUENCED=") {
        options.is_bisulfite_sequenced = picard_analysis::args::is_true(&value);
    }
    if let Some(value) = arg(&args, "COLLECT_ALIGNMENT_INFORMATION=") {
        options.collect_alignment_information = picard_analysis::args::is_true(&value);
    }
    if let Some(value) = arg(&args, "MAX_INSERT_SIZE=") {
        options.max_insert_size = value.parse()?;
    }
    // Appended to the default like every collection argument here; `null` empties it first.
    for value in args_all(&args, "ADAPTER_SEQUENCE=") {
        if value == "null" {
            options.adapter_sequence.clear();
        } else {
            options.adapter_sequence.push(value);
        }
    }
    let stop_after: i64 = arg(&args, "STOP_AFTER=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(0);
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

    // Picard streams, and so does this: the input is checked whole first, so that a damaged
    // file is refused as reading it whole refused it, and its records are then decoded as the
    // collectors take them.
    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();
    t = phase("check_input", t);

    let contigs = match &reference {
        Some(path) => read_fasta_file(path).map_err(|e| format!("{e:?}"))?,
        None => Vec::new(),
    };
    t = phase("read_reference", t);

    let mut collectors = LevelCollectors::new(&options, &levels, &header);

    // `makeItSo` stops once `STOP_AFTER` records have been handed over; zero or less reads all.
    let limit = match usize::try_from(stop_after) {
        Ok(n) if n > 0 => n,
        _ => usize::MAX,
    };
    for rec in records.take(limit) {
        let rec = &rec;
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
        collectors.accept(rec, bases);
    }
    let file = collectors.metrics_file(&picard_analysis::command_line::default_headers(
        "CollectAlignmentSummaryMetrics",
    ))?;
    t = phase("decode_and_collect", t);

    std::fs::write(&output, file.write())?;
    phase("write_metrics", t);
    Ok(())
}
