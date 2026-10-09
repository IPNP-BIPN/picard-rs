//! `AccumulateVariantCallingMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.AccumulateVariantCallingMetrics.doWork` at tag 3.4.0. Every `INPUT` and the
//! `OUTPUT` are PREFIXES: the tool reads `<prefix>.variant_calling_detail_metrics` and
//! `<prefix>.variant_calling_summary_metrics` of each input, merges them and writes the same two
//! files under the output prefix. The merge itself is `picard_analysis::
//! accumulate_variant_calling_metrics::tool`.

use htsjdk_metrics::file::MetricsFile;
use picard_analysis::accumulate_variant_calling_metrics::tool::{run, InputText};
use picard_analysis::accumulate_variant_calling_metrics::{DETAIL_EXTENSION, SUMMARY_EXTENSION};
use picard_analysis::metrics_cli::{absolute, thrown, Args};

fn read_text(path: &str) -> String {
    if !std::path::Path::new(path).is_file() {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{path}"
        ));
    }
    std::fs::read_to_string(path).unwrap_or_else(|e| {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read file {path}: {e}"
        ))
    })
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT")]);
    let prefixes = args.all("INPUT");
    if prefixes.is_empty() {
        eprintln!("Argument 'INPUT' is required");
        std::process::exit(1);
    }
    let output_prefix = format!("{}.", absolute(&args.required("OUTPUT")));
    let detail_path = format!("{output_prefix}{DETAIL_EXTENSION}");
    let summary_path = format!("{output_prefix}{SUMMARY_EXTENSION}");

    let mut inputs = Vec::new();
    for prefix in &prefixes {
        let stem = format!("{}.", absolute(prefix));
        let detail = read_text(&format!("{stem}{DETAIL_EXTENSION}"));
        let summary = read_text(&format!("{stem}{SUMMARY_EXTENSION}"));
        inputs.push(InputText { detail, summary });
    }
    let (details, summary) = run(&inputs).unwrap_or_else(|message| thrown(&message));

    let mut detail_file = MetricsFile::new();
    detail_file.add_header("AccumulateVariantCallingMetrics <command line>");
    detail_file.add_header("Started on: <timestamp>");
    for row in &details {
        detail_file.add_metric(row);
    }
    let mut summary_file = MetricsFile::new();
    summary_file.add_header("AccumulateVariantCallingMetrics <command line>");
    summary_file.add_header("Started on: <timestamp>");
    summary_file.add_metric(&summary);
    for (path, text) in [
        (&detail_path, detail_file.write()),
        (&summary_path, summary_file.write()),
    ] {
        if let Err(e) = std::fs::write(path, text) {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}
