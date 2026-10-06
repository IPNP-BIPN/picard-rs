//! `CollectQualityYieldMetricsSNVQ` as a runnable binary: the covering array's port side.
//!
//! The tool is `CollectQualityYieldMetrics` with a second tally: at every base it also counts the
//! qualities of the three bases the read does NOT have, which it reads from the tags `qa`, `qc`,
//! `qg` and `qt`. `ALTERNATE_QUALITY_ATTRIBUTE` names a tag of FASTQ-encoded qualities to count in
//! place of the base qualities, and `INCLUDE_BQ_HISTOGRAM` adds the quality histograms and the
//! per-position means.

use htsjdk_metrics::file::MetricsFile;
use picard_analysis::metrics_cli::{read_input, Args};
use picard_analysis::single_pass_driver;
use picard_analysis::snvq::SnvqCollector;

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("R", "REFERENCE_SEQUENCE"),
        ("AQA", "ALTERNATE_QUALITY_ATTRIBUTE"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let (header, records) = read_input(&input);

    let mut collector = SnvqCollector::with_options(
        args.get("ALTERNATE_QUALITY_ATTRIBUTE").map(str::to_string),
        args.bool("INCLUDE_SECONDARY_ALIGNMENTS", false),
        args.bool("INCLUDE_SUPPLEMENTAL_ALIGNMENTS", false),
        args.bool("INCLUDE_BQ_HISTOGRAM", false),
    );
    single_pass_driver::run(&args, &input, &header, &records, |_, record| {
        collector.try_accept(record)
    });
    let (metrics, histograms) = collector.finish_with_histograms();

    let mut file = MetricsFile::new();
    file.add_header("CollectQualityYieldMetricsSNVQ <command line>");
    file.add_header("Started on: <timestamp>");
    file.add_metric(&metrics);
    file.histograms = histograms;
    if let Err(e) = std::fs::write(&output, file.write()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
