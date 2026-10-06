//! `CollectQualityYieldMetricsFlow` as a runnable binary: the covering array's port side.
//!
//! The tool counts FLOWS, not bases: every read that passes the vendor filter is converted to its
//! flow key and flow matrix (`picard_analysis::flow_based`), each flow gets a quality from the
//! probability its called homopolymer length was wrong, and the metrics are the counts of those
//! at 20 and 30. A read whose read group is not a flow platform is refused, which is every read
//! of an Illumina file.
//!
//! The arguments are the three the tool declares, the two flow arguments of its collection
//! (`--flow-ignore-t0-tag`, `--flow-fill-empty-bins-value`) and `SinglePassSamProgram`'s
//! `ASSUME_SORTED` and `STOP_AFTER`.

use htsjdk_metrics::file::MetricsFile;
use picard_analysis::flow_based::FlowArguments;
use picard_analysis::metrics_cli::{read_input, Args};
use picard_analysis::quality_yield_flow::FlowCollector;
use picard_analysis::single_pass_driver;

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("R", "REFERENCE_SEQUENCE")]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let (header, records) = read_input(&input);

    let mut collector = FlowCollector::new(
        args.bool("INCLUDE_SECONDARY_ALIGNMENTS", false),
        args.bool("INCLUDE_SUPPLEMENTAL_ALIGNMENTS", false),
        args.bool("INCLUDE_BQ_HISTOGRAM", false),
        FlowArguments {
            ignore_t0_tag: args.bool("flow-ignore-t0-tag", false),
            filling_value: args.double("flow-fill-empty-bins-value", 0.0),
        },
    );
    single_pass_driver::run(&args, &input, &header, &records, |header, record| {
        collector.accept(header, record)
    });
    let (metrics, histograms) = collector.finish();

    let mut file = MetricsFile::new();
    file.add_header("CollectQualityYieldMetricsFlow <command line>");
    file.add_header("Started on: <timestamp>");
    file.add_metric(&metrics);
    file.histograms = histograms;
    if let Err(e) = std::fs::write(&output, file.write()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
