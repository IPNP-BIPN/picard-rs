//! `CollectTargetedPcrMetrics` as a runnable binary: the covering array's port side.
//!
//! The command line and the collector are shared with the other targeted tool; see
//! `picard_analysis::targeted_cli` and `picard_analysis::targeted_metrics`.

fn main() {
    picard_analysis::targeted_cli::main_for(picard_analysis::targeted_cli::Tool::TargetedPcr);
}
