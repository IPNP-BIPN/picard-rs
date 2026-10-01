//! `CollectWgsMetrics` as a runnable binary: the covering array's port side.
//!
//! The command line, the walk and the metrics are shared with the other two WGS tools; see
//! `picard_analysis::wgs_cli` and `picard_analysis::wgs_walk`.

fn main() {
    picard_analysis::wgs_cli::main_for(picard_analysis::wgs_cli::Tool::Wgs);
}
