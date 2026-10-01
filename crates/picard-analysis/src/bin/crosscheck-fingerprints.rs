//! `CrosscheckFingerprints` as a runnable binary: the covering array's port side.
//!
//! The tool lives in `picard_analysis::crosscheck_cli`, which its deprecated wrapper
//! `CrosscheckReadGroupFingerprints` shares.

fn main() {
    picard_analysis::crosscheck_cli::run(false);
}
