//! `CrosscheckReadGroupFingerprints` as a runnable binary: the covering array's port side.
//!
//! The deprecated wrapper around `CrosscheckFingerprints`, in `picard_analysis::crosscheck_cli`.

fn main() {
    picard_analysis::crosscheck_cli::run(true);
}
