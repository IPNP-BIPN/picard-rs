//! `ExtractFingerprint` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.ExtractFingerprint` at tag 3.4.0; the run is
//! `picard_analysis::extract_run`.

use picard_analysis::extract_run::run;
use picard_analysis::metrics_cli::Args;

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("R", "REFERENCE_SEQUENCE"),
        ("C", "CONTAMINATION"),
    ]);
    run(
        &args,
        args.bool("EXTRACT_CONTAMINATION", false),
        args.int("LOCUS_MAX_READS", 50),
    );
}
