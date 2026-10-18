//! `IdentifyContaminant` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.IdentifyContaminant` at tag 3.4.0, which builds an
//! `ExtractFingerprint` and runs it with `EXTRACT_CONTAMINATION = !EXTRACT_CONTAMINATED` and a
//! `LOCUS_MAX_READS` of 200 by default; the run is `picard_analysis::extract_run`.

use picard_analysis::extract_run::run;
use picard_analysis::metrics_cli::Args;

fn main() {
    picard_analysis::argument_check::check("IdentifyContaminant");
    picard_analysis::output::install();
    // The `ExtractFingerprint` it delegates to is handed the argument collection and not the
    // `REFERENCE_SEQUENCE` file its checker reads with, so a CRAM has no reference.
    picard_analysis::cram_input::set_reader_reference(None);
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("R", "REFERENCE_SEQUENCE"),
        ("C", "CONTAMINATION"),
    ]);
    run(
        &args,
        !args.bool("EXTRACT_CONTAMINATED", false),
        args.int("LOCUS_MAX_READS", 200),
    );
}
