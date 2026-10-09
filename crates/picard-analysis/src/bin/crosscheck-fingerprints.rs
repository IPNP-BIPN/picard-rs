//! `CrosscheckFingerprints` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CrosscheckFingerprints.doWork` at tag 3.4.0; the run is
//! `picard_analysis::crosscheck_run`, and this is its argument surface and its writer.

use picard_analysis::crosscheck_run::{run, write_outcome, DataType, Options};
use picard_analysis::metrics_cli::{thrown, Args};

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("SI", "SECOND_INPUT"),
        ("O", "OUTPUT"),
        ("MO", "MATRIX_OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("LOD", "LOD_THRESHOLD"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let options = Options {
        tool: "CrosscheckFingerprints".to_string(),
        inputs: args.collection("INPUT", &[]),
        second_inputs: args.collection("SECOND_INPUT", &[]),
        haplotype_map: args.required("HAPLOTYPE_MAP"),
        output: args.get("OUTPUT").map(str::to_string),
        matrix_output: args.get("MATRIX_OUTPUT").map(str::to_string),
        sample_individual_map: args.get("SAMPLE_INDIVIDUAL_MAP").map(str::to_string),
        check_all_others: args.get("CROSSCHECK_MODE") == Some("CHECK_ALL_OTHERS"),
        lod_threshold: args.double("LOD_THRESHOLD", 0.0),
        crosscheck_by: DataType::parse(args.get("CROSSCHECK_BY").unwrap_or("READGROUP"))
            .unwrap_or(DataType::ReadGroup),
        tumor_aware: args.bool("CALCULATE_TUMOR_AWARE_RESULTS", true),
        allow_duplicate_reads: args.bool("ALLOW_DUPLICATE_READS", false),
        output_errors_only: args.bool("OUTPUT_ERRORS_ONLY", false),
        loss_of_het_rate: args.double("LOSS_OF_HET_RATE", 0.5),
        expect_all_groups_to_match: args.bool("EXPECT_ALL_GROUPS_TO_MATCH", false),
        exit_code_when_mismatch: args.int("EXIT_CODE_WHEN_MISMATCH", 1) as i32,
        exit_code_when_no_valid_checks: args.int("EXIT_CODE_WHEN_NO_VALID_CHECKS", 1) as i32,
        max_effect_of_each_haplotype_block: args.double("MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK", 3.0),
        require_index_files: args.bool("REQUIRE_INDEX_FILES", false),
        strict: !matches!(
            args.get("VALIDATION_STRINGENCY"),
            Some("LENIENT") | Some("SILENT")
        ),
    };
    let outcome = run(&options).unwrap_or_else(|e| thrown(&e));
    std::process::exit(write_outcome(&options, &outcome));
}
