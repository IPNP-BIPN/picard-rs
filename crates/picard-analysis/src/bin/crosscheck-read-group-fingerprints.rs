//! `CrosscheckReadGroupFingerprints` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CrosscheckReadGroupFingerprints` at tag 3.4.0: its parent with the
//! roll-up moved into two booleans. Rolling up to samples or libraries writes the MATRIX to the
//! file named by `OUTPUT` and sends the table to `/dev/null`
//! (`picard_analysis::crosscheck_read_group_fingerprints::destination`).

use picard_analysis::crosscheck_read_group_fingerprints::{destination, Options as Rollup};
use picard_analysis::crosscheck_run::{run, write_outcome, DataType, Options};
use picard_analysis::metrics_cli::{thrown, Args};

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("LOD", "LOD_THRESHOLD"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let rollup = Rollup {
        crosscheck_samples: args.bool("CROSSCHECK_SAMPLES", false),
        crosscheck_libraries: args.bool("CROSSCHECK_LIBRARIES", false),
        expect_all_read_groups_to_match: args.bool("EXPECT_ALL_READ_GROUPS_TO_MATCH", false),
        ..Rollup::default()
    };
    let output = args.get("OUTPUT").unwrap_or("/dev/stdout").to_string();
    let to = destination(&output, &rollup);
    let crosscheck_by = match to.crosscheck_by {
        picard_analysis::crosscheck_fingerprints::DataType::Library => DataType::Library,
        picard_analysis::crosscheck_fingerprints::DataType::Sample => DataType::Sample,
        picard_analysis::crosscheck_fingerprints::DataType::File => DataType::File,
        _ => DataType::ReadGroup,
    };
    let options = Options {
        tool: "CrosscheckReadGroupFingerprints".to_string(),
        inputs: args.collection("INPUT", &[]),
        second_inputs: Vec::new(),
        haplotype_map: args.required("HAPLOTYPE_MAP"),
        output: Some(to.output),
        matrix_output: to.matrix_output,
        sample_individual_map: args.get("SAMPLE_INDIVIDUAL_MAP").map(str::to_string),
        check_all_others: args.get("CROSSCHECK_MODE") == Some("CHECK_ALL_OTHERS"),
        lod_threshold: args.double("LOD_THRESHOLD", 0.0),
        crosscheck_by,
        tumor_aware: args.bool("CALCULATE_TUMOR_AWARE_RESULTS", true),
        allow_duplicate_reads: args.bool("ALLOW_DUPLICATE_READS", false),
        output_errors_only: args.bool("OUTPUT_ERRORS_ONLY", false),
        loss_of_het_rate: args.double("LOSS_OF_HET_RATE", 0.5),
        expect_all_groups_to_match: args.bool("EXPECT_ALL_READ_GROUPS_TO_MATCH", false)
            || args.bool("EXPECT_ALL_GROUPS_TO_MATCH", false),
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
