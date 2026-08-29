//! Conformance for `CollectArraysVariantCallingMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/arraysmetrics-conformance/CollectArraysVariantCallingMetricsDump.java`,
//! nineteen runs over an Illumina array VCF that is mostly header.
//!
//! # What this suite is for
//!
//!  * **which header lines are required and which are not**;
//!  * **the three fates a filter can have, of which one is counted as a pass**;
//!  * **the autocall test, whose two spellings do opposite things**;
//!  * **the twenty-three control codes as a file of their own**;
//!  * **and the derived fields, each over its own denominator.**

use std::io::Read;

use picard_analysis::collect_arrays_variant_calling_metrics::{
    assay_counts, call_rate, file_names, het_pct, is_autocall, is_counted_as_non_filtered,
    is_zcalled, missing_header_line_message, parse_control_header, passes_call_rate,
    sex_concordance, ControlCode, CONTROL_CODES, DEFAULT_CALL_RATE_PF_THRESHOLD,
    OPTIONAL_HEADER_LINES, REQUIRED_HEADER_LINES, REQUIRED_HEADER_LINES_ALSO,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/collect_arrays_variant_calling_metrics.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn table(text: &str, kind: &str, case: &str) -> Vec<Vec<(String, String)>> {
    let prefix = format!("{kind}\t{case}\t");
    let body = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
        .unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = body.split('\n');
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .map(str::to_string)
        .collect();
    lines
        .filter(|line| !line.is_empty())
        .map(|line| {
            header
                .iter()
                .cloned()
                .zip(line.split('\t').map(str::to_string))
                .collect()
        })
        .collect()
}

fn value(text: &str, kind: &str, case: &str, column: &str) -> String {
    table(text, kind, case)
        .into_iter()
        .next()
        .unwrap_or_else(|| panic!("{kind}/{case} has a row"))
        .into_iter()
        .find(|(name, _)| name == column)
        .map(|(_, value)| value)
        .unwrap_or_else(|| panic!("{kind}/{case}/{column}"))
}

fn count(text: &str, case: &str, column: &str) -> i64 {
    value(text, "summary", case, column)
        .parse()
        .unwrap_or_else(|_| panic!("{case}/{column} is a number"))
}

fn refusal(text: &str, case: &str) -> Option<String> {
    let prefix = format!("error\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
}

/// A missing header line is a refusal naming the line, and only some lines are required.
#[test]
fn the_header_is_a_contract() {
    let text = corpus();
    let cases = [
        ("without-the-autocall-version", "autocallVersion"),
        ("without-the-cluster-file", "clusterFile"),
        ("without-p95-red", "p95Red"),
    ];
    for (case, line) in cases {
        let written = refusal(&text, case).unwrap_or_else(|| panic!("{case}"));
        assert!(
            written.ends_with(&missing_header_line_message(line)),
            "{written}"
        );
        assert!(
            REQUIRED_HEADER_LINES.contains(&line) || REQUIRED_HEADER_LINES_ALSO.contains(&line),
            "{line}"
        );
    }
    // And the optional ones are absences the run survives, with the same numbers.
    for case in [
        "without-zcall",
        "without-the-pipeline-version",
        "without-the-genders",
        "without-the-call-rate",
    ] {
        assert_eq!(refusal(&text, case), None, "{case}");
        assert_eq!(
            count(&text, case, "NUM_ASSAYS"),
            count(&text, "plain", "NUM_ASSAYS")
        );
    }
    for line in OPTIONAL_HEADER_LINES {
        assert!(!REQUIRED_HEADER_LINES.contains(&line), "{line}");
    }
    // A gender that is not there is NotReported rather than an error, and the concordance is a
    // vote of three: two unknowns and one female still pass.
    assert_eq!(
        value(&text, "detail", "without-the-genders", "REPORTED_GENDER"),
        "N"
    );
    assert_eq!(
        value(
            &text,
            "detail",
            "without-the-genders",
            "GENDER_CONCORDANCE_PF"
        ),
        "N"
    );
    assert!(!sex_concordance("N", "N", "F"));
    assert!(sex_concordance("F", "F", "U"));
    assert!(!sex_concordance("F", "M", "F"));
    assert!(!sex_concordance("U", "U", "U"));
    assert_eq!(
        value(&text, "detail", "plain", "GENDER_CONCORDANCE_PF"),
        "Y"
    );
    assert!(sex_concordance("F", "F", "F"));
}

/// A filter has three fates, and one of them is counted as a pass.
#[test]
fn a_duplicate_assay_is_counted_as_a_pass() {
    let text = corpus();
    assert!(is_counted_as_non_filtered(&[]));
    assert!(is_counted_as_non_filtered(&["DUPE"]));
    assert!(!is_counted_as_non_filtered(&["TRIALLELIC"]));
    let plain = assay_counts(&[], true, None, "0/1");
    assert_eq!((plain.non_filtered_assays, plain.filtered_assays), (1, 0));
    let duplicate = assay_counts(&["DUPE"], true, None, "0/1");
    assert_eq!(
        (duplicate.non_filtered_assays, duplicate.filtered_assays),
        (1, 0)
    );
    let filtered = assay_counts(&["TRIALLELIC"], true, None, "0/1");
    assert_eq!(
        (filtered.non_filtered_assays, filtered.filtered_assays),
        (0, 1)
    );
    let zeroed = assay_counts(&["ZEROED_OUT_ASSAY"], true, None, "0/1");
    assert_eq!(
        (zeroed.filtered_assays, zeroed.zeroed_out_assays),
        (1, 1),
        "a zeroed-out assay is filtered AND zeroed out"
    );
    // Which is what the golden's three cases say, over two assays each.
    assert_eq!(
        count(&text, "a-duplicate-assay", "NUM_NON_FILTERED_ASSAYS"),
        2
    );
    assert_eq!(count(&text, "a-duplicate-assay", "NUM_FILTERED_ASSAYS"), 0);
    assert_eq!(count(&text, "a-filtered-assay", "NUM_FILTERED_ASSAYS"), 1);
    assert_eq!(count(&text, "a-filtered-assay", "NUM_ZEROED_OUT_ASSAYS"), 0);
    assert_eq!(count(&text, "a-zeroed-out-assay", "NUM_FILTERED_ASSAYS"), 1);
    assert_eq!(
        count(&text, "a-zeroed-out-assay", "NUM_ZEROED_OUT_ASSAYS"),
        1
    );
}

/// The autocall test's two spellings do opposite things.
#[test]
fn the_autocall_test_compares_against_a_no_call_genotype() {
    let text = corpus();
    // No GTA at all: the default is the genotype's own string, which is a call.
    assert!(is_autocall(None, "0/1"));
    assert_eq!(count(&text, "an-autocall", "NUM_AUTOCALL_CALLS"), 1);
    // A single dot is a MISSING attribute, so the default applies and it is an autocall.
    assert_eq!(
        count(
            &text,
            "a-call-whose-autocall-is-empty",
            "NUM_AUTOCALL_CALLS"
        ),
        1
    );
    // And `./.` is the value the test refuses, so it is not.
    assert!(!is_autocall(Some("./."), "0/1"));
    assert_eq!(
        count(
            &text,
            "a-call-whose-autocall-is-a-no-call-genotype",
            "NUM_AUTOCALL_CALLS"
        ),
        0
    );
    assert_eq!(
        count(
            &text,
            "a-call-whose-autocall-is-a-no-call-genotype",
            "NUM_CALLS"
        ),
        1,
        "it is still a call"
    );
}

/// The control codes are a file of their own, parsed out of the header.
#[test]
fn the_control_codes_are_a_file_of_their_own() {
    let text = corpus();
    let rows = table(&text, "controls", "plain");
    assert_eq!(rows.len(), CONTROL_CODES.len());
    for (row, (control, category)) in rows.iter().zip(CONTROL_CODES) {
        let field = |name: &str| {
            row.iter()
                .find(|(column, _)| column == name)
                .map(|(_, value)| value.clone())
                .unwrap_or_else(|| panic!("{name}"))
        };
        assert_eq!(field("CONTROL"), control);
        assert_eq!(field("CATEGORY"), category);
    }
    // The value the header carries is four fields separated by a pipe.
    assert_eq!(
        parse_control_header("DNP(High)|Staining|10|20"),
        Some(ControlCode {
            control: "DNP(High)".to_string(),
            category: "Staining".to_string(),
            red: 10,
            green: 20,
        })
    );
    assert_eq!(parse_control_header("DNP(High)|Staining|10"), None);
    // And the three files are one prefix, a dot and their own extensions.
    let mut written: Vec<String> = text
        .lines()
        .find(|line| line.starts_with("files\tplain\t"))
        .expect("the files")
        .split('\t')
        .nth(2)
        .expect("the names")
        .split(' ')
        .map(str::to_string)
        .collect();
    written.sort();
    let mut expected = file_names("m").to_vec();
    expected.sort();
    assert_eq!(written, expected);
}

/// The derived fields, each over its own denominator.
#[test]
fn each_rate_is_over_its_own_denominator() {
    let text = corpus();
    // Four assays, none filtered, three called and one het.
    assert_eq!(count(&text, "plain", "NUM_NON_FILTERED_ASSAYS"), 4);
    assert_eq!(count(&text, "plain", "NUM_CALLS"), 3);
    let written: f64 = value(&text, "detail", "plain", "CALL_RATE")
        .parse()
        .expect("a number");
    assert_eq!(call_rate(3, 4), written);
    let het: f64 = value(&text, "detail", "plain", "HET_PCT")
        .parse()
        .expect("a number");
    // The file rounds to six places, so the comparison is against the rounding and not the double.
    assert!((het_pct(1, 3) - het).abs() < 1e-6);
    // The threshold is against the AUTOCALL call rate, which is not the rate the header reported.
    assert_eq!(value(&text, "detail", "plain", "GTC_CALL_RATE"), "0.995");
    assert_eq!(
        value(&text, "detail", "plain", "AUTOCALL_CALL_RATE"),
        "0.75"
    );
    assert_eq!(value(&text, "detail", "plain", "AUTOCALL_PF"), "N");
    assert!(!passes_call_rate(0.75, DEFAULT_CALL_RATE_PF_THRESHOLD));
    assert_eq!(
        value(
            &text,
            "detail",
            "over-the-call-rate-threshold",
            "AUTOCALL_PF"
        ),
        "Y"
    );
    assert!(passes_call_rate(0.75, 0.5));
    assert_eq!(
        value(
            &text,
            "detail",
            "under-the-call-rate-threshold",
            "AUTOCALL_PF"
        ),
        "N"
    );
    assert!(!passes_call_rate(0.75, 0.999));
    // And the thresholds file, not the version, is what makes a run zcalled.
    assert_eq!(value(&text, "detail", "plain", "IS_ZCALLED"), "Y");
    assert!(is_zcalled(Some("thresholds.7.txt")));
    assert!(!is_zcalled(None));
    assert_eq!(value(&text, "detail", "without-zcall", "IS_ZCALLED"), "N");
}

/// Two processors do not change a number, which is what makes the multithreading safe to golden.
#[test]
fn the_number_of_processors_changes_nothing() {
    let text = corpus();
    for column in [
        "NUM_ASSAYS",
        "NUM_CALLS",
        "NUM_AUTOCALL_CALLS",
        "NUM_IN_DB_SNP",
        "NOVEL_SNPS",
    ] {
        assert_eq!(
            count(&text, "two-processors", column),
            count(&text, "plain", column),
            "{column}"
        );
    }
    // And the dbSNP file is what NUM_IN_DB_SNP counts against, so a dbSNP that knows nothing
    // moves every SNP into NOVEL_SNPS.
    assert_eq!(count(&text, "plain", "NUM_IN_DB_SNP"), 2);
    assert_eq!(count(&text, "nothing-in-dbsnp", "NUM_IN_DB_SNP"), 0);
    assert_eq!(
        count(&text, "nothing-in-dbsnp", "NOVEL_SNPS"),
        count(&text, "nothing-in-dbsnp", "NUM_SNPS")
    );
}
