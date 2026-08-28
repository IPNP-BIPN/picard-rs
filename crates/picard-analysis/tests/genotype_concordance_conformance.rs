//! Conformance for `GenotypeConcordance` against Picard 3.4.0.
//!
//! Golden from `tools/genotypeconcordance-conformance`. Each case carries the two VCFs, the three
//! files' names, the summary and contingency tables, and the detail rows that counted something.
//!
//! # What this suite is for
//!
//!  * **the output argument being a basename for three files**;
//!  * **a low quality, a low depth, a filter and a missing site being STATES**;
//!  * **`--IGNORE_FILTER_STATUS` reading a filtered call as it stands**;
//!  * **`--MISSING_SITES_HOM_REF` needing an interval list, and changing the scheme**;
//!  * **every detail row's contingency being the scheme's cell for that pair**;
//!  * **and `--OUTPUT_ALL_ROWS` keeping the pairs nothing was seen for.**

use std::io::Read;

use picard_analysis::genotype_concordance::{
    call_state, contingency, contingency_values, file_names, truth_state, CallState, Cell,
    ContingencyState, TruthState, CONTINGENCY_METRICS_FILE_EXTENSION,
    DETAILED_METRICS_FILE_EXTENSION, GA4GH, GA4GH_MISSING_AS_HOM_REF,
    SUMMARY_METRICS_FILE_EXTENSION, TRUTH_ORDER,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("genotype_concordance.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn table(text: &str, kind: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let body = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = body.lines().filter(|line| !line.is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

fn truth(name: &str) -> TruthState {
    truth_state(name).unwrap_or_else(|| panic!("{name}"))
}

fn call(name: &str) -> CallState {
    call_state(name).unwrap_or_else(|| panic!("{name}"))
}

/// The output argument is a basename for three files.
#[test]
fn the_output_is_a_basename_for_three_files() {
    let text = corpus();
    let written = field(&text, "files", "plain").expect("the file names");
    let names: std::collections::BTreeSet<String> =
        written.split(',').map(str::to_string).collect();
    assert_eq!(names.len(), 3);
    assert_eq!(
        file_names("out")
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>(),
        names
    );
    assert!(SUMMARY_METRICS_FILE_EXTENSION.ends_with("summary_metrics"));
    assert!(DETAILED_METRICS_FILE_EXTENSION.ends_with("detail_metrics"));
    assert!(CONTINGENCY_METRICS_FILE_EXTENSION.ends_with("contingency_metrics"));
}

/// Every detail row's contingency is the scheme's cell for that pair.
#[test]
fn every_row_carries_the_schemes_cell() {
    let text = corpus();
    let mut checked = 0;
    for case in [
        "plain",
        "min-gq",
        "min-dp",
        "ignore-filters",
        "a-missing-site",
    ] {
        for row in table(&text, "detail", case) {
            let truth = truth(&row["TRUTH_STATE"]);
            let call = call(&row["CALL_STATE"]);
            let cell = contingency(&GA4GH, call, truth).expect("a cell");
            assert_eq!(
                contingency_values(cell).as_deref(),
                Some(row["CONTINGENCY_VALUES"].as_str()),
                "{case}: {truth:?} against {call:?}"
            );
            checked += 1;
        }
    }
    assert!(checked >= 20, "{checked} rows");
    // A pair may contribute to more than one counter, and the golden shows three kinds.
    assert_eq!(
        contingency_values(
            contingency(&GA4GH, CallState::HetRefVar1, TruthState::HetRefVar1).expect("a cell")
        ),
        Some("TP,TN".to_string())
    );
    assert_eq!(
        contingency_values(
            contingency(&GA4GH, CallState::HomVar1, TruthState::HetRefVar1).expect("a cell")
        ),
        Some("TP,FP".to_string())
    );
    assert_eq!(
        contingency_values(
            contingency(&GA4GH, CallState::VcFiltered, TruthState::HetRefVar1).expect("a cell")
        ),
        Some("TN,FN".to_string())
    );
    // And some pairs the reference says its own code cannot reach.
    assert_eq!(
        contingency(&GA4GH, CallState::HetRefVar2, TruthState::Missing),
        Some(Cell::Unreachable)
    );
    assert_eq!(contingency_values(Cell::Unreachable), None);
}

/// A low quality and a low depth are states, so the row count goes up and not down.
#[test]
fn the_floors_move_a_row_rather_than_dropping_it() {
    let text = corpus();
    let rows = |case: &str| -> usize {
        field(&text, "rows", case)
            .unwrap_or_else(|| panic!("{case}"))
            .parse()
            .expect("a count")
    };
    assert_eq!(rows("plain"), 4);
    assert_eq!(rows("min-gq"), 5);
    assert_eq!(rows("min-dp"), 5);
    // The state each floor writes is its own.
    let gq = table(&text, "detail", "min-gq");
    assert!(gq.iter().any(|row| row["CALL_STATE"] == "LOW_GQ"));
    let dp = table(&text, "detail", "min-dp");
    assert!(dp.iter().any(|row| row["CALL_STATE"] == "LOW_DP"));
    // Both are `TN,FN` against a het truth, which is what the filtered call is too.
    for state in [CallState::LowGq, CallState::LowDp, CallState::VcFiltered] {
        assert_eq!(
            contingency_values(contingency(&GA4GH, state, TruthState::HetRefVar1).expect("a cell")),
            Some("TN,FN".to_string()),
            "{state:?}"
        );
    }
}

/// `--IGNORE_FILTER_STATUS` reads the filtered call as it stands.
#[test]
fn ignoring_the_filter_reads_the_call() {
    let text = corpus();
    let plain = table(&text, "detail", "plain");
    let ignored = table(&text, "detail", "ignore-filters");
    assert!(plain.iter().any(|row| row["CALL_STATE"] == "VC_FILTERED"));
    assert!(!ignored.iter().any(|row| row["CALL_STATE"] == "VC_FILTERED"));
    // The site is not lost: it joins the agreeing het pair, which is why the row count falls.
    let agreeing = |rows: &[std::collections::HashMap<String, String>]| -> i64 {
        rows.iter()
            .find(|row| row["TRUTH_STATE"] == "HET_REF_VAR1" && row["CALL_STATE"] == "HET_REF_VAR1")
            .map(|row| row["COUNT"].parse().expect("a count"))
            .unwrap_or(0)
    };
    assert_eq!(agreeing(&plain), 3);
    assert_eq!(agreeing(&ignored), 4);
}

/// `--MISSING_SITES_HOM_REF` needs an interval list, and changes the scheme when it has one.
#[test]
fn the_hom_ref_flag_needs_an_interval_list() {
    let text = corpus();
    let refusal =
        field(&text, "refusal", "missing-as-hom-ref-without-intervals").expect("a refusal");
    assert!(
        refusal.starts_with(
            "You cannot use the MISSING_HOM option without also supplying an interval list"
        ),
        "{refusal}"
    );
    assert_eq!(
        field(&text, "error", "missing-as-hom-ref-without-intervals").as_deref(),
        Some("exit 1")
    );
    // With one, every position of the interval list that neither side called becomes a
    // MISSING/MISSING pair worth a true negative, which the plain scheme calls unreachable.
    let with = table(&text, "detail", "missing-as-hom-ref");
    let both_missing = with
        .iter()
        .find(|row| row["TRUTH_STATE"] == "MISSING" && row["CALL_STATE"] == "MISSING")
        .expect("a MISSING/MISSING row");
    assert_eq!(both_missing["CONTINGENCY_VALUES"], "TN");
    assert!(both_missing["COUNT"].parse::<i64>().expect("a count") > 6000);
    assert!(!table(&text, "detail", "a-missing-site")
        .iter()
        .any(|row| row["TRUTH_STATE"] == "MISSING"));
    assert_eq!(
        contingency(
            &GA4GH_MISSING_AS_HOM_REF,
            CallState::Missing,
            TruthState::Missing
        ),
        Some(Cell::Values(&[ContingencyState::Tn]))
    );
    // The two schemes differ exactly where a missing call meets a truth: MISSING is unreachable in
    // one and a true negative in the other.
    assert_eq!(
        contingency(&GA4GH, CallState::Missing, TruthState::Missing),
        Some(Cell::Unreachable)
    );
    assert_eq!(
        contingency(
            &GA4GH_MISSING_AS_HOM_REF,
            CallState::LowGq,
            TruthState::Missing
        ),
        Some(Cell::Values(&[ContingencyState::Tn]))
    );
    assert_eq!(
        contingency(&GA4GH, CallState::LowGq, TruthState::Missing),
        Some(Cell::Values(&[]))
    );
}

/// `--OUTPUT_ALL_ROWS` keeps the pairs nothing was seen for.
#[test]
fn every_row_may_be_kept() {
    let text = corpus();
    let all: usize = field(&text, "rows", "all-rows")
        .expect("a count")
        .parse()
        .expect("a number");
    assert_eq!(all, 374);
    // Which is more than the seventeen call states times the eleven truth ones, because the file
    // carries a row per variant type as well.
    assert_eq!(GA4GH.len(), 17);
    assert_eq!(TRUTH_ORDER.len(), 11);
    assert!(all > GA4GH.len() * TRUTH_ORDER.len());
    // The summary is per variant type either way.
    let summary = table(&text, "summary", "plain");
    assert_eq!(summary.len(), 2);
    assert_eq!(summary[0]["VARIANT_TYPE"], "SNP");
    assert_eq!(summary[1]["VARIANT_TYPE"], "INDEL");
    assert_eq!(summary[0]["HET_SENSITIVITY"], "0.75");
    assert_eq!(summary[0]["HET_SPECIFICITY"], "?");
}
