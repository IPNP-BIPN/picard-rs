//! Conformance for `ConvertSequencingArtifactToOxoG` against Picard 3.4.0.
//!
//! Each case carries the OxoG table the tool wrote, and the two input tables are the fixture the
//! dump built. There is no sequence data anywhere: the port is given the same two tables and must
//! produce the same rows.
//!
//! # What this suite is for
//!
//!  * **the output reporting the C contexts only**;
//!  * **the pre-adapter figures coming from the reverse-complement context**;
//!  * **the bait-bias figures coming from the context itself**;
//!  * **only C>A and G>T being read**;
//!  * **TOTAL_SITES always being nought**;
//!  * **the oxidation rate having a floor of one base**;
//!  * **the two bait-bias rates having a floor of 1e-10, which reads as a Q of a hundred**;
//!  * **the file names being derived from a basename**;
//!  * **and naming neither being refused.**

use std::io::Read;

use picard_analysis::convert_artifact_to_oxog::{
    convert, derived_names, is_oxo_g, oxog_contexts, oxog_libraries, reverse_complement,
    BaitBiasDetail, PreAdapterDetail, BAIT_BIAS_DETAILS_EXT, NO_BAIT_BIAS_MESSAGE,
    NO_OXOG_OUT_MESSAGE, NO_PRE_ADAPTER_MESSAGE, OXOG_METRICS_EXT, PRE_ADAPTER_DETAILS_EXT,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("convert_artifact_to_oxog.txt.gz");
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

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// One case's OxoG table, as its rows keyed by column name.
fn written(text: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let table = field(text, "metrics", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = table.lines().filter(|line| !line.is_empty());
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

/// The pre-adapter rows the dump wrote, read back from the golden's own copy of the file.
fn pre_adapter(text: &str) -> Vec<PreAdapterDetail> {
    rows(text, "pre-adapter")
        .into_iter()
        .map(|c| PreAdapterDetail {
            sample_alias: c[0].clone(),
            library: c[1].clone(),
            reference_base: c[2].chars().next().expect("a base"),
            alternate_base: c[3].chars().next().expect("a base"),
            context: c[4].clone(),
            pro_ref_bases: c[5].parse().expect("a count"),
            pro_alt_bases: c[6].parse().expect("a count"),
            con_ref_bases: c[7].parse().expect("a count"),
            con_alt_bases: c[8].parse().expect("a count"),
        })
        .collect()
}

fn bait_bias(text: &str) -> Vec<BaitBiasDetail> {
    rows(text, "bait-bias")
        .into_iter()
        .map(|c| BaitBiasDetail {
            library: c[1].clone(),
            reference_base: c[2].chars().next().expect("a base"),
            alternate_base: c[3].chars().next().expect("a base"),
            context: c[4].clone(),
            forward_ref_bases: c[5].parse().expect("a count"),
            forward_alt_bases: c[6].parse().expect("a count"),
            reverse_ref_bases: c[7].parse().expect("a count"),
            reverse_alt_bases: c[8].parse().expect("a count"),
        })
        .collect()
}

/// The data rows of one of the golden's two input files.
fn rows(text: &str, name: &str) -> Vec<Vec<String>> {
    field(text, "in", name)
        .unwrap_or_else(|| panic!("the golden carries in/{name}"))
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .skip(1)
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

fn round(value: f64) -> String {
    let rendered = format!("{value:.6}");
    let trimmed = rendered.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

/// The `one-context` case, which is the fixture the golden prints, converts to the row it wrote.
#[test]
fn the_printed_fixture_converts_to_the_written_row() {
    let text = corpus();
    let ours = convert(&pre_adapter(&text), &bait_bias(&text)).expect("converted");
    let theirs = written(&text, "one-context");
    assert_eq!(ours.len(), 1);
    assert_eq!(theirs.len(), 1);
    let (row, expected) = (&ours[0], &theirs[0]);
    assert_eq!(row.context, expected["CONTEXT"]);
    assert_eq!(row.library, expected["LIBRARY"]);
    assert_eq!(row.total_sites.to_string(), expected["TOTAL_SITES"]);
    assert_eq!(row.total_bases.to_string(), expected["TOTAL_BASES"]);
    assert_eq!(row.ref_oxo_bases.to_string(), expected["REF_OXO_BASES"]);
    assert_eq!(
        row.ref_nonoxo_bases.to_string(),
        expected["REF_NONOXO_BASES"]
    );
    assert_eq!(row.alt_oxo_bases.to_string(), expected["ALT_OXO_BASES"]);
    assert_eq!(
        row.alt_nonoxo_bases.to_string(),
        expected["ALT_NONOXO_BASES"]
    );
    assert_eq!(
        round(row.oxidation_error_rate),
        expected["OXIDATION_ERROR_RATE"]
    );
    assert_eq!(round(row.oxidation_q), expected["OXIDATION_Q"]);
    assert_eq!(round(row.c_ref_oxo_q), expected["C_REF_OXO_Q"]);
    assert_eq!(round(row.g_ref_oxo_q), expected["G_REF_OXO_Q"]);
}

/// The pre-adapter figures come from the reverse-complement context and the bait-bias ones from
/// the context itself, so one output row draws on two different input rows.
#[test]
fn the_two_inputs_are_read_at_different_contexts() {
    let text = corpus();
    assert_eq!(reverse_complement("ACA"), "TGT");
    assert_eq!(reverse_complement("TCT"), "AGA");
    let row = &written(&text, "one-context")[0];
    assert_eq!(row["CONTEXT"], "ACA");
    // The pre-adapter row for TGT carries 1000 and 40 and 2000 and 10, which is 3050 together.
    assert_eq!(row["TOTAL_BASES"], "3050");
    assert_eq!(row["REF_OXO_BASES"], "1000");
    assert_eq!(row["ALT_OXO_BASES"], "40");
    // Where ACA's own pre-adapter row carries 100 and 5 and 200 and 3, and is not read.
    let aca = pre_adapter(&text)
        .into_iter()
        .find(|r| r.context == "ACA")
        .expect("the ACA row");
    assert_eq!(aca.pro_ref_bases, 100);
    // The bait-bias figures come from ACA, not TGT.
    assert_eq!(row["C_REF_REF_BASES"], "900");
    assert_eq!(row["G_REF_REF_BASES"], "800");
}

/// Only the C>A and G>T transitions are read, and only the C contexts are reported.
#[test]
fn only_the_oxog_transitions_and_the_c_contexts() {
    let text = corpus();
    assert!(is_oxo_g('C', 'A'));
    assert!(is_oxo_g('G', 'T'));
    assert!(!is_oxo_g('A', 'G'));
    assert!(!is_oxo_g('T', 'C'));
    // The fixture carries an A>G row in both inputs, and no output row names AAA.
    let rows = pre_adapter(&text);
    assert!(rows.iter().any(|r| r.reference_base == 'A'));
    assert_eq!(oxog_contexts(&rows), vec!["ACA".to_string()]);
    assert_eq!(oxog_libraries(&rows), vec!["lib1".to_string()]);
    assert!(written(&text, "one-context")
        .iter()
        .all(|row| row["CONTEXT"] == "ACA"));
}

/// TOTAL_SITES is always nought, the input not carrying it.
#[test]
fn total_sites_is_always_nought() {
    let text = corpus();
    for case in [
        "one-context",
        "oxidation-floor",
        "bait-bias-floor",
        "bait-bias-reversed",
        "two-libraries-two-contexts",
    ] {
        for row in written(&text, case) {
            assert_eq!(row["TOTAL_SITES"], "0", "{case}");
        }
    }
}

/// The oxidation rate has a floor of ONE BASE, so a context with fewer oxidised alternates than
/// unoxidised ones reports `1 / TOTAL_BASES`.
#[test]
fn the_oxidation_floor_is_one_base() {
    let text = corpus();
    let row = &written(&text, "oxidation-floor")[0];
    assert_eq!(row["ALT_OXO_BASES"], "2");
    assert_eq!(row["ALT_NONOXO_BASES"], "50");
    assert_eq!(row["TOTAL_BASES"], "3052");
    // Not (2 - 50) / 3052, which is negative, but 1 / 3052.
    assert_eq!(row["OXIDATION_ERROR_RATE"], round(1.0 / 3052.0));
    let ours = convert(
        &[
            PreAdapterDetail {
                sample_alias: "sample1".to_string(),
                library: "lib1".to_string(),
                reference_base: 'G',
                alternate_base: 'T',
                context: "TGT".to_string(),
                pro_ref_bases: 1000,
                pro_alt_bases: 2,
                con_ref_bases: 2000,
                con_alt_bases: 50,
            },
            PreAdapterDetail {
                sample_alias: "sample1".to_string(),
                library: "lib1".to_string(),
                reference_base: 'C',
                alternate_base: 'A',
                context: "ACA".to_string(),
                pro_ref_bases: 0,
                pro_alt_bases: 0,
                con_ref_bases: 0,
                con_alt_bases: 0,
            },
        ],
        &[BaitBiasDetail {
            library: "lib1".to_string(),
            reference_base: 'C',
            alternate_base: 'A',
            context: "ACA".to_string(),
            forward_ref_bases: 900,
            forward_alt_bases: 30,
            reverse_ref_bases: 800,
            reverse_alt_bases: 10,
        }],
    )
    .expect("converted");
    assert_eq!(
        round(ours[0].oxidation_error_rate),
        row["OXIDATION_ERROR_RATE"]
    );
}

/// The two bait-bias rates are floored at 1e-10 and are opposite differences of the same two
/// numbers, so at most one is above the floor and the other is a Q of exactly a hundred.
#[test]
fn the_bait_bias_floor_is_a_hundred_on_the_q() {
    let text = corpus();
    // Equal rates: both differences are nought, so both are floored and both Qs are a hundred.
    let both = &written(&text, "bait-bias-floor")[0];
    assert_eq!(both["C_REF_OXO_Q"], "100");
    assert_eq!(both["G_REF_OXO_Q"], "100");
    // The forward rate the larger.
    let forward = &written(&text, "one-context")[0];
    assert_eq!(forward["G_REF_OXO_Q"], "100");
    assert_ne!(forward["C_REF_OXO_Q"], "100");
    // And the reverse rate the larger, which swaps which of the two is floored.
    let reverse = &written(&text, "bait-bias-reversed")[0];
    assert_eq!(reverse["C_REF_OXO_Q"], "100");
    assert_ne!(reverse["G_REF_OXO_Q"], "100");
    // -10 * log10(1e-10) is a hundred exactly.
    assert_eq!(round(-10.0 * (1e-10f64).log10()), "100");
}

/// Two libraries by two contexts is four rows, one per pair.
#[test]
fn every_library_meets_every_context() {
    let text = corpus();
    let theirs = written(&text, "two-libraries-two-contexts");
    assert_eq!(theirs.len(), 4);
    let mut pairs: Vec<(String, String)> = theirs
        .iter()
        .map(|row| (row["LIBRARY"].clone(), row["CONTEXT"].clone()))
        .collect();
    pairs.sort();
    assert_eq!(
        pairs,
        vec![
            ("lib1".to_string(), "ACA".to_string()),
            ("lib1".to_string(), "TCT".to_string()),
            ("lib2".to_string(), "ACA".to_string()),
            ("lib2".to_string(), "TCT".to_string()),
        ]
    );
}

/// The three file names are derived from a basename, and `--OUTPUT_BASE` defaults to the input's.
#[test]
fn the_names_are_derived_from_a_basename() {
    let text = corpus();
    let (pre, bait, out) = derived_names(Some("in"), None).expect("derived");
    assert_eq!(pre, format!("in{PRE_ADAPTER_DETAILS_EXT}"));
    assert_eq!(bait, format!("in{BAIT_BIAS_DETAILS_EXT}"));
    assert_eq!(out, format!("in{OXOG_METRICS_EXT}"));
    assert_eq!(
        field(&text, "name", "one-context").as_deref(),
        Some("in.oxog_metrics")
    );
    // A separate output basename.
    let (_, _, other) = derived_names(Some("in"), Some("other")).expect("derived");
    assert_eq!(other, "other.oxog_metrics");
    assert_eq!(
        field(&text, "name", "output-base").as_deref(),
        Some("other.oxog_metrics")
    );
    // And the file named outright, which bypasses the derivation altogether.
    assert_eq!(
        field(&text, "name", "oxog-out").as_deref(),
        Some("named.txt")
    );
}

/// Naming neither a basename nor the files it would derive is refused, by both messages at once.
#[test]
fn naming_neither_is_refused() {
    let text = corpus();
    let errors = derived_names(None, None).expect_err("refused");
    assert_eq!(
        errors,
        vec![
            NO_PRE_ADAPTER_MESSAGE.to_string(),
            NO_BAIT_BIAS_MESSAGE.to_string(),
            NO_OXOG_OUT_MESSAGE.to_string(),
        ]
    );
    // The tool answers a failing exit code rather than an exception.
    assert_eq!(
        field(&text, "error", "no-basename").as_deref(),
        Some("exit 1")
    );
}
