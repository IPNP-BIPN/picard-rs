//! Conformance for `CollectOxoGMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/oxog-conformance`. Each case carries the input as SAM, the context of every
//! row in the file's own order, and the rows that had a site in them.
//!
//! # What this suite is for
//!
//!  * **the row order being a HashMap's over the contexts and not a sorted one**;
//!  * **a G site folded into the reverse complement of its context**;
//!  * **which end of the pair, on which strand, reaching which counter**;
//!  * **the oxidation rate floored at one base rather than at nought**;
//!  * **the reference-bias rates floored at 1e-10, capping their Q at a hundred**;
//!  * **a context no read covered writing `?`**;
//!  * **the base-quality floor reading the OQ tag unless it is turned off**;
//!  * **each whole-read filter dropping its read and leaving the clean pair**;
//!  * **a site inside the contig's own context never being assayed**;
//!  * **and the two contexts the validation refuses, by name.**

use std::io::Read;

use picard_analysis::collect_oxo_g_metrics::{
    accept, context_at, contexts, finish, reverse_complement, row_order, validate_context, Counts,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_oxo_g_metrics.txt.gz");
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

fn rows(text: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let body = field(text, "rows", case).unwrap_or_else(|| panic!("rows/{case}"));
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

/// The one context every case's reads cover, which is where the counters land.
fn assayed(text: &str, case: &str) -> std::collections::HashMap<String, String> {
    rows(text, case)
        .into_iter()
        .find(|row| row["CONTEXT"] == "ACG")
        .unwrap_or_else(|| panic!("an ACG row in {case}"))
}

fn number(row: &std::collections::HashMap<String, String>, name: &str) -> f64 {
    let value = row.get(name).unwrap_or_else(|| panic!("{name}"));
    if value == "?" {
        f64::NAN
    } else {
        value.parse().unwrap_or_else(|_| panic!("{name}={value}"))
    }
}

/// The rows come out in a HashMap's order over the contexts, which is neither sorted nor the
/// order they were built in.
#[test]
fn the_row_order_is_the_tables_and_not_a_sorted_one() {
    let text = corpus();
    let theirs: Vec<String> = field(&text, "order", "reference-only")
        .expect("an order")
        .split(',')
        .map(str::to_string)
        .collect();
    let all = contexts(1);
    assert_eq!(all.len(), 16);
    let ours: Vec<String> = row_order(&all, &["lib1".to_string()])
        .into_iter()
        .map(|(context, _)| context)
        .collect();
    assert_eq!(ours, theirs);
    // Which is not the sorted order the contexts were generated in.
    assert_ne!(theirs, all);
    assert_eq!(theirs[0], "CCA");
    // Two libraries put the libraries of one context together, not the contexts of one library.
    let two: Vec<String> = field(&text, "order", "two-libraries")
        .expect("an order")
        .split(',')
        .map(str::to_string)
        .collect();
    assert_eq!(two.len(), 32);
    assert_eq!(two[0], two[1]);
    let ours_two: Vec<String> = row_order(&all, &["lib1".to_string(), "lib2".to_string()])
        .into_iter()
        .map(|(context, _)| context)
        .collect();
    assert_eq!(ours_two, two);
}

/// A G site is folded into the reverse complement of its own context.
#[test]
fn a_g_site_is_folded_into_the_reverse_complement() {
    let text = corpus();
    let reference = field(&text, "ref", "chr1").expect("the reference");
    let bases = reference.as_bytes();
    // Position 22 is a G whose window is CGT, filed under ACG; position 21 is the C of that same
    // window, filed under its own.
    assert_eq!(bases[21] as char, 'G');
    assert_eq!(context_at(bases, 22, 1).as_deref(), Some("ACG"));
    assert_eq!(context_at(bases, 21, 1).as_deref(), Some("ACG"));
    assert_eq!(reverse_complement("CGT"), "ACG");
    // So one row carries both strands, told apart only by the C_REF and G_REF columns.
    let row = assayed(&text, "alt-read-one-forward");
    assert_ne!(row["C_REF_REF_BASES"], "0");
    assert_ne!(row["G_REF_REF_BASES"], "0");
}

/// Which end of the pair, on which strand, decides the counter.
#[test]
fn the_end_and_the_strand_decide_the_counter() {
    let text = corpus();
    // Read one forward and read two reverse are the oxidised state; the other two are the control.
    for case in ["alt-read-one-forward", "alt-read-two-reverse"] {
        let row = assayed(&text, case);
        assert_eq!(row["ALT_OXO_BASES"], "1", "{case}");
        assert_eq!(row["ALT_NONOXO_BASES"], "0", "{case}");
    }
    for case in ["alt-read-one-reverse", "alt-read-two-forward"] {
        let row = assayed(&text, case);
        assert_eq!(row["ALT_OXO_BASES"], "0", "{case}");
        assert_eq!(row["ALT_NONOXO_BASES"], "1", "{case}");
    }
    // At a C site the alternate is an A and the two ends swap roles.
    assert_eq!(
        assayed(&text, "c-site-read-two-forward")["ALT_OXO_BASES"],
        "1"
    );
    assert_eq!(
        assayed(&text, "c-site-read-one-forward")["ALT_NONOXO_BASES"],
        "1"
    );
    // Which is what the port's own counters do with the same four reads.
    let mut oxidised = Counts::default();
    accept(&mut oxidised, b'G', b'T', 1, false);
    assert_eq!(oxidised.ref_g_oxidated_a, 1);
    let mut also = Counts::default();
    accept(&mut also, b'G', b'T', 2, true);
    assert_eq!(also.ref_g_oxidated_a, 1);
    let mut control = Counts::default();
    accept(&mut control, b'G', b'T', 1, true);
    assert_eq!(control.ref_g_control_a, 1);
    let mut c_site = Counts::default();
    accept(&mut c_site, b'C', b'A', 2, false);
    assert_eq!(c_site.ref_c_oxidated_a, 1);
}

/// The oxidation rate is floored at one base, so a row with no alternate still has a finite Q.
#[test]
fn the_oxidation_rate_is_floored_at_one_base() {
    let text = corpus();
    let clean = assayed(&text, "reference-only");
    assert_eq!(clean["ALT_OXO_BASES"], "0");
    assert_eq!(clean["ALT_NONOXO_BASES"], "0");
    let total = number(&clean, "TOTAL_BASES");
    // The rate is written with six fraction digits, so 1/12 comes back as 0.083333.
    assert!((number(&clean, "OXIDATION_ERROR_RATE") - 1.0 / total).abs() < 1e-6);
    assert!(number(&clean, "OXIDATION_Q").is_finite());
    // Which is the port's own arithmetic: nought minus nought floored at one.
    let counts = Counts {
        ref_c_control_c: 6,
        ref_g_control_c: 6,
        ..Counts::default()
    };
    let metrics = finish(&counts);
    assert_eq!(metrics.total_bases, 12);
    assert_eq!(metrics.oxidation_error_rate, 1.0 / 12.0);
}

/// The reference-bias rates are floored at 1e-10, which caps their Q at a hundred.
#[test]
fn the_reference_bias_rates_are_floored() {
    let text = corpus();
    let clean = assayed(&text, "reference-only");
    assert_eq!(clean["C_REF_OXO_Q"], "100");
    assert_eq!(clean["G_REF_OXO_Q"], "100");
    // The side that carries the alternate drops off the cap; the other stays on it.
    let one = assayed(&text, "alt-read-one-forward");
    assert_eq!(one["C_REF_OXO_Q"], "100");
    assert!(number(&one, "G_REF_OXO_Q") < 100.0);
    let counts = Counts {
        ref_c_control_c: 6,
        ref_g_control_c: 5,
        ref_g_oxidated_a: 1,
        ..Counts::default()
    };
    let metrics = finish(&counts);
    assert_eq!(metrics.c_ref_oxo_error_rate, 1e-10);
    assert!((metrics.c_ref_oxo_q - 100.0).abs() < 1e-9);
    assert!(metrics.g_ref_oxo_q < 100.0);
}

/// A context no read covered divides nought by nought, which the writer renders as `?`.
#[test]
fn a_context_no_read_covered_writes_a_question_mark() {
    let text = corpus();
    let row = rows(&text, "contig-start")
        .into_iter()
        .find(|row| row["CONTEXT"] == "TCA")
        .expect("a TCA row");
    assert_eq!(row["C_REF_OXO_Q"], "?");
    assert_eq!(row["G_REF_OXO_Q"], "?");
    // The site at position one is never assayed, whatever covers it: it lies inside the contig's
    // own context.
    let reference = field(&text, "ref", "chr1").expect("the reference");
    assert_eq!(reference.as_bytes()[0] as char, 'C');
    assert_eq!(context_at(reference.as_bytes(), 1, 1), None);
    assert_eq!(context_at(reference.as_bytes(), reference.len(), 1), None);
    // And a base that is neither a C nor a G is not a site at all.
    assert_eq!(context_at(b"AAAAA", 3, 1), None);
}

/// The base-quality floor reads the OQ tag unless it is turned off.
#[test]
fn the_quality_floor_reads_the_oq_tag() {
    let text = corpus();
    let honoured = assayed(&text, "original-qualities");
    let ignored = assayed(&text, "original-qualities-ignored");
    // The same read, counted or not according to a tag the alignment did not have to keep.
    assert_eq!(honoured["ALT_OXO_BASES"], "0");
    assert_eq!(ignored["ALT_OXO_BASES"], "1");
    assert_eq!(
        number(&honoured, "TOTAL_BASES") + 1.0,
        number(&ignored, "TOTAL_BASES")
    );
    // The plain quality floor drops the same base without any tag at all.
    assert_eq!(assayed(&text, "low-base-quality")["ALT_OXO_BASES"], "0");
}

/// Each whole-read filter drops its read and leaves the clean pair standing.
#[test]
fn each_filter_drops_its_read_and_leaves_the_pair() {
    let text = corpus();
    let counted = number(&assayed(&text, "original-qualities-ignored"), "TOTAL_BASES");
    for case in [
        "low-mapping-quality",
        "insert-too-small",
        "insert-too-large",
        "duplicate",
        "secondary",
    ] {
        let row = assayed(&text, case);
        assert_eq!(row["ALT_OXO_BASES"], "0", "{case}");
        // Six fewer bases than the run that counted the same read: the clean pair alone.
        assert_eq!(number(&row, "TOTAL_BASES"), counted - 6.0, "{case}");
    }
}

/// The two contexts the validation refuses, by name.
#[test]
fn the_validation_refuses_two_kinds_of_context() {
    let text = corpus();
    assert_eq!(
        field(&text, "refusal", "context-without-a-c").as_deref(),
        Some("Middle base of context sequence TTT must be C")
    );
    assert_eq!(
        field(&text, "refusal", "context-of-the-wrong-length").as_deref(),
        Some("Context ACGTA is not 3 long as implied by CONTEXT_SIZE=1")
    );
    assert_eq!(
        field(&text, "error", "context-without-a-c").as_deref(),
        Some("exit 1")
    );
    assert_eq!(
        validate_context("TTT", 1),
        Err("Middle base of context sequence TTT must be C".to_string())
    );
    assert_eq!(
        validate_context("ACGTA", 1),
        Err("Context ACGTA is not 3 long as implied by CONTEXT_SIZE=1".to_string())
    );
    assert_eq!(validate_context("ACG", 1), Ok(()));
    // A context size of nought leaves one context, `C`.
    assert_eq!(contexts(0), vec!["C".to_string()]);
    assert_eq!(
        field(&text, "order", "context-size-zero").as_deref(),
        Some("C")
    );
}
