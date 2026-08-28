//! Conformance for `AccumulateVariantCallingMetrics` against Picard 3.4.0.
//!
//! Each case carries the merged detail and summary tables. The golden prints one input pair, so
//! the port is driven by the reference's own files for that case and by the fixture's own numbers
//! for the rest.
//!
//! # What this suite is for
//!
//!  * **the arguments being prefixes with two fixed extensions**;
//!  * **the merge being lossy, visible on a single input**;
//!  * **`invertFromRatio` rounding, and NaN reconstructing as nought**;
//!  * **the counts adding while the ratios are recomputed**;
//!  * **the merge being per SAMPLE_ALIAS**;
//!  * **the summary's bias resting on the detail file beside it**;
//!  * **and a summary of more than one row being refused.**

use std::io::Read;

use picard_analysis::accumulate_variant_calling_metrics::{
    accumulate, file_names, invert_from_ratio, wrong_summary_row_count_message, DetailMetrics,
    Input, SummaryMetrics, DETAIL_EXTENSION, SUMMARY_EXTENSION,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("accumulate_variant_calling_metrics.txt.gz");
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

/// One table's rows, keyed by column name.
fn rows(text: &str, kind: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let table = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
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

/// A detail row of the fixture, before any reconstruction.
fn detail(
    sample: &str,
    total_snps: i64,
    in_db_snp: i64,
    dbsnp_titv: f64,
    bias: f64,
    het_homvar: f64,
    het_depth: i64,
) -> DetailMetrics {
    DetailMetrics {
        sample_alias: sample.to_string(),
        summary: SummaryMetrics {
            total_snps,
            num_in_db_snp: in_db_snp,
            novel_snps: 0,
            dbsnp_titv,
            novel_titv: f64::NAN,
            snp_reference_bias: bias,
            ..Default::default()
        },
        total_het_depth: het_depth,
        het_homvar_ratio: het_homvar,
        ..Default::default()
    }
}

fn summary(total_snps: i64, in_db_snp: i64, dbsnp_titv: f64, bias: f64) -> SummaryMetrics {
    SummaryMetrics {
        total_snps,
        num_in_db_snp: in_db_snp,
        novel_snps: 0,
        dbsnp_titv,
        novel_titv: f64::NAN,
        snp_reference_bias: bias,
        ..Default::default()
    }
}

fn even() -> Input {
    Input {
        detail: vec![detail("s1", 300, 300, 2.0, 0.5, 1.0, 100)],
        summary: vec![summary(300, 300, 2.0, 0.5)],
    }
}

fn number(value: f64) -> String {
    if value.is_nan() {
        return "?".to_string();
    }
    let rendered = format!("{value:.6}");
    let trimmed = rendered.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

/// The arguments are prefixes, not files.
#[test]
fn the_arguments_are_prefixes() {
    let (d, s) = file_names("out");
    assert_eq!(d, "out.variant_calling_detail_metrics");
    assert_eq!(s, "out.variant_calling_summary_metrics");
    assert_eq!(DETAIL_EXTENSION, "variant_calling_detail_metrics");
    assert_eq!(SUMMARY_EXTENSION, "variant_calling_summary_metrics");
}

/// The merge is lossy, and one input shows it: 301 SNPs at a TI/TV of 2.0 come back at 2.01.
#[test]
fn a_single_input_is_already_lossy() {
    let text = corpus();
    let odd = rows(&text, "detail", "one-input-odd");
    assert_eq!(odd[0]["TOTAL_SNPS"], "301");
    assert_eq!(odd[0]["DBSNP_TITV"], "2.01");
    // Where 300 divides evenly and comes back unchanged.
    let even_rows = rows(&text, "detail", "one-input-even");
    assert_eq!(even_rows[0]["DBSNP_TITV"], "2");
    // Which the port reaches by the same route.
    let (merged, _) = accumulate(&[even()]).expect("merged");
    assert_eq!(number(merged[0].summary.dbsnp_titv), "2");
    let odd_input = Input {
        detail: vec![detail("s1", 301, 301, 2.0, 0.5, 1.0, 100)],
        summary: vec![summary(301, 301, 2.0, 0.5)],
    };
    let (merged, _) = accumulate(&[odd_input]).expect("merged");
    assert_eq!(number(merged[0].summary.dbsnp_titv), "2.01");
}

/// `invertFromRatio` rounds, and NaN reconstructs as nought.
#[test]
fn the_reconstruction_rounds_and_nan_is_nought() {
    // 301 / (2 + 1) is 100.33, which rounds to 100, leaving 201 transitions: 2.01.
    assert_eq!(invert_from_ratio(301, 2.0), 100);
    assert_eq!(invert_from_ratio(300, 2.0), 100);
    assert_eq!(invert_from_ratio(0, 2.0), 0);
    assert_eq!(invert_from_ratio(300, f64::NAN), 0);
    let text = corpus();
    // And the NaN comes back out as nought, not as NaN.
    let nan = rows(&text, "detail", "nan-ratio");
    assert_eq!(nan[0]["DBSNP_TITV"], "0");
    assert_eq!(nan[0]["SNP_REFERENCE_BIAS"], "0");
}

/// The counts add while the ratios are recomputed.
#[test]
fn the_counts_add() {
    let text = corpus();
    let one = rows(&text, "detail", "one-input-even");
    let two = rows(&text, "detail", "two-inputs-one-sample");
    assert_eq!(one[0]["TOTAL_SNPS"], "300");
    assert_eq!(two[0]["TOTAL_SNPS"], "600");
    assert_eq!(two[0]["TOTAL_HET_DEPTH"], "200");
    // The ratio is unchanged, both halves dividing evenly.
    assert_eq!(two[0]["DBSNP_TITV"], "2");
    let (merged, _) = accumulate(&[even(), even()]).expect("merged");
    assert_eq!(merged[0].summary.total_snps, 600);
    assert_eq!(merged[0].total_het_depth, 200);
    assert_eq!(number(merged[0].summary.dbsnp_titv), "2");
}

/// The merge is per SAMPLE_ALIAS, whether the samples arrive in one file or two.
#[test]
fn the_merge_is_per_sample() {
    let text = corpus();
    for case in ["two-samples", "two-samples-one-file"] {
        let detail_rows = rows(&text, "detail", case);
        assert_eq!(detail_rows.len(), 2, "{case}");
        let mut names: Vec<&String> = detail_rows.iter().map(|r| &r["SAMPLE_ALIAS"]).collect();
        names.sort();
        assert_eq!(names, vec!["s1", "s2"], "{case}");
        // s1's counts are its own, not the pair's.
        let s1 = detail_rows
            .iter()
            .find(|r| r["SAMPLE_ALIAS"] == "s1")
            .expect("s1");
        assert_eq!(s1["TOTAL_SNPS"], "300", "{case}");
    }
    // The port keeps them apart the same way.
    let other = Input {
        detail: vec![detail("s2", 100, 100, 1.0, 0.25, 2.0, 40)],
        summary: vec![summary(100, 100, 1.0, 0.25)],
    };
    let (merged, _) = accumulate(&[even(), other]).expect("merged");
    assert_eq!(merged.len(), 2);
    assert_eq!(merged[0].sample_alias, "s1");
    assert_eq!(merged[1].sample_alias, "s2");
}

/// The summary is one row however many inputs there were, and its bias rests on the detail file
/// beside it.
#[test]
fn the_summary_is_one_row() {
    let text = corpus();
    for case in ["one-input-even", "two-inputs-one-sample", "two-samples"] {
        assert_eq!(rows(&text, "summary", case).len(), 1, "{case}");
    }
    let two = &rows(&text, "summary", "two-inputs-one-sample")[0];
    assert_eq!(two["TOTAL_SNPS"], "600");
    let (_, summary_row) = accumulate(&[even(), even()]).expect("merged");
    assert_eq!(summary_row.total_snps, 600);
    // The bias survives, both halves carrying the same one.
    assert_eq!(number(summary_row.snp_reference_bias), "0.5");
    assert_eq!(two["SNP_REFERENCE_BIAS"], "0.5");
}

/// A summary file of more than one row is refused by a message counting them.
#[test]
fn a_summary_of_two_rows_is_refused() {
    let text = corpus();
    let error = field(&text, "error", "two-summary-rows").expect("its refusal");
    assert_eq!(
        error,
        format!(
            "picard.PicardException:{}",
            wrong_summary_row_count_message(2)
        )
    );
    let two_rows = Input {
        detail: even().detail,
        summary: vec![summary(300, 300, 2.0, 0.5), summary(1, 1, 2.0, 0.5)],
    };
    assert_eq!(
        accumulate(&[two_rows]),
        Err(wrong_summary_row_count_message(2))
    );
}

/// A missing input is refused by htsjdk rather than by the tool, so the message names the file.
#[test]
fn a_missing_input_names_the_file() {
    let text = corpus();
    let error = field(&text, "error", "missing-input").expect("its refusal");
    assert!(
        error.starts_with("htsjdk.samtools.SAMException:"),
        "{error}"
    );
    assert!(error.contains(DETAIL_EXTENSION), "{error}");
    // The tool's own catch would have named the prefix, and it is not what answered.
    assert!(
        !error.contains("Cannot read from metrics files with prefix"),
        "{error}"
    );
}
