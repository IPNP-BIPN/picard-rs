//! Conformance for `EstimateLibraryComplexity` against Picard 3.4.0.
//!
//! Golden from `tools/libcomplexity-conformance`. Each case carries the input as SAM, the metrics
//! table and the histogram.
//!
//! # What this suite is for
//!
//!  * **the grouping being on the first `MIN_IDENTICAL_BASES` of both ends**;
//!  * **the comparison skipping that prefix**;
//!  * **the diff rate being over both ends' compared length, floored**;
//!  * **the quality floor being an integer division**;
//!  * **an `N` in the seed and a short read dropping the pair**;
//!  * **a bin under `MIN_GROUP_COUNT` leaving the metrics and staying in the histogram**;
//!  * **an optical duplicate being subtracted before the estimate**;
//!  * **and the libraries being counted apart, when the read name allows it.**

use std::io::Read;

use picard_analysis::estimate_library_complexity::{
    matches, metrics, passes_quality_check, same_group, Metrics, DEFAULT_MAX_DIFF_RATE,
    DEFAULT_MIN_GROUP_COUNT, DEFAULT_MIN_IDENTICAL_BASES, DEFAULT_MIN_MEAN_QUALITY,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("estimate_library_complexity.txt.gz");
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
    let body = field(text, "metrics", case).unwrap_or_else(|| panic!("metrics/{case}"));
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

fn row(text: &str, case: &str) -> std::collections::HashMap<String, String> {
    rows(text, case).remove(0)
}

fn number(row: &std::collections::HashMap<String, String>, name: &str) -> f64 {
    let value = row.get(name).unwrap_or_else(|| panic!("{name}"));
    if value.is_empty() {
        f64::NAN
    } else {
        value.parse().unwrap_or_else(|_| panic!("{name}={value}"))
    }
}

/// The histogram of a case, as (set size, group count) pairs.
fn histogram(text: &str, case: &str) -> Vec<(i64, i64)> {
    let body = field(text, "histogram", case).unwrap_or_else(|| panic!("histogram/{case}"));
    body.lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut columns = line.split('\t');
            (
                columns.next().expect("a size").parse().expect("a number"),
                columns.next().expect("a count").parse().expect("a number"),
            )
        })
        .collect()
}

/// The grouping is on the first `MIN_IDENTICAL_BASES` of both ends.
#[test]
fn the_grouping_is_on_the_seed_of_both_ends() {
    let text = corpus();
    // Two pairs differing at the sixth base are one group at five and two at six.
    assert_eq!(histogram(&text, "differ-at-the-sixth-base"), vec![(2, 1)]);
    assert_eq!(histogram(&text, "seed-of-six"), vec![(1, 2)]);
    // A difference INSIDE the seed separates them whatever the rate allows.
    assert_eq!(histogram(&text, "differ-inside-the-seed"), vec![(1, 2)]);
    let one = b"ACGTACGTACGTACGTACGTACGTACGTAC";
    let two = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAA";
    let sixth = b"ACGTATGTACGTACGTACGTACGTACGTAC";
    let third = b"ACTTACGTACGTACGTACGTACGTACGTAC";
    assert!(same_group(
        (one, two),
        (sixth, two),
        DEFAULT_MIN_IDENTICAL_BASES
    ));
    assert!(!same_group((one, two), (sixth, two), 6));
    assert!(!same_group(
        (one, two),
        (third, two),
        DEFAULT_MIN_IDENTICAL_BASES
    ));
}

/// The comparison skips the seed, and the rate is over both ends together.
#[test]
fn the_rate_is_over_both_ends_and_floored() {
    let text = corpus();
    // One difference over sixty compared bases: inside 0.03 and outside 0.01.
    assert_eq!(histogram(&text, "one-difference"), vec![(2, 1)]);
    assert_eq!(histogram(&text, "one-difference-strict"), vec![(1, 2)]);
    let one = b"ACGTACGTACGTACGTACGTACGTACGTAC";
    let two = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAA";
    let sixth = b"ACGTATGTACGTACGTACGTACGTACGTAC";
    assert!(matches(
        (one, two),
        (sixth, two),
        DEFAULT_MIN_IDENTICAL_BASES,
        DEFAULT_MAX_DIFF_RATE,
        0
    ));
    assert!(!matches(
        (one, two),
        (sixth, two),
        DEFAULT_MIN_IDENTICAL_BASES,
        0.01,
        0
    ));
    // Sixty compared bases at 0.03 allow one error: the floor of 1.8.
    assert_eq!(((30 + 30) as f64 * DEFAULT_MAX_DIFF_RATE).floor() as i64, 1);
    // A difference past --MAX_READ_LENGTH is not compared at all.
    assert_eq!(
        histogram(&text, "difference-past-the-truncation"),
        vec![(2, 1)]
    );
    assert_eq!(
        histogram(&text, "difference-inside-the-window"),
        vec![(1, 2)]
    );
    let late = b"ACGTACGTACGTACGTACGTACGTAGGTAC";
    assert!(matches(
        (one, two),
        (late, two),
        DEFAULT_MIN_IDENTICAL_BASES,
        0.0,
        20
    ));
    assert!(!matches(
        (one, two),
        (late, two),
        DEFAULT_MIN_IDENTICAL_BASES,
        0.0,
        0
    ));
}

/// The quality floor is an integer division, and the seed must be callable.
#[test]
fn the_quality_floor_and_the_seed_admit_the_pair() {
    let text = corpus();
    // Qualities of nineteen are dropped at twenty and kept at one.
    assert_eq!(histogram(&text, "low-mean-quality"), vec![(1, 1)]);
    assert_eq!(histogram(&text, "low-mean-quality-allowed"), vec![(2, 1)]);
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTAC";
    let nineteen = [19u8; 30];
    let twenty = [20u8; 30];
    assert!(!passes_quality_check(
        bases,
        &nineteen,
        DEFAULT_MIN_IDENTICAL_BASES,
        DEFAULT_MIN_MEAN_QUALITY,
        0
    ));
    assert!(passes_quality_check(
        bases,
        &twenty,
        DEFAULT_MIN_IDENTICAL_BASES,
        DEFAULT_MIN_MEAN_QUALITY,
        0
    ));
    // The division is integer: a mean of 19.9 is nineteen.
    let mut almost = [20u8; 30];
    almost[0] = 17;
    assert!(!passes_quality_check(
        bases,
        &almost,
        DEFAULT_MIN_IDENTICAL_BASES,
        DEFAULT_MIN_MEAN_QUALITY,
        0
    ));
    // An N in the seed, and a read shorter than it, each drop the pair.
    assert_eq!(histogram(&text, "an-n-in-the-seed"), vec![(1, 1)]);
    assert_eq!(histogram(&text, "shorter-than-the-seed"), vec![(1, 1)]);
    assert!(!passes_quality_check(
        b"ACNTACGTAC",
        &[40u8; 10],
        DEFAULT_MIN_IDENTICAL_BASES,
        DEFAULT_MIN_MEAN_QUALITY,
        0
    ));
    assert!(!passes_quality_check(
        b"ACG",
        &[40u8; 3],
        DEFAULT_MIN_IDENTICAL_BASES,
        DEFAULT_MIN_MEAN_QUALITY,
        0
    ));
}

/// A bin under `MIN_GROUP_COUNT` leaves the metrics and stays in the histogram.
#[test]
fn a_lonely_bin_leaves_the_metrics_alone() {
    let text = corpus();
    let lonely = row(&text, "two-identical");
    assert_eq!(number(&lonely, "READ_PAIRS_EXAMINED"), 0.0);
    assert_eq!(histogram(&text, "two-identical"), vec![(2, 1)]);
    // The same file with the floor at one reports the pair, and an estimate with it.
    let counted = row(&text, "one-group-counted");
    assert_eq!(number(&counted, "READ_PAIRS_EXAMINED"), 2.0);
    assert_eq!(number(&counted, "READ_PAIR_DUPLICATES"), 1.0);
    assert_eq!(number(&counted, "PERCENT_DUPLICATION"), 0.5);
    assert_eq!(counted["ESTIMATED_LIBRARY_SIZE"], "1");
    // Two groups of two clear the default floor on their own.
    let two = row(&text, "two-duplicate-groups");
    assert_eq!(number(&two, "READ_PAIRS_EXAMINED"), 4.0);
    assert_eq!(number(&two, "READ_PAIR_DUPLICATES"), 2.0);
    // Which is what the port makes of the same bins.
    assert_eq!(
        metrics(&[(2, 1, 0)], DEFAULT_MIN_GROUP_COUNT),
        Metrics::default()
    );
    let ours = metrics(&[(2, 1, 0)], 1);
    assert_eq!(ours.read_pairs_examined, 2);
    assert_eq!(ours.read_pair_duplicates, 1);
    assert_eq!(ours.percent_duplication, 0.5);
    assert_eq!(ours.estimated_library_size, Some(1));
    let both = metrics(&[(2, 2, 0)], DEFAULT_MIN_GROUP_COUNT);
    assert_eq!(
        (both.read_pairs_examined, both.read_pair_duplicates),
        (4, 2)
    );
}

/// An optical duplicate is subtracted before the estimate runs.
#[test]
fn an_optical_duplicate_leaves_nothing_to_estimate_from() {
    let text = corpus();
    let optical = row(&text, "optical-duplicates");
    assert_eq!(number(&optical, "READ_PAIR_DUPLICATES"), 1.0);
    assert_eq!(number(&optical, "READ_PAIR_OPTICAL_DUPLICATES"), 1.0);
    assert_eq!(optical["ESTIMATED_LIBRARY_SIZE"], "");
    // The deeper fixture, whose reads sit far apart on the flowcell, does answer.
    let deep = row(&text, "an-estimate");
    assert_eq!(number(&deep, "READ_PAIRS_EXAMINED"), 10.0);
    assert_eq!(number(&deep, "READ_PAIR_DUPLICATES"), 3.0);
    assert_eq!(number(&deep, "READ_PAIR_OPTICAL_DUPLICATES"), 0.0);
    assert_eq!(deep["ESTIMATED_LIBRARY_SIZE"], "13");
    // Both are the port's own arithmetic over the same bins.
    assert_eq!(metrics(&[(2, 1, 1)], 1).estimated_library_size, None);
    let ours = metrics(&[(1, 7, 0), (2, 2, 0), (3, 0, 0)], 1);
    assert_eq!(ours.read_pairs_examined, 11);
}

/// The libraries are counted apart, and the read name is what makes that possible.
#[test]
fn the_libraries_are_counted_apart() {
    let text = corpus();
    let two = rows(&text, "two-libraries");
    assert_eq!(two.len(), 2);
    assert_eq!(two[0]["LIBRARY"], "lib1");
    assert_eq!(two[1]["LIBRARY"], "lib2");
    // The histogram carries a column per library.
    let body = field(&text, "histogram", "two-libraries").expect("a histogram");
    assert!(body
        .lines()
        .next()
        .expect("a header")
        .contains("lib1\tlib2"));
    // Every case here names its reads by flowcell position, which is what records the read group.
    let sam = field(&text, "sam", "two-libraries").expect("the input");
    assert!(sam.contains("H0164ALXX140820:2:1101:"));
    assert!(sam.contains("RG:Z:rg-lib2"));
}
