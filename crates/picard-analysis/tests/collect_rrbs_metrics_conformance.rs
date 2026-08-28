//! Conformance for `CollectRrbsMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/rrbsmetrics-conformance/CollectRrbsMetricsDump.java`, which ran the tool
//! over a reference of `ACGTTCAACGTA` repeated: two CpG sites and one isolated cytosine per twelve
//! bases, so both branches of the collector run on the same read.
//!
//! # What this suite is for
//!
//!  * **a CpG being looked for in the reference and not in the read**;
//!  * **the last base of a block never being one**;
//!  * **the two quality thresholds being two different bars, on both branches**;
//!  * **the mismatch bound being rounded under a strictly greater test**;
//!  * **and the non-CpG counts reaching the totals even when the read carries no CpG, which is
//!    what the comment above that branch says they do not.**

use std::io::Read;

use picard_analysis::collect_rrbs_metrics::{
    block_counts, file_names, is_too_mismatched, is_too_short, mean_cpg_coverage, mismatch_bound,
    DEFAULT_C_QUALITY_THRESHOLD, DEFAULT_MAX_MISMATCH_RATE, DEFAULT_MINIMUM_READ_LENGTH,
    DEFAULT_NEXT_BASE_QUALITY_THRESHOLD, PLOT_BYTES_ARE_REPRODUCIBLE,
};

const MOTIF: &[u8] = b"ACGTTCAACGTA";

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/collect_rrbs_metrics.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

/// A table's one data row, by column name.
fn row(text: &str, kind: &str, case: &str) -> Vec<(String, String)> {
    let table = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = table.split('\n');
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let values: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    header
        .into_iter()
        .zip(values)
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

fn value(text: &str, case: &str, column: &str) -> String {
    row(text, "summary", case)
        .into_iter()
        .find(|(name, _)| name == column)
        .map(|(_, value)| value)
        .unwrap_or_else(|| panic!("{case}/{column}"))
}

fn count(text: &str, case: &str, column: &str) -> i64 {
    value(text, case, column).parse().expect("a number")
}

/// The detail table's positions, in the order it wrote them.
fn positions(text: &str, case: &str) -> Vec<usize> {
    let table = field(text, "detail", case).unwrap_or_else(|| panic!("detail/{case}"));
    table
        .split('\n')
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            line.split('\t')
                .nth(1)
                .expect("a position")
                .parse()
                .expect("a number")
        })
        .collect()
}

/// The reference over a one-based span, as the fixture built it.
fn reference(start: usize, length: usize) -> Vec<u8> {
    (start - 1..start - 1 + length)
        .map(|index| MOTIF[index % MOTIF.len()])
        .collect()
}

/// One read over the reference, with every quality at `I` unless another string is given.
fn counts_for(
    start: usize,
    read: &[u8],
    qualities: Option<&[u8]>,
) -> picard_analysis::collect_rrbs_metrics::Counts {
    let reference = reference(start, read.len());
    let default = vec![b'I'; read.len()];
    let qualities: Vec<u8> = qualities
        .map(|written| written.iter().map(|q| q - 33).collect())
        .unwrap_or_else(|| default.iter().map(|q| q - 33).collect());
    block_counts(
        &reference,
        read,
        &qualities,
        &qualities,
        start - 1,
        false,
        DEFAULT_C_QUALITY_THRESHOLD,
        DEFAULT_NEXT_BASE_QUALITY_THRESHOLD,
    )
}

/// A read that reads the reference back sees every CpG and converts none.
#[test]
fn a_cpg_is_looked_for_in_the_reference() {
    let text = corpus();
    let read = reference(1, 24);
    let counts = counts_for(1, &read, None);
    assert_eq!(
        counts.cpg_sites.len() as i64,
        count(&text, "unconverted", "CPG_BASES_SEEN")
    );
    assert_eq!(counts.converted_cpg_sites.len(), 0);
    assert_eq!(
        counts.non_cpg_total,
        count(&text, "unconverted", "NON_CPG_BASES")
    );
    assert_eq!(counts.cpg_sites, positions(&text, "unconverted"));
    // A read that reads TG over a reference CG is a CONVERTED site and not a mismatch, and the
    // isolated cytosine is a separate branch with a separate counter.
    let converted: Vec<u8> = String::from_utf8(read.clone())
        .expect("ascii")
        .replace("CG", "TG")
        .into_bytes();
    let counts = counts_for(1, &converted, None);
    assert_eq!(
        counts.converted_cpg_sites.len() as i64,
        count(&text, "cpg-converted", "CPG_BASES_CONVERTED")
    );
    assert_eq!(
        counts.non_cpg_converted,
        count(&text, "cpg-converted", "NON_CPG_CONVERTED_BASES")
    );
    let non_cpg: Vec<u8> = String::from_utf8(read.clone())
        .expect("ascii")
        .replace("CA", "TA")
        .into_bytes();
    let counts = counts_for(1, &non_cpg, None);
    assert_eq!(
        counts.non_cpg_converted,
        count(&text, "non-cpg-converted", "NON_CPG_CONVERTED_BASES")
    );
    assert_eq!(counts.converted_cpg_sites.len(), 0);
}

/// The last base of a block is never the C of a pair.
#[test]
fn the_last_base_is_never_a_cpg() {
    let text = corpus();
    let short = counts_for(1, &reference(1, 9), None);
    let long = counts_for(1, &reference(1, 10), None);
    assert_eq!(
        short.cpg_sites.len() as i64,
        count(&text, "cpg-at-the-last-base", "CPG_BASES_SEEN")
    );
    assert_eq!(
        long.cpg_sites.len() as i64,
        count(&text, "cpg-one-base-longer", "CPG_BASES_SEEN")
    );
    assert_eq!(long.cpg_sites.len(), short.cpg_sites.len() + 1);
}

/// The two thresholds are two different bars, and they gate both branches.
#[test]
fn the_two_thresholds_are_two_different_bars() {
    let text = corpus();
    // The neighbour's bar is the lower of the two, which is why a base that fails one may pass
    // the other: the golden shows a site kept at quality ten beside one dropped at nineteen.
    assert_eq!(
        DEFAULT_NEXT_BASE_QUALITY_THRESHOLD.min(DEFAULT_C_QUALITY_THRESHOLD),
        DEFAULT_NEXT_BASE_QUALITY_THRESHOLD
    );
    let read = reference(1, 12);
    // The cytosine of the first CpG, at and one under its own threshold.
    let at = counts_for(1, &read, Some(b"I5IIIIIIIIII"));
    let under = counts_for(1, &read, Some(b"I4IIIIIIIIII"));
    assert_eq!(
        at.cpg_sites.len() as i64,
        count(&text, "c-at-the-threshold", "CPG_BASES_SEEN")
    );
    assert_eq!(
        under.cpg_sites.len() as i64,
        count(&text, "c-under-the-threshold", "CPG_BASES_SEEN")
    );
    // Its neighbour, against the lower bar.
    let at = counts_for(1, &read, Some(b"II+IIIIIIIII"));
    let under = counts_for(1, &read, Some(b"II*IIIIIIIII"));
    assert_eq!(
        at.cpg_sites.len() as i64,
        count(&text, "neighbour-at-the-threshold", "CPG_BASES_SEEN")
    );
    assert_eq!(
        under.cpg_sites.len() as i64,
        count(&text, "neighbour-under-the-threshold", "CPG_BASES_SEEN")
    );
    // And the same two bars on the isolated cytosine, which is the other branch.
    let under = counts_for(1, &read, Some(b"IIIII4IIIIII"));
    assert_eq!(
        under.non_cpg_total,
        count(&text, "isolated-c-under-the-threshold", "NON_CPG_BASES")
    );
    let under = counts_for(1, &read, Some(b"IIIIII*IIIII"));
    assert_eq!(
        under.non_cpg_total,
        count(
            &text,
            "isolated-c-neighbour-under-the-threshold",
            "NON_CPG_BASES"
        )
    );
}

/// The mismatch bound is rounded and the test is strictly greater.
#[test]
fn the_mismatch_bound_is_rounded() {
    let text = corpus();
    assert_eq!(mismatch_bound(20, DEFAULT_MAX_MISMATCH_RATE), 2);
    assert!(!is_too_mismatched(2, 20, DEFAULT_MAX_MISMATCH_RATE));
    assert!(is_too_mismatched(3, 20, DEFAULT_MAX_MISMATCH_RATE));
    assert_eq!(
        count(&text, "two-mismatches", "READS_IGNORED_MISMATCHES"),
        0
    );
    assert_eq!(
        count(&text, "three-mismatches", "READS_IGNORED_MISMATCHES"),
        1
    );
    // Raising the rate keeps the same read, which is what says the bound and not the count moved.
    assert_eq!(
        count(
            &text,
            "three-mismatches-allowed",
            "READS_IGNORED_MISMATCHES"
        ),
        0
    );
    assert_eq!(mismatch_bound(20, 0.2), 4);
    // A conversion is not a mismatch: four of them on one read cost nothing.
    assert_eq!(
        count(
            &text,
            "four-conversions-are-not-mismatches",
            "READS_IGNORED_MISMATCHES"
        ),
        0
    );
    // And the length filter is its own, counted separately.
    assert!(is_too_short(4, DEFAULT_MINIMUM_READ_LENGTH));
    assert_eq!(count(&text, "short-read", "READS_IGNORED_SHORT"), 1);
    assert_eq!(count(&text, "short-read-allowed", "READS_IGNORED_SHORT"), 0);
}

/// The non-CpG counts reach the totals even when the read carries no CpG at all.
#[test]
fn a_read_with_no_cpg_still_counts_its_cytosines() {
    let text = corpus();
    let counts = counts_for(4, &reference(4, 6), None);
    assert!(counts.cpg_sites.is_empty());
    assert_eq!(
        counts.non_cpg_total,
        count(&text, "no-cpg", "NON_CPG_BASES")
    );
    assert_eq!(count(&text, "no-cpg", "READS_WITH_NO_CPG"), 1);
    assert_eq!(count(&text, "no-cpg", "CPG_BASES_SEEN"), 0);
    // With nothing in the histogram the mean is a NaN, which the file writes as a question mark.
    assert!(mean_cpg_coverage(&[]).is_nan());
    assert_eq!(value(&text, "no-cpg", "MEAN_CPG_COVERAGE"), "?");
    assert_eq!(value(&text, "no-cpg", "MEDIAN_CPG_COVERAGE"), "0");
}

/// The prefix gains a dot if it has none, and the plot is not bytes a golden can hold.
#[test]
fn the_prefix_gains_a_dot() {
    let text = corpus();
    let names = file_names("m");
    let mut written = field(&text, "files", "unconverted")
        .expect("the files")
        .split(' ')
        .map(str::to_string)
        .collect::<Vec<_>>();
    written.sort();
    let mut expected = names.to_vec();
    expected.sort();
    assert_eq!(written, expected);
    // A prefix that already ends in a dot does not gain a second one.
    assert_eq!(file_names("m."), file_names("m"));
    assert_eq!(
        field(&text, "files", "prefix-with-a-dot"),
        field(&text, "files", "unconverted")
    );
    // The plot is R's, and the golden holds its name and nothing about its bytes: no row of it
    // describes the file's contents, which is what the constant says.
    assert_eq!(
        PLOT_BYTES_ARE_REPRODUCIBLE,
        text.lines().any(|line| line.starts_with("plot\t"))
    );
    assert!(written.iter().any(|name| name.ends_with(".pdf")));
}
