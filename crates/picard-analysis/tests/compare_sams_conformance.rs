//! Conformance for `CompareSAMs` (queryname-sorted, strict) against Picard 3.4.0.
//!
//! Each case carries two SAM inputs, the stdout verdict Picard printed (`SAM files match.` /
//! `SAM files differ.` plus its return code), and the metrics report with the command-line and
//! timestamp banner stripped and the two file-path columns canonicalized to `LEFT` / `RIGHT`. The
//! port runs `compare_sams` on the same inputs and must reproduce the verdict and, after the same
//! banner stripping, the report byte-for-byte.

use std::io::Read;

use picard_analysis::compare_sams::{compare_sams, verdict, write_report};

fn corpus(name: &str) -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// The same banner stripping `CmpCorpus.java` applies: drop the `## htsjdk...` and `# ...` lines.
fn strip_banner(report: &str) -> String {
    let mut out = String::new();
    for line in report.split('\n') {
        if line.starts_with("## htsjdk") || line.starts_with("# ") {
            continue;
        }
        out.push_str(line);
        out.push('\n');
    }
    out
}

/// Group the flat corpus rows into `(input1, input2, verdict, report)` per case, in file order.
fn cases(name: &str) -> Vec<(String, String, String, String)> {
    let text = corpus(name);
    let mut it = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let mut p = l.splitn(3, '\t');
            let kind = p.next().unwrap().to_string();
            let _case = p.next().unwrap();
            (kind, unescape(p.next().unwrap_or("")))
        });
    let mut out = Vec::new();
    while let Some((k1, input1)) = it.next() {
        let (k2, input2) = it.next().expect("input2 row");
        let (k3, verdict) = it.next().expect("verdict row");
        let (k4, report) = it.next().expect("report row");
        assert_eq!(
            (k1.as_str(), k2.as_str(), k3.as_str(), k4.as_str()),
            ("input1", "input2", "verdict", "report")
        );
        out.push((input1, input2, verdict, report));
    }
    out
}

fn check(name: &str, expected: usize) {
    let cases = cases(name);
    assert_eq!(cases.len(), expected, "case count for {name}");
    for (input1, input2, want_verdict, want_report) in &cases {
        let metric = compare_sams(input1, input2, "LEFT", "RIGHT").expect("inputs parse");
        // The oracle records "SAM files match. rc=0" (or "differ. rc=1"): rc is 0 when equal, else 1.
        let rc = if metric.are_equal { 0 } else { 1 };
        assert_eq!(&format!("{} rc={rc}", verdict(&metric)), want_verdict);
        assert_eq!(&strip_banner(&write_report(&metric)), want_report);
    }
}

#[test]
fn every_queryname_case_is_byte_identical() {
    check("compare_sams.txt.gz", 7);
}

#[test]
fn every_coordinate_case_is_byte_identical() {
    check("compare_sams_coord.txt.gz", 6);
}
