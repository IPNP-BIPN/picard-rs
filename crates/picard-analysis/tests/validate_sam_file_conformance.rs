//! Conformance for `ValidateSamFile` (VERBOSE mode, SAM input, no reference) against Picard 3.4.0.
//!
//! The corpus carries, per case, a hand-built SAM input and the exact stdout of `ValidateSamFile
//! MODE=VERBOSE`. The cases exercise the header and per-record checks that need neither the
//! reference nor cross-record (mate/sort) state: a clean file, a mapped read missing `NM`, an empty
//! read-group header with a record missing its read group, a `QUAL == *` read, missing / invalid
//! read-group platform, an unacceptable header version, an empty dictionary with only an unmapped
//! read, and a mixed file that pins the interleaving order. The verbose output has no timestamp and
//! no banner, so each case is compared byte-for-byte.

use std::io::Read;

use picard_analysis::validate_sam_file::validate_sam_file_verbose;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/validate_sam_file.txt.gz");
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

/// Every case's input, keyed by case name, in file order: `(case, kind) -> payload`.
fn rows() -> Vec<(String, String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let mut it = l.splitn(3, '\t');
            let kind = it.next().unwrap().to_string();
            let case = it.next().unwrap().to_string();
            let payload = unescape(it.next().unwrap_or(""));
            (case, kind, payload)
        })
        .collect()
}

#[test]
fn every_case_is_byte_identical() {
    let rows = rows();
    let mut checked = 0;
    // Each case has an `input` row followed by an `output` row.
    for window in rows.chunks(2) {
        let (case, kind0, input) = &window[0];
        let (_, kind1, expected) = &window[1];
        assert_eq!(kind0, "input", "case {case} out of shape");
        assert_eq!(kind1, "output", "case {case} out of shape");

        let ours = validate_sam_file_verbose(input).expect("input parses");
        assert_eq!(&ours, expected, "case {case} diverged");
        checked += 1;
    }
    assert_eq!(checked, 9, "expected 9 conformance cases");
}
