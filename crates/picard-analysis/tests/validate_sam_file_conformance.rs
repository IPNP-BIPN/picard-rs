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

use picard_analysis::validate_sam_file::{
    validate_sam_file_verbose, validate_sam_file_verbose_with_reference,
};

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

/// Every case's input, keyed by case name, in file order: `(case, kind) -> payload`.
fn rows(name: &str) -> Vec<(String, String, String)> {
    corpus(name)
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

/// Runs every `(input, output)` pair in a corpus through the verbose validator and asserts
/// byte-identity, returning the number of cases checked.
fn check(name: &str) -> usize {
    let rows = rows(name);
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
    checked
}

/// Like [`check`], but each case is a `fasta`, `input`, `output` triple run through the
/// reference-taking validator.
fn check_with_ref(name: &str) -> usize {
    let rows = rows(name);
    let mut checked = 0;
    for window in rows.chunks(3) {
        let (case, kind0, fasta) = &window[0];
        let (_, kind1, input) = &window[1];
        let (_, kind2, expected) = &window[2];
        assert_eq!(kind0, "fasta", "case {case} out of shape");
        assert_eq!(kind1, "input", "case {case} out of shape");
        assert_eq!(kind2, "output", "case {case} out of shape");

        let ours = validate_sam_file_verbose_with_reference(input, fasta.as_bytes())
            .expect("input parses");
        assert_eq!(&ours, expected, "case {case} diverged");
        checked += 1;
    }
    checked
}

#[test]
fn every_header_and_record_case_is_byte_identical() {
    assert_eq!(
        check("validate_sam_file.txt.gz"),
        9,
        "expected 9 header/record cases"
    );
}

#[test]
fn every_nm_value_case_is_byte_identical() {
    assert_eq!(
        check_with_ref("validate_sam_file_nmvalue.txt.gz"),
        5,
        "expected 5 NM-value cases"
    );
}

#[test]
fn every_mate_case_is_byte_identical() {
    assert_eq!(
        check("validate_sam_file_mates.txt.gz"),
        6,
        "expected 6 mate cases"
    );
}

#[test]
fn every_isvalid_case_is_byte_identical() {
    assert_eq!(
        check("validate_sam_file_isvalid.txt.gz"),
        10,
        "expected 10 isValid() cases"
    );
}
