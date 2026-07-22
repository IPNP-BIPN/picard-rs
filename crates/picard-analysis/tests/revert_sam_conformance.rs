//! Conformance for `RevertSam` (SAM I/O, default option path) against Picard 3.4.0.
//!
//! The corpus carries a coordinate-sorted SAM whose queryname order differs from its coordinate order
//! (a duplicate read with an OQ to restore and NM/MD/AS to clear, a negative-strand read, and a
//! proper pair with an MC tag), and the output after RevertSam ran with default options. The port
//! reverts the same input and must reproduce the SAM byte-for-byte: the bare @HD/@RG header, the
//! queryname re-order, and every record stripped to an unmapped read. RevertSam adds no @PG and no
//! timestamp, so the whole file is compared raw.

use std::io::Read;

use picard_analysis::revert_sam::revert_sam;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/revert_sam.txt.gz");
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

fn payload(kind: &str) -> String {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .find_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let _case = it.next()?;
            let p = it.next().unwrap_or("");
            (k == kind).then(|| unescape(p))
        })
        .unwrap_or_else(|| panic!("no {kind} row"))
}

#[test]
fn the_reverted_sam_is_byte_identical() {
    let ours = revert_sam(&payload("input")).unwrap();
    let theirs = payload("output");
    if ours != theirs {
        let at = ours
            .lines()
            .zip(theirs.lines())
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        panic!(
            "first difference at line {at}\n  picard: {:?}\n  ours  : {:?}",
            theirs.lines().nth(at),
            ours.lines().nth(at)
        );
    }
}
