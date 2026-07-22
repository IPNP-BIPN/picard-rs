//! Conformance for `FixMateInformation` (single SAM input, default options) against Picard 3.4.0.
//!
//! The corpus carries a raw coordinate-sorted SAM whose pairs have absent mate info, and the cleaned
//! output. The port fixes the same input and must reproduce the SAM byte-for-byte. FixMateInformation
//! adds no @PG and no timestamp, so the whole file is compared raw. The cases include two FR pairs
//! and a pair split across contigs (which gets an explicit mate reference and a zero insert size).

use std::io::Read;

use picard_analysis::fix_mate_information::fix_mate_information;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fix_mate.txt.gz");
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
fn the_fixed_sam_is_byte_identical() {
    let ours = fix_mate_information(&payload("input"), None).unwrap();
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
