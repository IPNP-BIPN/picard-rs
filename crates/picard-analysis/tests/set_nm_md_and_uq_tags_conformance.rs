//! Conformance for `SetNmMdAndUqTags` (SAM I/O, default path) against Picard 3.4.0.
//!
//! The corpus carries a reference FASTA, a coordinate-sorted SAM (a perfect read, a low-quality
//! single mismatch, and a 2bp deletion), and the SAM after SetNmMdAndUqTags recomputed MD/NM/UQ. The
//! port runs the same input against the same reference and must reproduce the SAM byte-for-byte. The
//! tool adds no @PG and no timestamp, so the whole file is compared raw.

use std::io::Read;

use picard_analysis::set_nm_md_and_uq_tags::set_nm_md_and_uq_tags;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/set_nm_md_and_uq_tags.txt.gz");
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
fn the_tagged_sam_is_byte_identical() {
    let ours = set_nm_md_and_uq_tags(&payload("input"), payload("fasta").as_bytes()).unwrap();
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
