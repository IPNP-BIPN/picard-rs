//! Conformance for `RevertOriginalBaseQualitiesAndAddMateCigar` (SAM I/O, default path) against
//! Picard 3.4.0.
//!
//! The corpus carries a coordinate-sorted proper pair, both mapped on-reference, each with an OQ to
//! restore and no MC, and the output after the tool ran. The port runs the same input and must
//! reproduce the SAM byte-for-byte: QUAL restored from OQ (OQ dropped), the mate cigar `MC` and mate
//! mapping quality `MQ` stamped on each end, and the coordinate order kept. The tool adds no @PG and
//! no timestamp, so the whole file is compared raw.

use std::io::Read;

use picard_analysis::revert_original_quals_add_mate_cigar::revert_original_base_qualities_and_add_mate_cigar;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/revert_original_quals_add_mate_cigar.txt.gz");
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
fn the_reverted_and_mate_cigared_sam_is_byte_identical() {
    let ours = revert_original_base_qualities_and_add_mate_cigar(&payload("input")).unwrap();
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
