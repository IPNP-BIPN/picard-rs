//! Conformance for `SortSam` (SAM I/O, coordinate and queryname orders) against Picard 3.4.0.
//!
//! The corpus carries an htsjdk-generated unsorted SAM and the two SortSam outputs. The port sorts
//! the same input each way and must reproduce the output byte-for-byte. SortSam changes only the
//! header's `SO` field and adds no `@PG` and no timestamp, so the whole SAM is compared raw. The
//! input has reads out of both orders across two contigs, with a reverse-strand read sharing a
//! position, which pins the coordinate comparator's forward-before-reverse tie-break in a real file.

use std::io::Read;

use picard_analysis::sort_sam::{sort_sam, SortOrder};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sort_sam.txt.gz");
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

fn assert_sorted(order: SortOrder, kind: &str) {
    let ours = sort_sam(&payload("input"), order).unwrap();
    let theirs = payload(kind);
    if ours != theirs {
        let at = ours
            .lines()
            .zip(theirs.lines())
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        panic!(
            "{kind}: first difference at line {at}\n  picard: {:?}\n  ours  : {:?}",
            theirs.lines().nth(at),
            ours.lines().nth(at)
        );
    }
}

#[test]
fn coordinate_sorted_output_is_byte_identical() {
    assert_sorted(SortOrder::Coordinate, "coordinate");
}

#[test]
fn queryname_sorted_output_is_byte_identical() {
    assert_sorted(SortOrder::Queryname, "queryname");
}
