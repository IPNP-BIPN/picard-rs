//! Conformance for `FilterSamReads` (SAM I/O, read-list filters) against Picard 3.4.0.
//!
//! The corpus carries an input SAM, a `READ_LIST_FILE`, and the outputs of `FILTER=includeReadList`
//! and `FILTER=excludeReadList`. The port filters the same input with the same list and must reproduce
//! each SAM byte-for-byte. FilterSamReads keeps the input sort order and adds no @PG and no timestamp,
//! so the whole file is compared raw.

use std::io::Read;

use picard_analysis::filter_sam_reads::{filter_sam_reads, Filter};

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/filter_sam_reads.txt.gz");
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

fn assert_byte_identical(kind: &str, filter: Filter) {
    let ours = filter_sam_reads(&payload("input"), &payload("list"), filter).unwrap();
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
fn the_include_filter_is_byte_identical() {
    assert_byte_identical("include", Filter::IncludeReadList);
}

#[test]
fn the_exclude_filter_is_byte_identical() {
    assert_byte_identical("exclude", Filter::ExcludeReadList);
}
