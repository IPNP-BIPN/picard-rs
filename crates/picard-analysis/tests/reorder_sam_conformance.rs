//! Conformance for `ReorderSam` (SAM I/O, unindexed path) against Picard 3.4.0.
//!
//! The corpus carries a coordinate-sorted SAM over a read dictionary chr1, chr2, chr3, a
//! `SEQUENCE_DICTIONARY` (.dict) that swaps chr1/chr2 and drops chr3, and the output after ReorderSam
//! ran with `ALLOW_INCOMPLETE_DICT_CONCORDANCE`. The port reorders the same input against the same
//! dictionary and must reproduce the SAM byte-for-byte: the header is cloned with the new @SQ block,
//! the chr3 read and the chr3-mate become unmapped (MC removed), and the writer re-sorts into the new
//! coordinate order. ReorderSam adds no @PG and no timestamp, so the whole file is compared raw.

use std::io::Read;

use picard_analysis::reorder_sam::{reorder_sam, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reorder_sam.txt.gz");
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
fn the_reordered_sam_is_byte_identical() {
    // The golden was generated with ALLOW_INCOMPLETE_DICT_CONCORDANCE=true.
    let opts = Options {
        allow_incomplete_dict_concordance: true,
        allow_contig_length_discordance: false,
    };
    let ours = reorder_sam(&payload("input"), &payload("dict"), &opts).unwrap();
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
