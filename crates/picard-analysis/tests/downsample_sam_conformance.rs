//! Conformance for `DownsampleSam` (SAM I/O, ConstantMemory strategy) against Picard 3.4.0.
//!
//! The corpus carries a coordinate-sorted SAM of 12 distinctly-named reads and the output after
//! DownsampleSam ran with PROBABILITY=0.5 and RANDOM_SEED=1. The port downsamples the same input and
//! must reproduce the survivors byte-for-byte. DownsampleSam adds a `@PG` provenance record whose
//! `CL:` is the command line; that record is canonicalized away, so both sides drop `@PG` lines and
//! the surviving records plus the rest of the header are compared raw.

use std::io::Read;

use picard_analysis::downsample_sam::{downsample_sam, DEFAULT_SEED};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/downsample_sam.txt.gz");
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

/// Drop `@PG` lines, the canonicalized command-line provenance record.
fn strip_pg(sam: &str) -> String {
    sam.lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(|l| format!("{l}\n"))
        .collect()
}

#[test]
fn the_downsampled_sam_is_byte_identical_apart_from_the_pg_record() {
    let ours = downsample_sam(&payload("input"), 0.5, DEFAULT_SEED).unwrap();
    let theirs = strip_pg(&payload("output"));
    let ours = strip_pg(&ours);
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
