//! Conformance for `ViewSam` (SAM I/O, default no-interval path) against Picard 3.4.0.
//!
//! The corpus carries one input SAM (mapped, vendor-fail, and unmapped reads) plus the stdout ViewSam
//! produced under each status filter: default (identity), Aligned, Unaligned, PF, and NonPF. The port
//! runs the same input under each filter and must reproduce the printed SAM byte-for-byte. ViewSam
//! prints the header verbatim and adds no @PG and no timestamp, so the whole output is compared raw.

use std::io::Read;

use picard_analysis::view_sam::{view_sam, AlignmentStatus, Options, PfStatus};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/view_sam.txt.gz");
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

fn assert_byte_identical(kind: &str, opts: &Options) {
    let input = payload("input");
    let ours = view_sam(&input, opts).unwrap();
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
fn the_default_view_is_byte_identical() {
    assert_byte_identical("default", &Options::default());
}

#[test]
fn the_aligned_view_is_byte_identical() {
    assert_byte_identical(
        "aligned",
        &Options {
            alignment_status: AlignmentStatus::Aligned,
            ..Default::default()
        },
    );
}

#[test]
fn the_unaligned_view_is_byte_identical() {
    assert_byte_identical(
        "unaligned",
        &Options {
            alignment_status: AlignmentStatus::Unaligned,
            ..Default::default()
        },
    );
}

#[test]
fn the_pf_view_is_byte_identical() {
    assert_byte_identical(
        "pf",
        &Options {
            pf_status: PfStatus::Pf,
            ..Default::default()
        },
    );
}

#[test]
fn the_nonpf_view_is_byte_identical() {
    assert_byte_identical(
        "nonpf",
        &Options {
            pf_status: PfStatus::NonPf,
            ..Default::default()
        },
    );
}
