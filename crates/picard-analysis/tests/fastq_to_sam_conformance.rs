//! Conformance for `FastqToSam` (unpaired, default options, SAM output) against Picard 3.4.0.
//!
//! The corpus carries the input FASTQ and the SAM file Picard produced. The port converts the same
//! FASTQ and must reproduce the SAM byte-for-byte. FastqToSam writes no `@PG` and no timestamp, so
//! the whole file is compared raw. The input reads are in non-queryname order (r2, r10, r1, r3) so
//! the output must be re-sorted, and one carries a `/1` suffix to exercise the read-name cleanup.

use std::io::Read;

use picard_analysis::fastq_to_sam::{fastq_to_sam_unpaired, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fastq_to_sam.txt.gz");
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
fn the_sam_output_is_byte_identical() {
    let ours = fastq_to_sam_unpaired(&payload("fastq"), &Options::new("s1")).unwrap();
    let theirs = payload("sam");
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
