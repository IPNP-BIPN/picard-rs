//! Conformance for `SamToFastq` (unpaired reads, default options) against Picard 3.4.0.
//!
//! The corpus carries the input BAM (hex) and the FASTQ file Picard produced. The port decodes the
//! BAM, runs the unpaired path, and must reproduce the FASTQ bytes. The reads exercise a forward
//! read, a negative-strand read (reverse-complemented with reversed qualities), an N base, and
//! secondary / supplementary / vendor-fail reads dropped by default.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use picard_analysis::sam_to_fastq::{sam_to_fastq_unpaired, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sam_to_fastq.txt.gz");
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

fn from_hex(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
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
            (k == kind).then(|| p.to_string())
        })
        .unwrap_or_else(|| panic!("no {kind} row"))
}

#[test]
fn the_fastq_output_is_byte_identical() {
    let bam = htsjdk_bgzf::decompress_all(&from_hex(&payload("bam"))).unwrap();
    let reader = BamReader::new(&bam).expect("bam header");
    let records: Vec<_> = reader.map(|r| r.expect("record")).collect();

    let ours = sam_to_fastq_unpaired(&records, &Options::default());
    let theirs = unescape(&payload("fastq"));

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
