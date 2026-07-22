//! Conformance for `SamFormatConverter` (BAM -> SAM) against Picard 3.4.0.
//!
//! The corpus carries a small BAM (as hex) and the SAM Picard converted it to. The port decompresses
//! the BAM and runs the same conversion, which must reproduce the SAM byte-for-byte. SamFormatConverter
//! adds no @PG and no timestamp, so the whole file is compared raw. This is the `BamReader` ->
//! `write_sam` pipeline the benchmark exercises at scale, checked here at the tool level.

use std::io::Read;

use picard_analysis::sam_format_converter::bam_to_sam;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/sam_format_converter.txt.gz");
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
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
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
fn the_converted_sam_is_byte_identical() {
    let bam = from_hex(&payload("bam"));
    let plain = htsjdk_bgzf::decompress_all(&bam).expect("bam decompresses");
    let ours = bam_to_sam(&plain).unwrap();
    let theirs = unescape(&payload("sam"));
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
