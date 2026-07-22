//! Conformance for `SamFormatConverter` (SAM -> BAM) against Picard 3.4.0.
//!
//! The corpus carries a small SAM and the BAM Picard produced from it with `USE_JDK_DEFLATER=true`.
//! The port converts the same SAM with htsjdk-rs's `BamWriter` and must reproduce the BAM byte-for-byte
//! (BGZF framing, header, records, terminator block). The JDK deflater is what makes the BGZF blocks
//! reproducible; Picard's default GKL deflater is a separate surface.

use std::io::Read;

use picard_analysis::sam_format_converter::{bam_to_sam, sam_to_bam};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sam_to_bam.txt.gz");
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
fn the_converted_bam_is_byte_identical() {
    let ours = sam_to_bam(&unescape(&payload("sam"))).unwrap();
    let theirs = from_hex(&payload("bam"));
    assert_eq!(
        ours.len(),
        theirs.len(),
        "byte length differs: ours={} picard={}",
        ours.len(),
        theirs.len()
    );
    if let Some(at) = ours.iter().zip(&theirs).position(|(a, b)| a != b) {
        panic!(
            "first byte differs at offset {at}: ours={:#04x} picard={:#04x}",
            ours[at], theirs[at]
        );
    }
}

/// The BAM the port writes reads back to the same SAM it came from.
#[test]
fn the_bam_round_trips_back_to_the_input_sam() {
    let sam = unescape(&payload("sam"));
    let bam = sam_to_bam(&sam).unwrap();
    let plain = htsjdk_bgzf::decompress_all(&bam).expect("bam decompresses");
    assert_eq!(bam_to_sam(&plain).unwrap(), sam);
}
