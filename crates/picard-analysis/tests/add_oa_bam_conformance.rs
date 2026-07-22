//! Conformance for `AddOATag` (SAM -> BAM) against Picard 3.4.0.
//!
//! The corpus carries a SAM and the BAM Picard produced by tagging it with `USE_JDK_DEFLATER=true`.
//! The port stamps the same input and writes a BAM through htsjdk-rs's `BamWriter`, which must
//! reproduce the BAM byte-for-byte. This joins the multicore per-record transform to a byte-identical
//! BAM write.

use std::io::Read;

use picard_analysis::add_oa_tag::add_oa_to_bam;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/add_oa_bam.txt.gz");
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
fn the_oa_tagged_bam_is_byte_identical() {
    let ours = add_oa_to_bam(&unescape(&payload("sam"))).unwrap();
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
