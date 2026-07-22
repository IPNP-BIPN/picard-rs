//! Conformance for `CalculateReadGroupChecksum` (SAM I/O) against Picard 3.4.0.
//!
//! The corpus carries an input SAM with two read groups and the MD5 digest Picard wrote to the
//! `.read_group_md5` file. The port computes the digest over the same header and must reproduce the
//! 32-character hex string exactly (no trailing newline), so it is compared raw.

use std::io::Read;

use picard_analysis::calculate_read_group_checksum::calculate_read_group_checksum;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/calculate_read_group_checksum.txt.gz");
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
fn the_read_group_checksum_is_byte_identical() {
    let ours = calculate_read_group_checksum(&payload("input")).unwrap();
    assert_eq!(ours, payload("output"));
}
