//! Conformance for `AddOrReplaceReadGroups` (SAM I/O, default sort order) against Picard 3.4.0.
//!
//! The corpus carries an htsjdk-generated coordinate-sorted SAM (already carrying an @RG ID:1 and
//! RG:Z:1 tags) and the output after replacing the read group with ID:2. The port replaces the read
//! group and the RG tags on the same input and must reproduce the SAM byte-for-byte. The tool adds
//! no @PG and no timestamp and does not re-sort, so the whole file is compared raw.

use std::io::Read;

use picard_analysis::add_or_replace_read_groups::{add_or_replace_read_groups, Options};

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/add_replace_rg.txt.gz");
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
fn the_reheadered_sam_is_byte_identical() {
    let opts = Options {
        rgid: "2".to_string(),
        ..Options::new("lib1", "ILLUMINA", "unit1", "sample1")
    };
    let ours = add_or_replace_read_groups(&payload("input"), &opts).unwrap();
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
