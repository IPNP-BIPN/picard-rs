//! Conformance for `CreateSequenceDictionary` against Picard 3.4.0.
//!
//! The corpus carries a reference FASTA and the `.dict` Picard wrote for it. The port builds the
//! dictionary from the same FASTA and must reproduce it, apart from the `UR` field (the reference's
//! `file:` URI, which is path-dependent). Both sides strip `UR` and compare `@HD` and each
//! `SN`/`LN`/`M5` raw.

use std::io::Read;

use picard_analysis::create_sequence_dictionary::create_sequence_dictionary;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/create_sequence_dictionary.txt.gz");
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

/// Drop the `UR:` field (the path-dependent reference URI) from each line.
fn strip_ur(dict: &str) -> String {
    dict.lines()
        .map(|l| {
            let kept: Vec<&str> = l.split('\t').filter(|f| !f.starts_with("UR:")).collect();
            format!("{}\n", kept.join("\t"))
        })
        .collect()
}

#[test]
fn the_dictionary_is_byte_identical_apart_from_the_ur_field() {
    let ours =
        create_sequence_dictionary(payload("fasta").as_bytes(), "file:///placeholder").unwrap();
    let ours = strip_ur(&ours);
    let theirs = strip_ur(&payload("dict"));
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
