//! Conformance for `NormalizeFasta` against Picard 3.4.0.
//!
//! Each case carries the input FASTA, the `LINE_LENGTH` / `TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE`
//! options, and the normalized output Picard produced. The port runs `normalize_fasta` on the same
//! input and options and must reproduce the output byte-for-byte (plain text, compared raw).

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::normalize_fasta::{normalize_fasta, Options};

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/normalize_fasta.txt.gz");
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

#[derive(Default)]
struct Case {
    input: String,
    line_length: usize,
    truncate: bool,
    output: String,
}

fn cases() -> Vec<(String, Case)> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, Case> = HashMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let name = it.next().unwrap().to_string();
        let payload = unescape(it.next().unwrap_or(""));
        let case = map.entry(name.clone()).or_insert_with(|| {
            order.push(name.clone());
            Case::default()
        });
        match kind {
            "input" => case.input = payload,
            "line_length" => case.line_length = payload.parse().unwrap(),
            "truncate" => case.truncate = payload == "true",
            "output" => case.output = payload,
            "rc" => {}
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|n| (n.clone(), map.remove(&n).unwrap()))
        .collect()
}

#[test]
fn every_normalized_fasta_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 3, "case count");
    for (name, case) in &cases {
        let opts = Options {
            line_length: case.line_length,
            truncate_names_at_whitespace: case.truncate,
        };
        let got = normalize_fasta(&case.input, &opts).expect("normalize");
        assert_eq!(got, case.output, "{name}");
    }
}
