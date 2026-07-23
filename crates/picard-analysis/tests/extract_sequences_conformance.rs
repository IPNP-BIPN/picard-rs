//! Conformance for `ExtractSequences` against Picard 3.4.0.
//!
//! Each case carries the interval list, the reference FASTA, the `LINE_LENGTH`, and the FASTA Picard
//! extracted. The port runs `extract_sequences` on the same interval list, reference, and line length
//! and must reproduce the output byte-for-byte. (The interval list's `@SQ` header carries an absolute
//! `UR:` path from the oracle run; the port skips header lines, so that path never reaches the
//! output and does not affect the comparison.)

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::extract_sequences::{extract_sequences, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/extract_sequences.txt.gz");
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
    interval_list: String,
    reference: String,
    line_length: usize,
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
            "interval_list" => case.interval_list = payload,
            "reference" => case.reference = payload,
            "line_length" => case.line_length = payload.parse().unwrap(),
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
fn every_extracted_fasta_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 2, "case count");
    for (name, case) in &cases {
        let opts = Options {
            line_length: case.line_length,
        };
        let got = extract_sequences(&case.interval_list, &case.reference, &opts).expect("extract");
        assert_eq!(got, case.output, "{name}");
    }
}
