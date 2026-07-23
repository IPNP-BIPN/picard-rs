//! Conformance for `NonNFastaSize` (whole-genome, no intervals) against Picard 3.4.0.
//!
//! Each case carries the input FASTA and the count Picard wrote. The port runs `non_n_fasta_size` on
//! the same input and must reproduce the output (the count plus a newline) byte-for-byte.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::non_n_fasta_size::non_n_fasta_size;

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/non_n_fasta_size.txt.gz");
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

fn cases() -> Vec<(String, String, String)> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, (String, String)> = HashMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let name = it.next().unwrap().to_string();
        let payload = unescape(it.next().unwrap_or(""));
        let entry = map.entry(name.clone()).or_insert_with(|| {
            order.push(name.clone());
            (String::new(), String::new())
        });
        match kind {
            "input" => entry.0 = payload,
            "output" => entry.1 = payload,
            "rc" => {}
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|n| {
            let (i, o) = map.remove(&n).unwrap();
            (n, i, o)
        })
        .collect()
}

#[test]
fn every_count_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 2, "case count");
    for (name, input, output) in &cases {
        assert_eq!(&non_n_fasta_size(input).unwrap(), output, "{name}");
    }
}
