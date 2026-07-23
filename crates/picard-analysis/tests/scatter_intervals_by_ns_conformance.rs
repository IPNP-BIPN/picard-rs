//! Conformance for `ScatterIntervalsByNs` against Picard 3.4.0.
//!
//! Each case carries the reference FASTA, its sequence dictionary, the `OUTPUT_TYPE`/`MAX_TO_MERGE`
//! options, and the interval list Picard produced. The port runs `scatter_intervals_by_ns` on the same
//! reference, dictionary, and options and must reproduce the output byte-for-byte, `@SQ` header
//! included: the header's `@SQ` lines are the dictionary's verbatim, and the committed dictionary is
//! the one the oracle ran on, so the absolute `UR:` path matches.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::scatter_intervals_by_ns::{scatter_intervals_by_ns, Options, OutputType};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/scatter_intervals_by_ns.txt.gz");
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
    reference: String,
    dict: String,
    output_type: String,
    max_to_merge: i32,
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
            "reference" => case.reference = payload,
            "dict" => case.dict = payload,
            "output_type" => case.output_type = payload,
            "max_to_merge" => case.max_to_merge = payload.parse().unwrap(),
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

fn output_type(name: &str) -> OutputType {
    match name {
        "N" => OutputType::N,
        "ACGT" => OutputType::Acgt,
        "BOTH" => OutputType::Both,
        other => panic!("unexpected OUTPUT_TYPE {other}"),
    }
}

#[test]
fn every_interval_list_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 3, "case count");
    for (name, case) in &cases {
        let opts = Options {
            output_type: output_type(&case.output_type),
            max_to_merge: case.max_to_merge,
        };
        let got = scatter_intervals_by_ns(&case.reference, &case.dict, &opts).expect("scatter");
        assert_eq!(got, case.output, "{name}");
    }
}
