//! Conformance for `AccumulateQualityYieldMetrics` against Picard 3.4.0.
//!
//! The case carries the two input `QualityYieldMetrics` files and the combined file Picard wrote. The
//! port runs `accumulate_quality_yield_metrics` on the same inputs and must reproduce the output
//! byte-for-byte. The tool writes a bare `MetricsFile` with no command-line or start-time header
//! comments, so the whole file is compared raw with no canonicalization.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::accumulate_quality_yield_metrics::accumulate_quality_yield_metrics;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/accumulate_quality_yield.txt.gz");
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

fn rows() -> HashMap<String, String> {
    let mut map = HashMap::new();
    for line in corpus().lines() {
        if line.is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap().to_string();
        let _case = it.next().unwrap();
        let payload = unescape(it.next().unwrap_or(""));
        map.insert(kind, payload);
    }
    map
}

#[test]
fn the_combined_metrics_file_is_byte_identical() {
    let map = rows();
    let input1 = map.get("input1").expect("input1");
    let input2 = map.get("input2").expect("input2");
    let expected = map.get("output").expect("output");
    let got = accumulate_quality_yield_metrics(&[input1, input2]).expect("accumulate");
    assert_eq!(&got, expected);
}
