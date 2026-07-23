//! Conformance for `BedToIntervalList` against Picard 3.4.0.
//!
//! Each case carries the sequence dictionary, the BED, the `SORT`/`UNIQUE`/`KEEP_LENGTH_ZERO_INTERVALS`
//! flags, and the interval list Picard produced. The port runs `bed_to_interval_list` on the same
//! dictionary, BED, and flags and must reproduce the output byte-for-byte, `@SQ` header included: the
//! output's `@SQ` lines are the dictionary's `@SQ` lines verbatim, and the committed dictionary is the
//! one the oracle ran on, so the absolute `UR:` path matches.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::bed_to_interval_list::{bed_to_interval_list, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/bed_to_interval_list.txt.gz");
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
    dict: String,
    bed: String,
    sort: bool,
    unique: bool,
    keep_zero: bool,
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
            "dict" => case.dict = payload,
            "bed" => case.bed = payload,
            "sort" => case.sort = payload == "true",
            "unique" => case.unique = payload == "true",
            "keep_zero" => case.keep_zero = payload == "true",
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
fn every_interval_list_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 3, "case count");
    for (name, case) in &cases {
        let opts = Options {
            sort: case.sort,
            unique: case.unique,
            keep_length_zero_intervals: case.keep_zero,
            drop_missing_contigs: false,
        };
        let got = bed_to_interval_list(&case.dict, &case.bed, &opts).expect("convert");
        assert_eq!(got, case.output, "{name}");
    }
}
