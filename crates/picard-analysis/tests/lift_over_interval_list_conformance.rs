//! Conformance for `LiftOverIntervalList` against Picard 3.4.0.
//!
//! Each case carries the input interval list, the target `SEQUENCE_DICTIONARY` `.dict`, the UCSC
//! chain file, `MIN_LIFTOVER_PCT`, the process return code, and the output interval list Picard
//! wrote. The port runs `lift_over_interval_list` on the same inputs and must reproduce the output
//! byte-for-byte and the return code. (The `@SQ` line of the `.dict` carries an absolute `UR:` path
//! from the oracle run; it appears verbatim in both the committed dictionary and Picard's output, so
//! the port's verbatim passthrough reproduces it regardless of the path.)

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::lift_over_interval_list::lift_over_interval_list;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/liftover_interval_list.txt.gz");
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
    sd: String,
    chain: String,
    min_pct: f64,
    rc: i32,
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
            "sd" => case.sd = payload,
            "chain" => case.chain = payload,
            "min_pct" => case.min_pct = payload.parse().unwrap(),
            "rc" => case.rc = payload.parse().unwrap(),
            "output" => case.output = payload,
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|n| (n.clone(), map.remove(&n).unwrap()))
        .collect()
}

#[test]
fn every_lifted_interval_list_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 3, "case count");
    for (name, case) in &cases {
        let got = lift_over_interval_list(&case.input, &case.sd, &case.chain, case.min_pct)
            .expect("liftover");
        assert_eq!(got.output, case.output, "output for {name}");
        assert_eq!(got.return_code, case.rc, "rc for {name}");
    }
}
