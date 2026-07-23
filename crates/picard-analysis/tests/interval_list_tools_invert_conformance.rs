//! Conformance for `IntervalListTools` invert-dependent paths (SUBTRACT / SYMDIFF / the INVERT
//! option) against Picard 3.4.0.
//!
//! These all route through `IntervalList.invert`. As with the other slices, the tool adds a `@PG`
//! (stripped before comparing). The `@HD` differs by action and input count: SUBTRACT and the
//! single-`INPUT` INVERT case are emitted verbatim (no `SO`), while SYMDIFF's `union` yields
//! `SO:unsorted`.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::interval_list_tools::{interval_list_tools, Action, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/interval_list_tools_invert.txt.gz");
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

fn strip_pg(text: &str) -> String {
    text.lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(|l| format!("{l}\n"))
        .collect()
}

#[derive(Default)]
struct Case {
    input1: String,
    input2: String,
    action: String,
    invert: bool,
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
            "input1" => case.input1 = payload,
            "input2" => case.input2 = payload,
            "action" => case.action = payload,
            "invert" => case.invert = payload == "true",
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
fn every_invert_path_is_byte_identical_after_stripping_pg() {
    let cases = cases();
    assert_eq!(cases.len(), 3, "case count");
    for (name, case) in &cases {
        let opts = Options {
            action: match case.action.as_str() {
                "CONCAT" => Action::Concat,
                "SUBTRACT" => Action::Subtract,
                "SYMDIFF" => Action::Symdiff,
                other => panic!("unexpected action {other}"),
            },
            invert: case.invert,
            ..Options::default()
        };
        // SUBTRACT/SYMDIFF take input2 as SECOND_INPUT; the INVERT/CONCAT case has no second input.
        let second: Vec<&str> = if case.input2.is_empty() {
            vec![]
        } else {
            vec![case.input2.as_str()]
        };
        let got = interval_list_tools(&[case.input1.as_str()], &second, &opts).expect("tool");
        assert_eq!(got, strip_pg(&case.output), "{name}");
    }
}
