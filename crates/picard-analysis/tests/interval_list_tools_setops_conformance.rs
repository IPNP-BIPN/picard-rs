//! Conformance for `IntervalListTools` set-ops (INTERSECT / OVERLAPS) against Picard 3.4.0.
//!
//! INTERSECT reduces two `INPUT`s by `IntervalList.intersection`; OVERLAPS keeps whole `INPUT`
//! intervals overlapping any `SECOND_INPUT` interval. As with the CONCAT/UNION slice, the tool adds a
//! `@PG` whose `CL` is the command line, so the comparison strips `@PG` from Picard's output. The
//! `@HD` differs by action (INTERSECT is emitted verbatim without `SO`; OVERLAPS is `SO:unsorted`),
//! which is part of the byte-identity check.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::interval_list_tools::{interval_list_tools, Action, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/interval_list_tools_setops.txt.gz");
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
fn every_set_op_is_byte_identical_after_stripping_pg() {
    let cases = cases();
    assert_eq!(cases.len(), 2, "case count");
    for (name, case) in &cases {
        let opts = Options {
            action: match case.action.as_str() {
                "INTERSECT" => Action::Intersect,
                "OVERLAPS" => Action::Overlaps,
                other => panic!("unexpected action {other}"),
            },
            ..Options::default()
        };
        // INTERSECT takes both files as INPUT; OVERLAPS takes input1 as INPUT and input2 as
        // SECOND_INPUT.
        let got = match opts.action {
            Action::Intersect => {
                interval_list_tools(&[&case.input1, &case.input2], &[], &opts).expect("tool")
            }
            Action::Overlaps => {
                interval_list_tools(&[&case.input1], &[&case.input2], &opts).expect("tool")
            }
            _ => unreachable!(),
        };
        assert_eq!(got, strip_pg(&case.output), "{name}");
    }
}
