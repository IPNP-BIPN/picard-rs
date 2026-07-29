//! Conformance for `IntervalListTools` (CONCAT / UNION slice) against Picard 3.4.0.
//!
//! Each case carries two input interval lists, the ACTION/SORT/UNIQUE/DONT_MERGE_ABUTTING options,
//! and the output Picard wrote. IntervalListTools always adds a `@PG` whose `CL` is the command line
//! (non-reproducible), so the comparison strips `@PG` lines from Picard's output and the port emits
//! none; everything else (the `@HD VN:1.6 SO:unsorted` line, the `@SQ` lines, and the interval body)
//! is compared byte-for-byte.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::interval_list_tools::{interval_list_tools, Action, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/interval_list_tools.txt.gz");
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

/// Drop `@PG` lines, which carry the non-reproducible command line.
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
    sort: bool,
    unique: bool,
    dont_merge_abutting: bool,
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
            "sort" => case.sort = payload == "true",
            "unique" => case.unique = payload == "true",
            "dont_merge_abutting" => case.dont_merge_abutting = payload == "true",
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
fn every_interval_list_is_byte_identical_after_stripping_pg() {
    let cases = cases();
    assert_eq!(cases.len(), 4, "case count");
    for (name, case) in &cases {
        let action = match case.action.as_str() {
            "CONCAT" => Action::Concat,
            "UNION" => Action::Union,
            other => panic!("unexpected action {other}"),
        };
        let opts = Options {
            action,
            sort: case.sort,
            unique: case.unique,
            dont_merge_abutting: case.dont_merge_abutting,
            invert: false,
            // This suite predates PADDING and BREAK_BANDS and exercises neither; the two are
            // measured by their own suite, which has no golden yet.
            ..Options::default()
        };
        let got = interval_list_tools(&[&case.input1, &case.input2], &[], &opts).expect("tool");
        assert_eq!(got, strip_pg(&case.output), "{name}");
    }
}
