//! Conformance for `MergeSamFiles` (shared-header inputs) against Picard 3.4.0.
//!
//! Each case carries the `SORT_ORDER`, the input SAMs, and the merged output. The tool adds no `@PG`
//! and the only header change is the group order set to `none`, so the output compares raw. The port
//! runs `merge_sam_files` on the same inputs and must reproduce the output byte-for-byte.

use std::io::Read;

use picard_analysis::merge_sam_files::merge_sam_files;
use picard_analysis::sort_sam::SortOrder;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("merge_sam_files.txt.gz");
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

struct Case {
    name: String,
    order: SortOrder,
    merge_dictionaries: bool,
    inputs: Vec<String>,
    output: String,
}

#[derive(Default)]
struct Raw {
    so: String,
    msd: bool,
    inputs: Vec<String>,
    output: String,
}

fn cases() -> Vec<Case> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut map: std::collections::HashMap<String, Raw> = std::collections::HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let name = it.next().unwrap().to_string();
        let payload = unescape(it.next().unwrap_or(""));
        let entry = map.entry(name.clone()).or_insert_with(|| {
            order.push(name.clone());
            Raw::default()
        });
        match kind {
            "so" => entry.so = payload,
            "msd" => entry.msd = payload == "true",
            "input" => entry.inputs.push(payload),
            "rc" => {} // all corpus cases return 0
            "output" => entry.output = payload,
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|name| {
            let raw = map.remove(&name).unwrap();
            let order = match raw.so.as_str() {
                "coordinate" => SortOrder::Coordinate,
                "queryname" => SortOrder::Queryname,
                other => panic!("unhandled sort order {other}"),
            };
            Case {
                name,
                order,
                merge_dictionaries: raw.msd,
                inputs: raw.inputs,
                output: raw.output,
            }
        })
        .collect()
}

#[test]
fn every_merge_case_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 12, "case count");
    for case in &cases {
        let refs: Vec<&str> = case.inputs.iter().map(|s| s.as_str()).collect();
        let got = merge_sam_files(&refs, case.order, case.merge_dictionaries).expect("merge");
        assert_eq!(got, case.output, "{}", case.name);
    }
}
