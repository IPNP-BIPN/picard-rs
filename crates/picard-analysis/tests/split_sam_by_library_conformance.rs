//! Conformance for `SplitSamByLibrary` against Picard 3.4.0.
//!
//! Each case carries the input SAM and every output file as `file:<base name>` rows in output-name
//! order (`<library>` / `unknown`). The tool adds no `@PG` and writes each output's header with the
//! `@RG` block filtered to that library, so every output is compared raw. The port runs
//! `split_sam_by_library` and must reproduce the file base names and their SAM bytes.

use std::io::Read;

use picard_analysis::split_sam_by_library::split_sam_by_library;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("split_sam_by_library.txt.gz");
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
    input: String,
    files: Vec<(String, String)>,
}

fn cases() -> Vec<Case> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut map: std::collections::HashMap<String, Case> = std::collections::HashMap::new();
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
            Case {
                name: name.clone(),
                input: String::new(),
                files: Vec::new(),
            }
        });
        if kind == "input" {
            entry.input = payload;
        } else if kind == "rc" {
            // all corpus cases return 0
        } else if let Some(base) = kind.strip_prefix("file:") {
            entry.files.push((base.to_string(), payload));
        } else {
            panic!("unexpected row kind {kind}");
        }
    }
    order.into_iter().map(|n| map.remove(&n).unwrap()).collect()
}

#[test]
fn every_split_by_library_case_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 5, "case count");
    for case in &cases {
        let out = split_sam_by_library(&case.input).expect("split");
        assert_eq!(
            out.len(),
            case.files.len(),
            "{}: file count {} vs {}",
            case.name,
            out.len(),
            case.files.len()
        );
        for ((got_name, got_sam), (want_name, want_sam)) in out.iter().zip(&case.files) {
            assert_eq!(got_name, want_name, "{}: file name", case.name);
            assert_eq!(got_sam, want_sam, "{} file {got_name}", case.name);
        }
    }
}
