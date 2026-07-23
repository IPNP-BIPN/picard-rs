//! Full-file conformance for `MergeBamAlignment` (unpaired and paired, coordinate output) against
//! Picard 3.4.0.
//!
//! Each case carries the reference `.dict`, the unmapped and aligned SAM, and the complete output SAM
//! Picard produced (header included). The port runs `merge_bam_alignment` on the same inputs and must
//! reproduce the whole file byte-for-byte. Because the port reads the committed `.dict`, the `@SQ`
//! line (with its `M5` and absolute `UR` path) matches without any canonicalization.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::merge_bam_alignment::merge_bam_alignment;

const REF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/merge_bam_alignment_full.txt.gz");
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
    unmapped: String,
    aligned: String,
    full: String,
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
            "unmapped" => case.unmapped = payload,
            "aligned" => case.aligned = payload,
            "full" => case.full = payload,
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
fn the_whole_output_file_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 3, "case count");
    let reference_bases = HashMap::from([("chr1".to_string(), REF.to_vec())]);
    for (name, case) in &cases {
        let got = merge_bam_alignment(&case.dict, &case.unmapped, &case.aligned, &reference_bases)
            .unwrap();
        assert_eq!(got, case.full, "{name}");
    }
}
