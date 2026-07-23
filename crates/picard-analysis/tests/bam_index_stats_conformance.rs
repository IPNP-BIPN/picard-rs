//! Conformance for `BamIndexStats` against Picard 3.4.0.
//!
//! Each case carries a BAM (hex) and the exact stdout `BamIndexStats` printed for it. The port runs
//! `bam_index_stats` on the same bytes and must reproduce the text byte-for-byte.

use std::io::Read;

use picard_analysis::bam_index_stats::bam_index_stats;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("bam_index_stats.txt.gz");
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

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len() / 2)
        .map(|i| u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
}

struct Case {
    name: String,
    bam: Vec<u8>,
    stats: String,
}

fn cases() -> Vec<Case> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut bams = std::collections::HashMap::new();
    let mut stats = std::collections::HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let name = it.next().unwrap().to_string();
        let payload = it.next().unwrap_or("");
        match kind {
            "bam" => {
                if !order.contains(&name) {
                    order.push(name.clone());
                }
                bams.insert(name, unhex(payload));
            }
            "stats" => {
                stats.insert(name, unescape(payload));
            }
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|name| Case {
            bam: bams[&name].clone(),
            stats: stats[&name].clone(),
            name,
        })
        .collect()
}

#[test]
fn every_bam_index_stats_case_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 2, "case count");
    for case in &cases {
        let got = bam_index_stats(&case.bam).expect("stats");
        assert_eq!(got, case.stats, "{}", case.name);
    }
}
