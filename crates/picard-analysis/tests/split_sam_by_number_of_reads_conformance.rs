//! Conformance for `SplitSamByNumberOfReads` against Picard 3.4.0.
//!
//! Each case carries the input SAM, the `KEY=VALUE` options, Picard's return code, and each output
//! shard's SAM in order (`shard_0001.sam`, `shard_0002.sam`, ...). The tool adds no `@PG` and writes
//! the input header verbatim, so every shard is compared raw. The port runs
//! `split_sam_by_number_of_reads` on the same input and options and must reproduce the shards
//! byte-for-byte.

use std::io::Read;

use picard_analysis::split_sam_by_number_of_reads::{split_sam_by_number_of_reads, SplitOptions};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("split_sam_by_number_of_reads.txt.gz");
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

fn parse_opts(opts: &str) -> SplitOptions {
    let mut o = SplitOptions::default();
    for tok in opts.split_whitespace() {
        let (k, v) = tok.split_once('=').expect("KEY=VALUE");
        let n: i64 = v.parse().expect("int");
        match k {
            "SPLIT_TO_N_FILES" => o.split_to_n_files = n,
            "SPLIT_TO_N_READS" => o.split_to_n_reads = n,
            "TOTAL_READS" => o.total_reads_in_input = n,
            other => panic!("unhandled option {other}"),
        }
    }
    o
}

struct Case {
    name: String,
    opts: String,
    input: String,
    shards: Vec<String>,
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
                opts: String::new(),
                input: String::new(),
                shards: Vec::new(),
            }
        });
        match kind {
            "opts" => entry.opts = payload,
            "input" => entry.input = payload,
            "rc" => {} // all corpus cases return 0
            "shard" => entry.shards.push(payload),
            other => panic!("unexpected row kind {other}"),
        }
    }
    order.into_iter().map(|n| map.remove(&n).unwrap()).collect()
}

#[test]
fn every_split_case_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 5, "split case count");
    for case in &cases {
        let opts = parse_opts(&case.opts);
        let shards = split_sam_by_number_of_reads(&case.input, &opts).expect("split");
        assert_eq!(
            shards.len(),
            case.shards.len(),
            "{}: shard count {} vs {}",
            case.name,
            shards.len(),
            case.shards.len()
        );
        for (i, (got, want)) in shards.iter().zip(&case.shards).enumerate() {
            assert_eq!(got, want, "{} shard {i}", case.name);
        }
    }
}
