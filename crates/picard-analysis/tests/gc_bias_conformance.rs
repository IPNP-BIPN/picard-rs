//! Conformance for `CollectGcBiasMetrics` against Picard 3.4.0.
//!
//! Each case carries the input BAM, the reference FASTA, and **both** metrics files Picard
//! produced: the detail file of 101 GC bins and the summary. Two outputs from one run, so both
//! are compared rather than the convenient one.
//!
//! The cases exist to hold the findings in `src/gc.rs` in place. `forward_and_reverse_same_bases`
//! pins that the two strands are charged to different windows; `past_last_window_start` pins that
//! a read whose window was never computed lands in GC bin 0; `reverse_near_contig_start` pins
//! that such a read is dropped from the bins while still counting as aligned.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::gc::{GcBiasMetricsCollector, DEFAULT_SCAN_WINDOW_SIZE};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gc_bias.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn rows(kind: &str) -> Vec<(String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let name = it.next()?.to_string();
            let payload = it.next()?.to_string();
            (k == kind).then_some((name, payload))
        })
        .collect()
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

fn from_hex(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
}

const CANONICALIZED: &[&str] = &["# CollectGcBiasMetrics ", "# Started on:"];

fn canonicalize(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        if CANONICALIZED.iter().any(|p| line.starts_with(*p)) {
            out.push_str("<canonicalized>\n");
        } else {
            out.push_str(line);
        }
    }
    out
}

#[test]
fn every_detail_and_summary_file_matches_picards() {
    let bams = rows("bam");
    let fastas = rows("fasta");
    let details = rows("detail");
    let summaries = rows("summary");
    assert_eq!(bams.len(), details.len());
    assert_eq!(bams.len(), summaries.len());
    assert!(
        bams.len() >= 9,
        "expected at least 9 cases, got {}",
        bams.len()
    );

    let mut failures = Vec::new();
    for (i, (name, bam_hex)) in bams.iter().enumerate() {
        let reference = &fastas[i].1;
        let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
        let reader = BamReader::new(&plain).unwrap_or_else(|e| panic!("{name}: {e:?}"));

        let contigs = vec![reference.as_bytes().to_vec()];
        let mut collector = GcBiasMetricsCollector::new(&contigs, DEFAULT_SCAN_WINDOW_SIZE, false);
        for record in reader {
            let rec = record.unwrap_or_else(|e| panic!("{name}: {e:?}"));
            collector.accept(&rec, reference.as_bytes());
        }

        let (our_detail, our_summary) = match collector.rows() {
            Some(v) => v,
            None => {
                failures.push(format!("{name}: no rows, but Picard wrote a file"));
                continue;
            }
        };

        let mut detail_file = MetricsFile::new();
        detail_file.add_header("CollectGcBiasMetrics <command line>");
        detail_file.add_header("Started on: <timestamp>");
        for row in &our_detail {
            detail_file.add_metric(row);
        }
        let mut summary_file = MetricsFile::new();
        summary_file.add_header("CollectGcBiasMetrics <command line>");
        summary_file.add_header("Started on: <timestamp>");
        summary_file.add_metric(&our_summary);

        for (label, ours, theirs_raw) in [
            ("detail", detail_file.write(), &details[i].1),
            ("summary", summary_file.write(), &summaries[i].1),
        ] {
            let ours = canonicalize(&ours);
            let theirs = canonicalize(&unescape(theirs_raw));
            if ours != theirs {
                let at = ours
                    .lines()
                    .zip(theirs.lines())
                    .position(|(a, b)| a != b)
                    .unwrap_or(0);
                failures.push(format!(
                    "{name} [{label}]: first difference at line {at}\n  picard: {}\n  ours  : {}",
                    theirs.lines().nth(at).unwrap_or("<end>"),
                    ours.lines().nth(at).unwrap_or("<end>")
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} divergences over {} cases:\n{}",
        failures.len(),
        bams.len(),
        failures.join("\n")
    );
}
