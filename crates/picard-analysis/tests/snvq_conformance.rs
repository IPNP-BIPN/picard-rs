//! Conformance for `CollectQualityYieldMetricsSNVQ` against Picard 3.4.0.
//!
//! One case carries the exact SAM (reads tagged with the per-alternate-base qualities qa/qc/qg/qt)
//! and the metrics file Picard produced. The port parses the same SAM, runs the collector, and must
//! reproduce the metrics row after the two header lines are canonicalized. The case exercises a PF
//! read spanning the quality thresholds, a vendor-fail read (counted in TOTAL_* but not PF_*), an
//! `N` base (unequal to all four alternates, so four SNVQ observations), and secondary/supplementary
//! reads (excluded by default).

use std::io::Read;

use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::snvq::SnvqCollector;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/snvq.txt.gz");
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

fn payload(kind: &str) -> String {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .find_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let _case = it.next()?;
            let p = it.next().unwrap_or("");
            (k == kind).then(|| unescape(p))
        })
        .unwrap_or_else(|| panic!("no {kind} row"))
}

const CANONICALIZED: &[&str] = &["# CollectQualityYieldMetricsSNVQ ", "# Started on:"];

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
fn the_metrics_row_is_byte_identical() {
    let (_, records) = read_sam(&payload("sam")).expect("parse sam");
    let mut collector = SnvqCollector::new(false, false);
    for rec in &records {
        collector.accept(rec);
    }
    let metrics = collector.finish();

    let mut file = MetricsFile::new();
    file.add_header("CollectQualityYieldMetricsSNVQ <command line>");
    file.add_header("Started on: <timestamp>");
    file.add_metric(&metrics);
    let ours = canonicalize(&file.write());
    let theirs = canonicalize(&payload("metrics"));

    if ours != theirs {
        let at = ours
            .lines()
            .zip(theirs.lines())
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        panic!(
            "first difference at line {at}\n  picard: {}\n  ours  : {}",
            theirs.lines().nth(at).unwrap_or("<end>"),
            ours.lines().nth(at).unwrap_or("<end>")
        );
    }
}
