//! Conformance for `CollectRnaSeqMetrics` against Picard 3.4.0.
//!
//! Each case carries its exact inputs (the SAM, the refFlat, and the ribosomal interval_list when
//! there is one) and the `.rna_metrics` file Picard produced from them. The port parses the same
//! bytes, runs the collector, and must reproduce the metrics file after the two header lines
//! (command line and timestamp) are canonicalized.
//!
//! Two cases:
//!   `basic`    the reference's own testBasic: its 451bp transcript is below MINIMUM_LENGTH, so the
//!              histogram is empty and the MEDIAN_* metrics are 0. Pins the whole metrics row.
//!   `coverage` three long transcripts of deliberately different depth, so the normalized_coverage
//!              histogram is non-empty and its floating-point fold over the transcripts is
//!              order-sensitive in principle. Measurement showed the fold order is unobservable at
//!              printed precision (decision 0005), so the port folds in a deterministic content
//!              order; this case pins the whole histogram, the MEDIAN_* metrics, and the empty
//!              RIBOSOMAL_BASES column of a no-ribosomal run.

use std::io::Read;

use htsjdk_bam::interval::{Interval, IntervalList};
use htsjdk_bam::overlap::OverlapDetector;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::refflat;
use picard_analysis::rnaseq_metrics::{
    RnaSeqMetricsCollector, StrandSpecificity, DEFAULT_END_BIAS_BASES, DEFAULT_MINIMUM_LENGTH,
    DEFAULT_RRNA_FRAGMENT_PERCENTAGE,
};

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/rnaseq_metrics.txt.gz");
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

/// Every payload for one case, keyed by kind.
fn case(name: &str) -> std::collections::HashMap<String, String> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let kind = it.next()?;
            let c = it.next()?;
            let payload = it.next().unwrap_or("");
            (c == name).then(|| (kind.to_string(), unescape(payload)))
        })
        .collect()
}

const CANONICALIZED: &[&str] = &["# CollectRnaSeqMetrics ", "# Started on:"];

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

fn run_case(name: &str, strand: StrandSpecificity, ignore_sequence: &[&str]) -> String {
    let c = case(name);
    let sam_text = c.get("sam").expect("sam");
    let (header, records) = read_sam(sam_text).expect("parse sam");
    let seq_names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();

    let recognized: std::collections::HashSet<&str> =
        seq_names.iter().map(|s| s.as_str()).collect();
    let gene_overlap = refflat::load(c.get("refflat").expect("refflat"), |ctg| {
        recognized.contains(ctg)
    })
    .expect("refflat load");

    // Ribosomal detector: present only when the case carries an interval_list.
    let rib_text = c.get("ribosomal").map(|s| s.as_str()).unwrap_or("");
    let ribosomal_present = !rib_text.trim().is_empty();
    let mut ribosomal_overlap: OverlapDetector<Interval> = OverlapDetector::new(0, 0);
    if ribosomal_present {
        let list = IntervalList::parse_body(seq_names.clone(), rib_text).expect("interval_list");
        for iv in list.uniqued(true).intervals {
            let (contig, start, end) = (iv.contig.clone(), iv.start, iv.end);
            ribosomal_overlap.add(&contig, start, end, iv);
        }
    }
    let ribosomal_initial = if ribosomal_present { Some(0i64) } else { None };

    let ignored_indices: Vec<i32> = ignore_sequence
        .iter()
        .map(|name| {
            seq_names
                .iter()
                .position(|s| s == name)
                .map(|i| i as i32)
                .unwrap_or_else(|| panic!("unrecognized IGNORE_SEQUENCE {name}"))
        })
        .collect();

    let mut collector = RnaSeqMetricsCollector::new(
        &seq_names,
        &gene_overlap,
        &ribosomal_overlap,
        ribosomal_initial,
        &ignored_indices,
        DEFAULT_MINIMUM_LENGTH,
        strand,
        DEFAULT_RRNA_FRAGMENT_PERCENTAGE,
        DEFAULT_END_BIAS_BASES,
    );
    for rec in &records {
        collector.accept(rec);
    }
    let (metrics, histogram) = collector.finish();

    let mut file = MetricsFile::new();
    file.add_header("CollectRnaSeqMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    file.add_metric(&metrics);
    file.histograms.push(histogram);
    file.write()
}

fn assert_matches(name: &str, strand: StrandSpecificity, ignore_sequence: &[&str]) {
    let ours = canonicalize(&run_case(name, strand, ignore_sequence));
    let theirs = canonicalize(case(name).get("metrics").expect("golden"));
    if ours != theirs {
        let at = ours
            .lines()
            .zip(theirs.lines())
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        panic!(
            "{name}: first difference at line {at}\n  picard: {}\n  ours  : {}",
            theirs.lines().nth(at).unwrap_or("<end>"),
            ours.lines().nth(at).unwrap_or("<end>")
        );
    }
}

#[test]
fn basic_metrics_row_is_byte_identical() {
    assert_matches(
        "basic",
        StrandSpecificity::SecondReadTranscriptionStrand,
        &["chrM"],
    );
}

#[test]
fn coverage_histogram_and_medians_are_byte_identical() {
    assert_matches("coverage", StrandSpecificity::None, &[]);
}
