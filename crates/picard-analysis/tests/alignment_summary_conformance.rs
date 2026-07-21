//! Conformance for `CollectAlignmentSummaryMetrics` against Picard 3.4.0.
//!
//! Same discipline as the other three suites: each case carries the input BAM, the reference
//! FASTA and the metrics file Picard produced from them, and the two run-time header lines are
//! canonicalized by name rather than skipped silently.
//!
//! This is the second member of the `histogram` + `multi_level` + `single_pass` stratum, after
//! `CollectInsertSizeMetrics`. Its cost is what the archetype delta at the large end is measured
//! from; the number is recorded in `docs/decisions/0003`.
//!
//! Several cases exist for one reason: to hold the `BAD_CYCLES` divergence in place. `collide`
//! and `distinct` are the probe from `tools/asm-conformance/badcycle_probe.sh`, and
//! `three_blocks`, `insertion` and `split_match` reach the same code by three other routes. A
//! port that indexed the cycle by the read position - which is what Picard's own parameter name
//! asks for - fails all five.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricsFile};
use picard_analysis::alignment_summary::{GroupCollector, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/alignment_summary.txt.gz");
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

/// The only two lines that carry run-time state. Everything else is compared raw.
const CANONICALIZED: &[&str] = &["# CollectAlignmentSummaryMetrics ", "# Started on:"];

fn canonicalize(text: &str) -> (String, usize) {
    let mut fired = std::collections::BTreeSet::new();
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        match CANONICALIZED.iter().find(|p| line.starts_with(**p)) {
            Some(prefix) => {
                fired.insert(*prefix);
                out.push_str("<canonicalized>\n");
            }
            None => out.push_str(line),
        }
    }
    (out, fired.len())
}

/// The gate: every metrics file matches Picard's, under the two declared rules.
#[test]
fn every_metrics_file_matches_picards() {
    let bams = rows("bam");
    let fastas = rows("fasta");
    let metrics = rows("metrics");
    assert_eq!(bams.len(), metrics.len(), "corpus is inconsistent");
    assert_eq!(bams.len(), fastas.len(), "corpus is inconsistent");
    assert!(
        bams.len() >= 20,
        "expected at least 20 cases, got {}",
        bams.len()
    );

    let mut failures = Vec::new();
    for (i, (name, bam_hex)) in bams.iter().enumerate() {
        let (fasta_name, reference) = &fastas[i];
        let (metrics_name, expected_raw) = &metrics[i];
        assert_eq!(name, fasta_name, "case lists out of step");
        assert_eq!(name, metrics_name, "case lists out of step");

        let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
        let reader = BamReader::new(&plain).unwrap_or_else(|e| panic!("{name}: {e:?}"));

        let mut collector = GroupCollector::new(Options::default());
        for record in reader {
            let rec = record.unwrap_or_else(|e| panic!("{name}: {e:?}"));
            collector.accept(&rec, Some(reference.as_bytes()));
        }
        collector.finish().unwrap_or_else(|e| panic!("{name}: {e}"));

        let mut file = MetricsFile::new();
        file.add_header("CollectAlignmentSummaryMetrics <command line>");
        file.add_header("Started on: <timestamp>");
        for row in collector.rows() {
            file.add_metric(&row);
        }
        for (label, histogram) in collector.read_length_histograms() {
            file.histograms.push(OutHistogram {
                bin_label: "READ_LENGTH".to_string(),
                value_label: label.to_string(),
                // The read length is an int in Java, so the histogram's key class is Integer
                // even though the port carries the keys as f64.
                key_class: "java.lang.Integer".to_string(),
                bins: histogram
                    .bins()
                    .map(|(id, count)| (format!("{}", id as i64), count))
                    .collect(),
            });
        }

        let (ours, _) = canonicalize(&file.write());
        let (theirs, fired) = canonicalize(&unescape(expected_raw));
        assert_eq!(
            fired,
            CANONICALIZED.len(),
            "{name}: a canonicalized line went missing"
        );

        if ours != theirs {
            let at = ours
                .lines()
                .zip(theirs.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            let ours_line = ours.lines().nth(at).unwrap_or("<end>");
            let theirs_line = theirs.lines().nth(at).unwrap_or("<end>");
            failures.push(format!(
                "{name}: first difference at line {at}\n  picard: {theirs_line}\n  ours  : {ours_line}"
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "{} of {} cases diverge:\n{}",
        failures.len(),
        bams.len(),
        failures.join("\n")
    );
}

/// The divergence, asserted directly on the goldens rather than only through the whole file, so
/// that its removal is loud. `collide` and `distinct` are the same read against references that
/// differ in one base position; if the cycle were indexed by the read they would agree.
#[test]
fn the_bad_cycles_divergence_is_still_present_in_the_goldens() {
    let metrics = rows("metrics");
    let field = |case: &str, column: &str| -> String {
        let (_, raw) = metrics.iter().find(|(n, _)| n == case).expect(case);
        let text = unescape(raw);
        let mut lines = text.lines().skip_while(|l| !l.starts_with("CATEGORY"));
        let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
        let values: Vec<&str> = lines.next().expect("row").split('\t').collect();
        let i = header.iter().position(|h| *h == column).expect(column);
        values[i].to_string()
    };

    assert_eq!(field("collide", "BAD_CYCLES"), "1");
    assert_eq!(field("distinct", "BAD_CYCLES"), "2");
    assert_eq!(
        field("collide", "PF_HQ_MEDIAN_MISMATCHES"),
        field("distinct", "PF_HQ_MEDIAN_MISMATCHES"),
        "the control: both reads carry the same number of mismatches"
    );
}
