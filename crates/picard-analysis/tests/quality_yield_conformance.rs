//! Conformance for `CollectQualityYieldMetrics` against Picard 3.4.0.
//!
//! Goldens from `tools/metrics-conformance/QualityYieldDump.java` in the pinned picard-rs
//! oracle. Each case carries the input BAM **and** the metrics file Picard produced from it, so
//! the port is measured on the same bytes rather than on a reconstruction.
//!
//! Two header lines are canonicalized, and named here rather than skipped silently:
//!
//! - `# CollectQualityYieldMetrics INPUT=/tmp/...` records absolute temp paths chosen by the
//!   JVM at run time, which cannot be reproduced by construction.
//! - `# Started on: ...` is wall-clock time.
//!
//! Everything else, including both `## htsjdk.samtools.metrics.StringHeader` class lines, the
//! metrics class line, the column header and every data row, is compared raw. The test reports
//! which rules fired, so "identical" is never claimed without saying what was excused.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::{Options, QualityYieldMetricsCollector};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/quality_yield.txt.gz");
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

/// The two lines whose content cannot be reproduced by construction, with the reason.
const CANONICALIZED: &[(&str, &str)] = &[
    (
        "# CollectQualityYieldMetrics ",
        "records absolute temp paths chosen by the JVM at run time",
    ),
    ("# Started on:", "wall-clock time of the run"),
];

/// Replaces the canonicalized lines and returns which rules fired.
fn canonicalize(text: &str) -> (String, Vec<&'static str>) {
    let mut fired = Vec::new();
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        match CANONICALIZED.iter().find(|(p, _)| line.starts_with(p)) {
            Some((prefix, _)) => {
                if !fired.contains(prefix) {
                    fired.push(*prefix);
                }
                out.push_str("<canonicalized>\n");
            }
            None => out.push_str(line),
        }
    }
    (out, fired)
}

fn options_for(case: &str) -> Options {
    match case {
        "secondary_included" => Options {
            include_secondary_alignments: true,
            ..Options::default()
        },
        "supplemental_included" => Options {
            include_supplemental_alignments: true,
            ..Options::default()
        },
        "original_qualities_off" => Options {
            use_original_qualities: false,
            ..Options::default()
        },
        _ => Options::default(),
    }
}

/// The gate: every metrics file matches Picard's, under the two declared rules.
#[test]
fn every_metrics_file_matches_picards() {
    let bams = rows("bam");
    let metrics = rows("metrics");
    assert_eq!(bams.len(), metrics.len(), "corpus is inconsistent");
    assert!(bams.len() >= 10, "expected at least 10 cases");

    let mut failures = Vec::new();
    for ((name, bam_hex), (metrics_name, expected_raw)) in bams.iter().zip(&metrics) {
        assert_eq!(name, metrics_name, "case lists out of step");

        let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
        let reader = BamReader::new(&plain).unwrap_or_else(|e| panic!("{name}: {e:?}"));

        let mut collector = QualityYieldMetricsCollector::new(options_for(name));
        for record in reader {
            collector.accept(&record.unwrap_or_else(|e| panic!("{name}: {e:?}")));
        }
        collector.finish();

        let mut file = MetricsFile::new();
        // The two header lines Picard writes; both are canonicalized away, but they must be
        // present so the file has the right shape.
        file.add_header("CollectQualityYieldMetrics <command line>");
        file.add_header("Started on: <timestamp>");
        file.add_metric(collector.metrics());

        let (ours, _) = canonicalize(&file.write());
        let (theirs, fired) = canonicalize(&unescape(expected_raw));

        if ours != theirs {
            let at = ours
                .lines()
                .zip(theirs.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            failures.push(format!(
                "{name}: first differs at line {}\n  ours:   {:?}\n  picard: {:?}",
                at + 1,
                ours.lines().nth(at).unwrap_or("<eof>"),
                theirs.lines().nth(at).unwrap_or("<eof>"),
            ));
        } else {
            assert_eq!(
                fired.len(),
                2,
                "{name}: both canonicalization rules must fire, or one is dead"
            );
        }
    }

    assert!(
        failures.is_empty(),
        "{} of {} cases diverge from Picard:\n{}",
        failures.len(),
        bams.len(),
        failures.join("\n")
    );
}

/// Canonicalization must not be able to hide a real difference. If a data row changes, the
/// comparison must still fail.
#[test]
fn canonicalization_does_not_hide_a_data_difference() {
    let (_, raw) = rows("metrics")
        .into_iter()
        .find(|(n, _)| n == "one_read")
        .expect("the one_read case");
    let text = unescape(&raw);
    let (a, _) = canonicalize(&text);
    let (b, _) = canonicalize(&text.replace("\t50\t", "\t51\t"));
    assert_ne!(a, b, "a changed data row must survive canonicalization");
}

/// The corpus must exercise the branches the port claims to handle.
#[test]
fn the_corpus_covers_the_branches_that_matter() {
    let names: Vec<String> = rows("metrics").into_iter().map(|(n, _)| n).collect();
    for expected in [
        "empty",
        "vendor_failed",
        "secondary_supplementary",
        "secondary_included",
        "supplemental_included",
        "original_qualities",
        "original_qualities_off",
        "varied_lengths",
    ] {
        assert!(
            names.iter().any(|n| n == expected),
            "the corpus must contain the {expected} case"
        );
    }
}
