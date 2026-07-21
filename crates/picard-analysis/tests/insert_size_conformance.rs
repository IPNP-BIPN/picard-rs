//! Conformance for `CollectInsertSizeMetrics` against Picard 3.4.0.
//!
//! Same discipline as the quality-yield suite: each case carries the input BAM **and** the
//! metrics file Picard produced from it, and the two run-time header lines are canonicalized
//! by name rather than skipped silently.
//!
//! This is the first member of the calibration triple. Its cost is recorded in
//! `docs/decisions/0001` so the second and third can be measured against it.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::insert_size::{InsertSizeMetricsCollector, Options};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/insert_size.txt.gz");
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

const CANONICALIZED: &[&str] = &["# CollectInsertSizeMetrics ", "# Started on:"];

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

fn options_for(case: &str) -> Options {
    match case {
        "rare_orientation_pct0" => Options {
            minimum_pct: 0.0,
            ..Options::default()
        },
        "duplicates_included" => Options {
            include_duplicates: true,
            ..Options::default()
        },
        "long_tail_fixed_width" => Options {
            histogram_width: Some(1000),
            ..Options::default()
        },
        "long_tail_deviations" => Options {
            deviations: 2.0,
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

        let mut collector = InsertSizeMetricsCollector::new(options_for(name));
        for record in reader {
            collector.accept(&record.unwrap_or_else(|e| panic!("{name}: {e:?}")));
        }

        let mut file = MetricsFile::new();
        file.add_header("CollectInsertSizeMetrics <command line>");
        file.add_header("Started on: <timestamp>");
        for (metric, histogram) in collector.finish() {
            file.add_metric(&metric);
            file.histograms.push(histogram);
        }

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
            assert_eq!(fired, 2, "{name}: both rules must fire, or one is dead");
        }
    }

    assert!(
        failures.is_empty(),
        "{} of {} cases diverge from Picard:\n{}",
        failures.len(),
        bams.len(),
        failures.join("\n\n")
    );
}

/// Canonicalization must not be able to hide a real difference.
#[test]
fn canonicalization_does_not_hide_a_data_difference() {
    let (_, raw) = rows("metrics")
        .into_iter()
        .find(|(n, _)| n == "one_pair")
        .expect("the one_pair case");
    let text = unescape(&raw);
    let (a, _) = canonicalize(&text);
    let (b, _) = canonicalize(&text.replace("\t300\t", "\t301\t"));
    assert_ne!(a, b, "a changed data row must survive canonicalization");
}

/// The corpus must exercise the branches the port claims to handle.
#[test]
fn the_corpus_covers_the_branches_that_matter() {
    let names: Vec<String> = rows("metrics").into_iter().map(|(n, _)| n).collect();
    for expected in [
        "one_pair",
        "normal",
        "mixed_orientations",
        "rare_orientation",
        "rare_orientation_pct0",
        "duplicates",
        "duplicates_included",
        "long_tail",
        "long_tail_fixed_width",
        "long_tail_deviations",
    ] {
        assert!(
            names.iter().any(|n| n == expected),
            "the corpus must contain the {expected} case"
        );
    }
}
