//! Conformance for `QualityScoreDistribution` against Picard 3.4.0.
//!
//! The only suite here whose files have **no metric rows** — the whole body is one or two
//! histogram tables keyed on `java.lang.Byte`. Two things it is the first to exercise: the byte
//! key class, and a union of two histograms whose key sets interleave, which htsjdk re-sorts
//! rather than concatenates.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricsFile};
use picard_analysis::quality_score_distribution::{Options, QualityScoreDistribution};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/quality_score_distribution.txt.gz");
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

const CANONICALIZED: &[&str] = &["# QualityScoreDistribution ", "# Started on:"];

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

/// The arguments each case was generated with.
fn options_for(case: &str) -> Options {
    match case {
        "no_calls_included" => Options {
            include_no_calls: true,
            ..Options::default()
        },
        "pf_reads_only" => Options {
            pf_reads_only: true,
            ..Options::default()
        },
        "aligned_reads_only" => Options {
            aligned_reads_only: true,
            ..Options::default()
        },
        _ => Options::default(),
    }
}

#[test]
fn every_metrics_file_matches_picards() {
    let bams = rows("bam");
    let metrics = rows("metrics");
    assert_eq!(bams.len(), metrics.len());
    assert!(
        bams.len() >= 11,
        "expected at least 11 cases, got {}",
        bams.len()
    );

    let mut failures = Vec::new();
    for (i, (name, bam_hex)) in bams.iter().enumerate() {
        let (metrics_name, expected_raw) = &metrics[i];
        assert_eq!(name, metrics_name, "case lists out of step");

        let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
        let reader = BamReader::new(&plain).unwrap_or_else(|e| panic!("{name}: {e:?}"));

        let mut collector = QualityScoreDistribution::new(options_for(name));
        for record in reader {
            collector.accept(&record.unwrap_or_else(|e| panic!("{name}: {e:?}")));
        }
        let (q_bins, oq_bins) = collector.finish();

        let mut file = MetricsFile::new();
        file.add_header("QualityScoreDistribution <command line>");
        file.add_header("Started on: <timestamp>");
        let histogram = |label: &str, bins: &[(u8, f64)]| OutHistogram {
            bin_label: "QUALITY".to_string(),
            value_label: label.to_string(),
            key_class: "java.lang.Byte".to_string(),
            bins: bins.iter().map(|(k, v)| (k.to_string(), *v)).collect(),
        };
        file.histograms.push(histogram("COUNT_OF_Q", &q_bins));
        // The OQ histogram is added only when non-empty, mirroring `if (!oqHisto.isEmpty())`.
        // The guard is belt and braces: `MetricsFile::write` drops empty histograms anyway, so
        // removing it here changes no byte — verified by sabotage, and recorded so the next
        // reader does not mistake a redundant guard for a load-bearing one.
        if !oq_bins.is_empty() {
            file.histograms.push(histogram("COUNT_OF_OQ", &oq_bins));
        }

        let ours = canonicalize(&file.write());
        let theirs = canonicalize(&unescape(expected_raw));
        if ours != theirs {
            let at = ours
                .lines()
                .zip(theirs.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            failures.push(format!(
                "{name}: first difference at line {at}\n  picard: {:?}\n  ours  : {:?}",
                theirs.lines().nth(at).unwrap_or("<end>"),
                ours.lines().nth(at).unwrap_or("<end>")
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

/// The interleaved key union, asserted on the golden so its removal is loud. A concatenation of
/// the two histograms' key lists would put 30 first; htsjdk puts it where it sorts.
#[test]
fn the_two_histograms_key_sets_interleave_in_the_golden() {
    let metrics = rows("metrics");
    let (_, raw) = metrics
        .iter()
        .find(|(n, _)| n == "with_oq")
        .expect("with_oq");
    let text = unescape(raw);
    let keys: Vec<&str> = text
        .lines()
        .skip_while(|l| !l.starts_with("QUALITY"))
        .skip(1)
        .take_while(|l| !l.trim().is_empty())
        .filter_map(|l| l.split('\t').next())
        .collect();
    assert_eq!(keys, ["2", "3", "30", "40", "41"]);
}
