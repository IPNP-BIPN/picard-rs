//! Conformance for `MeanQualityByCycle` and `CollectBaseDistributionByCycle` against Picard.
//!
//! Both tools run over **identical inputs** in one harness. That is what makes the second
//! one's cost a delta rather than a separate measurement: they are stratum-mates in the
//! calibration gate, and any shared machinery is paid once here too.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::cycle::{CollectBaseDistributionByCycle, MeanQualityByCycle};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cycle.txt.gz");
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

fn canonicalize(text: &str, tool: &str) -> (String, usize) {
    let prefixes = [format!("# {tool} "), "# Started on:".to_string()];
    let mut fired = std::collections::BTreeSet::new();
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        match prefixes.iter().position(|p| line.starts_with(p.as_str())) {
            Some(i) => {
                fired.insert(i);
                out.push_str("<canonicalized>\n");
            }
            None => out.push_str(line),
        }
    }
    (out, fired.len())
}

fn records(bam_hex: &str) -> Vec<htsjdk_bam::record::BamRecord> {
    let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
    BamReader::new(&plain)
        .unwrap()
        .map(|r| r.unwrap())
        .collect()
}

/// The gate for both tools at once.
#[test]
fn both_cycle_tools_match_picard() {
    let bams = rows("bam");
    assert!(bams.len() >= 7, "expected at least 7 cases");

    let mean = rows("MeanQualityByCycle");
    let dist = rows("CollectBaseDistributionByCycle");
    assert_eq!(bams.len(), mean.len());
    assert_eq!(bams.len(), dist.len());

    let mut failures = Vec::new();
    for (i, (name, bam_hex)) in bams.iter().enumerate() {
        let recs = records(bam_hex);

        // MeanQualityByCycle: histograms only, no metric rows.
        let mut q = MeanQualityByCycle::default();
        for r in &recs {
            q.accept(r);
        }
        let mut file = MetricsFile::new();
        file.add_header("MeanQualityByCycle <command line>");
        file.add_header("Started on: <timestamp>");
        file.histograms = q.finish();
        let (ours, _) = canonicalize(&file.write(), "MeanQualityByCycle");
        let (theirs, fired) = canonicalize(&unescape(&mean[i].1), "MeanQualityByCycle");
        assert_eq!(name, &mean[i].0, "case lists out of step");
        if ours != theirs {
            let at = ours
                .lines()
                .zip(theirs.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            failures.push(format!(
                "MeanQualityByCycle/{name}: line {}\n  ours:   {:?}\n  picard: {:?}",
                at + 1,
                ours.lines().nth(at).unwrap_or("<eof>"),
                theirs.lines().nth(at).unwrap_or("<eof>")
            ));
        } else {
            assert_eq!(fired, 2, "{name}: both rules must fire");
        }

        // CollectBaseDistributionByCycle: metric rows only, no histogram.
        let mut d = CollectBaseDistributionByCycle::default();
        for r in &recs {
            d.accept(r);
        }
        let mut file = MetricsFile::new();
        file.add_header("CollectBaseDistributionByCycle <command line>");
        file.add_header("Started on: <timestamp>");
        for m in d.finish() {
            file.add_metric(&m);
        }
        let (ours, _) = canonicalize(&file.write(), "CollectBaseDistributionByCycle");
        let (theirs, _) = canonicalize(&unescape(&dist[i].1), "CollectBaseDistributionByCycle");
        if ours != theirs {
            let at = ours
                .lines()
                .zip(theirs.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            failures.push(format!(
                "CollectBaseDistributionByCycle/{name}: line {}\n  ours:   {:?}\n  picard: {:?}",
                at + 1,
                ours.lines().nth(at).unwrap_or("<eof>"),
                theirs.lines().nth(at).unwrap_or("<eof>")
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "{} divergences:\n{}",
        failures.len(),
        failures.join("\n\n")
    );
}

/// The corpus must exercise the branches both ports claim to handle.
#[test]
fn the_corpus_covers_the_branches_that_matter() {
    let names: Vec<String> = rows("bam").into_iter().map(|(n, _)| n).collect();
    for expected in [
        "single",
        "both_strands",
        "paired",
        "varied_lengths",
        "original_qualities",
        "flagged",
        "mixed_case_bases",
    ] {
        assert!(names.iter().any(|n| n == expected), "missing {expected}");
    }
}
