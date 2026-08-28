//! Conformance for `CollectHsMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/hsmetrics-conformance`. Each case carries the two interval lists, the input
//! as SAM, the metrics table and, where they were asked for, the per-target and per-base files.
//!
//! # What this suite is for
//!
//!  * **the baits and the targets being counted separately**;
//!  * **the three bait columns partitioning the aligned bases**;
//!  * **`--NEAR_DISTANCE` moving the middle one**;
//!  * **the two quality floors emptying the coverage and not the bait counts**;
//!  * **`--CLIP_OVERLAPPING_READS` changing nothing, because the coverage is per locus**;
//!  * **the per-target file's own columns**;
//!  * **the per-base file being a row per target base**;
//!  * **and the derived columns following from the counts.**

use std::io::Read;

use picard_analysis::collect_hs_metrics::{
    derived, placement, target_row, BaitPlacement, DEFAULT_MINIMUM_BASE_QUALITY,
    DEFAULT_MINIMUM_MAPPING_QUALITY,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_hs_metrics.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn table(text: &str, kind: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let body = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = body
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'));
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

fn metrics(text: &str, case: &str) -> std::collections::HashMap<String, String> {
    table(text, "metrics", case).remove(0)
}

fn number(row: &std::collections::HashMap<String, String>, name: &str) -> f64 {
    let value = row.get(name).unwrap_or_else(|| panic!("{name}"));
    if value == "?" || value.is_empty() {
        f64::NAN
    } else {
        value.parse().unwrap_or_else(|_| panic!("{name}={value}"))
    }
}

/// The baits and the targets are counted separately, and a read may be on one and not the other.
#[test]
fn the_baits_and_the_targets_are_counted_separately() {
    let text = corpus();
    let plain = metrics(&text, "plain");
    assert_eq!(number(&plain, "ON_BAIT_BASES"), 60.0);
    assert_eq!(number(&plain, "ON_TARGET_BASES"), 100.0);
    assert_eq!(number(&plain, "BAIT_TERRITORY"), 200.0);
    assert_eq!(number(&plain, "TARGET_TERRITORY"), 180.0);
    // One of the three targets lies outside every bait, which is where the difference comes from.
    let targets = table(&text, "per-target", "per-target");
    assert_eq!(targets.len(), 3);
    let orphan = targets
        .iter()
        .find(|row| row["name"] == "target-orphan")
        .expect("the orphan");
    assert_ne!(orphan["mean_coverage"], "0");
}

/// The three bait columns partition the aligned bases, and --NEAR_DISTANCE moves the middle one.
#[test]
fn the_three_bait_columns_partition_the_aligned_bases() {
    let text = corpus();
    for case in ["plain", "near-distance-zero"] {
        let row = metrics(&text, case);
        let sum = number(&row, "ON_BAIT_BASES")
            + number(&row, "NEAR_BAIT_BASES")
            + number(&row, "OFF_BAIT_BASES");
        assert_eq!(sum, number(&row, "PF_BASES_ALIGNED"), "{case}");
    }
    let plain = metrics(&text, "plain");
    let near_zero = metrics(&text, "near-distance-zero");
    assert_eq!(number(&plain, "NEAR_BAIT_BASES"), 120.0);
    assert_eq!(number(&near_zero, "NEAR_BAIT_BASES"), 0.0);
    assert_eq!(number(&near_zero, "OFF_BAIT_BASES"), 180.0);
    assert_eq!(
        number(&plain, "ON_BAIT_BASES"),
        number(&near_zero, "ON_BAIT_BASES")
    );
    // Which is the port's own placement of a base sixty past a bait's end.
    let baits = [(101, 200), (301, 400)];
    assert_eq!(placement(150, &baits, 250), BaitPlacement::On);
    assert_eq!(placement(260, &baits, 250), BaitPlacement::Near);
    assert_eq!(placement(260, &baits, 0), BaitPlacement::Off);
    assert_eq!(placement(700, &baits, 250), BaitPlacement::Off);
}

/// The derived columns follow from the counts.
#[test]
fn the_derived_columns_follow_from_the_counts() {
    let text = corpus();
    let row = metrics(&text, "plain");
    let ours = derived(
        number(&row, "PF_BASES") as i64,
        number(&row, "PF_BASES_ALIGNED") as i64,
        number(&row, "ON_BAIT_BASES") as i64,
        number(&row, "NEAR_BAIT_BASES") as i64,
        number(&row, "OFF_BAIT_BASES") as i64,
        number(&row, "ON_TARGET_BASES") as i64,
        number(&row, "BAIT_TERRITORY") as i64,
        number(&row, "TARGET_TERRITORY") as i64,
    );
    let close = |ours: f64, name: &str| {
        let theirs = number(&row, name);
        assert!(
            (ours - theirs).abs() < 1e-6,
            "{name}: {ours} against {theirs}"
        );
    };
    close(ours.pct_selected_bases, "PCT_SELECTED_BASES");
    close(ours.pct_off_bait, "PCT_OFF_BAIT");
    close(ours.on_bait_vs_selected, "ON_BAIT_VS_SELECTED");
    close(ours.mean_bait_coverage, "MEAN_BAIT_COVERAGE");
    close(ours.mean_target_coverage, "MEAN_TARGET_COVERAGE");
    close(ours.pct_usable_bases_on_bait, "PCT_USABLE_BASES_ON_BAIT");
    // The selected fraction counts the near bases, so it is not the on-bait fraction: 0.75 against
    // a third of that again.
    assert_eq!(number(&row, "PCT_SELECTED_BASES"), 0.75);
    assert_eq!(number(&row, "ON_BAIT_VS_SELECTED"), 0.333333);
}

/// The two quality floors empty the coverage and leave the bait counts alone.
#[test]
fn the_quality_floors_empty_the_coverage_alone() {
    let text = corpus();
    let plain = metrics(&text, "plain");
    for case in ["mapping-quality-floor", "base-quality-floor"] {
        let row = metrics(&text, case);
        assert_eq!(number(&row, "ON_TARGET_BASES"), 0.0, "{case}");
        assert_eq!(number(&row, "MEAN_TARGET_COVERAGE"), 0.0, "{case}");
        assert_eq!(
            number(&row, "ON_BAIT_BASES"),
            number(&plain, "ON_BAIT_BASES"),
            "{case}"
        );
    }
    // And the mapping-quality default is one, not nought, so a pair at quality nought is already
    // out of the coverage with no argument given.
    let zero = metrics(&text, "mapping-quality-zero");
    assert_eq!(DEFAULT_MINIMUM_MAPPING_QUALITY, 1);
    assert_eq!(number(&zero, "ON_BAIT_BASES"), 60.0);
    assert_eq!(number(&zero, "ON_TARGET_BASES"), 0.0);
    // The base-quality default is nought, which is what makes this tool count what a WGS one would
    // not.
    assert_eq!(DEFAULT_MINIMUM_BASE_QUALITY, 0);
}

/// The clipping argument changes nothing, because the coverage is counted per locus.
#[test]
fn clipping_the_overlap_changes_nothing() {
    let text = corpus();
    assert_eq!(
        field(&text, "metrics", "overlapping-pair"),
        field(&text, "metrics", "overlapping-pair-clipped")
    );
    assert_eq!(
        field(&text, "per-target", "overlapping-pair"),
        field(&text, "per-target", "overlapping-pair-clipped")
    );
    // The pair's two ends span thirty-five bases of target between them, and that is what is
    // counted: not the sixty their two lengths would give.
    let row = metrics(&text, "overlapping-pair");
    assert_eq!(number(&row, "ON_TARGET_BASES"), 35.0);
    assert_eq!(number(&row, "PF_BASES_ALIGNED"), 60.0);
}

/// The per-target file's columns, and the port's own row for the same coverage.
#[test]
fn the_per_target_file_is_a_row_per_target() {
    let text = corpus();
    let rows = table(&text, "per-target", "per-target");
    let uncovered = rows
        .iter()
        .find(|row| row["name"] == "target-b")
        .expect("target-b");
    assert_eq!(uncovered["mean_coverage"], "0");
    assert_eq!(uncovered["pct_0x"], "1");
    assert_eq!(uncovered["read_count"], "0");
    let covered = rows
        .iter()
        .find(|row| row["name"] == "target-a")
        .expect("target-a");
    assert_eq!(covered["read_count"], "2");
    // The port's row over the same sixty bases: ten covered once, fifty not.
    let mut coverage = vec![0i64; 60];
    for depth in coverage.iter_mut().take(50).skip(0) {
        *depth = 1;
    }
    let bases: Vec<u8> = (0..60).map(|i| b"ACGTGGCCATAT"[(i + 120) % 12]).collect();
    let run_mean = 0.555_556;
    let row = target_row("target-a", &coverage, &bases, 2, run_mean);
    assert_eq!(row.length, 60);
    assert!((row.mean_coverage - 50.0 / 60.0).abs() < 1e-9);
    assert!((row.pct_zero_coverage - 10.0 / 60.0).abs() < 1e-9);
    assert_eq!((row.minimum, row.maximum), (0, 1));
    assert!((row.normalized_coverage - row.mean_coverage / run_mean).abs() < 1e-9);
    // The GC of the golden's own target, which the fixture's repeating reference fixes at a half.
    assert_eq!(covered["%gc"], "0.5");
    assert!((row.gc - 0.5).abs() < 1e-9);
}

/// The per-base file is a row per target base.
#[test]
fn the_per_base_file_is_a_row_per_target_base() {
    let text = corpus();
    let body = field(&text, "per-base", "per-base").expect("the per-base file");
    assert!(body.contains("# 180 rows in all"));
    assert!(body.contains("chr1\t121\ttarget-a\t0"));
    // Three sixty-base targets, and the file is absent unless it is asked for.
    assert_eq!(field(&text, "per-base", "plain").as_deref(), Some("absent"));
    assert_eq!(
        field(&text, "per-target", "plain").as_deref(),
        Some("absent")
    );
}
