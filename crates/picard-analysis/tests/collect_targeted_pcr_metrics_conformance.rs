//! Conformance for `CollectTargetedPcrMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/pcrmetrics-conformance`, whose fixture is `CollectHsMetrics`' own, so the
//! two goldens can be read against each other.
//!
//! # What this suite is for
//!
//!  * **the columns being the amplicon's and not the bait's**;
//!  * **the amplicon arithmetic being the other tool's, number for number**;
//!  * **the target arithmetic NOT being, because of one line in a constructor**;
//!  * **`--NEAR_DISTANCE` moving the window here too**;
//!  * **`--CUSTOM_AMPLICON_SET_NAME` naming the set**;
//!  * **the two quality floors emptying the coverage**;
//!  * **and the same two coverage files being written.**

use std::io::Read;

use picard_analysis::collect_hs_metrics::{BaitPlacement, Counts};
use picard_analysis::collect_targeted_pcr_metrics::{
    bait_column, derived, placement, AMPLICON_COLUMNS, HS_METRICS_CLIP_OVERLAPPING_READS_DEFAULT,
    PCR_METRICS_CLIP_OVERLAPPING_READS_DEFAULT,
    SHARED_CLIP_OVERLAPPING_READS_DEFAULT as SHARED_DEFAULT,
};

fn read(name: &str) -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn corpus() -> String {
    read("collect_targeted_pcr_metrics.txt.gz")
}

fn hs_corpus() -> String {
    read("collect_hs_metrics.txt.gz")
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

fn metrics(text: &str, case: &str) -> std::collections::HashMap<String, String> {
    let body = field(text, "metrics", case).unwrap_or_else(|| panic!("metrics/{case}"));
    let mut lines = body.lines().filter(|line| !line.is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    header
        .iter()
        .zip(lines.next().expect("a row").split('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

fn number(row: &std::collections::HashMap<String, String>, name: &str) -> f64 {
    let value = row.get(name).unwrap_or_else(|| panic!("{name}"));
    if value == "?" || value.is_empty() {
        f64::NAN
    } else {
        value.parse().unwrap_or_else(|_| panic!("{name}={value}"))
    }
}

/// The columns are the amplicon's, and there is no bait column at all.
#[test]
fn the_columns_are_the_amplicons() {
    let text = corpus();
    let row = metrics(&text, "plain");
    for column in AMPLICON_COLUMNS {
        assert!(row.contains_key(column), "{column}");
        let bait = bait_column(column).unwrap_or_else(|| panic!("{column}"));
        assert!(!row.contains_key(bait), "{bait} is not written here");
    }
    assert_eq!(bait_column("NO_SUCH_COLUMN"), None);
}

/// The amplicon arithmetic is the other tool's, number for number.
#[test]
fn the_amplicon_numbers_are_the_bait_numbers() {
    let text = corpus();
    let hs = hs_corpus();
    let ours = metrics(&text, "plain");
    let theirs = metrics(&hs, "plain");
    for column in AMPLICON_COLUMNS {
        let bait = bait_column(column).expect("a bait column");
        if column == "CUSTOM_AMPLICON_SET" {
            // The name is the file's, and the two fixtures name their files differently.
            assert_eq!(ours[column], "amplicons");
            assert_eq!(theirs[bait], "baits");
            continue;
        }
        assert_eq!(ours[column], theirs[bait], "{column} against {bait}");
    }
    // And the port's own placement and arithmetic are the other module's.
    let amplicons = [(101, 200), (301, 400)];
    assert_eq!(placement(150, &amplicons, 250), BaitPlacement::On);
    assert_eq!(placement(260, &amplicons, 250), BaitPlacement::Near);
    assert_eq!(placement(260, &amplicons, 0), BaitPlacement::Off);
    let counts = Counts {
        pf_bases: number(&ours, "PF_BASES") as i64,
        pf_bases_aligned: number(&ours, "PF_BASES_ALIGNED") as i64,
        on_bait: number(&ours, "ON_AMPLICON_BASES") as i64,
        near_bait: number(&ours, "NEAR_AMPLICON_BASES") as i64,
        off_bait: number(&ours, "OFF_AMPLICON_BASES") as i64,
        on_target: number(&ours, "ON_TARGET_BASES") as i64,
        bait_territory: number(&ours, "AMPLICON_TERRITORY") as i64,
        target_territory: number(&ours, "TARGET_TERRITORY") as i64,
    };
    let computed = derived(&counts);
    assert!((computed.pct_selected_bases - number(&ours, "PCT_AMPLIFIED_BASES")).abs() < 1e-6);
    assert!((computed.pct_off_bait - number(&ours, "PCT_OFF_AMPLICON")).abs() < 1e-6);
    assert!((computed.on_bait_vs_selected - number(&ours, "ON_AMPLICON_VS_SELECTED")).abs() < 1e-6);
    assert!((computed.mean_target_coverage - number(&ours, "MEAN_TARGET_COVERAGE")).abs() < 1e-6);
}

/// The target arithmetic is NOT the other tool's, because of one line in a constructor.
#[test]
fn the_target_numbers_differ_by_the_clipping_default() {
    let text = corpus();
    let hs = hs_corpus();
    let ours = metrics(&text, "plain");
    let theirs = metrics(&hs, "plain");
    // The overlap of a pair is counted twice here and once there.
    assert_eq!(number(&ours, "ON_TARGET_BASES"), 120.0);
    assert_eq!(number(&theirs, "ON_TARGET_BASES"), 100.0);
    assert_eq!(number(&ours, "MEAN_TARGET_COVERAGE"), 0.666667);
    assert_eq!(number(&theirs, "MEAN_TARGET_COVERAGE"), 0.555556);
    // Which is the difference between the two defaults, and nothing else: turning the other tool's
    // clipping OFF gives this tool's numbers.
    let unclipped = metrics(&hs, "clip-overlapping-off");
    assert_eq!(
        number(&unclipped, "ON_TARGET_BASES"),
        number(&ours, "ON_TARGET_BASES")
    );
    assert_eq!(
        number(&unclipped, "MEAN_TARGET_COVERAGE"),
        number(&ours, "MEAN_TARGET_COVERAGE")
    );
    // The two defaults the port carries, which is the fact the two goldens just demonstrated.
    assert_ne!(
        HS_METRICS_CLIP_OVERLAPPING_READS_DEFAULT,
        PCR_METRICS_CLIP_OVERLAPPING_READS_DEFAULT
    );
    assert_eq!(PCR_METRICS_CLIP_OVERLAPPING_READS_DEFAULT, SHARED_DEFAULT);
}

/// The near window, the set name and the two quality floors behave as they do next door.
#[test]
fn the_arguments_behave_as_they_do_next_door() {
    let text = corpus();
    let plain = metrics(&text, "plain");
    let near_zero = metrics(&text, "near-distance-zero");
    assert_eq!(number(&plain, "NEAR_AMPLICON_BASES"), 120.0);
    assert_eq!(number(&near_zero, "NEAR_AMPLICON_BASES"), 0.0);
    assert_eq!(number(&near_zero, "OFF_AMPLICON_BASES"), 180.0);
    assert_eq!(
        metrics(&text, "amplicon-set-name")["CUSTOM_AMPLICON_SET"],
        "my-amplicon-set"
    );
    for case in ["mapping-quality-floor", "base-quality-floor"] {
        let row = metrics(&text, case);
        assert_eq!(number(&row, "ON_TARGET_BASES"), 0.0, "{case}");
        assert_eq!(
            number(&row, "ON_AMPLICON_BASES"),
            number(&plain, "ON_AMPLICON_BASES"),
            "{case}"
        );
    }
    // A pair at mapping quality nought is out of the coverage with no argument given.
    assert_eq!(
        number(&metrics(&text, "mapping-quality-zero"), "ON_TARGET_BASES"),
        0.0
    );
    // And a file with no reads still reports both territories.
    let empty = metrics(&text, "no-reads");
    assert_eq!(number(&empty, "AMPLICON_TERRITORY"), 200.0);
    assert_eq!(number(&empty, "TARGET_TERRITORY"), 180.0);
}

/// The same two coverage files, written only when they are asked for.
#[test]
fn the_two_coverage_files_are_the_same_two() {
    let text = corpus();
    assert_eq!(
        field(&text, "per-target", "plain").as_deref(),
        Some("absent")
    );
    assert_eq!(field(&text, "per-base", "plain").as_deref(), Some("absent"));
    let per_target = field(&text, "per-target", "per-target").expect("the per-target file");
    assert!(per_target.starts_with("chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage"));
    let per_base = field(&text, "per-base", "per-base").expect("the per-base file");
    assert!(per_base.contains("# 180 rows in all"));
    // The uncovered target is a row here as well, at nought.
    assert!(per_target.contains("target-b"));
}
