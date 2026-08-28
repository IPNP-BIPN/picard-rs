//! Conformance for `CollectUmiPrevalenceMetrics` against Picard 3.4.0.
//!
//! Each case carries the file the tool read, as SAM without its header, and the histogram it
//! wrote. Grouping into duplicate sets is not ported, so the port is handed the reads grouped by
//! alignment start, which is how this fixture's sets are formed.
//!
//! # What this suite is for
//!
//!  * **the barcode-quality filter being inverted**;
//!  * **a well-formed file therefore reporting nothing, and a lower floor reporting less**;
//!  * **an absent quality tag being the only other way past it**;
//!  * **the histogram being sets by UMI count**;
//!  * **the UMIs of a set being a set**;
//!  * **the four other filters each taking their read**;
//!  * **`--FILTER_UNPAIRED_READS` being on by default**;
//!  * **and a file whose every read is filtered writing an empty histogram.**

use std::io::Read as _;

use picard_analysis::collect_umi_prevalence_metrics::{
    barcode_quality_filters_out, barcode_tag_filters_out, decode_barcode_qualities, filters_out,
    histogram, Arguments, Read, DEFAULT_BARCODE_QUALITY_TAG, DEFAULT_BARCODE_TAG,
    DEFAULT_MINIMUM_BARCODE_BASE_QUALITY, DEFAULT_MINIMUM_MAPPING_QUALITY,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_umi_prevalence_metrics.txt.gz");
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

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// One case's histogram, as UMI count to set count.
fn written(text: &str, case: &str) -> std::collections::BTreeMap<usize, i64> {
    let payload = field(text, "histogram", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = payload.lines().filter(|line| !line.is_empty());
    match lines.next() {
        None => std::collections::BTreeMap::new(),
        Some(header) => {
            assert_eq!(header, "numUmis\tduplicateSets", "{case}");
            lines
                .map(|line| {
                    let (k, v) = line.split_once('\t').expect("a pair");
                    (k.parse().expect("a count"), v.parse().expect("a number"))
                })
                .collect()
        }
    }
}

/// The reads of one case, grouped into duplicate sets by alignment start, which is how the
/// fixture's sets are formed: every set sits at its own position.
fn sets(text: &str, case: &str) -> Vec<Vec<Read>> {
    let mut by_start: std::collections::BTreeMap<i64, Vec<Read>> =
        std::collections::BTreeMap::new();
    for line in field(text, "sam", case)
        .unwrap_or_else(|| panic!("{case} has an input"))
        .lines()
        .filter(|line| !line.is_empty())
    {
        let columns: Vec<&str> = line.split('\t').collect();
        let flags: u32 = columns[1].parse().expect("a flag word");
        let tag = |name: &str| {
            columns
                .iter()
                .find_map(|column| column.strip_prefix(&format!("{name}:Z:")))
                .map(str::to_string)
        };
        let start: i64 = columns[3].parse().expect("a start");
        by_start.entry(start).or_default().push(Read {
            unmapped: flags & 0x4 != 0,
            mapping_quality: columns[4].parse().expect("a mapping quality"),
            secondary_or_supplementary: flags & 0x100 != 0 || flags & 0x800 != 0,
            paired: flags & 0x1 != 0,
            barcode: tag(DEFAULT_BARCODE_TAG),
            barcode_qualities: tag(DEFAULT_BARCODE_QUALITY_TAG)
                .as_deref()
                .map(decode_barcode_qualities),
        });
    }
    by_start.into_values().collect()
}

fn arguments(case: &str) -> Arguments {
    let mut arguments = Arguments::default();
    if case == "unpaired-kept" {
        arguments.filter_unpaired_reads = false;
    }
    if case == "low-barcode-floor" {
        arguments.minimum_barcode_base_quality = 0;
    }
    arguments
}

const CASES: &[&str] = &[
    "three-sets-one-umi",
    "one-set-three-umis",
    "one-set-two-umis",
    "no-umi-tag",
    "barcode-quality-good",
    "barcode-quality-all-good",
    "no-barcode-quality",
    "unaligned",
    "low-mapping-quality",
    "secondary",
    "unpaired-filtered",
    "unpaired-kept",
    "other-tags",
    "low-barcode-floor",
    "everything-filtered",
    "empty",
];

/// Every case's histogram is what the port reaches.
#[test]
fn every_case_writes_the_same_histogram() {
    let text = corpus();
    for case in CASES {
        assert_eq!(
            histogram(&sets(&text, case), &arguments(case)),
            written(&text, case),
            "{case}"
        );
    }
}

/// The barcode-quality filter is inverted: a barcode entirely above the floor is dropped and one
/// with a base under it is kept.
#[test]
fn the_barcode_quality_filter_is_inverted() {
    let good = Read {
        unmapped: false,
        mapping_quality: 60,
        secondary_or_supplementary: false,
        paired: true,
        barcode: Some("AAAAAA".to_string()),
        barcode_qualities: Some(decode_barcode_qualities("IIIIII")),
    };
    let bad = Read {
        barcode_qualities: Some(decode_barcode_qualities("II#III")),
        ..good.clone()
    };
    assert!(barcode_quality_filters_out(
        &good,
        DEFAULT_MINIMUM_BARCODE_BASE_QUALITY
    ));
    assert!(!barcode_quality_filters_out(
        &bad,
        DEFAULT_MINIMUM_BARCODE_BASE_QUALITY
    ));
    // Which is what the two goldens show: the well-formed read is the one that vanishes.
    let text = corpus();
    assert_eq!(written(&text, "barcode-quality-good")[&1], 1);
    assert!(written(&text, "barcode-quality-all-good").is_empty());
}

/// Lowering the floor makes the tool report less rather than more.
#[test]
fn a_lower_floor_reports_less() {
    let text = corpus();
    // At the default floor the fixture's one read survives.
    assert_eq!(written(&text, "other-tags")[&1], 1);
    // At a floor of nought no base is under it, so the inversion drops the read.
    assert!(written(&text, "low-barcode-floor").is_empty());
    let read = Read {
        unmapped: false,
        mapping_quality: 60,
        secondary_or_supplementary: false,
        paired: true,
        barcode: Some("AAAAAA".to_string()),
        barcode_qualities: Some(decode_barcode_qualities("II#III")),
    };
    assert!(!barcode_quality_filters_out(&read, 30));
    assert!(barcode_quality_filters_out(&read, 0));
}

/// An absent quality tag is the only other way past that filter.
#[test]
fn an_absent_quality_tag_is_the_other_way_past() {
    let text = corpus();
    // The read with no BQ tag survives, so the set holds two UMIs where the good-quality one
    // leaves only one.
    assert_eq!(written(&text, "no-barcode-quality")[&2], 1);
    assert_eq!(written(&text, "barcode-quality-good")[&1], 1);
    let none = Read {
        unmapped: false,
        mapping_quality: 60,
        secondary_or_supplementary: false,
        paired: true,
        barcode: Some("AAAAAA".to_string()),
        barcode_qualities: None,
    };
    assert!(!barcode_quality_filters_out(
        &none,
        DEFAULT_MINIMUM_BARCODE_BASE_QUALITY
    ));
}

/// The histogram is sets by UMI count, and the UMIs of a set are a set.
#[test]
fn the_histogram_is_sets_by_umi_count() {
    let text = corpus();
    // Three sets of one UMI each.
    assert_eq!(written(&text, "three-sets-one-umi")[&1], 3);
    assert_eq!(written(&text, "three-sets-one-umi").len(), 1);
    // One set of three.
    assert_eq!(written(&text, "one-set-three-umis")[&3], 1);
    // And one set whose three reads carry two distinct tags.
    assert_eq!(written(&text, "one-set-two-umis")[&2], 1);
    assert_eq!(sets(&text, "one-set-two-umis")[0].len(), 3);
}

/// The four other filters each take their read.
#[test]
fn the_other_filters_each_take_their_read() {
    let text = corpus();
    for case in [
        "no-umi-tag",
        "unaligned",
        "low-mapping-quality",
        "secondary",
    ] {
        // Two reads in, one UMI out: the second never reached a set. The unaligned one sits at
        // its own start, so the count is over every set rather than over the first.
        assert_eq!(written(&text, case)[&1], 1, "{case}");
        let reads: usize = sets(&text, case).iter().map(Vec::len).sum();
        assert_eq!(reads, 2, "{case}");
        assert_eq!(written(&text, case).values().sum::<i64>(), 1, "{case}");
    }
    let arguments = Arguments::default();
    let base = Read {
        unmapped: false,
        mapping_quality: 60,
        secondary_or_supplementary: false,
        paired: true,
        barcode: Some("AAAAAA".to_string()),
        barcode_qualities: Some(decode_barcode_qualities("II#III")),
    };
    assert!(!filters_out(&base, &arguments));
    assert!(filters_out(
        &Read {
            unmapped: true,
            ..base.clone()
        },
        &arguments
    ));
    assert!(filters_out(
        &Read {
            mapping_quality: 5,
            ..base.clone()
        },
        &arguments
    ));
    assert!(filters_out(
        &Read {
            secondary_or_supplementary: true,
            ..base.clone()
        },
        &arguments
    ));
    assert!(barcode_tag_filters_out(&Read {
        barcode: None,
        ..base.clone()
    }));
    assert_eq!(DEFAULT_MINIMUM_MAPPING_QUALITY, 30);
}

/// The unpaired filter is on by default and can be turned off.
#[test]
fn the_unpaired_filter_is_on_by_default() {
    let text = corpus();
    assert_eq!(written(&text, "unpaired-filtered")[&1], 1);
    assert_eq!(written(&text, "unpaired-kept")[&2], 1);
    assert!(Arguments::default().filter_unpaired_reads);
    let unpaired = Read {
        unmapped: false,
        mapping_quality: 60,
        secondary_or_supplementary: false,
        paired: false,
        barcode: Some("AAAAAA".to_string()),
        barcode_qualities: Some(decode_barcode_qualities("II#III")),
    };
    assert!(filters_out(&unpaired, &Arguments::default()));
    assert!(!filters_out(
        &unpaired,
        &Arguments {
            filter_unpaired_reads: false,
            ..Arguments::default()
        }
    ));
}

/// A file whose every read is filtered writes an empty histogram, as does a file with no reads.
#[test]
fn an_empty_result_is_still_a_histogram() {
    let text = corpus();
    assert!(written(&text, "everything-filtered").is_empty());
    assert!(written(&text, "empty").is_empty());
    assert!(sets(&text, "empty").is_empty());
    // The first has a read, and it is filtered rather than absent.
    assert_eq!(sets(&text, "everything-filtered")[0].len(), 1);
    assert!(filters_out(
        &sets(&text, "everything-filtered")[0][0],
        &Arguments::default()
    ));
}

/// The quality tag is FASTQ, decoded by subtracting thirty-three, with spaces dropped first.
#[test]
fn the_quality_tag_is_fastq() {
    assert_eq!(decode_barcode_qualities("I"), vec![40]);
    assert_eq!(decode_barcode_qualities("#"), vec![2]);
    assert_eq!(decode_barcode_qualities("!"), vec![0]);
    assert_eq!(decode_barcode_qualities("I I"), vec![40, 40]);
    assert_eq!(decode_barcode_qualities(""), Vec::<i32>::new());
}
