//! Conformance for `UmiAwareMarkDuplicatesWithMateCigar` against Picard 3.4.0.
//!
//! Golden from `tools/matecigardup-conformance/UmiDuplicatesDump.java`: eleven runs with the
//! marked output, the duplication metrics and the UMI metrics table.
//!
//! # What this suite is for
//!
//!  * **the UMI splitting a set, and the edit distance deciding when it does not**;
//!  * **`MAX_EDIT_DISTANCE_TO_JOIN` moved either way**;
//!  * **an `N` counting as a difference like any other base**;
//!  * **the molecular identifier written back, which carries the position as well as the UMI**;
//!  * **a missing UMI refused by name unless `ALLOW_MISSING_UMIS` says otherwise**;
//!  * **and the UMI metrics: the counts, the two entropies, and the base quality whose `-1` is
//!    `Math.round(Infinity)` cast to an `int`.**

use std::io::Read;

use htsjdk_bam::text_parse::parse_cigar;
use picard_analysis::mark_duplicates::Record;
use picard_analysis::mate_cigar_duplicates::SortOrder;
use picard_analysis::umi_duplicates::{
    cluster, mark_with_umis, phred_from_error_probability, within_hamming_distance, UmiOptions,
    UmiRefusal,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/umi_duplicates.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
                .replace("\\\\", "\\")
        })
}

fn record(line: &str) -> Record {
    let columns: Vec<&str> = line.split('\t').collect();
    let tag = |name: &str| -> Option<String> {
        columns
            .iter()
            .skip(11)
            .find(|column| column.starts_with(&format!("{name}:")))
            .map(|column| column.rsplit(':').next().expect("a tag value").to_string())
    };
    Record {
        name: columns[0].to_string(),
        flags: columns[1].parse().expect("the flags"),
        reference_index: 0,
        alignment_start: columns[3].parse().expect("the position"),
        cigar: parse_cigar(columns[5]).expect("the cigar"),
        qualities: columns[10].bytes().map(|byte| byte - 33).collect(),
        mate_reference_index: 0,
        library: "lib1".to_string(),
        read_group: 0,
        barcode: tag("RX"),
        existing_dt: None,
        mate_cigar: tag("MC").map(|text| parse_cigar(&text).expect("the mate cigar")),
        mate_alignment_start: columns[7].parse().unwrap_or(0),
    }
}

fn records(text: &str, case: &str) -> Vec<Record> {
    field(text, "sam", case)
        .unwrap_or_else(|| panic!("{case}"))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(record)
        .collect()
}

/// The output as the golden wrote it: the name, the flags and the two tags.
fn marked(text: &str, case: &str) -> Option<Vec<(String, u16, Option<String>)>> {
    field(text, "marked", case).map(|body| {
        body.lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                let identifier = columns
                    .iter()
                    .skip(11)
                    .find(|column| column.starts_with("MI:"))
                    .map(|column| column["MI:Z:".len()..].to_string());
                (
                    columns[0].to_string(),
                    columns[1].parse().expect("the flags"),
                    identifier,
                )
            })
            .collect()
    })
}

/// One row of a metrics table, by column name.
fn row(text: &str, kind: &str, case: &str) -> std::collections::HashMap<String, String> {
    let table = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = table.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let values: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    header
        .iter()
        .zip(values)
        .map(|(name, value)| ((*name).to_string(), value.to_string()))
        .collect()
}

fn options(case: &str) -> UmiOptions {
    let mut options = UmiOptions::default();
    match case {
        "umis-one-base-apart-with-no-joining" => options.max_edit_distance_to_join = 0,
        "two-umis-joined-at-four" => options.max_edit_distance_to_join = 4,
        "one-umi-with-an-assigned-tag" => options.molecular_identifier_tag = Some("MI".to_string()),
        "no-umis-allowed" => options.allow_missing_umis = true,
        _ => {}
    }
    options
}

const CASES: [&str; 9] = [
    "two-umis",
    "umis-one-base-apart",
    "umis-one-base-apart-with-no-joining",
    "two-umis-joined-at-four",
    "one-umi",
    "one-umi-with-an-assigned-tag",
    "a-umi-with-an-n",
    "no-umis-allowed",
    "a-different-umi-tag",
];

/// Every case's output, record for record.
#[test]
fn every_case_marks_what_the_reference_marked() {
    let text = corpus();
    for case in CASES {
        let input = records(&text, case);
        let run = mark_with_umis(&input, SortOrder::Coordinate, &options(case))
            .unwrap_or_else(|refusal| panic!("{case}: {}", refusal.message()));
        let produced: Vec<(String, u16, Option<String>)> = input
            .iter()
            .enumerate()
            .filter(|(index, _)| run.marking.written[*index])
            .map(|(index, record)| {
                let mut flags = record.flags & !0x400;
                if run.marking.duplicate[index] {
                    flags |= 0x400;
                }
                (
                    record.name.clone(),
                    flags,
                    run.molecular_identifiers[index].clone(),
                )
            })
            .collect();
        assert_eq!(produced, marked(&text, case).expect("the golden"), "{case}");
    }
}

/// The UMI metrics beside the duplication ones.
#[test]
fn the_umi_metrics_are_the_reference_ones() {
    let text = corpus();
    for case in CASES {
        let input = records(&text, case);
        let run = mark_with_umis(&input, SortOrder::Coordinate, &options(case))
            .unwrap_or_else(|refusal| panic!("{case}: {}", refusal.message()));
        let expected = row(&text, "umi", case);
        let metrics = &run.umi_metrics[0];
        assert_eq!(metrics.library, expected["LIBRARY"], "{case}");
        for (column, value) in [
            ("OBSERVED_UNIQUE_UMIS", metrics.observed_unique_umis),
            ("INFERRED_UNIQUE_UMIS", metrics.inferred_unique_umis),
            ("OBSERVED_BASE_ERRORS", metrics.observed_base_errors),
            (
                "DUPLICATE_SETS_IGNORING_UMI",
                metrics.duplicate_sets_ignoring_umi,
            ),
            ("DUPLICATE_SETS_WITH_UMI", metrics.duplicate_sets_with_umi),
        ] {
            assert_eq!(
                expected[column].parse::<i64>().expect("a count"),
                value,
                "{case}/{column}"
            );
        }
        for (column, value) in [
            ("MEAN_UMI_LENGTH", metrics.mean_umi_length),
            ("OBSERVED_UMI_ENTROPY", metrics.observed_umi_entropy),
            ("INFERRED_UMI_ENTROPY", metrics.inferred_umi_entropy),
            ("PCT_UMI_WITH_N", metrics.percent_umi_with_n),
        ] {
            // `?` is what the metrics file writes for a NaN, which is what a run with no UMIs at
            // all divides its way to.
            if expected[column] == "?" {
                assert!(value.is_nan(), "{case}/{column}: {value}");
                continue;
            }
            let recorded: f64 = expected[column].parse().expect("a number");
            assert!(
                (recorded - value).abs() < 1e-6,
                "{case}/{column}: {recorded} vs {value}"
            );
        }
        assert_eq!(
            expected["UMI_BASE_QUALITIES"]
                .parse::<i32>()
                .expect("a quality"),
            metrics.umi_base_qualities,
            "{case}"
        );
    }
}

/// A missing UMI is refused by name, and the wording is the tool's own.
#[test]
fn a_missing_umi_is_refused() {
    let text = corpus();
    let input = records(&text, "no-umis");
    let refusal = mark_with_umis(&input, SortOrder::Coordinate, &UmiOptions::default())
        .expect_err("the refusal");
    let recorded = field(&text, "error", "no-umis").expect("the golden's refusal");
    assert_eq!(
        recorded,
        format!("{}:{}", refusal.exception(), refusal.message())
    );
    assert!(matches!(refusal, UmiRefusal::MissingUmi { .. }));
    // With the argument that allows it, the same file runs and nothing carries a UMI.
    let options = UmiOptions {
        allow_missing_umis: true,
        ..UmiOptions::default()
    };
    assert!(mark_with_umis(&input, SortOrder::Coordinate, &options).is_ok());
}

/// The two pieces the split is built from, on their own.
#[test]
fn the_join_is_transitive_and_the_overflow_is_the_references() {
    // `ATCC`, `AACC` and `AACG` are one apart in a chain, so all three join at a distance of one
    // even though the first and the last are two apart.
    let umis = vec!["ATCC".to_string(), "AACC".to_string(), "AACG".to_string()];
    let clusters = cluster(&umis, 1);
    assert_eq!(clusters[0], clusters[1]);
    assert_eq!(clusters[1], clusters[2]);
    // At a distance of zero nothing joins.
    let clusters = cluster(&umis, 0);
    assert_ne!(clusters[0], clusters[1]);
    // An N is a difference like any other.
    assert!(within_hamming_distance("AANA", "AAAA", 1));
    assert!(!within_hamming_distance("ANNA", "AAAA", 1));
    // A UMI of another length is never within any distance.
    assert!(!within_hamming_distance("AAA", "AAAA", 4));
    // And the quality of a run with no base errors is `-1`: `Math.round(Infinity)` is
    // `Long.MAX_VALUE`, whose low thirty-two bits are all ones.
    assert_eq!(phred_from_error_probability(0.0), -1);
    assert_eq!(phred_from_error_probability(0.1), 10);
}
