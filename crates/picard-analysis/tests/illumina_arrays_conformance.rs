//! Conformance for the genotyping-array formats and two tools against Picard 3.4.0.
//!
//! Goldens from `tools/arrays-conformance/`: `BpmToNormalizationManifestCsv` and
//! `CompareGtcFiles`, over manifests and call files written byte by byte.
//!
//! # What this suite is for
//!
//!  * **the varint string length the three formats share**;
//!  * **the normalization id being the file's id plus a hundred times the assay type**;
//!  * **the two cross-checks a manifest is refused by, in the parser's own order**;
//!  * **the CSV's rows, including the score's four decimal places**;
//!  * **and a different sample name not being a difference.**

use std::io::Read;

use picard_analysis::illumina_arrays::{
    compare, normalization_line, parse_string, validate, Comparison, Locus, ManifestRefusal,
    NormalizationRow, EXCLUDED_FROM_COMPARISON, NORMALIZATION_HEADER,
};

fn corpus(name: &str) -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.txt.gz"));
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
        })
}

/// The fixture's four loci, as the writer wrote them.
fn loci() -> Vec<Locus> {
    vec![
        Locus {
            name: "rs1".into(),
            index: 0,
            chrom: "1".into(),
            position: 1000,
            address_a: 11,
            address_b: 0,
            assay_type: 0,
            raw_normalization_id: 1,
        },
        Locus {
            name: "rs2".into(),
            index: 1,
            chrom: "1".into(),
            position: 2000,
            address_a: 12,
            address_b: 13,
            assay_type: 1,
            raw_normalization_id: 1,
        },
        Locus {
            name: "rs3".into(),
            index: 2,
            chrom: "2".into(),
            position: 3000,
            address_a: 14,
            address_b: 0,
            assay_type: 0,
            raw_normalization_id: 2,
        },
        Locus {
            name: "rs4".into(),
            index: 3,
            chrom: "2".into(),
            position: 4000,
            address_a: 15,
            address_b: 16,
            assay_type: 2,
            raw_normalization_id: 2,
        },
    ]
}

/// A varint length and then the bytes, which is how all three formats write a string.
#[test]
fn a_string_is_a_varint_and_its_bytes() {
    // A short string costs one byte of length.
    let mut bytes = vec![3u8];
    bytes.extend_from_slice(b"rs1");
    assert_eq!(parse_string(&bytes, 0), Some(("rs1".to_string(), 4)));
    // A string of two hundred bytes costs two, the first carrying its low seven bits.
    let long = "a".repeat(200);
    let mut bytes = vec![(200 & 0x7F) | 0x80, 200 >> 7];
    bytes.extend_from_slice(long.as_bytes());
    assert_eq!(parse_string(&bytes, 0), Some((long, 202)));
    // And an empty string is a single zero.
    assert_eq!(parse_string(&[0], 0), Some((String::new(), 1)));
}

/// The number the CSV reports is not the number in the file.
#[test]
fn the_normalization_id_folds_in_the_assay_type() {
    let text = corpus("bpm_to_csv");
    let recorded = field(&text, "csv", "four-loci").expect("the golden");
    let mut lines = recorded.lines();
    assert_eq!(lines.next().expect("a header"), NORMALIZATION_HEADER);

    let ids: Vec<i32> = loci().iter().map(Locus::normalization_id).collect();
    assert_eq!(ids, vec![1, 101, 2, 202]);
    for (line, id) in lines.zip(&ids) {
        assert_eq!(line.split(',').next_back().expect("the id"), id.to_string());
    }

    // One locus of each assay type, on its own, reports the same number.
    let one = field(&text, "csv", "one-locus-of-assay-type-one").expect("the golden");
    assert!(one.lines().nth(1).expect("a row").ends_with(",101"));
}

/// The CSV's rows, written the way the file writes them.
#[test]
fn the_rows_are_the_goldens() {
    let text = corpus("bpm_to_csv");
    let recorded = field(&text, "csv", "one-locus-of-assay-type-zero").expect("the golden");
    let row = NormalizationRow {
        index: 1,
        name: "rs1".to_string(),
        chromosome: "1".to_string(),
        position: 1000,
        // The score is the cluster file's, and the writer gives it four decimal places.
        gentrain_score: 0.5,
        snp: "[A/G]".to_string(),
        illumina_strand: "TOP".to_string(),
        customer_strand: "TOP".to_string(),
        normalization_id: 1,
    };
    assert_eq!(
        recorded.lines().nth(1).expect("a row"),
        normalization_line(&row)
    );
    assert!(normalization_line(&row).contains(",0.5000,"));
}

/// The two cross-checks, in the parser's own order and words.
#[test]
fn the_refusals_are_the_goldens() {
    let text = corpus("bpm_to_csv");
    for locus in loci() {
        assert_eq!(validate(&locus), Ok(()));
    }

    // An assay type of zero with a B address, and any other type without one.
    let with_b = Locus {
        address_b: 12,
        ..loci()[0].clone()
    };
    assert_eq!(
        validate(&with_b),
        Err(ManifestRefusal::AssayType {
            assay_type: 0,
            address_b: 12
        })
    );
    assert!(field(&text, "error", "assay-type-zero-with-a-b-address")
        .expect("the golden")
        .contains("Invalid assay_type '0' for address B '12'"));
    let without_b = Locus {
        assay_type: 1,
        address_b: 0,
        ..loci()[0].clone()
    };
    assert_eq!(
        validate(&without_b),
        Err(ManifestRefusal::AssayType {
            assay_type: 1,
            address_b: 0
        })
    );

    // A normalization id above a hundred is refused BEFORE the assay type is folded in, which is
    // what stops the addition from hiding it.
    let big = Locus {
        raw_normalization_id: 101,
        ..loci()[0].clone()
    };
    assert_eq!(
        validate(&big),
        Err(ManifestRefusal::NormalizationId {
            id: 101,
            name: "rs1".to_string()
        })
    );
    assert!(field(&text, "error", "a-normalization-id-above-a-hundred")
        .expect("the golden")
        .contains("Invalid normalization ID: 101 for name: rs1"));
}

/// A different sample name is not a difference; a different genotype is.
#[test]
fn the_comparison_excludes_what_two_runs_may_differ_in() {
    let text = corpus("compare_gtc");
    assert!(EXCLUDED_FROM_COMPARISON.contains(&"getSampleName"));

    let same = vec![
        (
            "getGenotypes",
            vec!["1".into(), "2".into()],
            vec!["1".into(), "2".into()],
        ),
        (
            "getSampleName",
            vec!["sample1".into()],
            vec!["sample2".into()],
        ),
    ];
    assert_eq!(compare(&same), Comparison::Same);
    assert_eq!(
        field(&text, "code", "two-identical-files").as_deref(),
        Some("0")
    );
    assert_eq!(
        field(&text, "code", "a-different-sample-name").as_deref(),
        Some("0")
    );

    let different = vec![(
        "getGenotypes",
        vec!["1".into(), "2".into()],
        vec!["1".into(), "3".into()],
    )];
    assert_eq!(compare(&different), Comparison::Different);
    assert_eq!(
        field(&text, "code", "a-different-genotype").as_deref(),
        Some("1")
    );

    // An array of a different length is a difference like any other.
    let shorter = vec![(
        "getGenotypes",
        vec!["1".into(), "2".into(), "3".into()],
        vec!["1".into(), "2".into()],
    )];
    assert_eq!(compare(&shorter), Comparison::Different);
    assert_eq!(field(&text, "code", "a-shorter-file").as_deref(), Some("1"));
}
