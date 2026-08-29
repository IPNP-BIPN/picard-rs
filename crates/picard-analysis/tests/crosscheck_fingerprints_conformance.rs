//! Conformance for `CrosscheckFingerprints` against Picard 3.4.0.
//!
//! Golden from `tools/checkfingerprint-conformance/CrosscheckFingerprintsDump.java`: thirteen runs
//! over the same three-site fixture `CheckFingerprint` uses.
//!
//! The LODs are the golden's, for the reason `check_fingerprint_conformance` gives: the pileup and
//! the haplotype likelihoods are the fingerprint's, and what this suite is about is what the tool
//! does with a number once it has one.
//!
//! # What this suite is for
//!
//!  * **all four verdicts, each earned by a fixture rather than asserted**;
//!  * **the threshold used with both signs, so a negative one makes everything inconclusive**;
//!  * **`--OUTPUT_ERRORS_ONLY` dropping agreement and keeping uncertainty**;
//!  * **the three exit codes, and an inconclusive row not counting as a mismatch**;
//!  * **and the matrix carrying the LODs and no verdict at all.**

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::crosscheck_fingerprints::{
    exit_code, matrix, rows, verdict, DataType, Fingerprint, Mode, Options, Verdict,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/crosscheck_fingerprints.txt.gz");
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

/// One case's rows: the two group values, the verdict and the LOD.
fn recorded(text: &str, case: &str) -> Vec<(String, String, String, f64)> {
    match field(text, "metrics", case) {
        None => Vec::new(),
        Some(table) => table
            .lines()
            .skip(1)
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                (
                    columns[0].to_string(),
                    columns[1].to_string(),
                    columns[2].to_string(),
                    columns[4].parse().expect("a LOD"),
                )
            })
            .collect(),
    }
}

/// The fingerprints a case compared, and the LOD it recorded for each pair.
fn fingerprints(text: &str, case: &str) -> (Vec<Fingerprint>, HashMap<(String, String), f64>) {
    let mut names: Vec<String> = Vec::new();
    let mut lods = HashMap::new();
    for (left, right, _, lod) in recorded(text, case) {
        for name in [&left, &right] {
            if !names.contains(name) {
                names.push(name.clone());
            }
        }
        lods.insert((left, right), lod);
    }
    let fingerprints = names
        .iter()
        .map(|name| Fingerprint {
            group: name.clone(),
            // The fixture's groups are named after their file, and the sample is the file's:
            // `f0unit1` is sample1 in both files, and `two-samples` is the case where it is not.
            sample: "sample1".to_string(),
            library: name.replace("unit", "lib"),
            file: name[..2].to_string(),
        })
        .collect();
    (fingerprints, lods)
}

/// Every case's verdicts, from the golden's own LODs.
#[test]
fn the_verdicts_are_the_goldens() {
    let text = corpus();
    for case in [
        "one-sample-agreeing",
        "one-sample-disagreeing",
        "expect-all-to-match",
        "a-lod-threshold",
    ] {
        let (prints, lods) = fingerprints(&text, case);
        let mut options = Options::default();
        if case == "a-lod-threshold" {
            options.lod_threshold = -100.0;
        }
        let produced = rows(
            &prints,
            &|left, right| {
                *lods
                    .get(&(left.group.clone(), right.group.clone()))
                    .unwrap_or(&0.0)
            },
            &options,
        );
        let expected = recorded(&text, case);
        assert_eq!(produced.len(), expected.len(), "{case}");
        for (row, (left, right, name, lod)) in produced.iter().zip(expected) {
            assert_eq!(row.left, left, "{case}");
            assert_eq!(row.right, right, "{case}");
            assert_eq!(row.verdict.name(), name, "{case}/{left}/{right}");
            assert_eq!(row.lod, lod, "{case}");
        }
    }
}

/// Two samples: what is expected comes from the names and not from the numbers.
#[test]
fn a_match_between_two_samples_is_unexpected() {
    let text = corpus();
    let expected = recorded(&text, "two-samples-agreeing");
    // The fixture's two files carry two samples, so a pair across them is not expected to match
    // however good its LOD is.
    let prints: Vec<Fingerprint> = ["f0unit1", "f1unit1"]
        .iter()
        .enumerate()
        .map(|(index, name)| Fingerprint {
            group: (*name).to_string(),
            sample: format!("sample{}", index + 1),
            library: format!("f{index}lib1"),
            file: format!("f{index}"),
        })
        .collect();
    let produced = rows(
        &prints,
        &|_left, _right| 1.828358,
        &Options {
            mode: Mode::CheckAllOthers,
            ..Options::default()
        },
    );
    assert_eq!(produced.len(), expected.len());
    for (row, (_, _, name, _)) in produced.iter().zip(&expected) {
        assert_eq!(row.verdict.name(), name);
    }
    // Which is two of the four: a fingerprint against itself is expected to match.
    assert_eq!(
        produced
            .iter()
            .filter(|row| row.verdict == Verdict::UnexpectedMatch)
            .count(),
        2
    );
    // And the default mode compares no pair across two samples at all.
    let same_sample_only = rows(&prints, &|_left, _right| 1.828358, &Options::default());
    assert_eq!(same_sample_only.len(), 2);
}

/// The exit codes, and what counts as a mismatch for them.
#[test]
fn the_exit_codes_are_the_goldens() {
    let text = corpus();
    for (case, code) in [
        ("one-sample-agreeing", 0),
        ("one-sample-disagreeing", 1),
        ("a-lod-threshold", 0),
        ("expect-all-to-match", 1),
    ] {
        let (prints, lods) = fingerprints(&text, case);
        let mut options = Options::default();
        if case == "a-lod-threshold" {
            options.lod_threshold = -100.0;
        }
        let produced = rows(
            &prints,
            &|left, right| {
                *lods
                    .get(&(left.group.clone(), right.group.clone()))
                    .unwrap_or(&0.0)
            },
            &options,
        );
        assert_eq!(exit_code(&produced, &options), code, "{case}");
        assert_eq!(
            field(&text, "code", case).as_deref(),
            Some(code.to_string().as_str()),
            "{case}"
        );
    }
    // An INCONCLUSIVE row is not an unexpected one, which is why a run whose every verdict is
    // inconclusive still exits zero.
    let (prints, lods) = fingerprints(&text, "a-lod-threshold");
    let options = Options {
        lod_threshold: -100.0,
        ..Options::default()
    };
    let produced = rows(
        &prints,
        &|left, right| {
            *lods
                .get(&(left.group.clone(), right.group.clone()))
                .unwrap_or(&0.0)
        },
        &options,
    );
    assert!(produced
        .iter()
        .all(|row| row.verdict == Verdict::Inconclusive));
    assert_eq!(exit_code(&produced, &options), 0);
    // And the code a run with nothing to compare returns is the argument's, not a constant.
    let nothing: Vec<picard_analysis::crosscheck_fingerprints::Row> = Vec::new();
    assert_eq!(exit_code(&nothing, &options), 1);
    assert_eq!(
        exit_code(
            &nothing,
            &Options {
                exit_code_when_no_valid_checks: 9,
                ..Options::default()
            }
        ),
        9
    );
}

/// `--OUTPUT_ERRORS_ONLY` drops agreement and keeps uncertainty.
#[test]
fn errors_only_keeps_what_is_not_agreement() {
    let text = corpus();
    // The golden's errors-only run wrote a file with a header and no rows, because every pair
    // agreed: the case's code is the one for a run that compared nothing new.
    assert!(recorded(&text, "errors-only").is_empty());
    let (prints, lods) = fingerprints(&text, "one-sample-agreeing");
    let produced = rows(
        &prints,
        &|left, right| {
            *lods
                .get(&(left.group.clone(), right.group.clone()))
                .unwrap_or(&0.0)
        },
        &Options {
            output_errors_only: true,
            ..Options::default()
        },
    );
    assert!(produced.is_empty());
    // A disagreeing run keeps its two unexpected rows and drops its two expected ones.
    let (prints, lods) = fingerprints(&text, "one-sample-disagreeing");
    let produced = rows(
        &prints,
        &|left, right| {
            *lods
                .get(&(left.group.clone(), right.group.clone()))
                .unwrap_or(&0.0)
        },
        &Options {
            output_errors_only: true,
            ..Options::default()
        },
    );
    assert_eq!(produced.len(), 2);
    assert!(produced
        .iter()
        .all(|row| row.verdict == Verdict::UnexpectedMismatch));
}

/// The matrix is the LODs and no verdict at all.
#[test]
fn the_matrix_is_the_lods() {
    let text = corpus();
    let recorded = field(&text, "matrix", "a-matrix").expect("the golden's matrix");
    let lines: Vec<&str> = recorded.lines().collect();
    let columns: Vec<&str> = lines[0].split('\t').skip(1).collect();
    let (prints, lods) = fingerprints(&text, "one-sample-disagreeing");
    let produced = matrix(&prints, &|left, right| {
        *lods
            .get(&(left.group.clone(), right.group.clone()))
            .unwrap_or(&0.0)
    });
    for (index, row) in lines[1..].iter().enumerate() {
        let values: Vec<&str> = row.split('\t').collect();
        assert_eq!(values[0], prints[index].group);
        for (column, name) in columns.iter().enumerate() {
            assert_eq!(*name, prints[column].group);
            // The file writes four decimal places, so the comparison is at that width.
            let expected: f64 = values[column + 1].parse().expect("a LOD");
            assert!(
                (expected - produced[index][column]).abs() < 5e-5,
                "{expected} vs {}",
                produced[index][column]
            );
        }
    }
}

/// The threshold on its own, which is used with both signs.
#[test]
fn the_threshold_is_used_with_both_signs() {
    // With the default of zero the middle is empty: every pair gets an answer.
    assert_eq!(verdict(true, 0.1, 0.0), Verdict::ExpectedMatch);
    assert_eq!(verdict(true, -0.1, 0.0), Verdict::UnexpectedMismatch);
    assert_eq!(verdict(false, 0.1, 0.0), Verdict::UnexpectedMatch);
    assert_eq!(verdict(false, -0.1, 0.0), Verdict::ExpectedMismatch);
    // A LOD exactly at the threshold is inconclusive, because both comparisons are strict.
    assert_eq!(verdict(true, 0.0, 0.0), Verdict::Inconclusive);
    // The band is between the threshold and its negation, so a POSITIVE threshold has no band at
    // all: it moves the cut instead, and a LOD of one under a threshold of three is a mismatch
    // rather than an uncertainty. Only a negative threshold opens a band, and the golden's
    // `-100` opens one that swallows every LOD the fixture produces.
    assert_eq!(verdict(true, 1.0, 3.0), Verdict::UnexpectedMismatch);
    assert_eq!(verdict(true, 5.0, 3.0), Verdict::ExpectedMatch);
    assert_eq!(verdict(true, -1.0, -3.0), Verdict::Inconclusive);
    assert_eq!(verdict(true, 100.0, -100.0), Verdict::Inconclusive);
    assert_eq!(verdict(false, -100.0, -100.0), Verdict::Inconclusive);
    // And the data types name themselves the way the file writes them.
    assert_eq!(DataType::ReadGroup.name(), "READGROUP");
    assert_eq!(DataType::FileAndSample.name(), "FILE_AND_SAMPLE");
}
