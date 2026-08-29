//! `CrosscheckFingerprints`: every pair of fingerprints against every other.
//!
//! [`crate::check_fingerprint`] asks whether one file's sample is the sample a set of genotypes
//! says it is. This one has no genotypes to compare against: it compares the inputs with each
//! other, and what is EXPECTED comes from the sample names rather than from a truth set. Two read
//! groups of one sample are expected to match; a match between two samples is unexpected however
//! good the LOD is.
//!
//! # The verdict is two questions, not one
//!
//! `getMatchResults` asks whether the pair was expected to match and then where the LOD falls, and
//! the threshold is used with BOTH signs. A LOD below `LOD_THRESHOLD` is a mismatch, one above
//! `-LOD_THRESHOLD` is a match, and anything between them is `INCONCLUSIVE`. With the default
//! threshold of zero the middle is empty and every pair gets an answer; with a negative threshold
//! the middle is everything, which is why `--LOD_THRESHOLD -100` makes a whole run inconclusive.
//!
//! # What is not ported
//!
//! The LOD itself, for the same reason `check_fingerprint` does not compute it from reads: the
//! pileup and the haplotype likelihoods are the fingerprint's, and what this module is about is
//! what the tool does with a number once it has one. The golden's LODs are the input to the tests.
//!
//! Ported from `picard.fingerprint.CrosscheckFingerprints` and
//! `picard.fingerprint.CrosscheckMetric` in Picard 3.4.0.

/// `CrosscheckMetric.FingerprintResult`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    ExpectedMatch,
    ExpectedMismatch,
    UnexpectedMatch,
    UnexpectedMismatch,
    Inconclusive,
}

impl Verdict {
    /// The name the metrics file writes.
    pub fn name(self) -> &'static str {
        match self {
            Verdict::ExpectedMatch => "EXPECTED_MATCH",
            Verdict::ExpectedMismatch => "EXPECTED_MISMATCH",
            Verdict::UnexpectedMatch => "UNEXPECTED_MATCH",
            Verdict::UnexpectedMismatch => "UNEXPECTED_MISMATCH",
            Verdict::Inconclusive => "INCONCLUSIVE",
        }
    }

    /// `isExpected()`, which is `None` for the inconclusive verdict rather than false.
    pub fn is_expected(self) -> Option<bool> {
        match self {
            Verdict::ExpectedMatch | Verdict::ExpectedMismatch => Some(true),
            Verdict::UnexpectedMatch | Verdict::UnexpectedMismatch => Some(false),
            Verdict::Inconclusive => None,
        }
    }

    /// `isMatch()`, `None` for the same reason.
    pub fn is_match(self) -> Option<bool> {
        match self {
            Verdict::ExpectedMatch | Verdict::UnexpectedMatch => Some(true),
            Verdict::ExpectedMismatch | Verdict::UnexpectedMismatch => Some(false),
            Verdict::Inconclusive => None,
        }
    }
}

/// `CrosscheckMetric.DataType`: what a row is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    ReadGroup,
    Library,
    Sample,
    File,
    FileAndSample,
}

impl DataType {
    pub fn name(self) -> &'static str {
        match self {
            DataType::ReadGroup => "READGROUP",
            DataType::Library => "LIBRARY",
            DataType::Sample => "SAMPLE",
            DataType::File => "FILE",
            DataType::FileAndSample => "FILE_AND_SAMPLE",
        }
    }
}

/// `getMatchResults`: the verdict, from what was expected and where the LOD fell.
///
/// The comparisons are strict on both sides, so a LOD exactly at the threshold is inconclusive
/// rather than either answer.
pub fn verdict(expected_to_match: bool, lod: f64, threshold: f64) -> Verdict {
    if expected_to_match {
        if lod < threshold {
            Verdict::UnexpectedMismatch
        } else if lod > -threshold {
            Verdict::ExpectedMatch
        } else {
            Verdict::Inconclusive
        }
    } else if lod > -threshold {
        Verdict::UnexpectedMatch
    } else if lod < threshold {
        Verdict::ExpectedMismatch
    } else {
        Verdict::Inconclusive
    }
}

/// `CROSSCHECK_MODE`, which decides which pairs are compared at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// The default: only pairs whose samples agree.
    CheckSameSample,
    /// Every pair, whatever the samples say.
    CheckAllOthers,
}

/// One fingerprint, reduced to what a crosscheck reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint {
    /// The value the row is keyed by, which is what `--CROSSCHECK_BY` chose.
    pub group: String,
    pub sample: String,
    pub library: String,
    pub file: String,
}

/// Whether a pair is compared under a mode.
pub fn is_compared(left: &Fingerprint, right: &Fingerprint, mode: Mode) -> bool {
    match mode {
        Mode::CheckAllOthers => true,
        Mode::CheckSameSample => left.sample == right.sample,
    }
}

/// One row of the metrics file.
#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub left: String,
    pub right: String,
    pub verdict: Verdict,
    pub data_type: DataType,
    pub lod: f64,
}

/// The arguments a run's shape depends on.
#[derive(Debug, Clone)]
pub struct Options {
    pub lod_threshold: f64,
    pub crosscheck_by: DataType,
    pub mode: Mode,
    pub output_errors_only: bool,
    pub expect_all_groups_to_match: bool,
    pub exit_code_when_mismatch: i32,
    pub exit_code_when_no_valid_checks: i32,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            lod_threshold: 0.0,
            crosscheck_by: DataType::ReadGroup,
            mode: Mode::CheckSameSample,
            output_errors_only: false,
            expect_all_groups_to_match: false,
            exit_code_when_mismatch: 1,
            exit_code_when_no_valid_checks: 1,
        }
    }
}

/// Every pair's row, in the order the fingerprints were given.
///
/// A pair is expected to match when the two samples agree, or when
/// `--EXPECT_ALL_GROUPS_TO_MATCH` says every pair is, which is the argument that turns a file of
/// several samples into a file that had better be one.
pub fn rows(
    fingerprints: &[Fingerprint],
    lods: &dyn Fn(&Fingerprint, &Fingerprint) -> f64,
    options: &Options,
) -> Vec<Row> {
    let mut rows = Vec::new();
    for left in fingerprints {
        for right in fingerprints {
            if !is_compared(left, right, options.mode) {
                continue;
            }
            let expected = options.expect_all_groups_to_match || left.sample == right.sample;
            let lod = lods(left, right);
            let verdict = verdict(expected, lod, options.lod_threshold);
            // `OUTPUT_ERRORS_ONLY` keeps a row that is inconclusive as well as one that is
            // unexpected: what it drops is agreement, not uncertainty.
            if options.output_errors_only
                && verdict != Verdict::Inconclusive
                && verdict.is_expected() == Some(true)
            {
                continue;
            }
            rows.push(Row {
                left: left.group.clone(),
                right: right.group.clone(),
                verdict,
                data_type: options.crosscheck_by,
                lod,
            });
        }
    }
    rows
}

/// The status a run exits with, over the rows it produced.
///
/// Three answers and their order matters. A run whose every LOD is zero compared nothing, whatever
/// its rows say, and returns `EXIT_CODE_WHEN_NO_VALID_CHECKS`. Otherwise a run with an unexpected
/// verdict returns `EXIT_CODE_WHEN_MISMATCH`, and a run without one returns zero. An INCONCLUSIVE
/// row is not an unexpected one: `isExpected()` is neither true nor false for it.
pub fn exit_code(rows: &[Row], options: &Options) -> i32 {
    if rows.iter().all(|row| row.lod == 0.0) {
        return options.exit_code_when_no_valid_checks;
    }
    let unexpected = rows
        .iter()
        .filter(|row| row.verdict.is_expected() == Some(false))
        .count();
    if unexpected > 0 {
        options.exit_code_when_mismatch
    } else {
        0
    }
}

/// `--MATRIX_OUTPUT`: the same numbers as a square, one row and one column per fingerprint.
///
/// The cell is the LOD and nothing else, so the matrix carries no verdict at all: a reader has to
/// know the threshold to read it.
pub fn matrix(
    fingerprints: &[Fingerprint],
    lods: &dyn Fn(&Fingerprint, &Fingerprint) -> f64,
) -> Vec<Vec<f64>> {
    fingerprints
        .iter()
        .map(|left| fingerprints.iter().map(|right| lods(left, right)).collect())
        .collect()
}
