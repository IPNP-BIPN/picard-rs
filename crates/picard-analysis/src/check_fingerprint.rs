//! `CheckFingerprint`: whether a file's sample is the sample its genotypes claim.
//!
//! The answer is a log-odds score over haplotype blocks: how much more likely the evidence is if
//! the two fingerprints came from one sample than if they came from two. Everything else the tool
//! writes is bookkeeping around that number.
//!
//! Ported from `picard.fingerprint.CheckFingerprint`, `picard.fingerprint.FingerprintChecker` and
//! `picard.fingerprint.HaplotypeProbabilities`.

/// The two suffixes `--OUTPUT` builds its files from.
pub const SUMMARY_FILE_SUFFIX: &str = ".fingerprinting_summary_metrics";
pub const DETAIL_FILE_SUFFIX: &str = ".fingerprinting_detail_metrics";

/// The two files a run writes, from one basename.
pub fn file_names(basename: &str) -> [String; 2] {
    [
        format!("{basename}{SUMMARY_FILE_SUFFIX}"),
        format!("{basename}{DETAIL_FILE_SUFFIX}"),
    ]
}

/// `EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND`: the genotypes do not carry the sample, so nothing
/// is written and the code is the answer.
pub const EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND: i32 = 1;
/// `EXIT_CODE_WHEN_NO_VALID_CHECKS`: the run produced its files and every comparison in them was
/// inconclusive, so a caller that reads the metrics without reading the code sees numbers that
/// answer nothing.
pub const EXIT_CODE_WHEN_NO_VALID_CHECKS: i32 = 2;

/// `scaledEvidenceProbabilityUsingGenotypeFrequencies`, which is the dot product of the
/// likelihoods and the frequencies used as priors.
pub fn scaled_evidence_probability(likelihoods: [f64; 3], frequencies: [f64; 3]) -> f64 {
    likelihoods
        .iter()
        .zip(frequencies.iter())
        .map(|(likelihood, frequency)| likelihood * frequency)
        .sum()
}

/// `shiftedLogEvidenceProbabilityUsingGenotypeFrequencies`: its base-ten logarithm.
pub fn shifted_log_evidence_probability(likelihoods: [f64; 3], frequencies: [f64; 3]) -> f64 {
    scaled_evidence_probability(likelihoods, frequencies).log10()
}

/// `getPosteriorLikelihoods`: the likelihoods multiplied by the priors, NOT normalised.
///
/// The name says posterior and the value is not one: it is the numerator of Bayes' rule, and it is
/// used as the prior of the other fingerprint, which is what makes the two-sample model a product
/// rather than a comparison.
pub fn posterior_likelihoods(likelihoods: [f64; 3], priors: [f64; 3]) -> [f64; 3] {
    [
        likelihoods[0] * priors[0],
        likelihoods[1] * priors[1],
        likelihoods[2] * priors[2],
    ]
}

/// One block's contribution to the two models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlockContribution {
    /// The evidence under one sample: the observed likelihoods against the expected posterior.
    pub no_swap: f64,
    /// The evidence under two: each fingerprint against the population priors.
    pub swap: f64,
}

/// What one haplotype block contributes, given both fingerprints' likelihoods and the block's own
/// population frequencies.
///
/// A block the reads do not cover is NOT skipped: its observed likelihoods are the population
/// priors, so it contributes a small positive term and keeps its row in the detail file. A file
/// that covers one of two blocks therefore writes the same two rows as one that covers both, with
/// the uncovered row reading zero observations and a LOD near zero.
pub fn block_contribution(
    observed: [f64; 3],
    expected: [f64; 3],
    priors: [f64; 3],
) -> BlockContribution {
    BlockContribution {
        no_swap: shifted_log_evidence_probability(
            observed,
            posterior_likelihoods(expected, priors),
        ),
        swap: shifted_log_evidence_probability(observed, priors)
            + shifted_log_evidence_probability(expected, priors),
    }
}

/// The summary's three numbers, accumulated over the blocks that had evidence on both sides.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatchResult {
    pub ll_expected_sample: f64,
    pub ll_random_sample: f64,
    pub lod_expected_sample: f64,
}

/// `calculateMatchResults` as far as the summary goes.
pub fn match_result(contributions: &[BlockContribution]) -> MatchResult {
    let ll_expected_sample: f64 = contributions.iter().map(|block| block.no_swap).sum();
    let ll_random_sample: f64 = contributions.iter().map(|block| block.swap).sum();
    MatchResult {
        ll_expected_sample,
        ll_random_sample,
        // The LOD is the difference of the two logs, so a sample that matches is positive.
        lod_expected_sample: ll_expected_sample - ll_random_sample,
    }
}
