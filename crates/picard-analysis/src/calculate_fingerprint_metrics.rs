//! `CalculateFingerprintMetrics`: how good a fingerprint is, with nothing to compare it against.
//!
//! Every other fingerprint tool compares two samples. This one takes one, and asks how far its
//! genotypes are from what the panel's allele frequencies would predict, which is what catches a
//! contaminated or a mixed-up file before any comparison is made.
//!
//! Most of the row is arithmetic over the same fingerprint: how many haplotypes there were, how
//! many had evidence, how many were called definitely, and how many of those were homozygous for
//! each allele. Two columns are not: the chi-squared p-values come from the incomplete gamma
//! function, and `DISCRIMINATORY_POWER` comes from a hundred random permutations of the
//! fingerprint drawn from a Mersenne Twister seeded with a constant the tool does not expose.
//!
//! Ported from `picard.fingerprint.CalculateFingerprintMetrics`,
//! `picard.fingerprint.FingerprintMetrics`, `picard.fingerprint.HaplotypeProbabilities`,
//! `picard.fingerprint.HaplotypeProbabilitiesUsingLogLikelihoods`,
//! `picard.fingerprint.HaplotypeProbabilitiesFromSequence`,
//! `picard.fingerprint.FingerprintChecker` and `picard.util.MathUtil` in Picard 3.4.0.

use crate::check_fingerprint::{block_contribution, match_result};
use crate::math3::{chi_square_test, next_permutation, MersenneTwister};

/// The LOD above which a genotype counts as definite.
pub const GENOTYPE_LOD_THRESHOLD: f64 = 3.0;
/// How many permutations the discriminatory power is averaged over.
pub const NUMBER_OF_SAMPLING: usize = 100;
/// The seed the tool fixes so that its sampled column is reproducible.
pub const RANDOM_SEED: i32 = 42;

/// The largest probability the reference lets a normalisation return, and with it the smallest.
///
/// A probability is clamped rather than allowed to reach one, which is why a fingerprint that is
/// certain of every genotype still leaves a hundredth of a quadrillionth in the other two.
pub const MAX_PROB_BELOW_ONE: f64 = 0.9999999999999999;

/// `MathUtil.pNormalizeVector`: to probabilities, clamped at both ends.
pub fn p_normalize_vector(values: &[f64]) -> Vec<f64> {
    let total: f64 = values.iter().sum();
    let max_p = MAX_PROB_BELOW_ONE;
    let min_p = (1.0 - MAX_PROB_BELOW_ONE) / (values.len() - 1) as f64;
    values
        .iter()
        .map(|value| (value / total).clamp(min_p, max_p))
        .collect()
}

/// `MathUtil.pNormalizeLogProbability`: the same, from base-ten logs, bumped so the largest is
/// ten to the three hundredth before anything is unlogged.
pub fn p_normalize_log_probability(values: [f64; 3]) -> [f64; 3] {
    let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let bump = 300.0 - maximum;
    let unlogged: Vec<f64> = values
        .iter()
        .map(|value| 10f64.powf(value + bump))
        .collect();
    let normalized = p_normalize_vector(&unlogged);
    [normalized[0], normalized[1], normalized[2]]
}

/// `MathUtil.klDivergance`, which is the cross-entropy of the counts against the expectation.
pub fn kl_divergence(measured: &[f64], distribution: &[f64]) -> f64 {
    let measured = p_normalize_vector(measured);
    let distribution = p_normalize_vector(distribution);
    -measured
        .iter()
        .zip(distribution.iter())
        .map(|(observed, expected)| observed * (expected / observed).ln())
        .sum::<f64>()
}

/// One haplotype block's evidence, in the log-likelihood form every fingerprint uses internally.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Haplotype {
    /// Base-ten log-likelihoods of the three genotypes, in the order homozygous for the major
    /// allele, heterozygous, homozygous for the minor one.
    pub log_likelihoods: [f64; 3],
    /// The block's own Hardy-Weinberg frequencies.
    pub priors: [f64; 3],
}

impl Haplotype {
    /// A block with no evidence at all.
    pub fn new(priors: [f64; 3]) -> Haplotype {
        Haplotype {
            log_likelihoods: [0.0; 3],
            priors,
        }
    }

    /// `HaplotypeProbabilitiesFromSequence.addToProbs`: one base of one read.
    ///
    /// A base that is neither of the block's two alleles is counted and contributes nothing, which
    /// is how a sequencing error at a fingerprint site is kept out of the likelihoods.
    pub fn add_base(&mut self, base: u8, allele1: u8, allele2: u8, quality: u8) {
        let error = 10f64.powf(-f64::from(quality) / 10.0);
        if base == allele1 {
            for genotype in 0..3 {
                let p_alt = genotype as f64 / 2.0;
                self.log_likelihoods[genotype] +=
                    ((1.0 - p_alt) * (1.0 - error) + p_alt * error).log10();
            }
        } else if base == allele2 {
            for genotype in 0..3 {
                let p_alt = 1.0 - genotype as f64 / 2.0;
                self.log_likelihoods[genotype] +=
                    ((1.0 - p_alt) * (1.0 - error) + p_alt * error).log10();
            }
        }
    }

    /// Whether anything was ever added: a block the reads did not reach has no evidence, and a
    /// pair of fingerprints only contributes to a LOD where both sides have some.
    pub fn has_evidence(&self) -> bool {
        self.log_likelihoods.iter().any(|value| *value != 0.0)
    }

    /// `getLikelihoods`: the evidence as probabilities.
    pub fn likelihoods(&self) -> [f64; 3] {
        p_normalize_log_probability(self.log_likelihoods)
    }

    /// `getShiftedLogPosterior`: the likelihoods with the priors folded in, still in logs.
    pub fn shifted_log_posterior(&self) -> [f64; 3] {
        [
            self.log_likelihoods[0] + self.priors[0].log10(),
            self.log_likelihoods[1] + self.priors[1].log10(),
            self.log_likelihoods[2] + self.priors[2].log10(),
        ]
    }

    /// `getPosteriorProbabilities`.
    pub fn posterior_probabilities(&self) -> [f64; 3] {
        p_normalize_log_probability(self.shifted_log_posterior())
    }

    /// `getLodMostProbableGenotype`: how much better the best genotype is than the next.
    pub fn lod_most_probable_genotype(&self) -> f64 {
        let mut sorted = self.shifted_log_posterior();
        sorted.sort_by(|left, right| right.partial_cmp(left).expect("no NaN"));
        sorted[0] - sorted[1]
    }
}

/// One base of one read over one of the map's sites.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Observation {
    pub read_name: String,
    /// Which haplotype block the site belongs to.
    pub block: usize,
    pub base: u8,
    pub allele1: u8,
    pub allele2: u8,
    pub quality: u8,
}

/// Add a file's bases to a fingerprint, at most one base per read.
///
/// A read is used ONCE, across every site and every block: two sites of one block are in linkage
/// with each other, so a read that spans both would otherwise count its single molecule twice.
/// The fixture's reads span two sites of the first block and one of the second, and the first
/// block ends up with exactly as much evidence as the second.
pub fn accumulate(haplotypes: &mut [Haplotype], observations: &[Observation]) {
    let mut used: Vec<&str> = Vec::new();
    for observation in observations {
        if used.contains(&observation.read_name.as_str()) {
            continue;
        }
        haplotypes[observation.block].add_base(
            observation.base,
            observation.allele1,
            observation.allele2,
            observation.quality,
        );
        used.push(&observation.read_name);
    }
}

/// One row of the metrics file.
#[derive(Debug, Clone, PartialEq)]
pub struct FingerprintMetrics {
    pub sample_alias: String,
    pub source: String,
    pub info: String,
    pub haplotypes: u64,
    pub haplotypes_with_evidence: u64,
    pub definite_genotypes: u64,
    pub num_hom_allele1: i64,
    pub num_hom_allele2: i64,
    pub num_hom_any: i64,
    pub num_het: i64,
    pub expected_hom_allele1: f64,
    pub expected_hom_allele2: f64,
    pub expected_het: f64,
    pub chi_squared_pvalue: f64,
    pub log10_chi_squared_pvalue: f64,
    pub cross_entropy_lod: f64,
    pub het_chi_squared_pvalue: f64,
    pub log10_het_chi_squared_pvalue: f64,
    pub het_cross_entropy_lod: f64,
    pub hom_chi_squared_pvalue: f64,
    pub log10_hom_chi_squared_pvalue: f64,
    pub hom_cross_entropy_lod: f64,
    pub lod_self_check: f64,
    pub discriminatory_power: f64,
}

/// The LOD of one fingerprint against another, block by block.
///
/// A block contributes only where both sides have evidence, which is why a fingerprint compared
/// with itself over blocks the reads never reached scores nothing rather than scoring infinitely
/// well.
pub fn lod(observed: &[Haplotype], expected: &[Haplotype]) -> f64 {
    let contributions: Vec<_> = observed
        .iter()
        .zip(expected.iter())
        .filter(|(left, right)| left.has_evidence() && right.has_evidence())
        .map(|(left, right)| {
            block_contribution(left.likelihoods(), right.likelihoods(), left.priors)
        })
        .collect();
    match_result(&contributions).lod_expected_sample
}

/// One randomization: each block's log-likelihoods shuffled, in the order the fingerprint's own
/// sorted map walks them.
///
/// The permutations come one after another from a single generator, so the blocks' order is part
/// of the answer.
pub fn randomize(haplotypes: &[Haplotype], rng: &mut MersenneTwister) -> Vec<Haplotype> {
    haplotypes
        .iter()
        .map(|haplotype| {
            let permutation = next_permutation(rng, 3);
            let mut permuted = [0.0; 3];
            for (target, source) in permutation.iter().enumerate() {
                permuted[target] = haplotype.log_likelihoods[*source];
            }
            Haplotype {
                log_likelihoods: permuted,
                priors: haplotype.priors,
            }
        })
        .collect()
}

/// The LODs of the hundred trials, in the order they are drawn.
pub fn randomized_lods(haplotypes: &[Haplotype], trials: usize, seed: i32) -> Vec<f64> {
    let mut rng = MersenneTwister::new(seed);
    (0..trials)
        .map(|_| lod(haplotypes, &randomize(haplotypes, &mut rng)))
        .collect()
}

/// The whole row for one fingerprint.
pub fn fingerprint_metrics(
    sample_alias: &str,
    source: &str,
    info: &str,
    haplotypes: &[Haplotype],
) -> FingerprintMetrics {
    // The counts are the posterior probabilities added up, so a fingerprint that is sure of every
    // genotype has whole numbers in them and one that is not has fractions.
    let mut counts = [0.0; 3];
    let mut expected = [0.0; 3];
    for haplotype in haplotypes {
        let posterior = haplotype.posterior_probabilities();
        for genotype in 0..3 {
            counts[genotype] += posterior[genotype];
            expected[genotype] += haplotype.priors[genotype];
        }
    }

    let hom_vs_het_expected = [expected[0] + expected[2], expected[1]];
    let hom_vs_het_counts = [counts[0] + counts[2], counts[1]];
    let hom1_vs_hom2_expected = [expected[0], expected[2]];
    let hom1_vs_hom2_counts = [counts[0], counts[2]];

    let round =
        |values: &[f64]| -> Vec<i64> { values.iter().map(|value| value.round() as i64).collect() };
    let rounded_counts = round(&counts);
    let rounded_hom_vs_het = round(&hom_vs_het_counts);
    let rounded_hom1_vs_hom2 = round(&hom1_vs_hom2_counts);

    let chi_squared = chi_square_test(&expected, &rounded_counts);
    let het_chi_squared = chi_square_test(&hom_vs_het_expected, &rounded_hom_vs_het);
    let hom_chi_squared = chi_square_test(&hom1_vs_hom2_expected, &rounded_hom1_vs_hom2);

    let self_check = lod(haplotypes, haplotypes);
    let trials = randomized_lods(haplotypes, NUMBER_OF_SAMPLING, RANDOM_SEED);
    let mean: f64 = trials.iter().sum::<f64>() / trials.len() as f64;

    FingerprintMetrics {
        sample_alias: sample_alias.to_string(),
        source: source.to_string(),
        info: info.to_string(),
        haplotypes: haplotypes.len() as u64,
        haplotypes_with_evidence: haplotypes
            .iter()
            .filter(|haplotype| haplotype.has_evidence())
            .count() as u64,
        definite_genotypes: haplotypes
            .iter()
            .filter(|haplotype| haplotype.lod_most_probable_genotype() >= GENOTYPE_LOD_THRESHOLD)
            .count() as u64,
        num_hom_allele1: rounded_counts[0],
        num_hom_allele2: rounded_counts[2],
        // The homozygous total is the second of the hom-versus-het pair, which is the HET count:
        // the pair is built with the homozygotes first and read back the other way round.
        num_hom_any: rounded_hom_vs_het[1],
        num_het: rounded_counts[1],
        expected_hom_allele1: expected[0],
        expected_hom_allele2: expected[2],
        expected_het: expected[1],
        chi_squared_pvalue: chi_squared,
        log10_chi_squared_pvalue: chi_squared.log10(),
        cross_entropy_lod: kl_divergence(&counts, &expected),
        het_chi_squared_pvalue: het_chi_squared,
        log10_het_chi_squared_pvalue: het_chi_squared.log10(),
        het_cross_entropy_lod: kl_divergence(&hom_vs_het_counts, &hom_vs_het_expected),
        hom_chi_squared_pvalue: hom_chi_squared,
        log10_hom_chi_squared_pvalue: hom_chi_squared.log10(),
        hom_cross_entropy_lod: kl_divergence(&hom1_vs_hom2_counts, &hom1_vs_hom2_expected),
        lod_self_check: self_check,
        discriminatory_power: self_check - mean,
    }
}

/// The metrics file's header.
pub const HEADER: &str = "SAMPLE_ALIAS\tSOURCE\tINFO\tHAPLOTYPES\tHAPLOTYPES_WITH_EVIDENCE\t\
DEFINITE_GENOTYPES\tNUM_HOM_ALLELE1\tNUM_HOM_ALLELE2\tNUM_HOM_ANY\tNUM_HET\tEXPECTED_HOM_ALLELE1\t\
EXPECTED_HOM_ALLELE2\tEXPECTED_HET\tCHI_SQUARED_PVALUE\tLOG10_CHI_SQUARED_PVALUE\t\
CROSS_ENTROPY_LOD\tHET_CHI_SQUARED_PVALUE\tLOG10_HET_CHI_SQUARED_PVALUE\tHET_CROSS_ENTROPY_LOD\t\
HOM_CHI_SQUARED_PVALUE\tLOG10_HOM_CHI_SQUARED_PVALUE\tHOM_CROSS_ENTROPY_LOD\tLOD_SELF_CHECK\t\
DISCRIMINATORY_POWER";

/// `FormatUtil.format(double)`: at most six decimal places, and no trailing zeroes.
pub fn format_double(value: f64) -> String {
    let text = format!("{value:.6}");
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

/// One row, as the metrics file writes it.
pub fn render_row(row: &FingerprintMetrics) -> String {
    let numbers = [
        row.expected_hom_allele1,
        row.expected_hom_allele2,
        row.expected_het,
        row.chi_squared_pvalue,
        row.log10_chi_squared_pvalue,
        row.cross_entropy_lod,
        row.het_chi_squared_pvalue,
        row.log10_het_chi_squared_pvalue,
        row.het_cross_entropy_lod,
        row.hom_chi_squared_pvalue,
        row.log10_hom_chi_squared_pvalue,
        row.hom_cross_entropy_lod,
        row.lod_self_check,
        row.discriminatory_power,
    ];
    let mut text = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.sample_alias,
        row.source,
        row.info,
        row.haplotypes,
        row.haplotypes_with_evidence,
        row.definite_genotypes,
        row.num_hom_allele1,
        row.num_hom_allele2,
        row.num_hom_any,
        row.num_het
    );
    for number in numbers {
        text.push('\t');
        text.push_str(&format_double(number));
    }
    text
}

/// The table, header and rows.
pub fn render(rows: &[FingerprintMetrics]) -> String {
    let mut text = String::from(HEADER);
    for row in rows {
        text.push('\n');
        text.push_str(&render_row(row));
    }
    text
}
