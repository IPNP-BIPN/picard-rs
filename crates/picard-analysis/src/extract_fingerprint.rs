//! `ExtractFingerprint`: the contaminating sample's genotype likelihoods at a haplotype map's
//! sites.
//!
//! Reading the BAM and building the pileup are not ported. What is ported is the model: the 3x3
//! log-likelihood matrix a base contributes, the sum over the main sample's genotype that turns it
//! into three likelihoods, and the phred scaling that writes them.
//!
//! Ported from `picard.fingerprint.ExtractFingerprint`,
//! `picard.fingerprint.HaplotypeProbabilitiesFromContaminatorSequence` and
//! `picard.fingerprint.IdentifyContaminant` in Picard 3.4.0.

/// `ExtractFingerprint.LOCUS_MAX_READS`.
pub const DEFAULT_LOCUS_MAX_READS: i32 = 50;
/// `IdentifyContaminant.LOCUS_MAX_READS`, which is not the same number.
pub const IDENTIFY_CONTAMINANT_LOCUS_MAX_READS: i32 = 200;
/// `doWork`, on a file naming more than one sample.
pub fn wrong_fingerprint_count_message(count: usize) -> String {
    format!("Expected exactly 1 fingerprint in Input file, found {count}")
}

/// The three genotypes, whose value is the count of alternate alleles.
pub const GENOTYPES: [usize; 3] = [0, 1, 2];

/// `QualityUtil.getErrorProbabilityFromPhredScore`.
pub fn error_probability(phred: u8) -> f64 {
    10f64.powf(-f64::from(phred) / 10.0)
}

/// `HaplotypeProbabilitiesFromContaminatorSequence`: the nine models kept apart until every read
/// has been seen.
#[derive(Debug, Clone, PartialEq)]
pub struct ContaminatorProbabilities {
    /// Indexed by the contaminant's genotype and then the main sample's.
    pub log_likelihoods: [[f64; 3]; 3],
    pub observed_allele1: i64,
    pub observed_allele2: i64,
    /// A base matching neither allele, which reaches no likelihood at all.
    pub observed_other: i64,
    pub contamination: f64,
}

impl ContaminatorProbabilities {
    pub fn new(contamination: f64) -> ContaminatorProbabilities {
        ContaminatorProbabilities {
            log_likelihoods: [[0.0; 3]; 3],
            observed_allele1: 0,
            observed_allele2: 0,
            observed_other: 0,
            contamination,
        }
    }

    /// `addToProbs`: one base against one SNP.
    ///
    /// A base matching neither allele is counted and returns, so it moves no likelihood. The nine
    /// models are updated together because the main sample's genotype cannot be summed out until
    /// every read has been seen.
    pub fn add(&mut self, base: u8, allele1: u8, allele2: u8, quality: u8) {
        let alternate = if base == allele1 {
            self.observed_allele1 += 1;
            false
        } else if base == allele2 {
            self.observed_allele2 += 1;
            true
        } else {
            self.observed_other += 1;
            return;
        };
        let error = error_probability(quality);
        for contaminant in GENOTYPES {
            for main in GENOTYPES {
                // The expected frequency of the alternate allele under the two genotypes.
                let theta = 0.5
                    * ((1.0 - self.contamination) * main as f64
                        + self.contamination * contaminant as f64);
                let matching = if alternate { theta } else { 1.0 - theta };
                let opposing = if alternate { 1.0 - theta } else { theta };
                self.log_likelihoods[contaminant][main] +=
                    (matching * (1.0 - error) + opposing * error).log10();
            }
        }
    }

    /// `updateLikelihoods`: the main sample's genotype summed out under the block's priors.
    ///
    /// The priors are the haplotype block's own frequencies, which come from its minor-allele
    /// frequency, so two blocks of the same pileup do not report the same likelihoods.
    pub fn log_likelihoods(&self, priors: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for contaminant in GENOTYPES {
            let mut total = 0.0;
            for main in GENOTYPES {
                total += priors[main] * 10f64.powf(self.log_likelihoods[contaminant][main]);
            }
            out[contaminant] = total.log10();
        }
        out
    }
}

/// `HaplotypeBlock.getHaplotypeFrequencies`: the Hardy-Weinberg frequencies of the three
/// genotypes, from the minor-allele frequency.
pub fn haplotype_frequencies(minor_allele_frequency: f64) -> [f64; 3] {
    let q = minor_allele_frequency;
    let p = 1.0 - q;
    [p * p, 2.0 * p * q, q * q]
}

/// The PLs a VCF carries: the log likelihoods scaled to phred and shifted so the best is nought,
/// then rounded.
pub fn phred_likelihoods(log_likelihoods: [f64; 3]) -> [i32; 3] {
    let best = log_likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    log_likelihoods.map(|value| (-10.0 * (value - best)).round() as i32)
}

/// `ExtractFingerprint.doWork`: the contamination the model is given.
///
/// `--EXTRACT_CONTAMINATION` flips the ARGUMENT and not the output, so the same number means
/// opposite things under the two settings. `IdentifyContaminant` sets that flag from the negation
/// of its own `--EXTRACT_CONTAMINATED`, which is the whole of that tool.
pub fn contamination_to_use(contamination: f64, extract_contamination: bool) -> f64 {
    if extract_contamination {
        contamination
    } else {
        1.0 - contamination
    }
}

/// `IdentifyContaminant.doWork`, which is the whole of that tool: it sets the other's
/// `EXTRACT_CONTAMINATION` from the NEGATION of its own `EXTRACT_CONTAMINATED` and delegates.
///
/// So its default is the opposite one. `EXTRACT_CONTAMINATED` is false by default, which makes
/// `EXTRACT_CONTAMINATION` true, and a default run reports the contaminant where
/// `ExtractFingerprint`'s default reports the contaminated sample.
pub fn extract_contamination_for_identify(extract_contaminated: bool) -> bool {
    !extract_contaminated
}

/// `getSampleToUse`: the sample column's name.
///
/// Without an alias the header's sample name gets `-contaminant` appended, and only when the
/// contaminant is what is being extracted. An alias replaces the name outright rather than adding
/// to it.
pub fn sample_to_use(
    header_sample: &str,
    alias: Option<&str>,
    extract_contamination: bool,
) -> String {
    match alias {
        Some(alias) => alias.to_string(),
        None if extract_contamination => format!("{header_sample}-contaminant"),
        None => header_sample.to_string(),
    }
}
