//! `CollectVariantCallingMetrics`: a VCF counted against a dbSNP one.
//!
//! Reading the VCFs and the dbSNP bitset are not ported. What is ported is the tally: which
//! variant reaches which counter, and the derived columns the counts produce.
//!
//! Ported from `picard.vcf.CollectVariantCallingMetrics` and
//! `picard.vcf.CallingMetricAccumulator` in Picard 3.4.0.

use std::collections::BTreeMap;

use crate::accumulate_variant_calling_metrics::{DETAIL_EXTENSION, SUMMARY_EXTENSION};

/// The two file names the `--OUTPUT` prefix stands for, which are the accumulator's own.
pub fn file_names(prefix: &str) -> (String, String) {
    (
        format!("{prefix}.{DETAIL_EXTENSION}"),
        format!("{prefix}.{SUMMARY_EXTENSION}"),
    )
}

/// A base change, which decides the TI/TV columns.
///
/// A transition is a purine for a purine or a pyrimidine for a pyrimidine: `A<->G` and `C<->T`.
/// Everything else is a transversion.
pub fn is_transition(reference: u8, alternate: u8) -> bool {
    matches!(
        (
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase()
        ),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

/// One variant, reduced to what the tally reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub reference: String,
    pub alternates: Vec<String>,
    pub filtered: bool,
    pub in_db_snp: bool,
    /// One genotype per sample, as its two allele indices, `None` for a no-call.
    pub genotypes: Vec<[Option<usize>; 2]>,
}

impl Variant {
    /// A SNP is a single reference base against single alternate bases.
    pub fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.iter().all(|a| a.len() == 1)
    }

    /// Which is multiallelic once it carries more than one alternate.
    pub fn is_multiallelic(&self) -> bool {
        self.is_snp() && self.alternates.len() > 1
    }

    pub fn is_indel(&self) -> bool {
        !self.is_snp()
    }

    /// A singleton is a variant exactly one sample carries an alternate of.
    pub fn is_singleton(&self) -> bool {
        self.genotypes
            .iter()
            .filter(|genotype| {
                genotype
                    .iter()
                    .any(|allele| matches!(allele, Some(i) if *i > 0))
            })
            .count()
            == 1
    }
}

/// The counts one row of either table holds.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Counts {
    pub total_snps: i64,
    pub num_in_db_snp: i64,
    pub novel_snps: i64,
    pub filtered_snps: i64,
    pub db_snp_transitions: i64,
    pub db_snp_transversions: i64,
    pub novel_transitions: i64,
    pub novel_transversions: i64,
    pub total_indels: i64,
    pub novel_indels: i64,
    pub filtered_indels: i64,
    pub num_in_db_snp_indels: i64,
    pub total_multiallelic_snps: i64,
    pub num_in_db_snp_multiallelic: i64,
    pub num_singletons: i64,
}

impl Counts {
    /// `PCT_DBSNP`, a division of the known by the total, which is NaN on an empty file.
    pub fn pct_db_snp(&self) -> f64 {
        self.num_in_db_snp as f64 / self.total_snps as f64
    }

    /// `DBSNP_TITV`, which is NOUGHT and not NaN when nothing is known: nought over nought is
    /// NaN, but nought transitions over nought transversions is written as a plain zero here
    /// because the counts are integers and the division is guarded nowhere. The golden shows both
    /// kinds of empty side by side.
    pub fn db_snp_titv(&self) -> f64 {
        self.db_snp_transitions as f64 / self.db_snp_transversions as f64
    }

    pub fn novel_titv(&self) -> f64 {
        self.novel_transitions as f64 / self.novel_transversions as f64
    }
}

/// `CallingMetricAccumulator.accumulate`: one variant against the counters.
///
/// A FILTERED variant is counted as filtered and nowhere else, so it reaches neither the known nor
/// the novel tally. A multiallelic SNP is counted in its own column and not among the plain ones.
pub fn accumulate(counts: &mut Counts, variant: &Variant) {
    if variant.is_indel() {
        counts.total_indels += 1;
        if variant.filtered {
            counts.filtered_indels += 1;
            return;
        }
        if variant.in_db_snp {
            counts.num_in_db_snp_indels += 1;
        } else {
            counts.novel_indels += 1;
        }
        if variant.is_singleton() {
            counts.num_singletons += 1;
        }
        return;
    }
    if variant.is_multiallelic() {
        counts.total_multiallelic_snps += 1;
        if variant.in_db_snp {
            counts.num_in_db_snp_multiallelic += 1;
        }
        if variant.is_singleton() {
            counts.num_singletons += 1;
        }
        return;
    }
    counts.total_snps += 1;
    if variant.filtered {
        counts.filtered_snps += 1;
        return;
    }
    let transition = is_transition(
        variant.reference.as_bytes()[0],
        variant.alternates[0].as_bytes()[0],
    );
    if variant.in_db_snp {
        counts.num_in_db_snp += 1;
        if transition {
            counts.db_snp_transitions += 1;
        } else {
            counts.db_snp_transversions += 1;
        }
    } else {
        counts.novel_snps += 1;
        if transition {
            counts.novel_transitions += 1;
        } else {
            counts.novel_transversions += 1;
        }
    }
    if variant.is_singleton() {
        counts.num_singletons += 1;
    }
}

/// The whole walk: one summary row for the file and one detail row per sample.
///
/// The summary counts a variant ONCE however many samples carry it, while a detail row counts it
/// only for the sample that does. The summary is therefore not the detail rows' sum whenever a
/// sample is homozygous reference somewhere.
pub fn collect(variants: &[Variant], samples: &[String]) -> (BTreeMap<String, Counts>, Counts) {
    let mut summary = Counts::default();
    let mut details: BTreeMap<String, Counts> = samples
        .iter()
        .map(|sample| (sample.clone(), Counts::default()))
        .collect();
    for variant in variants {
        accumulate(&mut summary, variant);
        for (index, sample) in samples.iter().enumerate() {
            let carries = variant.genotypes.get(index).is_some_and(|genotype| {
                genotype
                    .iter()
                    .any(|allele| matches!(allele, Some(i) if *i > 0))
            });
            if carries {
                accumulate(details.get_mut(sample).expect("a row"), variant);
            }
        }
    }
    (details, summary)
}
