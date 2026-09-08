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

/// The counters the twenty columns need beyond [`Counts`], which the earlier slice did not model.
///
/// They are the reference's own hidden fields plus the two indel and two complex-indel columns:
/// the metrics file writes the RATIOS, and the ratios are computed from these.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct HiddenCounts {
    pub db_snp_insertions: i64,
    pub db_snp_deletions: i64,
    pub novel_insertions: i64,
    pub novel_deletions: i64,
    pub total_complex_indels: i64,
    pub num_in_db_snp_complex_indels: i64,
    /// The reference and alternate depths of the HETEROZYGOUS calls, which the reference bias is
    /// the ratio of. A summary row sums them over every sample; a detail row over its own.
    pub reference_allele_observations: i64,
    pub alternate_allele_observations: i64,
    /// Detail rows only: the calls whose genotype quality is zero, and the two zygosities.
    pub total_gq0_variants: i64,
    pub number_of_hets: i64,
    pub number_of_hom_var: i64,
    pub total_het_depth: i64,
}

/// One sample's call at a site, as the accumulator reads it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Call {
    /// `None` for a no-call; otherwise the allele indices, `0` being the reference.
    pub alleles: Option<Vec<usize>>,
    pub gq: i32,
    /// One depth per allele, reference first.
    pub allele_depths: Option<Vec<i32>>,
}

impl Call {
    pub fn is_called(&self) -> bool {
        self.alleles.is_some()
    }
    fn indices(&self) -> &[usize] {
        self.alleles.as_deref().unwrap_or(&[])
    }
    pub fn is_hom_ref(&self) -> bool {
        !self.indices().is_empty() && self.indices().iter().all(|allele| *allele == 0)
    }
    pub fn is_hom_var(&self) -> bool {
        !self.indices().is_empty()
            && self
                .indices()
                .iter()
                .all(|allele| *allele == self.indices()[0])
            && self.indices()[0] != 0
    }
    pub fn is_het(&self) -> bool {
        let indices = self.indices();
        !indices.is_empty() && indices.iter().any(|allele| *allele != indices[0])
    }
}

/// One site, as the accumulator reads it: the alleles, the filter, both dbSNP answers, and one
/// call per sample in the header's order.
///
/// The dbSNP membership is TWO answers and not one: the reference keeps a bitset of SNP sites and
/// another of indel sites, and asks the one that matches the site's own type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub reference: String,
    pub alternates: Vec<String>,
    pub filtered: bool,
    pub in_db_snp_snps: bool,
    pub in_db_snp_indels: bool,
    pub calls: Vec<Call>,
}

impl Site {
    /// `VariantContext.isSNP()`: one reference base against single-base alternates.
    pub fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.iter().all(|allele| allele.len() == 1)
    }
    pub fn is_biallelic(&self) -> bool {
        self.alternates.len() == 1
    }
    /// `isComplexIndel`: an indel whose reference and alternate share neither a prefix nor the
    /// simple insertion or deletion shape.
    pub fn is_complex_indel(&self) -> bool {
        if self.is_snp() || self.alternates.len() != 1 {
            return false;
        }
        let alternate = &self.alternates[0];
        !(alternate.starts_with(&self.reference) || self.reference.starts_with(alternate.as_str()))
    }
    pub fn is_simple_insertion(&self) -> bool {
        self.alternates.len() == 1 && self.alternates[0].len() > self.reference.len()
    }
    /// `isVariantExcluded`: a site with no alternate, or one where every call is homozygous
    /// reference, never reaches a counter at all.
    pub fn is_excluded(&self) -> bool {
        self.alternates.is_empty() || self.calls.iter().all(|call| call.is_hom_ref())
    }

    /// `getSingletonSample`: the one sample carrying exactly one variant chromosome.
    ///
    /// The reference sums over the FIRST TWO het-or-hom-var calls only, a het counting one and a
    /// hom-var two, and asks whether that sum is one. Two hets therefore answer no as soon as the
    /// second is seen, and a lone hom-var answers no as well.
    pub fn singleton_sample(&self) -> Option<usize> {
        let mut seen = Vec::new();
        for (index, call) in self.calls.iter().enumerate() {
            if call.is_het() || call.is_hom_var() {
                seen.push((index, if call.is_het() { 1 } else { 2 }));
                if seen.len() == 2 {
                    break;
                }
            }
        }
        let total: i32 = seen.iter().map(|(_, weight)| *weight).sum();
        if total == 1 {
            seen.last().map(|(index, _)| *index)
        } else {
            None
        }
    }
}

/// `updateSummaryMetric` for one site, against one row's counters.
///
/// `call` is the row's own genotype, and `None` for the summary row -- which is what stops the
/// reference bias from being counted twice for a one-sample file: the summary's numbers come from
/// the DETAIL rows' calls, added there.
pub fn accumulate_site(
    counts: &mut Counts,
    hidden: &mut HiddenCounts,
    summary_hidden: Option<&mut HiddenCounts>,
    site: &Site,
    call: Option<&Call>,
    has_singleton: bool,
) {
    // "If this sample's genotype doesn't have any variation, exclude it."
    if let Some(call) = call {
        if !call.is_called() {
            return;
        }
    }
    if site.filtered {
        if site.is_snp() {
            counts.filtered_snps += 1;
        } else {
            counts.filtered_indels += 1;
        }
        return;
    }
    if has_singleton {
        counts.num_singletons += 1;
    }

    if site.is_biallelic() && site.is_snp() {
        counts.total_snps += 1;
        let transition = is_transition(
            site.reference.as_bytes()[0],
            site.alternates[0].as_bytes()[0],
        );
        if site.in_db_snp_snps {
            counts.num_in_db_snp += 1;
            if transition {
                counts.db_snp_transitions += 1;
            } else {
                counts.db_snp_transversions += 1;
            }
        } else if transition {
            counts.novel_transitions += 1;
        } else {
            counts.novel_transversions += 1;
        }

        // The reference bias, which only a HETEROZYGOUS call with allele depths contributes to,
        // and which is added to the summary's counters at the same time.
        if let Some(call) = call {
            if call.is_het() {
                if let Some(depths) = &call.allele_depths {
                    if depths.len() >= 2 {
                        hidden.reference_allele_observations += i64::from(depths[0]);
                        hidden.alternate_allele_observations += i64::from(depths[1]);
                        if let Some(summary) = summary_hidden {
                            summary.reference_allele_observations += i64::from(depths[0]);
                            summary.alternate_allele_observations += i64::from(depths[1]);
                        }
                    }
                }
            }
        }
    } else if site.is_snp() && site.alternates.len() > 1 {
        counts.total_multiallelic_snps += 1;
        if site.in_db_snp_snps {
            counts.num_in_db_snp_multiallelic += 1;
        }
    } else if !site.is_snp() && !site.is_complex_indel() {
        counts.total_indels += 1;
        let insertion = site.is_simple_insertion();
        if site.in_db_snp_indels {
            counts.num_in_db_snp_indels += 1;
            if insertion {
                hidden.db_snp_insertions += 1;
            } else {
                hidden.db_snp_deletions += 1;
            }
        } else if insertion {
            hidden.novel_insertions += 1;
        } else {
            hidden.novel_deletions += 1;
        }
    } else if site.is_complex_indel() {
        hidden.total_complex_indels += 1;
        if site.in_db_snp_indels {
            hidden.num_in_db_snp_complex_indels += 1;
        }
    }
}

/// `updateDetailMetric`: the summary path, and then the three counters only a detail row has.
pub fn accumulate_detail(
    counts: &mut Counts,
    hidden: &mut HiddenCounts,
    summary_hidden: &mut HiddenCounts,
    site: &Site,
    call: &Call,
    has_singleton: bool,
) {
    accumulate_site(
        counts,
        hidden,
        Some(summary_hidden),
        site,
        Some(call),
        has_singleton,
    );
    if !site.filtered {
        if call.gq == 0 {
            hidden.total_gq0_variants += 1;
        }
        if call.is_het() {
            hidden.number_of_hets += 1;
        } else if call.is_hom_var() {
            hidden.number_of_hom_var += 1;
        }
    }
}
