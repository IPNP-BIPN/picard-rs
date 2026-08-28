//! `AccumulateVariantCallingMetrics`: several variant-calling metrics tables merged into one.
//!
//! There is no VCF in this tool. Reading and writing the metrics files are not ported; the merge
//! is, and the merge is not a plain sum.
//!
//! Ported from `picard.vcf.AccumulateVariantCallingMetrics` and
//! `picard.vcf.CollectVariantCallingMetrics` in Picard 3.4.0.

use std::collections::BTreeMap;

/// `VariantCallingDetailMetrics.getFileExtension`.
pub const DETAIL_EXTENSION: &str = "variant_calling_detail_metrics";
/// `VariantCallingSummaryMetrics.getFileExtension`.
pub const SUMMARY_EXTENSION: &str = "variant_calling_summary_metrics";

/// `doWork`, on a summary file that does not hold exactly one row.
pub fn wrong_summary_row_count_message(count: usize) -> String {
    format!("Expected 1 row in the summary metrics file but saw {count}")
}

/// The two file names a prefix stands for. The arguments are PREFIXES and not files.
pub fn file_names(prefix: &str) -> (String, String) {
    (
        format!("{prefix}.{DETAIL_EXTENSION}"),
        format!("{prefix}.{SUMMARY_EXTENSION}"),
    )
}

/// `invertFromRatio`: given `X/Y` and `X+Y`, answers `Y`, ROUNDED.
///
/// This is where the loss comes from. A sum that the ratio does not divide evenly rounds, and the
/// recomputed ratio afterwards is then not the one that was read: 301 at a ratio of 2.0 gives 100,
/// and 301 minus 100 over 100 is 2.01.
///
/// A ratio of NaN answers NOUGHT rather than propagating, which turns "no ratio" into a ratio of
/// zero once it is recomputed.
pub fn invert_from_ratio(sum: i64, ratio: f64) -> i64 {
    if ratio.is_nan() {
        0
    } else {
        (sum as f64 / (ratio + 1.0)).round() as i64
    }
}

/// One row of the summary table, reduced to what the merge reads and writes.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SummaryMetrics {
    pub total_snps: i64,
    pub num_in_db_snp: i64,
    pub novel_snps: i64,
    pub dbsnp_titv: f64,
    pub novel_titv: f64,
    pub snp_reference_bias: f64,
    // The hidden fields, reconstructed on read and recomputed on write.
    pub dbsnp_transitions: i64,
    pub dbsnp_transversions: i64,
    pub novel_transitions: i64,
    pub novel_transversions: i64,
    pub reference_allele_observations: i64,
    pub alternate_allele_observations: i64,
}

impl SummaryMetrics {
    /// `calculateFromDerivedFields`: the hidden counts rebuilt from the printed ratios.
    ///
    /// The total het depth comes from the DETAIL file beside this summary, summed over its rows,
    /// which is why a pair of files is read together and a summary alone would rebuild
    /// differently.
    pub fn from_derived_fields(&mut self, total_het_depth: i64) {
        self.dbsnp_transversions = invert_from_ratio(self.num_in_db_snp, self.dbsnp_titv);
        self.dbsnp_transitions = self.num_in_db_snp - self.dbsnp_transversions;
        self.novel_transversions = invert_from_ratio(self.novel_snps, self.novel_titv);
        self.novel_transitions = self.novel_snps - self.novel_transversions;
        self.reference_allele_observations = if self.snp_reference_bias.is_nan() {
            0
        } else {
            (total_het_depth as f64 * self.snp_reference_bias).round() as i64
        };
        self.alternate_allele_observations = total_het_depth - self.reference_allele_observations;
    }

    /// `merge`: the counts add and the hidden counts add; the ratios are not touched here.
    pub fn merge(&mut self, other: &SummaryMetrics) {
        self.total_snps += other.total_snps;
        self.num_in_db_snp += other.num_in_db_snp;
        self.novel_snps += other.novel_snps;
        self.dbsnp_transitions += other.dbsnp_transitions;
        self.dbsnp_transversions += other.dbsnp_transversions;
        self.novel_transitions += other.novel_transitions;
        self.novel_transversions += other.novel_transversions;
        self.reference_allele_observations += other.reference_allele_observations;
        self.alternate_allele_observations += other.alternate_allele_observations;
    }

    /// `calculateDerivedFields`: the ratios recomputed from the merged hidden counts.
    pub fn derived_fields(&mut self) {
        self.dbsnp_titv = self.dbsnp_transitions as f64 / self.dbsnp_transversions as f64;
        self.novel_titv = self.novel_transitions as f64 / self.novel_transversions as f64;
        let total = self.reference_allele_observations + self.alternate_allele_observations;
        self.snp_reference_bias = self.reference_allele_observations as f64 / total as f64;
    }
}

/// One row of the detail table.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DetailMetrics {
    pub sample_alias: String,
    pub summary: SummaryMetrics,
    pub total_het_depth: i64,
    pub het_homvar_ratio: f64,
    pub number_of_hets: i64,
    pub number_of_hom_var: i64,
}

impl DetailMetrics {
    pub fn from_derived_fields(&mut self) {
        self.number_of_hom_var = invert_from_ratio(self.summary.total_snps, self.het_homvar_ratio);
        self.number_of_hets = self.summary.total_snps - self.number_of_hom_var;
        self.summary.from_derived_fields(self.total_het_depth);
    }

    pub fn merge(&mut self, other: &DetailMetrics) {
        self.summary.merge(&other.summary);
        self.number_of_hets += other.number_of_hets;
        self.number_of_hom_var += other.number_of_hom_var;
        if self.sample_alias.is_empty() {
            self.sample_alias = other.sample_alias.clone();
        }
    }

    pub fn derived_fields(&mut self) {
        self.summary.derived_fields();
        self.het_homvar_ratio = self.number_of_hets as f64 / self.number_of_hom_var as f64;
        self.total_het_depth =
            self.summary.reference_allele_observations + self.summary.alternate_allele_observations;
    }
}

/// One input: a detail table and the single summary row beside it.
#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    pub detail: Vec<DetailMetrics>,
    pub summary: Vec<SummaryMetrics>,
}

/// `doWork`: every input read, merged per sample, and the derived fields recomputed.
///
/// The detail rows are merged by SAMPLE_ALIAS. The summary is merged into one row, and the total
/// het depth handed to each input's summary is summed over THAT input's detail rows.
///
/// The reference walks a HashMap to write the detail rows; this returns them sorted by sample.
pub fn accumulate(inputs: &[Input]) -> Result<(Vec<DetailMetrics>, SummaryMetrics), String> {
    let mut by_sample: BTreeMap<String, DetailMetrics> = BTreeMap::new();
    let mut summary = SummaryMetrics::default();
    for input in inputs {
        let mut total_het_depth = 0;
        for row in &input.detail {
            let mut row = row.clone();
            row.from_derived_fields();
            total_het_depth += row.total_het_depth;
            by_sample
                .entry(row.sample_alias.clone())
                .or_default()
                .merge(&row);
        }
        if input.summary.len() != 1 {
            return Err(wrong_summary_row_count_message(input.summary.len()));
        }
        let mut row = input.summary[0].clone();
        row.from_derived_fields(total_het_depth);
        summary.merge(&row);
    }
    let mut detail: Vec<DetailMetrics> = by_sample.into_values().collect();
    for row in &mut detail {
        row.derived_fields();
    }
    summary.derived_fields();
    Ok((detail, summary))
}
