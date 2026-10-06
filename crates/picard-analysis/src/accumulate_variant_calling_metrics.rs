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

/// The tool as it runs: every column of both tables, read from the files of each prefix, merged
/// the way `MergeableMetricBase.merge` merges, and written back.
///
/// [`accumulate`] above is the merge's arithmetic on the columns the conformance suite prints;
/// this carries all of them (the indel, multiallelic and complex-indel columns, the GQ0 counts),
/// because a metrics file that lost a column is not the reference's file.
pub mod tool {
    use htsjdk_metrics::file::{MetricBean, Value};

    use crate::java_hash_map::JavaHashMap;

    /// `VariantCallingSummaryMetrics`, public columns in declaration order, then the ten hidden
    /// fields the merge adds and the derivation reads.
    #[derive(Debug, Clone, Default)]
    pub struct Summary {
        pub total_snps: i64,
        pub num_in_db_snp: i64,
        pub novel_snps: i64,
        pub filtered_snps: i64,
        pub pct_dbsnp: f32,
        pub dbsnp_titv: f64,
        pub novel_titv: f64,
        pub total_indels: i64,
        pub novel_indels: i64,
        pub filtered_indels: i64,
        pub pct_dbsnp_indels: f32,
        pub num_in_db_snp_indels: i64,
        pub dbsnp_ins_del_ratio: f64,
        pub novel_ins_del_ratio: f64,
        pub total_multiallelic_snps: f64,
        pub num_in_db_snp_multiallelic: f64,
        pub total_complex_indels: f64,
        pub num_in_db_snp_complex_indels: f64,
        pub snp_reference_bias: f64,
        pub num_singletons: i64,
        ref_allele_obs: i64,
        alt_allele_obs: i64,
        novel_deletions: i64,
        novel_insertions: i64,
        novel_transitions: i64,
        novel_transversions: i64,
        db_snp_deletions: i64,
        db_snp_insertions: i64,
        db_snp_transitions: i64,
        db_snp_transversions: i64,
    }

    /// `VariantCallingDetailMetrics`: the sample, three derived columns and the GQ0 count, then
    /// the summary's columns.
    #[derive(Debug, Clone, Default)]
    pub struct Detail {
        pub sample_alias: Option<String>,
        pub het_homvar_ratio: f64,
        pub pct_gq0_variants: f64,
        pub total_gq0_variants: i64,
        pub total_het_depth: i64,
        pub summary: Summary,
        num_hets: i64,
        num_hom_var: i64,
    }

    /// `invertFromRatio`: `Math.round(sum / (ratio + 1.0))`, nought for a NaN ratio.
    fn invert_from_ratio(sum: i64, ratio: f64) -> i64 {
        if ratio.is_nan() {
            0
        } else {
            jmath::math::round(sum as f64 / (ratio + 1.0))
        }
    }

    impl Summary {
        /// `calculateDerivedFields`. A ratio keeps its value (zero) unless its denominator is
        /// positive; the others are plain divisions that put `?` in the file when they fail.
        pub fn calculate_derived_fields(&mut self) {
            self.pct_dbsnp = self.num_in_db_snp as f32 / self.total_snps as f32;
            self.novel_snps = self.total_snps - self.num_in_db_snp;
            self.snp_reference_bias =
                self.ref_allele_obs as f64 / (self.ref_allele_obs + self.alt_allele_obs) as f64;
            if self.db_snp_transversions > 0 {
                self.dbsnp_titv = self.db_snp_transitions as f64 / self.db_snp_transversions as f64;
            }
            if self.novel_transversions > 0 {
                self.novel_titv = self.novel_transitions as f64 / self.novel_transversions as f64;
            }
            self.pct_dbsnp_indels = self.num_in_db_snp_indels as f32 / self.total_indels as f32;
            self.novel_indels = self.total_indels - self.num_in_db_snp_indels;
            if self.db_snp_deletions > 0 {
                self.dbsnp_ins_del_ratio =
                    self.db_snp_insertions as f64 / self.db_snp_deletions as f64;
            }
            if self.novel_deletions > 0 {
                self.novel_ins_del_ratio =
                    self.novel_insertions as f64 / self.novel_deletions as f64;
            }
        }

        /// `calculateFromDerivedFields(totalHetDepth)`: the hidden fields rebuilt from what the
        /// file printed.
        pub fn calculate_from_derived_fields(&mut self, total_het_depth: i64) {
            self.db_snp_transversions = invert_from_ratio(self.num_in_db_snp, self.dbsnp_titv);
            self.db_snp_transitions = self.num_in_db_snp - self.db_snp_transversions;
            self.novel_transversions = invert_from_ratio(self.novel_snps, self.novel_titv);
            self.novel_transitions = self.novel_snps - self.novel_transversions;
            self.db_snp_deletions =
                invert_from_ratio(self.num_in_db_snp_indels, self.dbsnp_ins_del_ratio);
            self.db_snp_insertions = self.num_in_db_snp_indels - self.db_snp_deletions;
            self.novel_deletions = invert_from_ratio(self.novel_indels, self.novel_ins_del_ratio);
            self.novel_insertions = self.novel_indels - self.novel_deletions;
            self.ref_allele_obs = if self.snp_reference_bias.is_nan() {
                0
            } else {
                jmath::math::round(total_het_depth as f64 * self.snp_reference_bias)
            };
            self.alt_allele_obs = total_het_depth - self.ref_allele_obs;
        }

        /// `merge`: every `@MergeByAdding` field adds, the derived ones are left alone.
        pub fn merge(&mut self, other: &Summary) {
            self.total_snps += other.total_snps;
            self.num_in_db_snp += other.num_in_db_snp;
            self.novel_snps += other.novel_snps;
            self.filtered_snps += other.filtered_snps;
            self.total_indels += other.total_indels;
            self.novel_indels += other.novel_indels;
            self.filtered_indels += other.filtered_indels;
            self.num_in_db_snp_indels += other.num_in_db_snp_indels;
            self.total_multiallelic_snps += other.total_multiallelic_snps;
            self.num_in_db_snp_multiallelic += other.num_in_db_snp_multiallelic;
            self.total_complex_indels += other.total_complex_indels;
            self.num_in_db_snp_complex_indels += other.num_in_db_snp_complex_indels;
            self.num_singletons += other.num_singletons;
            self.ref_allele_obs += other.ref_allele_obs;
            self.alt_allele_obs += other.alt_allele_obs;
            self.novel_deletions += other.novel_deletions;
            self.novel_insertions += other.novel_insertions;
            self.novel_transitions += other.novel_transitions;
            self.novel_transversions += other.novel_transversions;
            self.db_snp_deletions += other.db_snp_deletions;
            self.db_snp_insertions += other.db_snp_insertions;
            self.db_snp_transitions += other.db_snp_transitions;
            self.db_snp_transversions += other.db_snp_transversions;
        }

        fn values(&self) -> Vec<Value> {
            vec![
                Value::Long(self.total_snps),
                Value::Long(self.num_in_db_snp),
                Value::Long(self.novel_snps),
                Value::Long(self.filtered_snps),
                Value::Double(f64::from(self.pct_dbsnp)),
                Value::Double(self.dbsnp_titv),
                Value::Double(self.novel_titv),
                Value::Long(self.total_indels),
                Value::Long(self.novel_indels),
                Value::Long(self.filtered_indels),
                Value::Double(f64::from(self.pct_dbsnp_indels)),
                Value::Long(self.num_in_db_snp_indels),
                Value::Double(self.dbsnp_ins_del_ratio),
                Value::Double(self.novel_ins_del_ratio),
                Value::Double(self.total_multiallelic_snps),
                Value::Double(self.num_in_db_snp_multiallelic),
                Value::Double(self.total_complex_indels),
                Value::Double(self.num_in_db_snp_complex_indels),
                Value::Double(self.snp_reference_bias),
                Value::Long(self.num_singletons),
            ]
        }
    }

    impl Detail {
        /// `calculateFromDerivedFields()`.
        pub fn calculate_from_derived_fields(&mut self) {
            self.num_hom_var = invert_from_ratio(self.summary.total_snps, self.het_homvar_ratio);
            self.num_hets = self.summary.total_snps - self.num_hom_var;
            self.summary
                .calculate_from_derived_fields(self.total_het_depth);
        }

        /// `calculateDerivedFields`: the summary's, then the detail's own three.
        pub fn calculate_derived_fields(&mut self) {
            self.summary.calculate_derived_fields();
            self.het_homvar_ratio = self.num_hets as f64 / self.num_hom_var as f64;
            self.pct_gq0_variants =
                self.total_gq0_variants as f64 / (self.num_hets + self.num_hom_var) as f64;
            self.total_het_depth = self.summary.ref_allele_obs + self.summary.alt_allele_obs;
        }

        /// `merge`: the sample alias asserts equal (a null one takes the other's), the GQ0 count
        /// and the two hidden counts add, and so does everything the summary adds.
        pub fn merge(&mut self, other: &Detail) -> Result<(), String> {
            self.summary.merge(&other.summary);
            self.total_gq0_variants += other.total_gq0_variants;
            self.num_hets += other.num_hets;
            self.num_hom_var += other.num_hom_var;
            match (&self.sample_alias, &other.sample_alias) {
                (None, _) => self.sample_alias = other.sample_alias.clone(),
                (Some(a), Some(b)) if a != b => {
                    return Err(format!(
                        "java.lang.IllegalStateException: Field SAMPLE_ALIAS is annotated as \
                         @MergeByAssertEquals, but found two different values: {a} & {b}"
                    ))
                }
                _ => {}
            }
            Ok(())
        }
    }

    const SUMMARY_COLUMNS: [&str; 20] = [
        "TOTAL_SNPS",
        "NUM_IN_DB_SNP",
        "NOVEL_SNPS",
        "FILTERED_SNPS",
        "PCT_DBSNP",
        "DBSNP_TITV",
        "NOVEL_TITV",
        "TOTAL_INDELS",
        "NOVEL_INDELS",
        "FILTERED_INDELS",
        "PCT_DBSNP_INDELS",
        "NUM_IN_DB_SNP_INDELS",
        "DBSNP_INS_DEL_RATIO",
        "NOVEL_INS_DEL_RATIO",
        "TOTAL_MULTIALLELIC_SNPS",
        "NUM_IN_DB_SNP_MULTIALLELIC",
        "TOTAL_COMPLEX_INDELS",
        "NUM_IN_DB_SNP_COMPLEX_INDELS",
        "SNP_REFERENCE_BIAS",
        "NUM_SINGLETONS",
    ];

    const DETAIL_COLUMNS: [&str; 25] = [
        "SAMPLE_ALIAS",
        "HET_HOMVAR_RATIO",
        "PCT_GQ0_VARIANTS",
        "TOTAL_GQ0_VARIANTS",
        "TOTAL_HET_DEPTH",
        "TOTAL_SNPS",
        "NUM_IN_DB_SNP",
        "NOVEL_SNPS",
        "FILTERED_SNPS",
        "PCT_DBSNP",
        "DBSNP_TITV",
        "NOVEL_TITV",
        "TOTAL_INDELS",
        "NOVEL_INDELS",
        "FILTERED_INDELS",
        "PCT_DBSNP_INDELS",
        "NUM_IN_DB_SNP_INDELS",
        "DBSNP_INS_DEL_RATIO",
        "NOVEL_INS_DEL_RATIO",
        "TOTAL_MULTIALLELIC_SNPS",
        "NUM_IN_DB_SNP_MULTIALLELIC",
        "TOTAL_COMPLEX_INDELS",
        "NUM_IN_DB_SNP_COMPLEX_INDELS",
        "SNP_REFERENCE_BIAS",
        "NUM_SINGLETONS",
    ];

    impl MetricBean for Summary {
        fn class_name(&self) -> &str {
            "picard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics"
        }
        fn columns(&self) -> &[&'static str] {
            &SUMMARY_COLUMNS
        }
        fn values(&self) -> Vec<Value> {
            Summary::values(self)
        }
    }

    impl MetricBean for Detail {
        fn class_name(&self) -> &str {
            "picard.vcf.CollectVariantCallingMetrics$VariantCallingDetailMetrics"
        }
        fn columns(&self) -> &[&'static str] {
            &DETAIL_COLUMNS
        }
        fn values(&self) -> Vec<Value> {
            let mut out = vec![
                self.sample_alias
                    .clone()
                    .map(Value::Str)
                    .unwrap_or(Value::Null),
                Value::Double(self.het_homvar_ratio),
                Value::Double(self.pct_gq0_variants),
                Value::Long(self.total_gq0_variants),
                Value::Long(self.total_het_depth),
            ];
            out.extend(self.summary.values());
            out
        }
    }

    /// The rows of the first metrics table in a metrics file, each as column name to text:
    /// `MetricsFile.read` up to the blank line that ends the table.
    pub fn read_metrics_table(text: &str) -> Vec<Vec<(String, String)>> {
        let mut lines = text.lines();
        for line in lines.by_ref() {
            if line.starts_with("## METRICS CLASS") {
                break;
            }
        }
        let columns: Vec<&str> = match lines.next() {
            Some(header) => header.split('\t').collect(),
            None => return Vec::new(),
        };
        let mut rows = Vec::new();
        for line in lines {
            if line.is_empty() {
                break;
            }
            let cells: Vec<&str> = line.split('\t').collect();
            rows.push(
                columns
                    .iter()
                    .enumerate()
                    .map(|(i, name)| {
                        (
                            name.to_string(),
                            cells.get(i).map(|c| c.to_string()).unwrap_or_default(),
                        )
                    })
                    .collect(),
            );
        }
        rows
    }

    /// `FormatUtil.parseDouble`: `?` and `-?` are NaN.
    fn parse_double(text: &str) -> Result<f64, String> {
        if text == "?" || text == "-?" {
            return Ok(f64::NAN);
        }
        text.parse::<f64>()
            .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{text}\""))
    }

    fn parse_float(text: &str) -> Result<f32, String> {
        if text == "?" || text == "-?" {
            return Ok(f32::NAN);
        }
        text.parse::<f32>()
            .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{text}\""))
    }

    fn parse_long(text: &str) -> Result<i64, String> {
        text.parse::<i64>()
            .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{text}\""))
    }

    fn summary_from(row: &[(String, String)]) -> Result<Summary, String> {
        let mut summary = Summary::default();
        for (name, text) in row {
            match name.as_str() {
                "TOTAL_SNPS" => summary.total_snps = parse_long(text)?,
                "NUM_IN_DB_SNP" => summary.num_in_db_snp = parse_long(text)?,
                "NOVEL_SNPS" => summary.novel_snps = parse_long(text)?,
                "FILTERED_SNPS" => summary.filtered_snps = parse_long(text)?,
                "PCT_DBSNP" => summary.pct_dbsnp = parse_float(text)?,
                "DBSNP_TITV" => summary.dbsnp_titv = parse_double(text)?,
                "NOVEL_TITV" => summary.novel_titv = parse_double(text)?,
                "TOTAL_INDELS" => summary.total_indels = parse_long(text)?,
                "NOVEL_INDELS" => summary.novel_indels = parse_long(text)?,
                "FILTERED_INDELS" => summary.filtered_indels = parse_long(text)?,
                "PCT_DBSNP_INDELS" => summary.pct_dbsnp_indels = parse_float(text)?,
                "NUM_IN_DB_SNP_INDELS" => summary.num_in_db_snp_indels = parse_long(text)?,
                "DBSNP_INS_DEL_RATIO" => summary.dbsnp_ins_del_ratio = parse_double(text)?,
                "NOVEL_INS_DEL_RATIO" => summary.novel_ins_del_ratio = parse_double(text)?,
                "TOTAL_MULTIALLELIC_SNPS" => summary.total_multiallelic_snps = parse_double(text)?,
                "NUM_IN_DB_SNP_MULTIALLELIC" => {
                    summary.num_in_db_snp_multiallelic = parse_double(text)?
                }
                "TOTAL_COMPLEX_INDELS" => summary.total_complex_indels = parse_double(text)?,
                "NUM_IN_DB_SNP_COMPLEX_INDELS" => {
                    summary.num_in_db_snp_complex_indels = parse_double(text)?
                }
                "SNP_REFERENCE_BIAS" => summary.snp_reference_bias = parse_double(text)?,
                "NUM_SINGLETONS" => summary.num_singletons = parse_long(text)?,
                _ => {}
            }
        }
        Ok(summary)
    }

    /// One detail row as `MetricsFile.read` builds the bean.
    pub fn detail_from(row: &[(String, String)]) -> Result<Detail, String> {
        let mut detail = Detail {
            summary: summary_from(row)?,
            ..Default::default()
        };
        for (name, text) in row {
            match name.as_str() {
                "SAMPLE_ALIAS" => detail.sample_alias = Some(text.clone()),
                "HET_HOMVAR_RATIO" => detail.het_homvar_ratio = parse_double(text)?,
                "PCT_GQ0_VARIANTS" => detail.pct_gq0_variants = parse_double(text)?,
                "TOTAL_GQ0_VARIANTS" => detail.total_gq0_variants = parse_long(text)?,
                "TOTAL_HET_DEPTH" => detail.total_het_depth = parse_long(text)?,
                _ => {}
            }
        }
        Ok(detail)
    }

    /// One input's two tables, as text.
    pub struct InputText {
        pub detail: String,
        pub summary: String,
    }

    /// `doWork` after the files are read: each input's detail rows merged into the one per
    /// sample (in the order a `HashMap` of them is walked), and its summary row into the total.
    pub fn run(inputs: &[InputText]) -> Result<(Vec<Detail>, Summary), String> {
        // The reference's `HashMap<String, Detail>` is filled by `computeIfAbsent` and walked in
        // bucket order, which the map of positions below reproduces, head insertion included; the details themselves are kept in a vector so that merging
        // into one never moves it.
        let mut order: JavaHashMap<usize> = JavaHashMap::new();
        let mut index_of: std::collections::HashMap<String, usize> = Default::default();
        let mut samples: Vec<Detail> = Vec::new();
        let mut collapsed = Summary::default();
        for input in inputs {
            let mut total_het_depth = 0i64;
            for row in read_metrics_table(&input.detail) {
                let mut detailed = detail_from(&row)?;
                detailed.calculate_from_derived_fields();
                total_het_depth += detailed.total_het_depth;
                let key = detailed.sample_alias.clone().unwrap_or_default();
                let at = *index_of.entry(key.clone()).or_insert_with(|| {
                    samples.push(Detail::default());
                    order.insert_front_if_absent(&key, samples.len() - 1);
                    samples.len() - 1
                });
                samples[at].merge(&detailed)?;
            }
            let rows = read_metrics_table(&input.summary);
            if rows.len() != 1 {
                return Err(format!(
                    "picard.PicardException: {}",
                    super::wrong_summary_row_count_message(rows.len())
                ));
            }
            let mut summary = summary_from(&rows[0])?;
            summary.calculate_from_derived_fields(total_het_depth);
            collapsed.merge(&summary);
        }
        let mut details: Vec<Detail> = order.iter().map(|(_, &at)| samples[at].clone()).collect();
        for detail in &mut details {
            detail.calculate_derived_fields();
        }
        collapsed.calculate_derived_fields();
        Ok((details, collapsed))
    }
}
