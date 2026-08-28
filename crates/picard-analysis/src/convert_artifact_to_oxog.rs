//! `ConvertSequencingArtifactToOxoG`: two artifact tables rewritten as the older OxoG one.
//!
//! There is no sequence data involved: the tool is arithmetic over two tables, so all of it is
//! here bar the metrics reader and writer.
//!
//! Ported from `picard.analysis.artifacts.ConvertSequencingArtifactToOxoG` in Picard 3.4.0.

use std::collections::BTreeMap;

/// `ConvertSequencingArtifactToOxoG.OXOG_METRICS_EXT`.
pub const OXOG_METRICS_EXT: &str = ".oxog_metrics";
/// `SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT`.
pub const PRE_ADAPTER_DETAILS_EXT: &str = ".pre_adapter_detail_metrics";
/// `SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT`.
pub const BAIT_BIAS_DETAILS_EXT: &str = ".bait_bias_detail_metrics";

/// `customCommandLineValidation`, where neither a basename nor the file it would derive is named.
pub const NO_PRE_ADAPTER_MESSAGE: &str = "Must specify either INPUT_BASE or PRE_ADAPTER_IN";
pub const NO_BAIT_BIAS_MESSAGE: &str = "Must specify either INPUT_BASE or BAIT_BIAS_IN";
pub const NO_OXOG_OUT_MESSAGE: &str = "Must specify either OUTPUT_BASE or OXOG_OUT";

/// One row of the pre-adapter detail table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreAdapterDetail {
    pub sample_alias: String,
    pub library: String,
    pub reference_base: char,
    pub alternate_base: char,
    pub context: String,
    pub pro_ref_bases: i64,
    pub pro_alt_bases: i64,
    pub con_ref_bases: i64,
    pub con_alt_bases: i64,
}

/// One row of the bait-bias detail table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaitBiasDetail {
    pub library: String,
    pub reference_base: char,
    pub alternate_base: char,
    pub context: String,
    pub forward_ref_bases: i64,
    pub forward_alt_bases: i64,
    pub reverse_ref_bases: i64,
    pub reverse_alt_bases: i64,
}

/// `isOxoG`: only the C>A and G>T transitions are read, every other row being ignored.
pub fn is_oxo_g(reference_base: char, alternate_base: char) -> bool {
    matches!((reference_base, alternate_base), ('C', 'A') | ('G', 'T'))
}

/// The reverse complement of a context, which is how the pre-adapter row for an output context is
/// found: OxoG reverse-complements its contexts, so a row for `ACA` reads the input's `TGT`.
pub fn reverse_complement(context: &str) -> String {
    context
        .chars()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// One row of the OxoG file.
#[derive(Debug, Clone, PartialEq)]
pub struct CpcgMetrics {
    pub sample_alias: String,
    pub library: String,
    pub context: String,
    /// Always nought: the input does not carry it.
    pub total_sites: i64,
    pub total_bases: i64,
    pub ref_total_bases: i64,
    pub ref_nonoxo_bases: i64,
    pub ref_oxo_bases: i64,
    pub alt_nonoxo_bases: i64,
    pub alt_oxo_bases: i64,
    pub oxidation_error_rate: f64,
    pub oxidation_q: f64,
    pub c_ref_ref_bases: i64,
    pub g_ref_ref_bases: i64,
    pub c_ref_alt_bases: i64,
    pub g_ref_alt_bases: i64,
    pub c_ref_oxo_error_rate: f64,
    pub g_ref_oxo_error_rate: f64,
    pub c_ref_oxo_q: f64,
    pub g_ref_oxo_q: f64,
}

/// The contexts the output reports: those whose pre-adapter row has a reference base of `C`.
pub fn oxog_contexts(pre_adapter: &[PreAdapterDetail]) -> Vec<String> {
    let mut contexts: Vec<String> = pre_adapter
        .iter()
        .filter(|row| row.reference_base == 'C')
        .map(|row| row.context.clone())
        .collect();
    contexts.sort();
    contexts.dedup();
    contexts
}

/// The libraries the output reports, which is every library the pre-adapter table names.
pub fn oxog_libraries(pre_adapter: &[PreAdapterDetail]) -> Vec<String> {
    let mut libraries: Vec<String> = pre_adapter.iter().map(|row| row.library.clone()).collect();
    libraries.sort();
    libraries.dedup();
    libraries
}

/// `doWork`'s conversion, for every library and every reported context.
///
/// The pre-adapter figures come from the REVERSE COMPLEMENT context and the bait-bias figures from
/// the context itself, so one output row draws on two different input rows.
///
/// The two floors are different, and that is the part worth reading twice. The oxidation rate is
/// floored at ONE BASE, so a context with fewer oxidised alternates than unoxidised ones reports
/// `1 / TOTAL_BASES` rather than a negative rate. The two bait-bias rates are floored at `1e-10`,
/// and they are opposite differences of the same two numbers, so at most one is above that floor
/// and the other reads as a Q of exactly one hundred.
///
/// The reference emits the rows in the iteration order of two HashSets. This returns them sorted,
/// by library and then by context, because that order is a hash and not a behaviour.
pub fn convert(
    pre_adapter: &[PreAdapterDetail],
    bait_bias: &[BaitBiasDetail],
) -> Result<Vec<CpcgMetrics>, String> {
    let sample_alias = pre_adapter
        .first()
        .ok_or("the pre-adapter table is empty")?
        .sample_alias
        .clone();
    let mut pre_by_library: BTreeMap<(&str, &str), &PreAdapterDetail> = BTreeMap::new();
    for row in pre_adapter {
        if is_oxo_g(row.reference_base, row.alternate_base) {
            pre_by_library.insert((&row.library, &row.context), row);
        }
    }
    let mut bait_by_library: BTreeMap<(&str, &str), &BaitBiasDetail> = BTreeMap::new();
    for row in bait_bias {
        if is_oxo_g(row.reference_base, row.alternate_base) {
            bait_by_library.insert((&row.library, &row.context), row);
        }
    }
    let mut rows = Vec::new();
    for library in oxog_libraries(pre_adapter) {
        for context in oxog_contexts(pre_adapter) {
            let complement = reverse_complement(&context);
            let pre = pre_by_library
                .get(&(library.as_str(), complement.as_str()))
                .ok_or_else(|| format!("no pre-adapter row for {library} {complement}"))?;
            let bait = bait_by_library
                .get(&(library.as_str(), context.as_str()))
                .ok_or_else(|| format!("no bait-bias row for {library} {context}"))?;

            let total_bases =
                pre.pro_ref_bases + pre.pro_alt_bases + pre.con_ref_bases + pre.con_alt_bases;
            let oxidation_error_rate =
                (pre.pro_alt_bases - pre.con_alt_bases).max(1) as f64 / total_bases as f64;
            let c_rate = pre_rate(bait.forward_alt_bases, bait.forward_ref_bases);
            let g_rate = pre_rate(bait.reverse_alt_bases, bait.reverse_ref_bases);
            let c_ref_oxo_error_rate = (c_rate - g_rate).max(1e-10);
            let g_ref_oxo_error_rate = (g_rate - c_rate).max(1e-10);
            rows.push(CpcgMetrics {
                sample_alias: sample_alias.clone(),
                library: library.clone(),
                context: context.clone(),
                total_sites: 0,
                total_bases,
                ref_total_bases: pre.pro_ref_bases + pre.con_ref_bases,
                ref_nonoxo_bases: pre.con_ref_bases,
                ref_oxo_bases: pre.pro_ref_bases,
                alt_nonoxo_bases: pre.con_alt_bases,
                alt_oxo_bases: pre.pro_alt_bases,
                oxidation_error_rate,
                oxidation_q: -10.0 * oxidation_error_rate.log10(),
                c_ref_ref_bases: bait.forward_ref_bases,
                g_ref_ref_bases: bait.reverse_ref_bases,
                c_ref_alt_bases: bait.forward_alt_bases,
                g_ref_alt_bases: bait.reverse_alt_bases,
                c_ref_oxo_error_rate,
                g_ref_oxo_error_rate,
                c_ref_oxo_q: -10.0 * c_ref_oxo_error_rate.log10(),
                g_ref_oxo_q: -10.0 * g_ref_oxo_error_rate.log10(),
            });
        }
    }
    Ok(rows)
}

/// One side's error rate: the alternates over the alternates and the references together.
fn pre_rate(alt: i64, reference: i64) -> f64 {
    alt as f64 / (alt + reference) as f64
}

/// `customCommandLineValidation`'s file names, each derived from a basename by a fixed extension.
///
/// `--OUTPUT_BASE` defaults to `--INPUT_BASE`, so naming the input alone names all three files.
pub fn derived_names(
    input_base: Option<&str>,
    output_base: Option<&str>,
) -> Result<(String, String, String), Vec<String>> {
    let output_base = output_base.or(input_base);
    let mut errors = Vec::new();
    if input_base.is_none() {
        errors.push(NO_PRE_ADAPTER_MESSAGE.to_string());
        errors.push(NO_BAIT_BIAS_MESSAGE.to_string());
    }
    if output_base.is_none() {
        errors.push(NO_OXOG_OUT_MESSAGE.to_string());
    }
    if !errors.is_empty() {
        return Err(errors);
    }
    let input_base = input_base.expect("checked");
    Ok((
        format!("{input_base}{PRE_ADAPTER_DETAILS_EXT}"),
        format!("{input_base}{BAIT_BIAS_DETAILS_EXT}"),
        format!("{}{OXOG_METRICS_EXT}", output_base.expect("checked")),
    ))
}
