//! `CollectSequencingArtifactMetrics`: which read reaches which artifact counter, and the rates
//! the counters give.
//!
//! Walking the alignments is not ported. What is ported is the split a base is filed under, the
//! two ways the four counters are folded together, the derived rates, and the five files the
//! output prefix stands for.
//!
//! Ported from `picard.analysis.artifacts.CollectSequencingArtifactMetrics`,
//! `picard.analysis.artifacts.ContextAccumulator`,
//! `picard.analysis.artifacts.SequencingArtifactMetrics` and
//! `htsjdk.samtools.util.QualityUtil` in Picard 3.4.0.

/// The floor under every rate the two detail files report, which is what keeps their Q finite.
pub const MIN_ERROR: f64 = 1e-10;

/// The five extensions the `--OUTPUT` prefix stands for, in the order the tool assigns them.
pub const PRE_ADAPTER_SUMMARY_EXT: &str = ".pre_adapter_summary_metrics";
pub const PRE_ADAPTER_DETAILS_EXT: &str = ".pre_adapter_detail_metrics";
pub const BAIT_BIAS_SUMMARY_EXT: &str = ".bait_bias_summary_metrics";
pub const BAIT_BIAS_DETAILS_EXT: &str = ".bait_bias_detail_metrics";
pub const ERROR_SUMMARY_EXT: &str = ".error_summary_metrics";

/// `setup`: the five names, with `--FILE_EXTENSION` appended to each rather than replacing any.
pub fn file_names(prefix: &str, extension: Option<&str>) -> Vec<String> {
    let suffix = extension.unwrap_or("");
    [
        PRE_ADAPTER_SUMMARY_EXT,
        PRE_ADAPTER_DETAILS_EXT,
        BAIT_BIAS_SUMMARY_EXT,
        BAIT_BIAS_DETAILS_EXT,
        ERROR_SUMMARY_EXT,
    ]
    .iter()
    .map(|ext| format!("{prefix}{ext}{suffix}"))
    .collect()
}

/// The four counters one context and one called base keep, split by end and by strand.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Alignment {
    pub r1_pos: i64,
    pub r1_neg: i64,
    pub r2_pos: i64,
    pub r2_neg: i64,
}

impl Alignment {
    /// `AlignmentAccumulator.countRecord`: an unpaired read counts as read ONE.
    pub fn count(&mut self, negative_strand: bool, paired: bool, second_of_pair: bool) {
        let read_two = paired && second_of_pair;
        match (read_two, negative_strand) {
            (true, true) => self.r2_neg += 1,
            (true, false) => self.r2_pos += 1,
            (false, true) => self.r1_neg += 1,
            (false, false) => self.r1_pos += 1,
        }
    }
}

/// The four numbers a pre-adapter detail row holds.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PreAdapterCounts {
    pub pro_ref: i64,
    pub pro_alt: i64,
    pub con_ref: i64,
    pub con_alt: i64,
}

/// The four numbers a bait-bias detail row holds.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BaitBiasCounts {
    pub fwd_ref: i64,
    pub fwd_alt: i64,
    pub rev_ref: i64,
    pub rev_alt: i64,
}

/// A pre-adapter artifact is one whose direction follows the READ.
///
/// The propitious and contrary sides are the two ways an end and a strand can agree, so read one
/// on the forward strand and read two on the reverse strand are the same side. `--TANDEM_READS`
/// says the two ends were sequenced from the SAME strand, which swaps read two's half of each
/// sum: the same file then answers the other way round.
pub fn pre_adapter(
    forward_reference: Alignment,
    forward_alternate: Alignment,
    reverse_reference: Alignment,
    reverse_alternate: Alignment,
    tandem: bool,
) -> PreAdapterCounts {
    if tandem {
        PreAdapterCounts {
            pro_ref: forward_reference.r1_pos
                + forward_reference.r2_pos
                + reverse_reference.r1_neg
                + reverse_reference.r2_neg,
            pro_alt: forward_alternate.r1_pos
                + forward_alternate.r2_pos
                + reverse_alternate.r1_neg
                + reverse_alternate.r2_neg,
            con_ref: forward_reference.r1_neg
                + forward_reference.r2_neg
                + reverse_reference.r1_pos
                + reverse_reference.r2_pos,
            con_alt: forward_alternate.r1_neg
                + forward_alternate.r2_neg
                + reverse_alternate.r1_pos
                + reverse_alternate.r2_pos,
        }
    } else {
        PreAdapterCounts {
            pro_ref: forward_reference.r1_pos
                + forward_reference.r2_neg
                + reverse_reference.r1_neg
                + reverse_reference.r2_pos,
            pro_alt: forward_alternate.r1_pos
                + forward_alternate.r2_neg
                + reverse_alternate.r1_neg
                + reverse_alternate.r2_pos,
            con_ref: forward_reference.r1_neg
                + forward_reference.r2_pos
                + reverse_reference.r1_pos
                + reverse_reference.r2_neg,
            con_alt: forward_alternate.r1_neg
                + forward_alternate.r2_pos
                + reverse_alternate.r1_pos
                + reverse_alternate.r2_neg,
        }
    }
}

/// A bait-bias artifact is one whose direction follows the REFERENCE STRAND, so the end and the
/// read's own strand are both summed away and `--TANDEM_READS` cannot reach it.
pub fn bait_bias(
    forward_reference: Alignment,
    forward_alternate: Alignment,
    reverse_reference: Alignment,
    reverse_alternate: Alignment,
) -> BaitBiasCounts {
    let total = |a: Alignment| a.r1_pos + a.r1_neg + a.r2_pos + a.r2_neg;
    BaitBiasCounts {
        fwd_ref: total(forward_reference),
        fwd_alt: total(forward_alternate),
        rev_ref: total(reverse_reference),
        rev_alt: total(reverse_alternate),
    }
}

/// `PreAdapterDetailMetrics.calculateDerivedStatistics`.
///
/// The contrary count is subtracted from the propitious one, on the argument that damage from
/// other causes falls evenly on the two sides, and a row nothing was seen for keeps the floor
/// rather than dividing by nought.
pub fn pre_adapter_error_rate(counts: &PreAdapterCounts) -> f64 {
    let total = counts.pro_ref + counts.pro_alt + counts.con_ref + counts.con_alt;
    if total == 0 {
        return MIN_ERROR;
    }
    let raw = (counts.pro_alt - counts.con_alt) as f64 / total as f64;
    raw.max(MIN_ERROR)
}

/// `BaitBiasDetailMetrics.calculateDerivedStatistics`: each strand's rate is floored on its own,
/// and then their difference is floored again.
pub fn bait_bias_error_rates(counts: &BaitBiasCounts) -> (f64, f64, f64) {
    let forward = counts.fwd_ref + counts.fwd_alt;
    let reverse = counts.rev_ref + counts.rev_alt;
    let forward_rate = if forward == 0 {
        MIN_ERROR
    } else {
        (counts.fwd_alt as f64 / forward as f64).max(MIN_ERROR)
    };
    let reverse_rate = if reverse == 0 {
        MIN_ERROR
    } else {
        (counts.rev_alt as f64 / reverse as f64).max(MIN_ERROR)
    };
    (
        forward_rate,
        reverse_rate,
        (forward_rate - reverse_rate).max(MIN_ERROR),
    )
}

/// `QualityUtil.getPhredScoreFromErrorProbability`, whose answer is an INTEGER: the Q columns of
/// both detail files are rounded, so a rate of 1e-10 reports exactly a hundred.
pub fn phred_from_error_probability(probability: f64) -> i32 {
    htsjdk_bam::quality_util::phred_score_from_error_probability(probability)
}

/// The transitions a detail file holds a row for: every reference base against every other base.
pub fn transitions() -> Vec<(u8, u8)> {
    let bases = *b"ACGT";
    let mut out = Vec::with_capacity(12);
    for reference in bases {
        for alternate in bases {
            if reference != alternate {
                out.push((reference, alternate));
            }
        }
    }
    out
}

/// How many detail rows one library's file holds: a row per transition per context.
pub fn detail_rows(context_size: usize) -> usize {
    transitions().len() * 4usize.pow(2 * context_size as u32)
}

// ---------------------------------------------------------------------------------------------
// The walk's accumulators: `ContextAccumulator`, `ArtifactCounter` and the metric beans.
// ---------------------------------------------------------------------------------------------

use std::collections::{BTreeMap, HashMap, HashSet};

use htsjdk_metrics::file::{MetricBean, Value};

/// `Transition.baseIndexMap`: A, C, G, T are 0..4, anything else has no index.
pub fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

const BASES: [u8; 4] = *b"ACGT";

/// `Transition.ALT_VALUES`, the twelve transitions a summary file holds a row for, in order.
pub const ALT_TRANSITIONS: [(u8, u8); 12] = [
    (b'A', b'C'),
    (b'A', b'G'),
    (b'A', b'T'),
    (b'C', b'A'),
    (b'C', b'G'),
    (b'C', b'T'),
    (b'G', b'A'),
    (b'G', b'C'),
    (b'G', b'T'),
    (b'T', b'A'),
    (b'T', b'C'),
    (b'T', b'G'),
];

/// `SequenceUtil.complement`, which leaves anything but a base (an `N` padding here) alone.
fn complement(base: u8) -> u8 {
    htsjdk_bam::sequence::complement(base)
}

fn reverse_complement(context: &str) -> String {
    context
        .bytes()
        .rev()
        .map(|b| complement(b) as char)
        .collect()
}

/// `PreAdapterDetailMetrics`.
#[derive(Debug, Clone)]
pub struct PreAdapterDetail {
    pub sample_alias: String,
    pub library: String,
    pub ref_base: u8,
    pub alt_base: u8,
    pub context: String,
    pub counts: PreAdapterCounts,
    pub error_rate: f64,
    pub qscore: f64,
}

/// `BaitBiasDetailMetrics`.
#[derive(Debug, Clone)]
pub struct BaitBiasDetail {
    pub sample_alias: String,
    pub library: String,
    pub ref_base: u8,
    pub alt_base: u8,
    pub context: String,
    pub counts: BaitBiasCounts,
    pub fwd_error_rate: f64,
    pub rev_error_rate: f64,
    pub error_rate: f64,
    pub qscore: f64,
}

/// `PreAdapterSummaryMetrics` and `BaitBiasSummaryMetrics`, which have the same fields.
#[derive(Debug, Clone)]
pub struct Summary {
    pub bait_bias: bool,
    pub sample_alias: String,
    pub library: String,
    pub ref_base: u8,
    pub alt_base: u8,
    pub total_qscore: f64,
    pub worst_cxt: String,
    pub worst_cxt_qscore: f64,
    pub worst_pre_cxt: String,
    pub worst_pre_cxt_qscore: f64,
    pub worst_post_cxt: String,
    pub worst_post_cxt_qscore: f64,
}

impl Summary {
    /// `inferArtifactName`, which names a different pair of transitions in each file.
    pub fn artifact_name(&self) -> &'static str {
        match (self.bait_bias, self.ref_base, self.alt_base) {
            (false, b'G', b'T') => "OxoG",
            (false, b'C', b'T') => "Deamination",
            (true, b'G', b'T') => "Gref",
            (true, b'C', b'A') => "Cref",
            _ => "NA",
        }
    }
}

fn ch(base: u8) -> Value {
    Value::Str((base as char).to_string())
}

impl MetricBean for PreAdapterDetail {
    fn class_name(&self) -> &str {
        "picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "SAMPLE_ALIAS",
            "LIBRARY",
            "REF_BASE",
            "ALT_BASE",
            "CONTEXT",
            "PRO_REF_BASES",
            "PRO_ALT_BASES",
            "CON_REF_BASES",
            "CON_ALT_BASES",
            "ERROR_RATE",
            "QSCORE",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.sample_alias.clone()),
            Value::Str(self.library.clone()),
            ch(self.ref_base),
            ch(self.alt_base),
            Value::Str(self.context.clone()),
            Value::Long(self.counts.pro_ref),
            Value::Long(self.counts.pro_alt),
            Value::Long(self.counts.con_ref),
            Value::Long(self.counts.con_alt),
            Value::Double(self.error_rate),
            Value::Double(self.qscore),
        ]
    }
}

impl MetricBean for BaitBiasDetail {
    fn class_name(&self) -> &str {
        "picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "SAMPLE_ALIAS",
            "LIBRARY",
            "REF_BASE",
            "ALT_BASE",
            "CONTEXT",
            "FWD_CXT_REF_BASES",
            "FWD_CXT_ALT_BASES",
            "REV_CXT_REF_BASES",
            "REV_CXT_ALT_BASES",
            "FWD_ERROR_RATE",
            "REV_ERROR_RATE",
            "ERROR_RATE",
            "QSCORE",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.sample_alias.clone()),
            Value::Str(self.library.clone()),
            ch(self.ref_base),
            ch(self.alt_base),
            Value::Str(self.context.clone()),
            Value::Long(self.counts.fwd_ref),
            Value::Long(self.counts.fwd_alt),
            Value::Long(self.counts.rev_ref),
            Value::Long(self.counts.rev_alt),
            Value::Double(self.fwd_error_rate),
            Value::Double(self.rev_error_rate),
            Value::Double(self.error_rate),
            Value::Double(self.qscore),
        ]
    }
}

impl MetricBean for Summary {
    fn class_name(&self) -> &str {
        if self.bait_bias {
            "picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasSummaryMetrics"
        } else {
            "picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterSummaryMetrics"
        }
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "SAMPLE_ALIAS",
            "LIBRARY",
            "REF_BASE",
            "ALT_BASE",
            "TOTAL_QSCORE",
            "WORST_CXT",
            "WORST_CXT_QSCORE",
            "WORST_PRE_CXT",
            "WORST_PRE_CXT_QSCORE",
            "WORST_POST_CXT",
            "WORST_POST_CXT_QSCORE",
            "ARTIFACT_NAME",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.sample_alias.clone()),
            Value::Str(self.library.clone()),
            ch(self.ref_base),
            ch(self.alt_base),
            Value::Double(self.total_qscore),
            Value::Str(self.worst_cxt.clone()),
            Value::Double(self.worst_cxt_qscore),
            Value::Str(self.worst_pre_cxt.clone()),
            Value::Double(self.worst_pre_cxt_qscore),
            Value::Str(self.worst_post_cxt.clone()),
            Value::Double(self.worst_post_cxt_qscore),
            Value::Str(self.artifact_name().to_string()),
        ]
    }
}

/// `ErrorSummaryMetrics`.
#[derive(Debug, Clone)]
pub struct ErrorSummary {
    pub ref_base: u8,
    pub alt_base: u8,
    pub substitution: String,
    pub ref_count: i64,
    pub alt_count: i64,
}

impl ErrorSummary {
    /// `calculateDerivedFields`.
    pub fn substitution_rate(&self) -> f64 {
        let total = (self.ref_count + self.alt_count) as f64;
        if total == 0.0 {
            0.0
        } else {
            self.alt_count as f64 / total
        }
    }
}

impl MetricBean for ErrorSummary {
    fn class_name(&self) -> &str {
        "picard.analysis.artifacts.ErrorSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "REF_BASE",
            "ALT_BASE",
            "SUBSTITUTION",
            "REF_COUNT",
            "ALT_COUNT",
            "SUBSTITUTION_RATE",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            ch(self.ref_base),
            ch(self.alt_base),
            Value::Str(self.substitution.clone()),
            Value::Long(self.ref_count),
            Value::Long(self.alt_count),
            Value::Double(self.substitution_rate()),
        ]
    }
}

/// `compareTo` on either detail class: Q, then the bases, then the context.
fn worse(q: f64, r: u8, a: u8, c: &str, than: (f64, u8, u8, &str)) -> bool {
    match q.partial_cmp(&than.0) {
        Some(std::cmp::Ordering::Less) => return true,
        Some(std::cmp::Ordering::Greater) => return false,
        _ => {}
    }
    (r, a, c) < (than.1, than.2, than.3)
}

type DetailPair = (PreAdapterDetail, BaitBiasDetail);

/// `ContextAccumulator`: the four counters for every called base of every context.
pub struct ContextAccumulator {
    tandem: bool,
    map: BTreeMap<String, [Alignment; 4]>,
}

impl ContextAccumulator {
    pub fn new<'a>(contexts: impl IntoIterator<Item = &'a String>, tandem: bool) -> Self {
        ContextAccumulator {
            tandem,
            map: contexts
                .into_iter()
                .map(|c| (c.clone(), [Alignment::default(); 4]))
                .collect(),
        }
    }

    /// `countRecord`: a context the accumulator does not know is ignored.
    pub fn count(&mut self, context: &str, called: u8, negative: bool, paired: bool, second: bool) {
        if let (Some(accumulators), Some(i)) = (self.map.get_mut(context), base_index(called)) {
            accumulators[i].count(negative, paired, second);
        }
    }

    fn merge(into: &mut [Alignment; 4], from: &[Alignment; 4]) {
        for (a, b) in into.iter_mut().zip(from) {
            a.r1_pos += b.r1_pos;
            a.r1_neg += b.r1_neg;
            a.r2_pos += b.r2_pos;
            a.r2_neg += b.r2_neg;
        }
    }

    /// `fillHalfRecords`. With a context size of nought the leading and trailing keys are the
    /// same key, and Picard merges into it twice; so does this.
    pub fn fill_half(&mut self, full: &ContextAccumulator, context_size: usize) {
        let padding = "N".repeat(context_size);
        for (key, accumulators) in &full.map {
            let central = &key[context_size..context_size + 1];
            let leading = format!("{}{central}{padding}", &key[..context_size]);
            let trailing = format!("{padding}{central}{}", &key[context_size + 1..]);
            Self::merge(self.map.get_mut(&trailing).expect("trailing"), accumulators);
            Self::merge(self.map.get_mut(&leading).expect("leading"), accumulators);
        }
    }

    /// `fillZeroRecords`.
    pub fn fill_zero(&mut self, full: &ContextAccumulator, context_size: usize) {
        let padding = "N".repeat(context_size);
        for (key, accumulators) in &full.map {
            let central = &key[context_size..context_size + 1];
            let zero = format!("{padding}{central}{padding}");
            Self::merge(self.map.get_mut(&zero).expect("zero"), accumulators);
        }
    }

    /// `calculateMetrics`: a detail pair per context (in sorted order) per called base.
    pub fn calculate(&self, sample: &str, library: &str) -> HashMap<(u8, u8), Vec<DetailPair>> {
        let mut out: HashMap<(u8, u8), Vec<DetailPair>> = HashMap::new();
        for (context, accumulators) in &self.map {
            let ref_base = context.as_bytes()[context.len() / 2];
            let reverse = &self.map[&reverse_complement(context)];
            for alt_base in BASES {
                let r = base_index(ref_base).expect("ref");
                let a = base_index(alt_base).expect("alt");
                let rc = base_index(complement(ref_base)).expect("ref");
                let ac = base_index(complement(alt_base)).expect("alt");
                let (fr, fa, rr, ra) = (accumulators[r], accumulators[a], reverse[rc], reverse[ac]);
                let pre_counts = pre_adapter(fr, fa, rr, ra, self.tandem);
                let bait_counts = bait_bias(fr, fa, rr, ra);
                let error_rate = pre_adapter_error_rate(&pre_counts);
                let (fwd, rev, bait_rate) = bait_bias_error_rates(&bait_counts);
                let pre = PreAdapterDetail {
                    sample_alias: sample.to_string(),
                    library: library.to_string(),
                    ref_base,
                    alt_base,
                    context: context.clone(),
                    counts: pre_counts,
                    error_rate,
                    qscore: phred_from_error_probability(error_rate) as f64,
                };
                let bait = BaitBiasDetail {
                    sample_alias: sample.to_string(),
                    library: library.to_string(),
                    ref_base,
                    alt_base,
                    context: context.clone(),
                    counts: bait_counts,
                    fwd_error_rate: fwd,
                    rev_error_rate: rev,
                    error_rate: bait_rate,
                    qscore: phred_from_error_probability(bait_rate) as f64,
                };
                out.entry((ref_base, alt_base))
                    .or_default()
                    .push((pre, bait));
            }
        }
        out
    }
}

/// `ArtifactCounter`: one library's accumulators, and the four lists its `finish` builds.
pub struct ArtifactCounter {
    sample_alias: String,
    library: String,
    context_size: usize,
    full: ContextAccumulator,
    half: ContextAccumulator,
    zero: ContextAccumulator,
    leading: HashSet<String>,
    trailing: HashSet<String>,
}

/// What `ArtifactCounter.finish` leaves behind, in the order its getters hand it out.
#[derive(Debug, Default)]
pub struct CounterMetrics {
    pub pre_adapter_summary: Vec<Summary>,
    pub pre_adapter_detail: Vec<PreAdapterDetail>,
    pub bait_bias_summary: Vec<Summary>,
    pub bait_bias_detail: Vec<BaitBiasDetail>,
}

impl ArtifactCounter {
    pub fn new(sample_alias: &str, library: &str, context_size: usize, tandem: bool) -> Self {
        let length = 2 * context_size + 1;
        let mut full: Vec<String> = vec![String::new()];
        for _ in 0..length {
            full = full
                .iter()
                .flat_map(|prefix| BASES.iter().map(move |&b| format!("{prefix}{}", b as char)))
                .collect();
        }
        let padding = "N".repeat(context_size);
        let mut leading = HashSet::new();
        let mut trailing = HashSet::new();
        let mut zero = HashSet::new();
        for context in &full {
            let central = &context[context_size..context_size + 1];
            leading.insert(format!("{}{central}{padding}", &context[..context_size]));
            trailing.insert(format!(
                "{padding}{central}{}",
                &context[context_size + 1..]
            ));
            zero.insert(format!("{padding}{central}{padding}"));
        }
        let half: HashSet<String> = leading.union(&trailing).cloned().collect();
        ArtifactCounter {
            sample_alias: sample_alias.to_string(),
            library: library.to_string(),
            context_size,
            full: ContextAccumulator::new(&full, tandem),
            half: ContextAccumulator::new(&half, tandem),
            zero: ContextAccumulator::new(&zero, tandem),
            leading,
            trailing,
        }
    }

    pub fn count(&mut self, context: &str, called: u8, negative: bool, paired: bool, second: bool) {
        self.full.count(context, called, negative, paired, second);
    }

    /// `getWorstMetrics`: the first strictly-worst of each kind, chosen independently.
    fn worst(pairs: &[&DetailPair]) -> (PreAdapterDetail, BaitBiasDetail) {
        let mut pre: Option<&PreAdapterDetail> = None;
        let mut bait: Option<&BaitBiasDetail> = None;
        for (p, b) in pairs.iter().map(|pair| (&pair.0, &pair.1)) {
            if pre.is_none_or(|w| {
                worse(
                    p.qscore,
                    p.ref_base,
                    p.alt_base,
                    &p.context,
                    (w.qscore, w.ref_base, w.alt_base, &w.context),
                )
            }) {
                pre = Some(p);
            }
            if bait.is_none_or(|w| {
                worse(
                    b.qscore,
                    b.ref_base,
                    b.alt_base,
                    &b.context,
                    (w.qscore, w.ref_base, w.alt_base, &w.context),
                )
            }) {
                bait = Some(b);
            }
        }
        (pre.expect("worst").clone(), bait.expect("worst").clone())
    }

    /// `finish`: the detail metrics of every alternate transition, and a summary row for each.
    pub fn finish(mut self) -> CounterMetrics {
        let details = self.full.calculate(&self.sample_alias, &self.library);
        // getSummaryMetrics recomputes the full metrics; they are the same rows.
        self.half.fill_half(&self.full, self.context_size);
        let half = self.half.calculate(&self.sample_alias, &self.library);
        self.zero.fill_zero(&self.full, self.context_size);
        let zero = self.zero.calculate(&self.sample_alias, &self.library);

        let mut out = CounterMetrics::default();
        for transition in ALT_TRANSITIONS {
            let full_rows = &details[&transition];
            let zero_rows = &zero[&transition];
            if zero_rows.len() != 1 {
                crate::metrics_cli::thrown(&format!(
                    "picard.PicardException: Should have exactly one context-free metric pair for transition: {}>{}",
                    transition.0 as char, transition.1 as char
                ));
            }
            let mut leading = Vec::new();
            let mut trailing = Vec::new();
            for pair in &half[&transition] {
                if self.leading.contains(&pair.0.context) {
                    leading.push(pair);
                }
                if self.trailing.contains(&pair.0.context) {
                    trailing.push(pair);
                }
            }
            let total = &zero_rows[0];
            let worst_full = Self::worst(&full_rows.iter().collect::<Vec<_>>());
            let worst_leading = Self::worst(&leading);
            let worst_trailing = Self::worst(&trailing);
            out.pre_adapter_summary.push(Summary {
                bait_bias: false,
                sample_alias: self.sample_alias.clone(),
                library: self.library.clone(),
                ref_base: transition.0,
                alt_base: transition.1,
                total_qscore: total.0.qscore,
                worst_cxt: worst_full.0.context.clone(),
                worst_cxt_qscore: worst_full.0.qscore,
                worst_pre_cxt: worst_leading.0.context.clone(),
                worst_pre_cxt_qscore: worst_leading.0.qscore,
                worst_post_cxt: worst_trailing.0.context.clone(),
                worst_post_cxt_qscore: worst_trailing.0.qscore,
            });
            out.bait_bias_summary.push(Summary {
                bait_bias: true,
                sample_alias: self.sample_alias.clone(),
                library: self.library.clone(),
                ref_base: transition.0,
                alt_base: transition.1,
                total_qscore: total.1.qscore,
                worst_cxt: worst_full.1.context.clone(),
                worst_cxt_qscore: worst_full.1.qscore,
                worst_pre_cxt: worst_leading.1.context.clone(),
                worst_pre_cxt_qscore: worst_leading.1.qscore,
                worst_post_cxt: worst_trailing.1.context.clone(),
                worst_post_cxt_qscore: worst_trailing.1.qscore,
            });
            for (pre, bait) in full_rows {
                out.pre_adapter_detail.push(pre.clone());
                out.bait_bias_detail.push(bait.clone());
            }
        }
        out
    }
}

/// `CollectSequencingArtifactMetrics.finish`'s last step: the pre-adapter details folded onto
/// the six substitutions read from the reference's A or C, in `TreeSet` order of `"R>A"`.
pub fn error_summaries(details: &[PreAdapterDetail]) -> Vec<ErrorSummary> {
    let mut by_error: BTreeMap<String, ErrorSummary> = BTreeMap::new();
    for m in details {
        let (r, a) = if m.ref_base == b'G' || m.ref_base == b'T' {
            (complement(m.ref_base), complement(m.alt_base))
        } else {
            (m.ref_base, m.alt_base)
        };
        let key = format!("{}>{}", r as char, a as char);
        let entry = by_error.entry(key.clone()).or_insert(ErrorSummary {
            ref_base: r,
            alt_base: a,
            substitution: key,
            ref_count: 0,
            alt_count: 0,
        });
        entry.ref_count += m.counts.pro_ref + m.counts.con_ref;
        entry.alt_count += m.counts.pro_alt + m.counts.con_alt;
    }
    by_error.into_values().collect()
}
