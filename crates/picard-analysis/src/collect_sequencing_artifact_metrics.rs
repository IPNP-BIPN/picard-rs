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
