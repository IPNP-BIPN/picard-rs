//! `CollectJumpingLibraryMetrics`: the pairs that point outwards, the pairs that do not, and the
//! ones the tool calls chimeric.
//!
//! What is ported is which pair lands in which bucket and the arithmetic on top of those counts:
//! the orientation test, the three chimera kinds in the order they are tried, the histogram trim,
//! and the library-size estimate.
//!
//! Ported from `picard.analysis.CollectJumpingLibraryMetrics`, `picard.sam.DuplicationMetrics` and
//! `htsjdk.samtools.util.Histogram` in Picard 3.4.0.

use std::collections::BTreeMap;

/// `SAMPLE_FOR_MODE`: how many outward pairs the first pass looks at to find the mode.
pub const SAMPLE_FOR_MODE: usize = 50000;

/// The default floor below which an insert is not oversized.
pub const DEFAULT_CHIMERA_KB_MIN: i64 = 100000;

/// The default tail limit the histogram is trimmed by.
pub const DEFAULT_TAIL_LIMIT: i32 = 10000;

/// One pair, reduced to what the tool reads off its first read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pair {
    pub reference_index: i32,
    pub mate_reference_index: i32,
    pub reverse: bool,
    pub mate_reverse: bool,
    pub insert_size: i64,
    pub duplicate: bool,
    /// The `MQ` tag, which is only consulted when it is there.
    pub mate_quality: Option<i32>,
    pub mapping_quality: i32,
    pub unmapped: bool,
    pub mate_unmapped: bool,
}

/// What one pair is counted as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bucket {
    /// One end mapped and the other not.
    Fragment,
    /// Outward-facing: a jump.
    Jump,
    /// Inward-facing: not a jump.
    NonJump,
    /// Too long, tandem, or across two chromosomes.
    Chimera,
    /// Below the mapping-quality floor, or otherwise not counted at all.
    Skipped,
    /// Both ends unmapped and unplaced: the traversal STOPS here.
    Terminator,
}

/// `SamPairUtil.getPairOrientation`, for a pair whose two ends are on the same contig.
///
/// The orientation is the two STRANDS read in position order, so the sign of the insert size has
/// nothing to do with it: a reverse read whose mate is forward and to its right is outward, and
/// the same two strands the other way round are inward.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Forward then reverse.
    Fr,
    /// Reverse then forward.
    Rf,
    /// Both on the same strand.
    Tandem,
}

pub fn orientation(reverse: bool, mate_reverse: bool, mate_is_to_the_right: bool) -> Orientation {
    if reverse == mate_reverse {
        return Orientation::Tandem;
    }
    // The leftmost read's own strand decides.
    let leftmost_is_reverse = if mate_is_to_the_right {
        reverse
    } else {
        mate_reverse
    };
    if leftmost_is_reverse {
        Orientation::Rf
    } else {
        Orientation::Fr
    }
}

/// The arguments that decide which bucket a pair lands in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arguments {
    pub minimum_mapping_quality: i32,
    pub tail_limit: i32,
    pub chimera_kb_min: i64,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            minimum_mapping_quality: 0,
            tail_limit: DEFAULT_TAIL_LIMIT,
            chimera_kb_min: DEFAULT_CHIMERA_KB_MIN,
        }
    }
}

/// Whether the quality floor lets a pair through.
///
/// The MQ tag is consulted ONLY when it is present, so a pair with no tag passes on its own
/// mapping quality alone rather than being refused for having no mate quality.
pub fn passes_quality(pair: &Pair, arguments: &Arguments) -> bool {
    if let Some(mate_quality) = pair.mate_quality {
        if mate_quality < arguments.minimum_mapping_quality {
            return false;
        }
    }
    pair.mapping_quality >= arguments.minimum_mapping_quality
}

/// Which bucket one pair lands in, given the chimera threshold the first pass settled.
///
/// The order of the three chimera tests is the behaviour: an oversized insert is counted as
/// oversized even when it is also tandem, and a tandem pair is counted as tandem even when it is
/// also across two chromosomes.
pub fn classify(
    pair: &Pair,
    chimera_threshold: f64,
    mate_is_to_the_right: bool,
    arguments: &Arguments,
) -> Bucket {
    if pair.unmapped {
        if !pair.mate_unmapped {
            return Bucket::Fragment;
        }
        // Both ends unmapped: an unplaced one ends the file, a placed one is passed over.
        if pair.reference_index < 0 {
            return Bucket::Terminator;
        }
        return Bucket::Skipped;
    }
    if pair.mate_unmapped {
        return Bucket::Fragment;
    }
    if !passes_quality(pair, arguments) {
        return Bucket::Skipped;
    }
    let absolute = pair.insert_size.abs() as f64;
    if absolute > chimera_threshold {
        return Bucket::Chimera;
    }
    if pair.mate_reverse == pair.reverse {
        return Bucket::Chimera;
    }
    if pair.mate_reference_index != pair.reference_index {
        return Bucket::Chimera;
    }
    match orientation(pair.reverse, pair.mate_reverse, mate_is_to_the_right) {
        Orientation::Rf => Bucket::Jump,
        Orientation::Fr => Bucket::NonJump,
        // The tandem case was answered above, so this is unreachable in the reference too: it
        // throws an IllegalStateException there.
        Orientation::Tandem => Bucket::Chimera,
    }
}

/// The mode of the outward inserts, taken in a first pass over the same file.
///
/// The pass stops after `SAMPLE_FOR_MODE` outward pairs, and it applies the quality floor and the
/// orientation test but NOT the chimera one: the mode is what the chimera test is built from.
pub fn outward_mode(pairs: &[(Pair, bool)], arguments: &Arguments) -> f64 {
    let mut histogram: BTreeMap<i64, i64> = BTreeMap::new();
    let mut sampled = 0usize;
    for (pair, mate_is_to_the_right) in pairs {
        if sampled >= SAMPLE_FOR_MODE {
            break;
        }
        if pair.unmapped && pair.reference_index < 0 {
            break;
        }
        if pair.unmapped || pair.mate_unmapped {
            continue;
        }
        if !passes_quality(pair, arguments)
            || pair.mate_reverse == pair.reverse
            || pair.mate_reference_index != pair.reference_index
            || orientation(pair.reverse, pair.mate_reverse, *mate_is_to_the_right)
                != Orientation::Rf
        {
            continue;
        }
        *histogram.entry(pair.insert_size.abs()).or_insert(0) += 1;
        sampled += 1;
    }
    mode(&histogram)
}

/// `Histogram.getMode`, which is the id of the bin with the largest count.
///
/// A tie goes to the FIRST such bin in key order, and an empty histogram has no mode at all: the
/// reference returns zero there, which is what makes the chimera floor fall back on the argument.
pub fn mode(histogram: &BTreeMap<i64, i64>) -> f64 {
    let mut best: Option<(i64, i64)> = None;
    for (id, count) in histogram {
        match best {
            Some((_, best_count)) if *count <= best_count => {}
            _ => best = Some((*id, *count)),
        }
    }
    best.map_or(0.0, |(id, _)| id as f64)
}

/// The threshold an insert has to EXCEED to be called oversized.
pub fn chimera_threshold(mode: f64, arguments: &Arguments) -> f64 {
    mode.max(arguments.chimera_kb_min as f64)
}

/// `Histogram.trimByTailLimit`.
///
/// It keeps every bin up to and including the mode, then walks forward only while each bin
/// follows the last by EXACTLY ONE and holds at least the mode's count over the limit. Bins that
/// are further apart than one are therefore cut whatever the limit says, which on a set of
/// inserts a hundred apart leaves the mode and everything before it and nothing after.
pub fn trim_by_tail_limit(histogram: &BTreeMap<i64, i64>, tail_limit: i32) -> BTreeMap<i64, i64> {
    if histogram.is_empty() {
        return BTreeMap::new();
    }
    let mode_id = mode(histogram) as i64;
    let mode_size = *histogram.get(&mode_id).unwrap_or(&0) as f64;
    let minimum = mode_size / tail_limit as f64;
    let mut kept = BTreeMap::new();
    let mut last: Option<i64> = None;
    for (id, count) in histogram {
        if *id <= mode_id {
            kept.insert(*id, *count);
        } else if last.is_some_and(|last| last != *id - 1) || (*count as f64) < minimum {
            break;
        } else {
            kept.insert(*id, *count);
        }
        last = Some(*id);
    }
    kept
}

/// `Histogram.getMean`, over the bins' ids weighted by their counts.
pub fn histogram_mean(histogram: &BTreeMap<i64, i64>) -> f64 {
    let total: i64 = histogram.values().sum();
    if total == 0 {
        return f64::NAN;
    }
    let sum: f64 = histogram
        .iter()
        .map(|(id, count)| *id as f64 * *count as f64)
        .sum();
    sum / total as f64
}

/// `Histogram.getStandardDeviation`, which divides by the count and not by the count less one.
pub fn histogram_standard_deviation(histogram: &BTreeMap<i64, i64>) -> f64 {
    let total: i64 = histogram.values().sum();
    if total == 0 {
        return f64::NAN;
    }
    let mean = histogram_mean(histogram);
    let mut sum = 0.0;
    for (id, count) in histogram {
        let difference = *id as f64 - mean;
        sum += difference * difference * *count as f64;
    }
    (sum / total as f64).sqrt()
}

/// `DuplicationMetrics.estimateLibrarySize`, by the bisection the reference uses.
///
/// It answers nothing at all when there are no duplicates, which the metric then writes as zero
/// rather than as a blank.
pub fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    let duplicates = read_pairs - unique_read_pairs;
    if read_pairs <= 0 || duplicates <= 0 {
        return None;
    }
    let f = |x: f64, c: f64, n: f64| c / x - 1.0 + (-n / x).exp();
    let c = unique_read_pairs as f64;
    let n = read_pairs as f64;
    if unique_read_pairs >= read_pairs || f(c, c, n) < 0.0 {
        return None;
    }
    let mut m = 1.0f64;
    let mut big = 100.0f64;
    while f(big * c, c, n) > 0.0 {
        big *= 10.0;
    }
    // Forty halvings, however close the answer already is.
    for _ in 0..40 {
        let r = (m + big) / 2.0;
        let u = f(r * c, c, n);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            m = r;
        } else {
            big = r;
        }
    }
    Some((c * (m + big) / 2.0) as i64)
}

/// The metrics one run produces.
#[derive(Debug, Clone, PartialEq)]
pub struct JumpingLibraryMetrics {
    pub jump_pairs: i64,
    pub jump_duplicate_pairs: i64,
    pub jump_duplicate_pct: f64,
    pub jump_library_size: i64,
    pub jump_mean_insert_size: f64,
    pub jump_stdev_insert_size: f64,
    pub nonjump_pairs: i64,
    pub nonjump_duplicate_pairs: i64,
    pub nonjump_duplicate_pct: f64,
    pub nonjump_library_size: i64,
    pub nonjump_mean_insert_size: f64,
    pub nonjump_stdev_insert_size: f64,
    pub chimeric_pairs: i64,
    pub fragments: i64,
    pub pct_jumps: f64,
    pub pct_nonjumps: f64,
    pub pct_chimeras: f64,
}

/// The column order the metrics file writes.
pub const COLUMNS: [&str; 17] = [
    "JUMP_PAIRS",
    "JUMP_DUPLICATE_PAIRS",
    "JUMP_DUPLICATE_PCT",
    "JUMP_LIBRARY_SIZE",
    "JUMP_MEAN_INSERT_SIZE",
    "JUMP_STDEV_INSERT_SIZE",
    "NONJUMP_PAIRS",
    "NONJUMP_DUPLICATE_PAIRS",
    "NONJUMP_DUPLICATE_PCT",
    "NONJUMP_LIBRARY_SIZE",
    "NONJUMP_MEAN_INSERT_SIZE",
    "NONJUMP_STDEV_INSERT_SIZE",
    "CHIMERIC_PAIRS",
    "FRAGMENTS",
    "PCT_JUMPS",
    "PCT_NONJUMPS",
    "PCT_CHIMERAS",
];

/// One whole run: the buckets, then the arithmetic on top of them.
///
/// Every ratio answers zero when its denominator is zero rather than a NaN, which is the one
/// place the metric departs from the histogram's own conventions.
pub fn collect(pairs: &[(Pair, bool)], arguments: &Arguments) -> JumpingLibraryMetrics {
    let threshold = chimera_threshold(outward_mode(pairs, arguments), arguments);
    let mut jump_histogram: BTreeMap<i64, i64> = BTreeMap::new();
    let mut nonjump_histogram: BTreeMap<i64, i64> = BTreeMap::new();
    let (mut jumps, mut nonjumps, mut chimeras, mut fragments) = (0i64, 0i64, 0i64, 0i64);
    let (mut jump_duplicates, mut nonjump_duplicates) = (0i64, 0i64);
    for (pair, mate_is_to_the_right) in pairs {
        match classify(pair, threshold, *mate_is_to_the_right, arguments) {
            Bucket::Terminator => break,
            Bucket::Fragment => fragments += 1,
            Bucket::Chimera => chimeras += 1,
            Bucket::Jump => {
                jumps += 1;
                *jump_histogram.entry(pair.insert_size.abs()).or_insert(0) += 1;
                if pair.duplicate {
                    jump_duplicates += 1;
                }
            }
            Bucket::NonJump => {
                nonjumps += 1;
                *nonjump_histogram.entry(pair.insert_size.abs()).or_insert(0) += 1;
                if pair.duplicate {
                    nonjump_duplicates += 1;
                }
            }
            Bucket::Skipped => {}
        }
    }
    let jump_trimmed = trim_by_tail_limit(&jump_histogram, arguments.tail_limit);
    let nonjump_trimmed = trim_by_tail_limit(&nonjump_histogram, arguments.tail_limit);
    let total = (jumps + nonjumps + chimeras) as f64;
    JumpingLibraryMetrics {
        jump_pairs: jumps,
        jump_duplicate_pairs: jump_duplicates,
        jump_duplicate_pct: if jumps != 0 {
            jump_duplicates as f64 / jumps as f64
        } else {
            0.0
        },
        jump_library_size: if jumps > 0 && jump_duplicates > 0 {
            estimate_library_size(jumps, jumps - jump_duplicates).unwrap_or(0)
        } else {
            0
        },
        jump_mean_insert_size: histogram_mean(&jump_trimmed),
        jump_stdev_insert_size: histogram_standard_deviation(&jump_trimmed),
        nonjump_pairs: nonjumps,
        nonjump_duplicate_pairs: nonjump_duplicates,
        nonjump_duplicate_pct: if nonjumps != 0 {
            nonjump_duplicates as f64 / nonjumps as f64
        } else {
            0.0
        },
        nonjump_library_size: if nonjumps > 0 && nonjump_duplicates > 0 {
            estimate_library_size(nonjumps, nonjumps - nonjump_duplicates).unwrap_or(0)
        } else {
            0
        },
        nonjump_mean_insert_size: histogram_mean(&nonjump_trimmed),
        nonjump_stdev_insert_size: histogram_standard_deviation(&nonjump_trimmed),
        chimeric_pairs: chimeras,
        fragments,
        pct_jumps: if total != 0.0 {
            jumps as f64 / total
        } else {
            0.0
        },
        pct_nonjumps: if total != 0.0 {
            nonjumps as f64 / total
        } else {
            0.0
        },
        pct_chimeras: if total != 0.0 {
            chimeras as f64 / total
        } else {
            0.0
        },
    }
}

/// The refusal a file that is not coordinate-sorted produces, misspelling and all.
pub fn unsorted_message(name: &str) -> String {
    format!("SAM file must {name} must be sorted in coordintate order")
}
