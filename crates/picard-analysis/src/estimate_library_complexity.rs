//! `EstimateLibraryComplexity`: which pairs a group calls duplicates, and what the histogram gives.
//!
//! Reading the file and sorting the pairs are not ported. What is ported is the quality check that
//! admits a pair, the rule that calls two pairs duplicates, the floor that keeps a bin out of the
//! metrics, and the way the estimate is asked for.
//!
//! Ported from `picard.sam.markduplicates.EstimateLibraryComplexity`,
//! `picard.sam.markduplicates.ElcIdenticalBasesDuplicatesFinder` and
//! `picard.sam.DuplicationMetrics` in Picard 3.4.0.

use crate::jumping_library::estimate_library_size;

/// `MIN_IDENTICAL_BASES`: the prefix of BOTH ends that a group agrees on.
pub const DEFAULT_MIN_IDENTICAL_BASES: usize = 5;
/// `MAX_DIFF_RATE`: how far the rest may differ.
pub const DEFAULT_MAX_DIFF_RATE: f64 = 0.03;
/// `MIN_MEAN_QUALITY`: the floor a pair'send must clear.
pub const DEFAULT_MIN_MEAN_QUALITY: i32 = 20;
/// `MIN_GROUP_COUNT`: how many groups a bin needs before the METRICS count it.
pub const DEFAULT_MIN_GROUP_COUNT: i64 = 2;

/// `passesQualityCheck`.
///
/// The mean is an INTEGER division over the read's length, so qualities averaging nineteen and a
/// half are dropped at twenty. A read shorter than the seed fails outright, and an `N` anywhere in
/// the seed fails whatever the qualities say.
pub fn passes_quality_check(
    bases: &[u8],
    qualities: &[u8],
    seed_length: usize,
    minimum_quality: i32,
    max_read_length: usize,
) -> bool {
    if bases.len() < seed_length {
        return false;
    }
    if bases[..seed_length]
        .iter()
        .any(|base| base.eq_ignore_ascii_case(&b'N'))
    {
        return false;
    }
    let read_length = if max_read_length == 0 {
        bases.len()
    } else {
        bases.len().min(max_read_length)
    };
    let total: i32 = qualities[..read_length].iter().map(|q| i32::from(*q)).sum();
    total / read_length as i32 >= minimum_quality
}

/// Whether two pairs land in the same group: the first `seed_length` bases of BOTH ends agree.
pub fn same_group(left: (&[u8], &[u8]), right: (&[u8], &[u8]), seed_length: usize) -> bool {
    left.0.len() >= seed_length
        && right.0.len() >= seed_length
        && left.1.len() >= seed_length
        && right.1.len() >= seed_length
        && left.0[..seed_length] == right.0[..seed_length]
        && left.1[..seed_length] == right.1[..seed_length]
}

/// `ElcIdenticalBasesDuplicatesFinder.matches`.
///
/// The comparison starts AT the seed, because the grouping has already settled it, so a difference
/// inside the prefix is not an error: it is what put the two pairs in different groups. The
/// allowance is a rate over the two ends' compared lengths together, floored, and the ends are
/// compared over the SHORTER of the two.
pub fn matches(
    left: (&[u8], &[u8]),
    right: (&[u8], &[u8]),
    seed_length: usize,
    max_diff_rate: f64,
    max_read_length: usize,
) -> bool {
    let truncate = |length: usize| {
        if max_read_length == 0 {
            length
        } else {
            length.min(max_read_length)
        }
    };
    let read_one = truncate(left.0.len().min(right.0.len()));
    let read_two = truncate(left.1.len().min(right.1.len()));
    let max_errors = ((read_one + read_two) as f64 * max_diff_rate).floor() as i64;
    let mut errors = 0;
    for (a, b) in [
        (&left.0[..read_one], &right.0[..read_one]),
        (&left.1[..read_two], &right.1[..read_two]),
    ] {
        for index in seed_length..a.len() {
            if a[index] != b[index] {
                errors += 1;
                if errors > max_errors {
                    return false;
                }
            }
        }
    }
    true
}

/// The metrics one library's histogram gives.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Metrics {
    pub read_pairs_examined: i64,
    pub read_pair_duplicates: i64,
    pub read_pair_optical_duplicates: i64,
    pub percent_duplication: f64,
    pub estimated_library_size: Option<i64>,
}

/// The metrics a library's duplicate-set histogram produces.
///
/// A bin holding fewer than `min_group_count` groups is dropped HERE and nowhere else: the
/// histogram file still carries it, which is why a single duplicate pair reports nothing examined
/// beside a histogram that says there were two.
///
/// `bins` maps a duplicate-set size to how many sets of that size a library has, and `optical` to
/// how many of those were optical duplicates.
pub fn metrics(bins: &[(i64, i64, i64)], min_group_count: i64) -> Metrics {
    let mut out = Metrics::default();
    for (size, groups, optical) in bins {
        if *groups >= min_group_count {
            out.read_pairs_examined += size * groups;
            out.read_pair_duplicates += (size - 1) * groups;
            out.read_pair_optical_duplicates += optical;
        }
    }
    out.percent_duplication = if out.read_pairs_examined == 0 {
        0.0
    } else {
        (out.read_pair_duplicates * 2) as f64 / (out.read_pairs_examined * 2) as f64
    };
    // `calculateDerivedFields` takes the OPTICAL duplicates off the pair count before it estimates,
    // so a library whose every duplicate is optical has nothing left to estimate from.
    out.estimated_library_size = estimate_library_size(
        out.read_pairs_examined - out.read_pair_optical_duplicates,
        out.read_pairs_examined - out.read_pair_duplicates,
    );
    out
}
