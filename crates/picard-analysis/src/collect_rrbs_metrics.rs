//! `CollectRrbsMetrics`: which cytosines an RRBS run counts, and which of them are CpG.
//!
//! The tool walks a bisulfite-converted alignment against its reference and counts cytosines twice
//! over: the ones in a CpG context, per site, and the ones outside it, in bulk. A CpG is looked
//! for in the REFERENCE and not in the read, which is what makes a read that reads `TG` over a
//! reference `CG` a converted site rather than a mismatch.
//!
//! The rules around that do not compose the way the code reads, and three of them are carried here
//! deliberately:
//!
//!  * the last base of an alignment block is never a CpG, the loop stopping one short so the pair
//!    has a second base to be;
//!  * the comment above the non-CpG branch says those counts are held back until the read is known
//!    to carry a CpG, and the code increments them straight into the totals;
//!  * and the non-CpG branch reads its quality off the WHOLE read with the block's own index,
//!    where the CpG branch reads the block's own qualities.
//!
//! Ported from `picard.analysis.CollectRrbsMetrics` and `picard.analysis.RrbsMetricsCollector`.

/// `MINIMUM_READ_LENGTH`, under which a read is counted and dropped.
pub const DEFAULT_MINIMUM_READ_LENGTH: usize = 5;
/// `C_QUALITY_THRESHOLD`, which the cytosine's own quality must reach.
pub const DEFAULT_C_QUALITY_THRESHOLD: i32 = 20;
/// `NEXT_BASE_QUALITY_THRESHOLD`, which its neighbour's must reach, and which is a lower bar.
pub const DEFAULT_NEXT_BASE_QUALITY_THRESHOLD: i32 = 10;
/// `MAX_MISMATCH_RATE`, whose bound is rounded before a strictly greater test.
pub const DEFAULT_MAX_MISMATCH_RATE: f64 = 0.1;

/// The three files the tool writes, appended to a prefix that has gained a dot if it had none.
pub const DETAIL_FILE_EXTENSION: &str = "rrbs_detail_metrics";
pub const SUMMARY_FILE_EXTENSION: &str = "rrbs_summary_metrics";
pub const PDF_FILE_EXTENSION: &str = "rrbs_qc.pdf";

/// `METRICS_FILE_PREFIX` as `doWork` leaves it: a dot is appended when it has none.
///
/// The three names are then a plain concatenation, which is why the prefix and not the tool owns
/// the separator.
pub fn file_names(prefix: &str) -> [String; 3] {
    let prefix = if prefix.ends_with('.') {
        prefix.to_string()
    } else {
        format!("{prefix}.")
    };
    [
        format!("{prefix}{SUMMARY_FILE_EXTENSION}"),
        format!("{prefix}{DETAIL_FILE_EXTENSION}"),
        format!("{prefix}{PDF_FILE_EXTENSION}"),
    ]
}

/// Whether a read is dropped for its length, before anything else is asked of it.
pub fn is_too_short(read_length: usize, minimum_read_length: usize) -> bool {
    read_length < minimum_read_length
}

/// `SequenceUtil.countMismatches(...) > Math.round(readLength * maxMismatchRate)`.
///
/// The bound is ROUNDED to a whole number and the test is strictly greater, so a twenty-base read
/// at the default rate allows two mismatches and is dropped at three.
pub fn mismatch_bound(read_length: usize, max_mismatch_rate: f64) -> i64 {
    (read_length as f64 * max_mismatch_rate).round() as i64
}

/// Whether a read is dropped for its mismatches.
pub fn is_too_mismatched(mismatches: i64, read_length: usize, max_mismatch_rate: f64) -> bool {
    mismatches > mismatch_bound(read_length, max_mismatch_rate)
}

/// `SequenceUtil.basesEqual`, which is case insensitive and lets `N` equal nothing.
fn bases_equal(left: u8, right: u8) -> bool {
    left.eq_ignore_ascii_case(&right) && !left.eq_ignore_ascii_case(&b'N')
}

/// `SequenceUtil.bisulfiteBasesEqual(false, read, reference)`: a read `T` over a reference `C` is
/// equal, and nothing else changes.
pub fn bisulfite_bases_equal(read: u8, reference: u8) -> bool {
    if reference.eq_ignore_ascii_case(&b'C') && read.eq_ignore_ascii_case(&b'T') {
        return true;
    }
    bases_equal(read, reference)
}

/// `isBisulfiteConverted`: the reference was a C and the read is a T.
pub fn is_bisulfite_converted(read: u8, reference: u8) -> bool {
    reference.eq_ignore_ascii_case(&b'C') && read.eq_ignore_ascii_case(&b'T')
}

/// `isC`: the reference is a cytosine and the read agrees with it under bisulfite rules.
pub fn is_c(reference: u8, read: u8) -> bool {
    bases_equal(reference, b'C') && bisulfite_bases_equal(read, reference)
}

/// `isAboveCytoQcThreshold`: the base's own quality and its NEIGHBOUR's, against two thresholds.
///
/// The index bound is what stops the last base of the array from ever qualifying, which is a
/// second reason a cytosine at the end of a read is not counted.
pub fn is_above_cyto_qc_threshold(
    qualities: &[u8],
    index: usize,
    c_quality_threshold: i32,
    next_base_quality_threshold: i32,
) -> bool {
    index + 1 < qualities.len()
        && i32::from(qualities[index]) >= c_quality_threshold
        && i32::from(qualities[index + 1]) >= next_base_quality_threshold
}

/// `isValidCpg`: the C matching under bisulfite rules, the G matching EXACTLY, and the qualities.
pub fn is_valid_cpg(
    reference: &[u8],
    read: &[u8],
    qualities: &[u8],
    index: usize,
    c_quality_threshold: i32,
    next_base_quality_threshold: i32,
) -> bool {
    is_c(reference[index], read[index])
        && bases_equal(reference[index + 1], read[index + 1])
        && is_above_cyto_qc_threshold(
            qualities,
            index,
            c_quality_threshold,
            next_base_quality_threshold,
        )
}

/// `getCurRefIndex`: where a site is reported, which is not the mirror of the index on the
/// negative strand.
///
/// A negative-strand block is reverse complemented and then treated as a positive one, so the
/// index has to be turned back. The reference turns it back one base FURTHER than the mirror,
/// which is what puts the site on the C of the pair rather than on the G.
pub fn current_reference_index(
    reference_start: usize,
    block_length: usize,
    index: usize,
    negative: bool,
) -> usize {
    if negative {
        reference_start + (block_length - 1) - index - 1
    } else {
        reference_start + index
    }
}

/// What one alignment block contributes.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Counts {
    /// The sites a CpG was seen at, as zero-based reference positions.
    pub cpg_sites: Vec<usize>,
    /// The sites of those where the read had converted the cytosine.
    pub converted_cpg_sites: Vec<usize>,
    /// Cytosines outside a CpG context.
    pub non_cpg_total: i64,
    /// Those of them the read had converted.
    pub non_cpg_converted: i64,
}

/// The loop over one alignment block, in the reference's own order.
///
/// `qualities` is the block's own, and `read_qualities` is the WHOLE read's: the CpG branch reads
/// the first and the non-CpG branch reads the second with the block's index, which is a difference
/// on any read whose alignment does not start at its first base.
#[allow(clippy::too_many_arguments)]
pub fn block_counts(
    reference: &[u8],
    read: &[u8],
    qualities: &[u8],
    read_qualities: &[u8],
    reference_start: usize,
    negative: bool,
    c_quality_threshold: i32,
    next_base_quality_threshold: i32,
) -> Counts {
    let mut counts = Counts::default();
    let block_length = reference.len();
    let mut index = 0;
    // The loop stops one short, so the last base of a block is never the C of a pair.
    while index + 1 < block_length {
        let site = current_reference_index(reference_start, block_length, index, negative);
        if bases_equal(reference[index], b'C') && bases_equal(reference[index + 1], b'G') {
            // A CpG in the reference is taken out of the running whether or not it is valid, which
            // is what stops its cytosine from being counted as one outside a CpG.
            if is_valid_cpg(
                reference,
                read,
                qualities,
                index,
                c_quality_threshold,
                next_base_quality_threshold,
            ) {
                counts.cpg_sites.push(site);
                if is_bisulfite_converted(read[index], reference[index]) {
                    counts.converted_cpg_sites.push(site);
                }
            }
            index += 2;
            continue;
        }
        if is_c(reference[index], read[index])
            && is_above_cyto_qc_threshold(
                read_qualities,
                index,
                c_quality_threshold,
                next_base_quality_threshold,
            )
            && bisulfite_bases_equal(read[index + 1], reference[index + 1])
        {
            counts.non_cpg_total += 1;
            if is_bisulfite_converted(read[index], reference[index]) {
                counts.non_cpg_converted += 1;
            }
        }
        index += 1;
    }
    counts
}

/// `finish`: the two rates, which are zero rather than a division by zero when nothing was seen.
pub fn conversion_rate(converted: i64, total: i64) -> f64 {
    if total == 0 {
        0.0
    } else {
        converted as f64 / total as f64
    }
}

/// `MEAN_CPG_COVERAGE`, which is a histogram's mean bin size and therefore NaN over an empty one.
///
/// The metrics file writes a NaN as `?`, which is how an empty run is told from a run that saw
/// sites and averaged one.
pub fn mean_cpg_coverage(bins: &[i64]) -> f64 {
    if bins.is_empty() {
        return f64::NAN;
    }
    bins.iter().sum::<i64>() as f64 / bins.len() as f64
}

/// `MEDIAN_CPG_COVERAGE`, taken as an integer, which is zero over an empty histogram.
pub fn median_cpg_coverage(bins: &[i64]) -> i64 {
    if bins.is_empty() {
        return 0;
    }
    let mut sorted = bins.to_vec();
    sorted.sort_unstable();
    let middle = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[middle]
    } else {
        (sorted[middle - 1] + sorted[middle]) / 2
    }
}

/// Whether the bytes of the QC plot are a claim this port makes, which they are not.
///
/// The plot is drawn by R from the two metrics files, and R's output is not byte stable.
pub const PLOT_BYTES_ARE_REPRODUCIBLE: bool = false;
