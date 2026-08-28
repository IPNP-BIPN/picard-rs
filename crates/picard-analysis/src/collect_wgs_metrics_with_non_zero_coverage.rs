//! `CollectWgsMetricsWithNonZeroCoverage`: one traversal reported twice.
//!
//! The second row is not a second walk. It is the same depth histogram with the depth-zero bin set
//! to zero, put back through the same arithmetic, which is why a fully covered reference makes the
//! two rows identical and why the exclusion columns do not move.
//!
//! The traversal itself is [`crate::collect_wgs_metrics`]. What is ported here is what the second
//! pass changes.
//!
//! Ported from `picard.analysis.CollectWgsMetricsWithNonZeroCoverage` in Picard 3.4.0.

/// The `CATEGORY` column, which the parent's own metrics have no trace of.
pub const WHOLE_GENOME: &str = "WHOLE_GENOME";
pub const NON_ZERO_REGIONS: &str = "NON_ZERO_REGIONS";

/// The two histogram columns, which are one table and not two sections.
pub const WHOLE_GENOME_COLUMN: &str = "count_WHOLE_GENOME";
pub const NON_ZERO_COLUMN: &str = "count_NON_ZERO_REGIONS";

/// `addToMetricsFile`: `highQualityDepthHistogramArray[0] = 0` before the second pass.
///
/// The uncovered loci do not leave the traversal, they leave the histogram, which is the whole of
/// the difference between the two rows.
pub fn drop_the_zero_bin(histogram: &[u64]) -> Vec<u64> {
    let mut bins = histogram.to_vec();
    if let Some(first) = bins.first_mut() {
        *first = 0;
    }
    bins
}

/// `getDepthHistogramNonZero`, which builds its own histogram from bin ONE upwards.
///
/// It never increments bin zero, so that bin has no entry rather than an entry of nought. The
/// metrics file prints the two histograms as one table keyed by depth, so the missing entry comes
/// out as a `0` and the two roads meet.
pub fn non_zero_column(histogram: &[u64]) -> Vec<u64> {
    drop_the_zero_bin(histogram)
}

/// `GENOME_TERRITORY`: the loci the histogram holds, which is its counts and not its bins.
pub fn territory(histogram: &[u64]) -> u64 {
    histogram.iter().sum()
}

/// `MEAN_COVERAGE`, over the territory the row has and not over the covered bases.
///
/// A second row whose territory is nought divides by nought, and the writer renders the NaN as
/// `?` rather than as a number.
pub fn mean_coverage(histogram: &[u64]) -> f64 {
    let total: u64 = histogram
        .iter()
        .enumerate()
        .map(|(depth, count)| depth as u64 * count)
        .sum();
    total as f64 / territory(histogram) as f64
}

/// The two rows of one traversal, in the order they are written.
pub fn rows(histogram: &[u64]) -> [(&'static str, Vec<u64>); 2] {
    [
        (WHOLE_GENOME, histogram.to_vec()),
        (NON_ZERO_REGIONS, drop_the_zero_bin(histogram)),
    ]
}
