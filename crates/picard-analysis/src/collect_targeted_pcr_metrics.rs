//! `CollectTargetedPcrMetrics`: [`crate::collect_hs_metrics`]' collector under the amplicon column
//! names.
//!
//! Nothing about the counting is different. What is different is one line in a constructor and the
//! names the columns are written under, and both are here so a caller can ask which tool it is
//! reproducing.
//!
//! Ported from `picard.analysis.directed.CollectTargetedPcrMetrics`,
//! `picard.analysis.directed.TargetedPcrMetricsCollector` and
//! `picard.analysis.directed.CollectTargetedMetrics` in Picard 3.4.0.

use crate::collect_hs_metrics::{BaitPlacement, Counts, Derived};

/// `CollectTargetedMetrics.CLIP_OVERLAPPING_READS`, which both tools inherit.
///
/// `CollectHsMetrics`'s constructor sets it to TRUE and this tool leaves it alone, which is the
/// whole of the difference between the two tools' target columns: the overlap of a pair is counted
/// once there and twice here.
pub const SHARED_CLIP_OVERLAPPING_READS_DEFAULT: bool = false;
pub const HS_METRICS_CLIP_OVERLAPPING_READS_DEFAULT: bool = true;
pub const PCR_METRICS_CLIP_OVERLAPPING_READS_DEFAULT: bool = SHARED_CLIP_OVERLAPPING_READS_DEFAULT;

/// The amplicon columns, in the order the metrics file writes them.
pub const AMPLICON_COLUMNS: [&str; 8] = [
    "CUSTOM_AMPLICON_SET",
    "AMPLICON_TERRITORY",
    "ON_AMPLICON_BASES",
    "NEAR_AMPLICON_BASES",
    "OFF_AMPLICON_BASES",
    "PCT_AMPLIFIED_BASES",
    "PCT_OFF_AMPLICON",
    "ON_AMPLICON_VS_SELECTED",
];

/// The bait column each amplicon column stands in for, which is what says the two tools write one
/// collector's numbers under two vocabularies.
pub fn bait_column(amplicon_column: &str) -> Option<&'static str> {
    Some(match amplicon_column {
        "CUSTOM_AMPLICON_SET" => "BAIT_SET",
        "AMPLICON_TERRITORY" => "BAIT_TERRITORY",
        "ON_AMPLICON_BASES" => "ON_BAIT_BASES",
        "NEAR_AMPLICON_BASES" => "NEAR_BAIT_BASES",
        "OFF_AMPLICON_BASES" => "OFF_BAIT_BASES",
        "PCT_AMPLIFIED_BASES" => "PCT_SELECTED_BASES",
        "PCT_OFF_AMPLICON" => "PCT_OFF_BAIT",
        "ON_AMPLICON_VS_SELECTED" => "ON_BAIT_VS_SELECTED",
        _ => return None,
    })
}

/// The placement rule, which is the other tool's: an amplicon is a bait by another name.
pub fn placement(position: i32, amplicons: &[(i32, i32)], near_distance: i32) -> BaitPlacement {
    crate::collect_hs_metrics::placement(position, amplicons, near_distance)
}

/// The derived columns, which are the other tool's arithmetic over the same counts.
pub fn derived(counts: &Counts) -> Derived {
    crate::collect_hs_metrics::derived(counts)
}
