//! `CollectHsMetrics`: a hybrid-selection experiment counted against its baits and its targets.
//!
//! Walking the alignments is not ported. What is ported is the partition a base falls into, the
//! arithmetic of the derived columns, and the per-target row.
//!
//! Ported from `picard.analysis.directed.CollectHsMetrics`,
//! `picard.analysis.directed.TargetMetricsCollector` and
//! `picard.analysis.directed.HsMetricCollector` in Picard 3.4.0.

/// `TargetedPcrMetricsCollector.NEAR_PROBE_DISTANCE_DEFAULT`, the window either side of a bait.
pub const DEFAULT_NEAR_DISTANCE: i32 = 250;
/// The mapping-quality floor, which is ONE and not nought: a read at quality nought never
/// contributes coverage, whatever else is asked for.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 1;
/// The base-quality floor, which IS nought: this tool counts bases `CollectWgsMetrics` would not.
pub const DEFAULT_MINIMUM_BASE_QUALITY: i32 = 0;

/// Where one aligned base falls with respect to the baits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaitPlacement {
    /// Inside a bait.
    On,
    /// Outside every bait but within `near_distance` of one.
    Near,
    /// Neither.
    Off,
}

/// The partition, which is the three columns that sum to `PF_BASES_ALIGNED`.
///
/// The near window is a distance either side of every bait, so a base sixty away from a bait's end
/// is near at the default and off at a distance of nought.
pub fn placement(position: i32, baits: &[(i32, i32)], near_distance: i32) -> BaitPlacement {
    if baits
        .iter()
        .any(|(start, end)| position >= *start && position <= *end)
    {
        return BaitPlacement::On;
    }
    if near_distance > 0
        && baits.iter().any(|(start, end)| {
            position >= start - near_distance && position <= end + near_distance
        })
    {
        return BaitPlacement::Near;
    }
    BaitPlacement::Off
}

/// The columns derived from the three bait counts and the two territories.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Derived {
    pub pct_selected_bases: f64,
    pub pct_off_bait: f64,
    pub on_bait_vs_selected: f64,
    pub mean_bait_coverage: f64,
    pub mean_target_coverage: f64,
    pub fold_enrichment: f64,
    pub pct_usable_bases_on_bait: f64,
}

/// The counts one run produced, which the derived columns are read off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Counts {
    pub pf_bases: i64,
    pub pf_bases_aligned: i64,
    pub on_bait: i64,
    pub near_bait: i64,
    pub off_bait: i64,
    pub on_target: i64,
    pub bait_territory: i64,
    pub target_territory: i64,
    /// The reference's own length, which only the enrichment column reads.
    pub genome_size: i64,
}

/// `TargetMetricsCollector.calculateDerivedMetrics`, on the counts a traversal produced.
///
/// `PCT_SELECTED_BASES` counts the NEAR bases as selected, so it is not `ON_BAIT / ALIGNED`, and
/// `ON_BAIT_VS_SELECTED` is the ratio of the two: a file whose reads are mostly near rather than on
/// its baits reads 0.75 and 0.333333 on the same run.
pub fn derived(counts: &Counts) -> Derived {
    let selected = counts.on_bait + counts.near_bait;
    let aligned = counts.on_bait + counts.near_bait + counts.off_bait;
    Derived {
        pct_selected_bases: selected as f64 / aligned as f64,
        pct_off_bait: counts.off_bait as f64 / aligned as f64,
        on_bait_vs_selected: counts.on_bait as f64 / selected as f64,
        mean_bait_coverage: counts.on_bait as f64 / counts.bait_territory as f64,
        mean_target_coverage: counts.on_target as f64 / counts.target_territory as f64,
        // `FOLD_ENRICHMENT` divides by the three-column denominator and not by PF_BASES_ALIGNED,
        // and its second half is the bait territory over the whole GENOME: a fifth of a thousand
        // bases baited, a quarter of the aligned bases on bait, and the answer is 1.25.
        fold_enrichment: (counts.on_bait as f64 / aligned as f64)
            / (counts.bait_territory as f64 / counts.genome_size as f64),
        pct_usable_bases_on_bait: counts.on_bait as f64 / counts.pf_bases as f64,
    }
}

/// One row of the `--PER_TARGET_COVERAGE` file.
#[derive(Debug, Clone, PartialEq)]
pub struct TargetRow {
    pub name: String,
    pub length: i64,
    pub gc: f64,
    pub mean_coverage: f64,
    pub normalized_coverage: f64,
    pub pct_zero_coverage: f64,
    pub minimum: i64,
    pub maximum: i64,
    pub read_count: i64,
}

/// The row a target's own coverage produces, given the run's mean over every target.
///
/// `normalized_coverage` is this target's mean over that one, so a target covered like the average
/// reads one and an uncovered one reads nought. `pct_0x` is the fraction of the target's bases no
/// read reached, which is why an uncovered target reads one there and a partly covered one reads
/// the rest.
pub fn target_row(
    name: &str,
    coverage: &[i64],
    bases: &[u8],
    read_count: i64,
    run_mean: f64,
) -> TargetRow {
    let length = coverage.len() as i64;
    let total: i64 = coverage.iter().sum();
    let mean = total as f64 / length as f64;
    let zero = coverage.iter().filter(|depth| **depth == 0).count() as f64;
    let gc = bases
        .iter()
        .filter(|base| matches!(base.to_ascii_uppercase(), b'G' | b'C'))
        .count() as f64
        / bases.len() as f64;
    TargetRow {
        name: name.to_string(),
        length,
        gc,
        mean_coverage: mean,
        normalized_coverage: mean / run_mean,
        pct_zero_coverage: zero / length as f64,
        minimum: *coverage.iter().min().unwrap_or(&0),
        maximum: *coverage.iter().max().unwrap_or(&0),
        read_count,
    }
}
