//! Two summaries of a run: `CollectIlluminaBasecallingMetrics` and `CollectHiSeqXPfFailMetrics`.
//!
//! The first counts a lane's clusters per barcode, splitting what passed the filter from the rest.
//! The second asks why a cluster failed, and classifies the failures rather than counting them
//! again: the classes partition the filter file's own answer.
//!
//! Ported from `picard.illumina.CollectIlluminaBasecallingMetrics` and
//! `picard.illumina.quality.CollectHiSeqXPfFailMetrics` in Picard 3.4.0.

use crate::illumina_basecalls::Cluster;
use crate::illumina_files::{Segment, SegmentKind};

/// One row of the basecalling metrics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BasecallingMetrics {
    pub lane: i32,
    pub barcode: String,
    pub total_bases: i64,
    pub pf_bases: i64,
    pub total_reads: i64,
    pub pf_reads: i64,
    pub total_clusters: i64,
    pub pf_clusters: i64,
}

/// How many of a cluster's cycles are BASES, which the read structure decides.
///
/// A barcode segment's cycles are not bases and neither are a skip's, so the same lane is sixteen
/// bases under `4T` and eight under `2T2B`.
pub fn base_cycles(structure: &[Segment]) -> usize {
    structure
        .iter()
        .filter(|segment| segment.kind == SegmentKind::Template)
        .map(|segment| segment.cycles)
        .sum()
}

/// The rows a lane makes: one per barcode in the order given, then the lane's own.
///
/// The lane row is the total over every cluster whatever its barcode, which is why a run with two
/// barcodes writes three rows and their base counts do not add up to the third: the barcode rows
/// count only the template cycles they were asked about.
pub fn basecalling_metrics(
    lane: i32,
    clusters: &[(Cluster, Option<String>)],
    structure: &[Segment],
    barcodes: &[String],
) -> Vec<BasecallingMetrics> {
    let bases = base_cycles(structure) as i64;
    let reads = structure
        .iter()
        .filter(|segment| segment.kind == SegmentKind::Template)
        .count() as i64;
    let mut rows: Vec<BasecallingMetrics> = barcodes
        .iter()
        .map(|barcode| BasecallingMetrics {
            lane,
            barcode: barcode.clone(),
            ..BasecallingMetrics::default()
        })
        .collect();
    rows.push(BasecallingMetrics {
        lane,
        barcode: String::new(),
        ..BasecallingMetrics::default()
    });
    let whole = rows.len() - 1;
    for (cluster, barcode) in clusters {
        let mut touched = vec![whole];
        if let Some(barcode) = barcode {
            if let Some(index) = barcodes.iter().position(|declared| declared == barcode) {
                touched.push(index);
            }
        }
        for row in touched {
            rows[row].total_clusters += 1;
            rows[row].total_reads += reads;
            rows[row].total_bases += bases;
            if cluster.passed_filter {
                rows[row].pf_clusters += 1;
                rows[row].pf_reads += reads;
                rows[row].pf_bases += bases;
            }
        }
    }
    rows
}

/// `CollectHiSeqXPfFailMetrics`' classification of a cluster that failed the filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PfFailure {
    /// No signal at all: the read is mostly no-calls.
    Empty,
    /// Signal, but from more than one molecule: the bases are called and poor.
    Polyclonal,
    /// Signal that could not be placed.
    Misaligned,
    /// Failed for none of the above.
    Unknown,
}

/// The number of cycles the tool judges by, which is NOT `--N_CYCLES`.
///
/// The read structure is a final field initialised from the argument's DEFAULT, and a field
/// initialiser runs before the parser assigns anything, so every run looks at twenty-four cycles
/// whatever the command line said. A run directory with fewer is refused outright.
pub const CYCLES_JUDGED: usize = 24;

/// `classifyCluster`: a cluster with too few called bases is empty, one with called bases and poor
/// qualities is polyclonal, and the rest are misaligned.
pub fn classify(cluster: &Cluster, cycles: usize) -> PfFailure {
    let looked_at = cluster.calls.iter().take(cycles);
    let no_calls = looked_at.clone().filter(|call| call.base == b'.').count();
    let above_two = looked_at.filter(|call| call.quality > 2).count();
    if no_calls * 2 > cycles {
        PfFailure::Empty
    } else if above_two >= cycles.min(cluster.calls.len()) {
        PfFailure::Polyclonal
    } else {
        PfFailure::Misaligned
    }
}

/// One tile's row of the PF-fail metrics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PfFailMetrics {
    pub tile: String,
    pub reads: i64,
    pub pf_fail_reads: i64,
    pub pf_fail_empty: i64,
    pub pf_fail_polyclonal: i64,
    pub pf_fail_misaligned: i64,
    pub pf_fail_unknown: i64,
}

/// The classes over one tile, which sum to the count the filter file gives.
pub fn pf_fail_metrics(tile: &str, clusters: &[Cluster], cycles: usize) -> PfFailMetrics {
    let mut row = PfFailMetrics {
        tile: tile.to_string(),
        ..PfFailMetrics::default()
    };
    for cluster in clusters {
        row.reads += 1;
        if cluster.passed_filter {
            continue;
        }
        row.pf_fail_reads += 1;
        match classify(cluster, cycles) {
            PfFailure::Empty => row.pf_fail_empty += 1,
            PfFailure::Polyclonal => row.pf_fail_polyclonal += 1,
            PfFailure::Misaligned => row.pf_fail_misaligned += 1,
            PfFailure::Unknown => row.pf_fail_unknown += 1,
        }
    }
    row
}
