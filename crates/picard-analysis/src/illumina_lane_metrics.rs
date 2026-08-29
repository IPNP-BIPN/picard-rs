//! `CollectIlluminaLaneMetrics`: what a lane's tile metrics say about it.
//!
//! The numbers come off the METRIC CODES of `InterOp/TileMetricsOut.bin` and off nothing else: 100
//! is the cluster count, 101 the count that passed the filter, 102 and 103 their densities, and
//! 200 and 201 the phasing and prephasing of a read, offset by twice the read DESCRIPTOR's index.
//!
//! Two rules the reference enforces and this reproduces. A file missing either the counts or the
//! densities is refused rather than reported with a gap, and a phasing PAIR has to be complete:
//! half a pair is refused by name. So the tool has no partial answer.
//!
//! What it reports is a MEAN over tiles rather than a sum: two tiles of a thousand and five
//! hundred clusters make a lane of seven hundred and fifty.
//!
//! Ported from `picard.illumina.CollectIlluminaLaneMetrics`,
//! `picard.illumina.parser.TileMetricsUtil` and `picard.illumina.parser.IlluminaMetricsCode` in
//! Picard 3.4.0.

use std::collections::BTreeMap;

use crate::illumina_files::TileMetric;

/// The codes the tool reads.
pub const CLUSTER_COUNT: u16 = 100;
pub const PF_CLUSTER_COUNT: u16 = 101;
pub const CLUSTER_DENSITY: u16 = 102;
pub const PF_CLUSTER_DENSITY: u16 = 103;
pub const PHASING_BASE: u16 = 200;
pub const PREPHASING_BASE: u16 = 201;

/// `IlluminaMetricsCode.getPhasingCode`: the base plus twice the descriptor's index.
pub fn phasing_code(descriptor: usize, base: u16) -> u16 {
    base + (descriptor as u16) * 2
}

/// One lane's row.
#[derive(Debug, Clone, PartialEq)]
pub struct LaneMetrics {
    pub lane: u16,
    pub cluster_density: f64,
}

/// One phasing row: a lane, which template read it is, and the two values.
#[derive(Debug, Clone, PartialEq)]
pub struct PhasingMetrics {
    pub lane: u16,
    pub read: usize,
    pub phasing_applied: f64,
    pub prephasing_applied: f64,
}

/// What a file can be wrong about, in the reference's own words.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// `Expected to find cluster and density record codes (102 and 100) in records read for tile
    /// location`.
    MissingCounts { lane: u16, tile: u16 },
    /// `Don't have both phasing and prephasing values for <which> read cycle <cycle>.`
    HalfAPhasingPair {
        which: &'static str,
        cycle: usize,
        phasing: u16,
        prephasing: u16,
    },
}

/// The name a read descriptor is reported by, which is not its index.
pub fn descriptor_name(index: usize) -> &'static str {
    match index {
        0 => "FIRST",
        1 => "SECOND",
        _ => "THIRD",
    }
}

/// The lane rows and the phasing rows a file makes, or the first refusal it earns.
///
/// `templates` is the descriptor index of each template read: `4T` is `[0]` and `4T8B4T` is
/// `[0, 2]`, which is what makes the phasing codes 200/201 and 204/205 rather than 202/203.
pub fn collect(
    metrics: &[TileMetric],
    templates: &[usize],
) -> Result<(Vec<LaneMetrics>, Vec<PhasingMetrics>), Refusal> {
    // Group by lane and tile, keeping the codes of each.
    let mut tiles: BTreeMap<(u16, u16), BTreeMap<u16, f32>> = BTreeMap::new();
    for metric in metrics {
        tiles
            .entry((metric.lane, metric.tile))
            .or_default()
            .insert(metric.code, metric.value);
    }

    let mut densities: BTreeMap<u16, Vec<f64>> = BTreeMap::new();
    let mut phasing: Vec<PhasingMetrics> = Vec::new();
    for ((lane, tile), codes) in &tiles {
        if !codes.contains_key(&CLUSTER_COUNT) || !codes.contains_key(&CLUSTER_DENSITY) {
            return Err(Refusal::MissingCounts {
                lane: *lane,
                tile: *tile,
            });
        }
        densities
            .entry(*lane)
            .or_default()
            .push(codes[&CLUSTER_COUNT] as f64);
        for (read, descriptor) in templates.iter().enumerate() {
            let phasing_code_ = phasing_code(*descriptor, PHASING_BASE);
            let prephasing_code = phasing_code(*descriptor, PREPHASING_BASE);
            match (codes.get(&phasing_code_), codes.get(&prephasing_code)) {
                (Some(first), Some(second)) => phasing.push(PhasingMetrics {
                    lane: *lane,
                    read,
                    // The file writes the values as percentages of themselves.
                    phasing_applied: f64::from(*first) * 100.0,
                    prephasing_applied: f64::from(*second) * 100.0,
                }),
                _ => {
                    return Err(Refusal::HalfAPhasingPair {
                        which: descriptor_name(read),
                        cycle: descriptor * 2 + 1,
                        phasing: phasing_code_,
                        prephasing: prephasing_code,
                    })
                }
            }
        }
    }

    let lanes = densities
        .into_iter()
        .map(|(lane, values)| LaneMetrics {
            lane,
            cluster_density: values.iter().sum::<f64>() / values.len() as f64,
        })
        .collect();
    Ok((lanes, phasing))
}
