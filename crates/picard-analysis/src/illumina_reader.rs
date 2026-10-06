//! The clusters of an Illumina run, as `picard.illumina.parser.IlluminaDataProvider` hands them
//! out for a per-tile `.bcl` run: per tile, per cluster, the bases and qualities of each output
//! read, the filter's verdict and the position.
//!
//! Ported from `IlluminaDataProviderFactory`, `IlluminaDataProvider`, `BclParser`,
//! `BaseBclReader`, `FilterParser`, `LocsFileReader` and `BclQualityEvaluationStrategy` at tag
//! 3.4.0. A basecall byte's low two bits are the base and its high six the quality, a byte of zero
//! is a no-call (`N`), and every quality is raised to at least `minimum_quality` (Illumina's
//! alleged minimum, two). A position is `Math.round(coordinate * 10 + 1000)` in float arithmetic.

use crate::illumina_dir::{long_lane, per_tile_files, Layout, PerTilePerCycle};

/// One cluster.
#[derive(Debug, Clone)]
pub struct Cluster {
    pub tile: i32,
    pub x: i32,
    pub y: i32,
    pub pf: bool,
    /// Per output read, its bases and qualities.
    pub reads: Vec<(Vec<u8>, Vec<u8>)>,
}

/// `AbstractIlluminaPositionFileReader.posToQSeqCoord`.
pub fn qseq_coordinate(pos: f32) -> i32 {
    let value = pos * 10.0 + 1000.0;
    (value + 0.5).floor() as i32
}

/// The `.locs` positions, as floats.
fn locs(bytes: &[u8]) -> Vec<(f32, f32)> {
    let mut out = Vec::new();
    if bytes.len() < 12 {
        return out;
    }
    let count = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize;
    for i in 0..count {
        let at = 12 + i * 8;
        if at + 8 > bytes.len() {
            break;
        }
        let x = f32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]);
        let y = f32::from_le_bytes([bytes[at + 4], bytes[at + 5], bytes[at + 6], bytes[at + 7]]);
        out.push((x, y));
    }
    out
}

/// A run directory's per-tile files for one lane.
pub struct Run {
    pub layout: Layout,
    pub lane: i32,
    pub bcl: PerTilePerCycle,
}

impl Run {
    pub fn new(basecalls: &std::path::Path, lane: i32) -> Run {
        let layout = Layout::new(basecalls);
        let bcl = PerTilePerCycle::new(".bcl", &basecalls.join(long_lane(lane)), lane);
        Run { layout, lane, bcl }
    }

    /// `getAvailableTiles`: the tiles the BCLs cover, in `TILE_NUMBER_COMPARATOR` order (by
    /// number).
    pub fn available_tiles(&self) -> Result<Vec<i32>, String> {
        let mut tiles: Vec<i32> = self
            .bcl
            .files
            .values()
            .flat_map(|m| m.keys().copied())
            .collect();
        tiles.sort_unstable();
        tiles.dedup();
        if tiles.is_empty() {
            return Err(format!(
                "picard.PicardException: No available tiles were found, make sure that {} has a lane {}",
                crate::fingerprint::reference_path(&self.layout.basecalls.display().to_string()),
                self.lane
            ));
        }
        Ok(tiles)
    }

    /// `makeDataProvider(tile)` drained: every cluster of the tile, `reads` being each output
    /// read's cycles.
    pub fn clusters(
        &self,
        tile: i32,
        reads: &[Vec<i32>],
        minimum_quality: u8,
    ) -> Result<Vec<Cluster>, String> {
        let cycles: Vec<i32> = reads.iter().flatten().copied().collect();
        let available = self.bcl.files.len();
        if cycles.iter().any(|c| !self.bcl.files.contains_key(c)) {
            return Err(format!(
                "picard.PicardException: Expected CycledIlluminaFileMap to contain {} cycles but only {available} were found!",
                cycles.len()
            ));
        }
        let read_bcl = |cycle: i32| -> Vec<u8> {
            self.bcl
                .files
                .get(&cycle)
                .and_then(|m| m.get(&tile))
                .and_then(|p| std::fs::read(p).ok())
                .unwrap_or_default()
        };
        let lane_dir = self.layout.basecalls.join(long_lane(self.lane));
        let filter = per_tile_files(&lane_dir, self.lane, ".filter")
            .get(&tile)
            .and_then(|p| std::fs::read(p).ok())
            .unwrap_or_default();
        let pf: Vec<bool> = filter.iter().skip(12).map(|b| b & 1 == 1).collect();
        let intensity_lane = self.layout.intensities.join(long_lane(self.lane));
        let positions = match per_tile_files(&intensity_lane, self.lane, ".locs").get(&tile) {
            Some(p) => locs(&std::fs::read(p).unwrap_or_default()),
            None => {
                locs(&std::fs::read(self.layout.intensities.join("s.locs")).unwrap_or_default())
            }
        };
        let data: Vec<Vec<u8>> = cycles.iter().map(|c| read_bcl(*c)).collect();
        let count = data
            .first()
            .filter(|d| d.len() >= 4)
            .map(|d| u32::from_le_bytes([d[0], d[1], d[2], d[3]]) as usize)
            .unwrap_or(0);
        let mut out = Vec::with_capacity(count);
        for i in 0..count {
            let mut at = 0;
            let mut read_data = Vec::with_capacity(reads.len());
            for read in reads {
                let mut bases = Vec::with_capacity(read.len());
                let mut quals = Vec::with_capacity(read.len());
                for _ in read {
                    let byte = data[at].get(4 + i).copied().unwrap_or(0);
                    at += 1;
                    if byte == 0 {
                        bases.push(b'N');
                        quals.push(minimum_quality);
                    } else {
                        bases.push(b"ACGT"[(byte & 3) as usize]);
                        quals.push((byte >> 2).max(minimum_quality));
                    }
                }
                read_data.push((bases, quals));
            }
            let (x, y) = positions.get(i).copied().unwrap_or((0.0, 0.0));
            out.push(Cluster {
                tile,
                x: qseq_coordinate(x),
                y: qseq_coordinate(y),
                pf: pf.get(i).copied().unwrap_or(false),
                reads: read_data,
            });
        }
        Ok(out)
    }
}
