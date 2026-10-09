//! An Illumina run directory as `picard.illumina.parser.IlluminaFileUtil` finds its way around it:
//! which lane and cycle directories exist, which tiles each per-tile file covers, and which tiles
//! the lane has according to `InterOp/TileMetricsOut.bin`.
//!
//! Ported from `picard.illumina.parser.IlluminaFileUtil`, `ParameterizedFileUtil`,
//! `PerTileFileUtil`, `PerTileOrPerRunFileUtil` and `PerTilePerCycleFileUtil` at tag 3.4.0, for
//! the per-tile formats (`.bcl`, `.filter`, `.locs`, `s.locs`). The multi-tile and `.cbcl`
//! layouts need a tile index or compressed cycle files and are reported as unavailable, as they
//! are for a directory without them.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// `IlluminaFileUtil.longLaneStr`: `L` and the lane, zero-padded to three digits.
pub fn long_lane(lane: i32) -> String {
    format!("L{lane:03}")
}

/// The places a run's files sit, from its basecalls directory.
pub struct Layout {
    pub basecalls: PathBuf,
    pub intensities: PathBuf,
    pub tile_metrics: PathBuf,
}

impl Layout {
    /// `new IlluminaFileUtil(basecallDir, lane)`: the intensities directory is the basecalls one's
    /// parent, and `InterOp` sits beside the parent of that.
    pub fn new(basecalls: &Path) -> Layout {
        let intensities = basecalls
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_default();
        let data = intensities
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_default();
        let run = data.parent().map(Path::to_path_buf).unwrap_or_default();
        Layout {
            basecalls: basecalls.to_path_buf(),
            intensities,
            tile_metrics: run.join("InterOp").join("TileMetricsOut.bin"),
        }
    }
}

/// The non-empty files in `dir` named `s_<lane>_<tile><extension>`, by tile: an `IlluminaFileMap`,
/// which is a `TreeMap`.
pub fn per_tile_files(dir: &Path, lane: i32, extension: &str) -> BTreeMap<i32, PathBuf> {
    let mut map = BTreeMap::new();
    let prefix = format!("s_{lane}_");
    let Ok(entries) = std::fs::read_dir(dir) else {
        return map;
    };
    for entry in entries.flatten() {
        let name = entry.file_name().to_string_lossy().into_owned();
        let Some(rest) = name.strip_prefix(&prefix) else {
            continue;
        };
        let Some(tile) = rest.strip_suffix(extension) else {
            continue;
        };
        if tile.is_empty() || tile.len() > 5 || !tile.bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }
        let path = entry.path();
        if std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0) > 0 {
            map.insert(tile.parse().unwrap_or(0), path);
        }
    }
    map
}

/// The cycle directories of a lane, `C<cycle>.1`, by cycle.
pub fn cycle_dirs(lane_dir: &Path) -> BTreeMap<i32, PathBuf> {
    let mut map = BTreeMap::new();
    let Ok(entries) = std::fs::read_dir(lane_dir) else {
        return map;
    };
    for entry in entries.flatten() {
        let name = entry.file_name().to_string_lossy().into_owned();
        if let Some(cycle) = name.strip_prefix('C').and_then(|r| r.strip_suffix(".1")) {
            if !cycle.is_empty() && cycle.bytes().all(|b| b.is_ascii_digit()) {
                map.insert(cycle.parse().unwrap_or(0), entry.path());
            }
        }
    }
    map
}

/// `IlluminaFileUtil.getExpectedTiles`: the tiles `TileMetricsOut.bin` names for the lane, sorted.
pub fn expected_tiles(layout: &Layout, lane: i32) -> Result<Vec<i32>, String> {
    let bytes = std::fs::read(&layout.tile_metrics).map_err(|_| {
        format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
            crate::fingerprint::reference_path(&layout.tile_metrics.display().to_string())
        )
    })?;
    let metrics = crate::illumina_files::parse_tile_metrics(&bytes).unwrap_or_default();
    let mut tiles: Vec<i32> = metrics
        .iter()
        .filter(|m| i32::from(m.lane) == lane)
        .map(|m| i32::from(m.tile))
        .collect();
    tiles.sort_unstable();
    tiles.dedup();
    Ok(tiles)
}

/// `PerTilePerCycleFileUtil` for one extension.
pub struct PerTilePerCycle {
    pub extension: String,
    pub base: PathBuf,
    /// Cycle to tile to file.
    pub files: BTreeMap<i32, BTreeMap<i32, PathBuf>>,
}

impl PerTilePerCycle {
    pub fn new(extension: &str, lane_dir: &Path, lane: i32) -> PerTilePerCycle {
        let files = cycle_dirs(lane_dir)
            .into_iter()
            .map(|(cycle, dir)| (cycle, per_tile_files(&dir, lane, extension)))
            .collect();
        PerTilePerCycle {
            extension: extension.to_string(),
            base: lane_dir.to_path_buf(),
            files,
        }
    }

    pub fn files_available(&self) -> bool {
        self.files.values().any(|m| !m.is_empty())
    }

    /// `verify`: per expected cycle, per expected tile, a file of the same length as the tile's
    /// first.
    pub fn verify(&self, tiles: &[i32], cycles: &[i32]) -> Vec<String> {
        let mut failures = Vec::new();
        let base = crate::fingerprint::reference_path(&self.base.display().to_string());
        if !self.base.exists() {
            failures.push(format!("Base directory({base}) does not exist!"));
            return failures;
        }
        let mut first_length: BTreeMap<i32, u64> = BTreeMap::new();
        for &cycle in cycles {
            let Some(map) = self.files.get(&cycle) else {
                failures.push(format!(
                    "Missing file for cycle {cycle} in directory {base} for file type {}",
                    self.extension
                ));
                continue;
            };
            for &tile in tiles {
                match map.get(&tile) {
                    Some(file) => {
                        let length = std::fs::metadata(file).map(|m| m.len()).unwrap_or(0);
                        match first_length.get(&tile) {
                            None => {
                                first_length.insert(tile, length);
                            }
                            Some(&first) if self.extension != ".bcl.gz" && first != length => {
                                failures.push(format!(
                                    "File type {} has cycles files of different length.  Current cycle ({cycle}) Length of first non-empty file ({first}) length of current cycle ({length}) File({})",
                                    self.extension,
                                    crate::fingerprint::reference_path(&file.display().to_string())
                                ));
                            }
                            _ => {}
                        }
                    }
                    None => failures.push(format!(
                        "File type {} is missing a file for cycle {cycle} and tile {tile}",
                        self.extension
                    )),
                }
            }
        }
        failures
    }
}

/// `PerTileFileUtil`, and `PerTileOrPerRunFileUtil` when `run_file` is set.
pub struct PerTile {
    pub extension: String,
    pub base: PathBuf,
    pub files: BTreeMap<i32, PathBuf>,
    pub run_file: Option<PathBuf>,
}

impl PerTile {
    pub fn new(extension: &str, base: &Path, lane: i32) -> PerTile {
        PerTile {
            extension: extension.to_string(),
            base: base.to_path_buf(),
            files: per_tile_files(base, lane, extension),
            run_file: None,
        }
    }

    /// `PerTileOrPerRunFileUtil`: also `s<extension>` in the base's parent.
    pub fn or_per_run(extension: &str, base: &Path, lane: i32) -> PerTile {
        let mut util = PerTile::new(extension, base, lane);
        let candidate = base
            .parent()
            .map(|p| p.join(format!("s{extension}")))
            .filter(|p| p.is_file());
        util.run_file = candidate;
        util
    }

    pub fn files_available(&self) -> bool {
        !self.files.is_empty() || self.run_file.is_some()
    }

    /// `setTilesForPerRunFile`, then `verify`.
    pub fn verify(&mut self, tiles: &[i32]) -> Vec<String> {
        let base = match &self.run_file {
            Some(run) => {
                for &t in tiles {
                    self.files.insert(t, run.clone());
                }
                run.parent().map(Path::to_path_buf).unwrap_or_default()
            }
            None => self.base.clone(),
        };
        let mut failures = Vec::new();
        if !base.exists() {
            failures.push(format!(
                "Base directory({}) does not exist!",
                crate::fingerprint::reference_path(&base.display().to_string())
            ));
        } else {
            let missing: Vec<String> = tiles
                .iter()
                .filter(|t| !self.files.contains_key(t))
                .map(|t| t.to_string())
                .collect();
            if !missing.is_empty() {
                failures.push(format!(
                    "Missing tile [{}] for file type {}.",
                    missing.join(", "),
                    self.extension
                ));
            }
        }
        failures
    }
}

/// An `htsjdk.samtools.util.Log` line at INFO, with a clock the harness drops.
pub fn log_info(class: &str, message: &str) {
    eprintln!("INFO\t1970-01-01 00:00:00\t{class}\t{message}");
}
