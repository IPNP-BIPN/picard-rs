//! The clusters of an Illumina run, as `picard.illumina.parser.IlluminaDataProvider` hands them
//! out for a per-tile `.bcl` run: per tile, per cluster, the bases and qualities of each output
//! read, the filter's verdict and the position.
//!
//! Ported from `picard.illumina.parser.IlluminaDataProviderFactory`, `IlluminaDataProvider`,
//! `BclParser`, `BaseBclReader`, `FilterParser`, `LocsFileReader` and
//! `BclQualityEvaluationStrategy` at tag 3.4.0. A basecall byte's low two bits are the base and
//! its high six the quality, a byte of zero is a no-call (`N`), and every quality is raised to at
//! least `minimum_quality` (Illumina's alleged minimum, two). A position is
//! `Math.round(coordinate * 10 + 1000)` in float arithmetic.

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

fn absolute(path: &std::path::Path) -> String {
    std::path::absolute(path)
        .unwrap_or_else(|_| path.to_path_buf())
        .display()
        .to_string()
}

fn u32_at(bytes: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]])
}

/// The file read whole, a missing or unreadable one being `IOUtil`'s refusal.
fn read(path: &std::path::Path) -> Result<Vec<u8>, String> {
    std::fs::read(path).map_err(|e| {
        format!(
            "picard.PicardException: IOException opening cluster binary file {}: {e}",
            path.display()
        )
    })
}

/// `MMapBackedIteratorFactory`'s header check and `assertTotalElementsEqual`.
fn elements(
    bytes: &[u8],
    header: usize,
    size: usize,
    expected: u64,
    path: &std::path::Path,
) -> Result<(), String> {
    if bytes.len() < header {
        return Err(format!(
            "picard.PicardException: Header size({header}) is greater than file size({}) for file {}",
            bytes.len(),
            absolute(path)
        ));
    }
    let found = ((bytes.len() - header) / size) as u64;
    if found != expected {
        return Err(format!(
            "picard.PicardException: Expected {expected} elements in file but found {found} elements! File({})",
            absolute(path)
        ));
    }
    if !(bytes.len() - header).is_multiple_of(size) {
        return Err(format!(
            "picard.PicardException: Malformed file, expected {} bytes in file, found {} bytes for file({})",
            header as u64 + expected * size as u64,
            bytes.len(),
            absolute(path)
        ));
    }
    Ok(())
}

/// `FilterFileReader`: four zero bytes, version three, a count, a byte per cluster.
fn filter(path: &std::path::Path) -> Result<Vec<bool>, String> {
    let bytes = read(path)?;
    if bytes.len() < 12 {
        return elements(&bytes, 12, 1, 0, path).map(|_| Vec::new());
    }
    if let Some(i) = bytes[..4].iter().position(|&b| b != 0) {
        return Err(format!(
            "picard.PicardException: The first four bytes of a Filter File should be 0 but byte {i} was {} in file {}",
            bytes[i] as i8,
            absolute(path)
        ));
    }
    let version = u32_at(&bytes, 4) as i32;
    if version != 3 {
        return Err(format!(
            "picard.PicardException: Expected version is 3 but version found was {version} in file {}",
            absolute(path)
        ));
    }
    let count = u32_at(&bytes, 8);
    elements(&bytes, 12, 1, u64::from(count), path)?;
    Ok(bytes[12..].iter().map(|b| b & 1 == 1).collect())
}

/// `LocsFileReader`: a one, a float 1.0, a count, and a pair of floats per cluster.
fn locs(path: &std::path::Path) -> Result<Vec<(f32, f32)>, String> {
    let bytes = read(path)?;
    if bytes.len() < 12 {
        return elements(&bytes, 12, 4, 0, path).map(|_| Vec::new());
    }
    let first = u32_at(&bytes, 0) as i32;
    if first != 1 || f32::from_le_bytes(bytes[4..8].try_into().expect("four bytes")) != 1.0 {
        let version = if first != 1 { "1" } else { "1.0" };
        return Err(format!(
            "picard.PicardException: First header byte of locs files should be {version} value found({first})"
        ));
    }
    let count = u32_at(&bytes, 8) as usize;
    elements(&bytes, 12, 4, 2 * count as u64, path)?;
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let at = 12 + i * 8;
        let x = f32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]);
        let y = f32::from_le_bytes([bytes[at + 4], bytes[at + 5], bytes[at + 6], bytes[at + 7]]);
        out.push((x, y));
    }
    Ok(out)
}

/// `PerTileParser`'s refusal of a tile it has no file for.
fn no_file(tile: i32, files: &std::collections::BTreeMap<i32, std::path::PathBuf>) -> String {
    let keys: Vec<String> = files.keys().map(i32::to_string).collect();
    format!(
        "picard.PicardException: PerTileParser does not contain key({tile}) keys available ({})",
        keys.join(",")
    )
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
        let tile = self.tile(tile, reads, minimum_quality)?;
        Ok((0..tile.len()).map(|i| tile.cluster(i)).collect())
    }

    /// `makeDataProvider(tile)`: the tile's files read, its clusters decoded one at a time by
    /// [`Tile::cluster`].
    pub fn tile(&self, tile: i32, reads: &[Vec<i32>], minimum_quality: u8) -> Result<Tile, String> {
        let cycles: Vec<i32> = reads.iter().flatten().copied().collect();
        let available = self.bcl.files.len();
        if cycles.iter().any(|c| !self.bcl.files.contains_key(c)) {
            return Err(format!(
                "picard.PicardException: Expected CycledIlluminaFileMap to contain {} cycles but only {available} were found!",
                cycles.len()
            ));
        }
        // `CycleIlluminaFileMap.assertValid`, which counts the tiles of the lowest cycle alone.
        if let Some(first) = cycles.iter().min().and_then(|c| self.bcl.files.get(c)) {
            let tiles = self.available_tiles()?;
            let found = tiles.iter().filter(|t| first.contains_key(t)).count();
            if found != tiles.len() {
                return Err(format!(
                    "picard.PicardException: Expected CycledIlluminaFileMap to contain {} tiles but only {found} were found!",
                    tiles.len()
                ));
            }
        }
        // `BclReader`: every cycle's file there, its header whole and its count its length.
        let read_bcl = |index: usize, cycle: i32| -> Result<Vec<u8>, String> {
            let path = self.bcl.files.get(&cycle).and_then(|m| m.get(&tile)).ok_or_else(|| {
                format!("htsjdk.samtools.util.RuntimeIOException: Could not find BCL file for cycle {index}")
            })?;
            let bytes = std::fs::read(path).map_err(|_| {
                format!(
                    "picard.PicardException: File not found: ({})",
                    path.display()
                )
            })?;
            if bytes.len() < 4 {
                return Err(format!(
                    "htsjdk.samtools.util.RuntimeIOException: BCL {} has invalid header structure.",
                    absolute(path)
                ));
            }
            let count = u32_at(&bytes, 0) as i32;
            if i64::from(count) != bytes.len() as i64 - 4 {
                return Err(format!(
                    "picard.PicardException: Expected {count} in file {} but found {}",
                    absolute(path),
                    bytes.len() - 4
                ));
            }
            Ok(bytes)
        };
        let lane_dir = self.layout.basecalls.join(long_lane(self.lane));
        let filters = per_tile_files(&lane_dir, self.lane, ".filter");
        let tile_locs = per_tile_files(
            &self.layout.intensities.join(long_lane(self.lane)),
            self.lane,
            ".locs",
        );
        let run_locs = self.layout.intensities.join("s.locs");
        // `IlluminaDataProviderFactory`: a data type with no file for any tile of the lane.
        // shortcut: `.clocs` and `_pos.txt` positions are not read, so a run with only those
        // is refused here too.
        let unmatched: Vec<&str> = [
            ("PF", filters.is_empty()),
            ("Position", tile_locs.is_empty() && !run_locs.exists()),
        ]
        .iter()
        .filter(|(_, missing)| *missing)
        .map(|(name, _)| *name)
        .collect();
        if !unmatched.is_empty() {
            return Err(format!(
                "picard.PicardException: Could not find a format with available files for the following data types: {}",
                unmatched.join(", ")
            ));
        }
        let pf = filter(filters.get(&tile).ok_or_else(|| no_file(tile, &filters))?)?;
        let positions = match tile_locs.get(&tile) {
            Some(p) => locs(p)?,
            None if run_locs.exists() => locs(&run_locs)?,
            None => return Err(no_file(tile, &tile_locs)),
        };
        // shortcut: a tile's BCLs are held whole, a byte per cluster per cycle (12 MB for 200,000
        // clusters of 58 cycles); read them a block of clusters at a time if tiles grow far past that.
        let data = cycles
            .iter()
            .enumerate()
            .map(|(index, &cycle)| read_bcl(index, cycle))
            .collect::<Result<Vec<_>, _>>()?;
        // The clusters are counted by the first BCL, or by the filter file when no cycle is read
        // (a provider of PF and positions alone).
        let count = match data.first() {
            Some(d) => u32_at(d, 0) as usize,
            None => pf.len(),
        };
        Ok(Tile {
            tile,
            count,
            lengths: reads.iter().map(Vec::len).collect(),
            minimum_quality,
            data,
            pf,
            positions,
        })
    }
}

/// One tile's files as read: a byte per cluster per cycle, the filter's verdicts and the positions.
pub struct Tile {
    pub tile: i32,
    count: usize,
    lengths: Vec<usize>,
    minimum_quality: u8,
    data: Vec<Vec<u8>>,
    pf: Vec<bool>,
    positions: Vec<(f32, f32)>,
}

impl Tile {
    /// The clusters the tile holds.
    pub fn len(&self) -> usize {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Whether cluster `i` passed the filter.
    pub fn pf(&self, i: usize) -> bool {
        self.pf.get(i).copied().unwrap_or(false)
    }

    /// Where cluster `i` sat, as a read name carries it.
    pub fn position(&self, i: usize) -> (i32, i32) {
        let (x, y) = self.positions.get(i).copied().unwrap_or((0.0, 0.0));
        (qseq_coordinate(x), qseq_coordinate(y))
    }

    /// Cluster `i`, decoded.
    pub fn cluster(&self, i: usize) -> Cluster {
        let minimum_quality = self.minimum_quality;
        let mut at = 0;
        let mut read_data = Vec::with_capacity(self.lengths.len());
        for &length in &self.lengths {
            let mut bases = Vec::with_capacity(length);
            let mut quals = Vec::with_capacity(length);
            for _ in 0..length {
                let byte = self.data[at].get(4 + i).copied().unwrap_or(0);
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
        let (x, y) = self.position(i);
        Cluster {
            tile: self.tile,
            x,
            y,
            pf: self.pf(i),
            reads: read_data,
        }
    }
}

/// `BclParser.runEamssForReadInPlace`: from the end of a read, a tally that falls by two at a
/// quality of thirty or more and rises by one under fifteen; where it peaks at one or more, every
/// quality from there on is masked to two, the mask first pulled back over a run of ten or more
/// `G`s (allowing one exception per ten).
pub fn eamss(bases: &[u8], qualities: &mut [u8]) {
    let mut tally: i32 = 0;
    let mut max_tally = i32::MIN;
    let mut index_of_max: i32 = -1;
    for i in (0..bases.len()).rev() {
        let quality = qualities[i];
        if quality >= 30 {
            tally -= 2;
        } else if quality < 15 {
            tally += 1;
        }
        if tally >= max_tally {
            index_of_max = i as i32;
            max_tally = tally;
        }
    }
    if max_tally < 1 {
        return;
    }
    let mut gs: i32 = 0;
    let mut exceptions: i32 = 0;
    let mut i = index_of_max;
    while i >= 0 {
        if bases[i as usize] == b'G' {
            gs += 1;
        } else {
            // `skipBy`: the nearest `G` behind, within the exceptions a run this long allows.
            let mut skip = None;
            for backup in 1..=i {
                let limit = ((gs + backup) / 10).max(1);
                if exceptions + backup > limit {
                    break;
                }
                if bases[(i - backup) as usize] == b'G' {
                    skip = Some(backup);
                    break;
                }
            }
            match skip {
                Some(s) => {
                    exceptions += s;
                    gs += s;
                    i -= s - 1;
                }
                None => break,
            }
        }
        i -= 1;
    }
    if gs >= 10 {
        index_of_max = (index_of_max + 1) - gs;
    }
    for q in qualities.iter_mut().skip(index_of_max.max(0) as usize) {
        *q = 2;
    }
}

/// `BarcodeFileReader`: per cluster of a tile, the barcode it matched (`Y` in the second column,
/// the third column then), from `s_<lane>_<tile>_barcode.txt` in `dir`. `None` when the tile has
/// no such file.
pub fn matched_barcodes(
    dir: &std::path::Path,
    lane: i32,
    tile: i32,
) -> Option<Vec<Option<String>>> {
    let files = per_tile_files(dir, lane, "_barcode.txt");
    let path = files.get(&tile)?;
    let text = std::fs::read_to_string(path).ok()?;
    Some(
        text.lines()
            .map(|line| {
                let fields: Vec<&str> = line.split_whitespace().collect();
                (fields.get(1) == Some(&"Y")).then(|| fields.get(2).unwrap_or(&"").to_string())
            })
            .collect(),
    )
}

/// Whether a lane has any barcode file in `dir`.
pub fn has_barcode_files(dir: &std::path::Path, lane: i32) -> bool {
    !per_tile_files(dir, lane, "_barcode.txt").is_empty()
}
