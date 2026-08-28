//! `PositionBasedDownsampleSam`: which reads a positional mask keeps.
//!
//! Reading the file and parsing the read names are not ported. What is ported is the selection,
//! which has no randomness in it at all: the same file downsampled twice keeps exactly the same
//! reads.
//!
//! Ported from `picard.sam.PositionBasedDownsampleSam` in Picard 3.4.0.

use std::collections::BTreeMap;

/// `PositionBasedDownsampleSam.ACCEPTABLE_FUDGE_FACTOR`, which only decides whether to warn.
pub const ACCEPTABLE_FUDGE_FACTOR: f64 = 0.2;
/// `PositionBasedDownsampleSam.PG_PROGRAM_NAME`, which the second-run guard looks for.
pub const PG_PROGRAM_NAME: &str = "PositionBasedDownsampleSam";
/// `customCommandLineValidation`, on a fraction outside the unit interval.
pub fn fraction_out_of_range_message(fraction: f64) -> String {
    format!("FRACTION must be a value between 0 and 1, found: {fraction}")
}

/// One read's place on the flowcell, which the read-name parser fills in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhysicalLocation {
    pub tile: i16,
    pub x: i32,
    pub y: i32,
}

/// `PositionBasedDownsampleSam.Coord`: one tile's extent and how many reads it holds.
///
/// It starts at ZERO on all four sides rather than at the first read's coordinates, so a tile
/// whose reads all sit past a thousand still has a minimum of nought and the mask is never
/// centred on the reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Coord {
    pub min_x: i32,
    pub min_y: i32,
    pub max_x: i32,
    pub max_y: i32,
    pub count: i32,
}

/// `fillTileMinMaxCoord`: the first pass, which reads every record to learn each tile's extent
/// and then widens it.
///
/// The widening is the span divided by the READ COUNT, added to each side, so a tile of few reads
/// has its boundary moved a long way and a tile of many barely at all. It is an integer division,
/// so a tile with more reads than span is not widened at all.
pub fn tile_extents(locations: &[PhysicalLocation]) -> BTreeMap<i16, Coord> {
    let mut tiles: BTreeMap<i16, Coord> = BTreeMap::new();
    for location in locations {
        let coord = tiles.entry(location.tile).or_default();
        coord.max_x = coord.max_x.max(location.x);
        coord.min_x = coord.min_x.min(location.x);
        coord.max_y = coord.max_y.max(location.y);
        coord.min_y = coord.min_y.min(location.y);
        coord.count += 1;
    }
    for coord in tiles.values_mut() {
        let span_x = coord.max_x - coord.min_x;
        let span_y = coord.max_y - coord.min_y;
        coord.max_x += span_x / coord.count;
        coord.min_x -= span_x / coord.count;
        coord.max_y += span_y / coord.count;
        coord.min_y -= span_y / coord.count;
    }
    tiles
}

/// `PositionBasedDownsampleSam.CircleSelector`, whose circle is chosen so that its area is the
/// fraction asked for and so that it overlaps each edge of the unit square by that fraction too.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircleSelector {
    pub radius_squared: f64,
    pub offset: f64,
    /// False above a half, where the circle is built for `1 - fraction` and its INSIDE is kept.
    pub positive_selection: bool,
}

impl CircleSelector {
    pub fn new(fraction: f64) -> CircleSelector {
        let (p, positive_selection) = if fraction > 0.5 {
            (1.0 - fraction, false)
        } else {
            (fraction, true)
        };
        let radius_squared = p / std::f64::consts::PI;
        CircleSelector {
            radius_squared,
            offset: (radius_squared - p * p / 4.0).sqrt(),
            positive_selection,
        }
    }

    /// `roundedPart`, the signed distance to the nearest whole number, which is what makes the
    /// mask repeat across the tile rather than sit once in the middle of it.
    ///
    /// `Math.round` is `floor(x + 0.5)`, so a value exactly halfway rounds UP even when negative.
    pub fn rounded_part(x: f64) -> f64 {
        x - (x + 0.5).floor()
    }

    /// `select`: whether the read's place in its tile falls outside the circle.
    ///
    /// The comparison is `>` and it is then exclusive-ored with the selection's sense, so above a
    /// half the reads INSIDE the circle are the ones kept.
    pub fn select(&self, location: &PhysicalLocation, tile: &Coord) -> bool {
        let x = Self::rounded_part(
            (f64::from(location.x - tile.min_x) / f64::from(tile.max_x - tile.min_x)) - self.offset,
        );
        let y = Self::rounded_part(
            (f64::from(location.y - tile.min_y) / f64::from(tile.max_y - tile.min_y)) - self.offset,
        );
        (x * x + y * y > self.radius_squared) != self.positive_selection
    }
}

/// `outputSamRecords`: the second pass, which keeps the reads the mask selects.
///
/// The two passes read the same records in the same order, so `STOP_AFTER` cuts both and moves
/// the tile extents as well as the read count.
pub fn keep(locations: &[PhysicalLocation], fraction: f64, stop_after: Option<usize>) -> Vec<bool> {
    let examined = match stop_after {
        Some(limit) => &locations[..limit.min(locations.len())],
        None => locations,
    };
    let tiles = tile_extents(examined);
    let selector = CircleSelector::new(fraction);
    let mut kept = vec![false; locations.len()];
    for (index, location) in examined.iter().enumerate() {
        let tile = tiles.get(&location.tile).copied().unwrap_or_default();
        kept[index] = selector.select(location, &tile);
    }
    kept
}

/// `doWork`'s closing check, which warns when the rate misses the fraction by more than a fifth
/// of the smaller of the two.
pub fn misses_the_fraction(kept: usize, total: usize, fraction: f64) -> bool {
    if total == 0 {
        return false;
    }
    let rate = kept as f64 / total as f64;
    (rate - fraction).abs() / (rate.min(fraction) + 1e-10) > ACCEPTABLE_FUDGE_FACTOR
}
