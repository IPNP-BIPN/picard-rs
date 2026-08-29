//! `CheckIlluminaDirectory`: whether a basecalls directory has the files a run would need.
//!
//! The cheapest of the Illumina tools, and the one that says what the others will refuse. It walks
//! the tiles a lane declares and asks, for each, whether every file the requested data types need
//! is there and readable.
//!
//! Two things decide the answer beyond the files themselves. The READ STRUCTURE says how many
//! cycles are wanted, so the same directory passes under `3T` and fails under `6T`. And the DATA
//! TYPES say which files count: the basecalls and the filter are asked for by default and the
//! POSITIONS are not, so a run with no `s.locs` passes until `--DATA_TYPES=Position` asks for it.
//!
//! Ported from `picard.illumina.CheckIlluminaDirectory` in Picard 3.4.0.

use crate::illumina_files::{total_cycles, Segment};

/// `IlluminaDataType`, of which three decide which files are wanted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    BaseCalls,
    QualityScores,
    PF,
    Position,
    Barcodes,
}

/// The default set: the basecalls, their qualities and the filter, and not the positions.
pub fn default_data_types() -> Vec<DataType> {
    vec![DataType::BaseCalls, DataType::QualityScores, DataType::PF]
}

/// One file a run needs, named the way the directory names it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Needed {
    /// `L001/C<cycle>.1/s_<lane>_<tile>.bcl`.
    BaseCall { cycle: usize },
    /// `L001/s_<lane>_<tile>.filter`.
    Filter,
    /// `s.locs`, which sits beside the basecalls directory rather than in it.
    Positions,
}

/// Which files one tile needs, under a structure and a set of data types.
///
/// The basecalls and their qualities are the same file, so asking for both asks for one; the
/// filter is the `PF` type's; and the positions are the lane's `s.locs` rather than a per-tile
/// file, which is why removing the tile's own `.locs` changes nothing.
pub fn needed(structure: &[Segment], types: &[DataType]) -> Vec<Needed> {
    let mut files = Vec::new();
    if types
        .iter()
        .any(|type_| matches!(type_, DataType::BaseCalls | DataType::QualityScores))
    {
        for cycle in 1..=total_cycles(structure) {
            files.push(Needed::BaseCall { cycle });
        }
    }
    if types.contains(&DataType::PF) {
        files.push(Needed::Filter);
    }
    if types.contains(&DataType::Position) {
        files.push(Needed::Positions);
    }
    files
}

/// The status a run exits with, which is the NUMBER of failures rather than one.
///
/// A directory with everything is zero. A tile the lane does not declare is not a failure but a
/// refusal, which the caller reports separately: the tool throws rather than counting it.
pub fn failures(present: &dyn Fn(&Needed) -> bool, wanted: &[Needed]) -> usize {
    wanted.iter().filter(|file| !present(file)).count()
}
