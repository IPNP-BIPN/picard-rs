//! The binary files an Illumina run directory carries, and the read structure that cuts them.
//!
//! Picard's Illumina tools read a DIRECTORY rather than a file, and four formats in it. None of
//! them is self-describing, so each is a fixed layout the reader knows by heart:
//!
//!  * a `.bcl` is an unsigned count of clusters and then one byte per cluster, whose low two bits
//!    are the base and whose high six are the quality. A byte of ZERO is a no-call whatever those
//!    bits say, which is why the quality cannot be read on its own;
//!  * a `.filter` is three unsigned ints (a zero, the version, which must be three, and the count)
//!    and then one byte per cluster;
//!  * a `.locs` is a little-endian one, a float version of 1.0, a count, and two floats per
//!    cluster;
//!  * a `TileMetricsOut.bin` is a version byte, a record size, and then that many bytes per
//!    record: a lane, a tile and a metric code as unsigned shorts, and a float.
//!
//! Everything is little-endian.
//!
//! Ported from `picard.illumina.parser.readers.BclReader`, `FilterFileReader`, `LocsFileReader`
//! and `TileMetricsOutReader` in Picard 3.4.0.

/// What a cluster's basecall byte says.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseCall {
    /// `A`, `C`, `G`, `T`, or `.` for a no-call.
    pub base: u8,
    /// The quality, which is zero for a no-call.
    pub quality: u8,
}

/// `BclReader`: the base is the low two bits and the quality is the rest.
///
/// A byte of zero is a no-call, and the reference writes it as a dot rather than as an `N`.
pub fn decode_basecall(byte: u8) -> BaseCall {
    if byte == 0 {
        return BaseCall {
            base: b'.',
            quality: 0,
        };
    }
    BaseCall {
        base: b"ACGT"[(byte & 0b11) as usize],
        quality: byte >> 2,
    }
}

/// What a `.bcl` carries: one call per cluster.
pub fn parse_bcl(bytes: &[u8]) -> Option<Vec<BaseCall>> {
    let count = read_u32(bytes, 0)? as usize;
    if bytes.len() < 4 + count {
        return None;
    }
    Some(
        bytes[4..4 + count]
            .iter()
            .map(|b| decode_basecall(*b))
            .collect(),
    )
}

/// `FilterFileReader`: the version must be three, and a cluster passed when its byte's low bit is
/// set.
pub fn parse_filter(bytes: &[u8]) -> Result<Vec<bool>, FilterError> {
    if bytes.len() < 12 {
        return Err(FilterError::TooShort);
    }
    let version = read_u32(bytes, 4).ok_or(FilterError::TooShort)?;
    if version != 3 {
        return Err(FilterError::Version(version));
    }
    let count = read_u32(bytes, 8).ok_or(FilterError::TooShort)? as usize;
    if bytes.len() < 12 + count {
        return Err(FilterError::TooShort);
    }
    Ok(bytes[12..12 + count].iter().map(|b| b & 1 == 1).collect())
}

/// What a filter file can be wrong about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterError {
    TooShort,
    /// The reference refuses any version but three, by name.
    Version(u32),
}

/// A cluster's position on the tile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Position {
    pub x: f32,
    pub y: f32,
}

/// `LocsFileReader`: a one, a float version of 1.0, a count, and a pair of floats per cluster.
pub fn parse_locs(bytes: &[u8]) -> Option<Vec<Position>> {
    if bytes.len() < 12 || read_u32(bytes, 0)? != 1 {
        return None;
    }
    if read_f32(bytes, 4)? != 1.0 {
        return None;
    }
    let count = read_u32(bytes, 8)? as usize;
    let mut positions = Vec::with_capacity(count);
    for index in 0..count {
        let at = 12 + index * 8;
        positions.push(Position {
            x: read_f32(bytes, at)?,
            y: read_f32(bytes, at + 4)?,
        });
    }
    Some(positions)
}

/// One record of a tile metrics file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TileMetric {
    pub lane: u16,
    pub tile: u16,
    pub code: u16,
    pub value: f32,
}

/// `TileMetricsOutReader` at version two: a version byte, a record size, and ten bytes a record.
pub fn parse_tile_metrics(bytes: &[u8]) -> Option<Vec<TileMetric>> {
    if bytes.len() < 2 || bytes[0] != 2 {
        return None;
    }
    let record = bytes[1] as usize;
    let mut metrics = Vec::new();
    let mut at = 2;
    while at + record <= bytes.len() {
        metrics.push(TileMetric {
            lane: read_u16(bytes, at)?,
            tile: read_u16(bytes, at + 2)?,
            code: read_u16(bytes, at + 4)?,
            value: read_f32(bytes, at + 6)?,
        });
        at += record;
    }
    Some(metrics)
}

fn read_u16(bytes: &[u8], at: usize) -> Option<u16> {
    Some(u16::from_le_bytes(bytes.get(at..at + 2)?.try_into().ok()?))
}

fn read_u32(bytes: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_le_bytes(bytes.get(at..at + 4)?.try_into().ok()?))
}

fn read_f32(bytes: &[u8], at: usize) -> Option<f32> {
    Some(f32::from_le_bytes(bytes.get(at..at + 4)?.try_into().ok()?))
}

/// One segment of a read structure: how many cycles, and what they are for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment {
    pub cycles: usize,
    pub kind: SegmentKind,
}

/// `ReadType`: what a segment's cycles become.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentKind {
    /// A read that is written out.
    Template,
    /// A barcode, which is read and not written as a read.
    Barcode,
    /// A molecular index.
    MolecularIndex,
    /// Cycles nobody reads.
    Skip,
}

/// `ReadStructure`: `4T`, `2T2B`, `4T8B4T`, and so on.
///
/// The letters are the reference's: `T`, `B`, `M` and `S`. A structure with a letter it does not
/// know, or a segment with no count, is not a structure.
pub fn parse_read_structure(text: &str) -> Option<Vec<Segment>> {
    let mut segments = Vec::new();
    let mut digits = String::new();
    for character in text.chars() {
        if character.is_ascii_digit() {
            digits.push(character);
            continue;
        }
        let kind = match character.to_ascii_uppercase() {
            'T' => SegmentKind::Template,
            'B' => SegmentKind::Barcode,
            'M' => SegmentKind::MolecularIndex,
            'S' => SegmentKind::Skip,
            _ => return None,
        };
        let cycles = digits.parse().ok()?;
        digits.clear();
        segments.push(Segment { cycles, kind });
    }
    if !digits.is_empty() || segments.is_empty() {
        return None;
    }
    Some(segments)
}

/// How many cycles a structure asks for, which is what decides whether a directory is complete.
pub fn total_cycles(structure: &[Segment]) -> usize {
    structure.iter().map(|segment| segment.cycles).sum()
}

/// The cycles of each segment, one-based, in order.
pub fn segment_cycles(structure: &[Segment]) -> Vec<Vec<usize>> {
    let mut cycle = 1;
    let mut out = Vec::new();
    for segment in structure {
        out.push((cycle..cycle + segment.cycles).collect());
        cycle += segment.cycles;
    }
    out
}
