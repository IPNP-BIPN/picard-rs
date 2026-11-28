//! The tile metrics file of an Illumina run directory, and the read structure that cuts its cycles.
//!
//! A `TileMetricsOut.bin` is a version byte, a record size, and then that many bytes per record: a
//! lane, a tile and a metric code as unsigned shorts, and a float, all little-endian. The
//! basecall, filter and locs files are read by `illumina_reader`.
//!
//! Ported from `picard.illumina.parser.readers.TileMetricsOutReader` and
//! `picard.illumina.parser.ReadStructure` in Picard 3.4.0.

/// One record of a tile metrics file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TileMetric {
    pub lane: u16,
    pub tile: u16,
    pub code: u16,
    pub value: f32,
}

/// `TileMetricsOutReader` at version two: a version byte, a record size, and ten bytes a record.
/// Java refuses any other record size; a 0 here would never advance.
pub fn parse_tile_metrics(bytes: &[u8]) -> Option<Vec<TileMetric>> {
    if bytes.len() < 2 || bytes[0] != 2 || bytes[1] != 10 {
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

#[cfg(test)]
mod tile_metrics_record_size {
    use super::parse_tile_metrics;

    /// A record size of 0 would never advance; Java refuses anything but 10.
    #[test]
    fn only_ten_byte_records_are_read() {
        let mut bytes = vec![2u8, 10];
        bytes.extend_from_slice(&[1, 0, 2, 0, 100, 0, 0, 0, 128, 63]);
        assert_eq!(parse_tile_metrics(&bytes).expect("ten").len(), 1);
        bytes[1] = 0;
        assert!(parse_tile_metrics(&bytes).is_none());
        bytes[1] = 6;
        assert!(parse_tile_metrics(&bytes).is_none());
    }
}
