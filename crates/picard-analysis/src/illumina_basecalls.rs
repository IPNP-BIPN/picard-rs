//! What a cluster becomes: `IlluminaBasecallsToFastq` and `IlluminaBasecallsToSam`.
//!
//! Both walk the same clusters and cut them the same way; what differs is the destination, and
//! what a BAM can say that a FASTQ cannot.
//!
//! A FASTQ writes one file per template read and puts the run's identity in every read NAME. A BAM
//! writes one record per template read and puts that identity once, in a read group; two template
//! reads become a PAIR, with the flags to prove it; the filter file's verdict becomes the
//! vendor-check flag on a record that is still there rather than a record that is gone; and the
//! barcode becomes a tag.
//!
//! Ported from `picard.illumina.IlluminaBasecallsToFastq` and
//! `picard.illumina.IlluminaBasecallsToSam` in Picard 3.4.0.

use crate::illumina_files::{BaseCall, Segment, SegmentKind};

/// One cluster: its calls in cycle order, whether it passed the filter, and where it sat.
#[derive(Debug, Clone, PartialEq)]
pub struct Cluster {
    pub calls: Vec<BaseCall>,
    pub passed_filter: bool,
    pub x: i32,
    pub y: i32,
}

/// Where a run came from, which is what a read name or a read group carries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Run {
    pub machine: String,
    pub run_barcode: String,
    pub flowcell: String,
    pub lane: i32,
    pub tile: i32,
}

/// The coordinates a read name carries, which are NOT the floats in the `.locs` file.
///
/// Illumina's own convention, which the reference reproduces: ten times the position, rounded, plus
/// a thousand. So a cluster written at `100.0, 200.0` is named `2000:3000`.
pub fn position_in_name(x: f32, y: f32) -> (i32, i32) {
    (
        (x * 10.0).round() as i32 + 1000,
        (y * 10.0).round() as i32 + 1000,
    )
}

/// `ReadNameFormat`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadNameFormat {
    /// `<machine>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>`, which is the default.
    Casava18,
    /// `<run>:<lane>:<tile>:<x>:<y>`.
    Illumina,
}

/// `IlluminaBasecallsToFastq`'s read name, which carries the position and the filter's verdict.
///
/// The CASAVA form is followed by a space, the read's ordinal where there is more than one, the
/// filter's verdict as `Y` or `N`, a zero, and the barcode. A single-read run leaves the ordinal
/// empty, which is why its names carry a bare space before the colon.
pub fn read_name(
    run: &Run,
    cluster: &Cluster,
    format: ReadNameFormat,
    ordinal: Option<usize>,
    barcode: Option<&str>,
) -> String {
    match format {
        ReadNameFormat::Illumina => format!(
            "{}:{}:{}:{}:{}",
            run.run_barcode, run.lane, run.tile, cluster.x, cluster.y
        ),
        ReadNameFormat::Casava18 => format!(
            "{}:{}:{}:{}:{}:{}:{} {}:{}:0:{}",
            run.machine,
            run.run_barcode,
            run.flowcell,
            run.lane,
            run.tile,
            cluster.x,
            cluster.y,
            ordinal.map(|value| value.to_string()).unwrap_or_default(),
            // The filter's verdict is inverted in the name: `Y` means the read was FILTERED.
            if cluster.passed_filter { "N" } else { "Y" },
            barcode.unwrap_or("")
        ),
    }
}

/// The bases and qualities of one segment of a cluster.
pub fn segment(cluster: &Cluster, structure: &[Segment], index: usize) -> (Vec<u8>, Vec<u8>) {
    let mut start = 0;
    for (position, part) in structure.iter().enumerate() {
        if position == index {
            let calls = &cluster.calls[start..start + part.cycles];
            return (
                calls.iter().map(|call| call.base).collect(),
                calls.iter().map(|call| call.quality).collect(),
            );
        }
        start += part.cycles;
    }
    (Vec::new(), Vec::new())
}

/// Which segments are written as reads, which is neither the barcodes nor the skips.
pub fn written_segments(structure: &[Segment]) -> Vec<usize> {
    structure
        .iter()
        .enumerate()
        .filter(|(_, part)| part.kind == SegmentKind::Template)
        .map(|(index, _)| index)
        .collect()
}

/// A quality as FASTQ writes it: the value plus thirty-three.
pub fn phred33(qualities: &[u8]) -> Vec<u8> {
    qualities.iter().map(|quality| quality + 33).collect()
}

/// The flags a BAM record carries, which is where the pairing lives.
///
/// A single template read is unmapped and nothing else: `0x4`. Two are a pair, so the first is
/// `0x1 | 0x4 | 0x8 | 0x40` and the second `0x1 | 0x4 | 0x8 | 0x80`. A cluster that failed the
/// filter adds `0x200` rather than being dropped.
pub fn sam_flags(reads: usize, ordinal: usize, passed_filter: bool) -> u16 {
    let mut flags: u16 = 0x4;
    if reads > 1 {
        flags |= 0x1 | 0x8;
        flags |= if ordinal == 0 { 0x40 } else { 0x80 };
    }
    if !passed_filter {
        flags |= 0x200;
    }
    flags
}

/// The read group's identifier and platform unit, which carry the run and, where asked, the
/// barcode.
pub fn read_group(run: &Run, barcode: Option<&str>, include_barcode: bool) -> (String, String) {
    let id = format!("{}.{}", run.run_barcode, run.lane);
    let unit = match (barcode, include_barcode) {
        (Some(barcode), _) => format!("{}.{}.{}", run.run_barcode, run.lane, barcode),
        (None, _) => id.clone(),
    };
    let _ = include_barcode;
    (id, unit)
}
