//! `BasecallsConverter`: the clusters `IlluminaBasecallsToFastq` and `IlluminaBasecallsToSam`
//! write, tile by tile, with the barcode each was assigned.
//!
//! Ported from `picard.illumina.BasecallsConverter`, `SortedBasecallsConverter`,
//! `UnsortedBasecallsConverter` and the per-tile `IlluminaDataProviderFactory` at tag 3.4.0. The
//! tiles of every lane are walked in `TILE_NUMBER_COMPARATOR` order, cut by `FIRST_TILE` and
//! `TILE_LIMIT`; each read is EAMSS-filtered when asked; a cluster failing the filter is dropped
//! unless non-PF reads are included; and a cluster's barcode is the per-tile barcode file's, or,
//! with an extractor, the inline match of its barcode reads, counted into the metrics the
//! extractor shares with the program.

use std::path::{Path, PathBuf};

use crate::barcode_extractor::{Extractor, Metric};
use crate::illumina_files::{Segment, SegmentKind};
use crate::illumina_reader::{eamss, has_barcode_files, matched_barcodes, Run};
use crate::metrics_cli::thrown;

/// A cluster as the converter hands it to a writer: every read the structure does not skip.
pub struct ClusterData {
    pub lane: i32,
    pub tile: i32,
    pub x: i32,
    pub y: i32,
    pub pf: bool,
    pub barcode: Option<String>,
    pub reads: Vec<(SegmentKind, Vec<u8>, Vec<u8>)>,
}

impl ClusterData {
    /// The reads of one kind, in order.
    pub fn of(&self, kind: SegmentKind) -> Vec<&(SegmentKind, Vec<u8>, Vec<u8>)> {
        self.reads.iter().filter(|r| r.0 == kind).collect()
    }
}

/// `TILE_NUMBER_COMPARATOR`: by the decimal string, except that a prefix sorts after.
pub fn tile_order(a: &i32, b: &i32) -> std::cmp::Ordering {
    let (s1, s2) = (a.to_string(), b.to_string());
    if s1.len() < s2.len() && s2.starts_with(&s1) {
        std::cmp::Ordering::Greater
    } else if s2.len() < s1.len() && s1.starts_with(&s2) {
        std::cmp::Ordering::Less
    } else {
        s1.cmp(&s2)
    }
}

pub struct Converter {
    runs: Vec<(i32, Run, Vec<i32>)>,
    pub tiles: Vec<i32>,
    cycles: Vec<Vec<i32>>,
    kinds: Vec<SegmentKind>,
    include_non_pf: bool,
    apply_eamss: bool,
    demultiplex: bool,
    barcodes_from: Option<PathBuf>,
    pub extractor: Option<Extractor>,
    /// The declared barcodes' metrics, which the extractor counts into.
    pub metrics: Vec<Metric>,
}

pub struct Options<'a> {
    pub basecalls: &'a Path,
    pub lanes: &'a [i32],
    pub structure: &'a [Segment],
    pub demultiplex: bool,
    /// Where the per-tile barcode files are, or `None` when the barcodes are matched inline.
    pub barcodes_dir: Option<PathBuf>,
    pub first_tile: Option<i32>,
    pub tile_limit: Option<usize>,
    pub include_non_pf: bool,
    pub apply_eamss: bool,
    pub extractor: Option<Extractor>,
    pub metrics: Vec<Metric>,
}

impl Converter {
    /// The constructor: a provider factory per lane, which is where a missing data type and a lane
    /// without tiles are refused, then the tile limits.
    pub fn new(o: Options) -> Converter {
        let has_sample_barcode = o.structure.iter().any(|s| s.kind == SegmentKind::Barcode);
        let wants_barcodes = has_sample_barcode && o.demultiplex && o.barcodes_dir.is_some();
        let mut runs = Vec::new();
        for &lane in o.lanes {
            let run = Run::new(o.basecalls, lane);
            if wants_barcodes
                && !has_barcode_files(o.barcodes_dir.as_deref().unwrap_or(o.basecalls), lane)
            {
                thrown("picard.PicardException: Could not find a format with available files for the following data types: Barcodes");
            }
            let tiles = run.available_tiles().unwrap_or_else(|e| thrown(&e));
            runs.push((lane, run, tiles));
        }
        let mut tiles: Vec<i32> = runs.iter().flat_map(|(_, _, t)| t.clone()).collect();
        tiles.sort_by(tile_order);
        tiles.dedup();
        if let Some(first) = o.first_tile {
            if let Some(i) = tiles.iter().position(|t| *t == first) {
                tiles = tiles.split_off(i);
            }
            if tiles.first() != Some(&first) {
                thrown(&format!(
                    "picard.PicardException: firstTile={first}, but that tile was not found."
                ));
            }
        }
        if let Some(limit) = o.tile_limit {
            tiles.truncate(limit);
        }
        let mut cycles = Vec::new();
        let mut kinds = Vec::new();
        let mut cycle = 1;
        for s in o.structure {
            let range: Vec<i32> = (cycle..cycle + s.cycles as i32).collect();
            cycle += s.cycles as i32;
            if s.kind != SegmentKind::Skip {
                cycles.push(range);
                kinds.push(s.kind);
            }
        }
        Converter {
            runs,
            tiles,
            cycles,
            kinds,
            include_non_pf: o.include_non_pf,
            apply_eamss: o.apply_eamss,
            demultiplex: o.demultiplex,
            barcodes_from: if wants_barcodes { o.barcodes_dir } else { None },
            extractor: o.extractor,
            metrics: o.metrics,
        }
    }

    /// One tile's clusters, per lane that has it, in the order the provider reads them.
    pub fn tile(&mut self, tile: i32) -> Vec<ClusterData> {
        let mut out = Vec::new();
        let sample: Vec<usize> = self
            .kinds
            .iter()
            .enumerate()
            .filter(|(_, k)| **k == SegmentKind::Barcode)
            .map(|(i, _)| i)
            .collect();
        for (lane, run, available) in &self.runs {
            if !available.contains(&tile) {
                continue;
            }
            let clusters = run
                .clusters(tile, &self.cycles, 2)
                .unwrap_or_else(|e| thrown(&e));
            let from_file = self
                .barcodes_from
                .as_deref()
                .and_then(|dir| matched_barcodes(dir, *lane, tile));
            for (i, c) in clusters.into_iter().enumerate() {
                if !(self.include_non_pf || c.pf) {
                    continue;
                }
                let mut data = ClusterData {
                    lane: *lane,
                    tile: c.tile,
                    x: c.x,
                    y: c.y,
                    pf: c.pf,
                    barcode: from_file.as_ref().and_then(|f| f.get(i).cloned().flatten()),
                    reads: c
                        .reads
                        .into_iter()
                        .zip(&self.kinds)
                        .map(|((bases, mut quals), kind)| {
                            if self.apply_eamss {
                                eamss(&bases, &mut quals);
                            }
                            (*kind, bases, quals)
                        })
                        .collect(),
                };
                if self.demultiplex {
                    if let Some(ex) = self.extractor.as_mut() {
                        // `maybeDemultiplex`: an inline match, counted per tile and merged into
                        // the declared barcodes' metrics. The unmatched count goes to the
                        // extractor's own row, which nothing writes.
                        let read: Vec<Vec<u8>> =
                            sample.iter().map(|r| data.reads[*r].1.clone()).collect();
                        let quals: Vec<Vec<u8>> =
                            sample.iter().map(|r| data.reads[*r].2.clone()).collect();
                        let m = ex.find(&read, Some(&quals), true);
                        if m.matched {
                            if let Some(t) = self
                                .metrics
                                .iter_mut()
                                .find(|x| x.barcode.replace('-', "") == m.barcode)
                            {
                                t.reads += 1;
                                if data.pf {
                                    t.pf_reads += 1;
                                }
                                if m.mismatches == 0 {
                                    t.perfect += 1;
                                    if data.pf {
                                        t.pf_perfect += 1;
                                    }
                                } else if m.mismatches == 1 {
                                    t.one += 1;
                                    if data.pf {
                                        t.pf_one += 1;
                                    }
                                }
                            }
                        }
                        data.barcode = m.matched.then_some(m.barcode);
                    }
                } else {
                    data.barcode = None;
                }
                out.push(data);
            }
        }
        out
    }
}

/// The record of a cluster whose barcode has no writer when unexpected barcodes are not ignored.
/// The sorted converter meets it on a tile thread, so it surfaces as that executor's failure.
pub fn unexpected(key: &Option<String>, sorted: bool) -> ! {
    let message = format!(
        "picard.PicardException: Read records with barcode {}, but this barcode was not expected.  (Is it referenced in the parameters file?)",
        key.as_deref().unwrap_or("null")
    );
    if sorted {
        thrown(&format!(
            "picard.PicardException: Exceptions in tile processing. There were 0 tasks still running or queued and they have been cancelled. Errors: {message}"
        ));
    }
    thrown(&message);
}
