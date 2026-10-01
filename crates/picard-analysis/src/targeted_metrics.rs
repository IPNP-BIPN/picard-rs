//! The engine behind `CollectHsMetrics` and `CollectTargetedPcrMetrics`: one pass over the
//! records, one `PerUnitTargetMetricCollector` per accumulation unit, and the `TargetMetrics` row
//! each converts into its tool's metric class.
//!
//! [`crate::collect_hs_metrics`] and [`crate::collect_targeted_pcr_metrics`] hold the partition
//! and the derived arithmetic in isolation. This applies them to a file, and reproduces two
//! things a tidier design would not:
//!
//! * **The depth and quality histograms belong to the collector, not to the unit.** Every unit's
//!   `finish` adds its targets to the same `highQualityDepthHistogram` and
//!   `unfilteredDepthHistogram`, so the second unit's median and sensitivity are computed over
//!   both units' bases, and its own consistency check -- "numbers of bases in the base quality
//!   histogram and the coverage histogram are not equal" -- compares its own quality counts
//!   against both. Any run with two units that saw on-target bases ends in that exception.
//! * **The metrics file is given the shared histograms once per unit**, so a run that survives
//!   with N units writes each histogram column N times.
//!
//! Ported from Picard 3.4.0 `TargetMetricsCollector` (and its `PerUnitTargetMetricCollector` and
//! `Coverage`), `MultiLevelCollector`, `HsMetricCollector.convertMetric`/`calculateHsPenalty`,
//! `TargetedPcrMetricsCollector.convertMetric`, and htsjdk 4.2.0
//! `SAMUtils.getNumOverlappingAlignedBasesToClip`/`clipOverlappingAlignedBases`.

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::cigar::{soft_clip_end_of_read, Cigar, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_metrics::histogram::Histogram;

use crate::adapter::AdapterUtility;
use crate::mark_duplicates::{estimate_library_size, estimate_roi};
use crate::theoretical_sensitivity::{het_snp_sensitivity, java_round, normalize};

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;
const SECONDARY: u16 = 0x100;
const VENDOR_FAILED: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

const LOG_ODDS_THRESHOLD: f64 = 3.0;

/// A Java exception, as `Exception in thread "main"` prints it.
pub type Thrown = String;

/// `MetricAccumulationLevel`, in the order `MultiLevelCollector.setup` lays the levels out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    AllReads,
    Sample,
    Library,
    ReadGroup,
}

/// One `@RG`, in the three attributes the levels key on.
#[derive(Debug, Clone)]
pub struct ReadGroup {
    pub id: String,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub platform_unit: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Options {
    pub levels: Vec<Level>,
    pub near_distance: i32,
    pub minimum_mapping_quality: i32,
    pub minimum_base_quality: i32,
    pub clip_overlapping_reads: bool,
    pub include_indels: bool,
    pub coverage_cap: i32,
    pub sample_size: i32,
    pub probe_set_name: String,
}

/// `(sequence, start, end)`, 1-based inclusive.
pub type Span = (i32, i32, i32);

/// The panel: uniqued targets and probes, the raw target count, the genome size, and each
/// target's GC when a reference was given.
#[derive(Debug, Clone)]
pub struct Panel {
    pub targets: Vec<Span>,
    pub probes: Vec<Span>,
    pub raw_target_count: usize,
    pub genome_size: i64,
    pub target_gc: Option<Vec<f64>>,
}

/// `TargetMetrics`, before its conversion.
#[derive(Debug, Clone, Default)]
pub struct TargetMetrics {
    pub probe_set: Option<String>,
    pub probe_territory: i64,
    pub on_probe_bases: i64,
    pub near_probe_bases: i64,
    pub off_probe_bases: i64,
    pub pct_selected_bases: f64,
    pub pct_off_probe: f64,
    pub on_probe_vs_selected: f64,
    pub mean_probe_coverage: f64,
    pub fold_enrichment: f64,
    // TargetMetricsBase
    pub pf_selected_pairs: i64,
    pub pf_selected_unique_pairs: i64,
    pub on_target_from_pair_bases: i64,
    // PanelMetricsBase
    pub target_territory: i64,
    pub genome_size: i64,
    pub total_reads: i64,
    pub pf_reads: i64,
    pub pf_bases: i64,
    pub pf_unique_reads: i64,
    pub pf_uq_reads_aligned: i64,
    pub pf_bases_aligned: i64,
    pub pf_uq_bases_aligned: i64,
    pub on_target_bases: i64,
    pub pct_pf_reads: f64,
    pub pct_pf_uq_reads: f64,
    pub pct_pf_uq_reads_aligned: f64,
    pub mean_target_coverage: f64,
    pub median_target_coverage: f64,
    pub max_target_coverage: i64,
    pub min_target_coverage: i64,
    pub zero_cvg_targets_pct: f64,
    pub pct_exc_dupe: f64,
    pub pct_exc_adapter: f64,
    pub pct_exc_mapq: f64,
    pub pct_exc_baseq: f64,
    pub pct_exc_overlap: f64,
    pub pct_exc_off_target: f64,
    pub fold_80_base_penalty: f64,
    /// `PCT_TARGET_BASES_1X` .. `PCT_TARGET_BASES_100000X`, seventeen of them.
    pub pct_target_bases: Vec<f64>,
    pub at_dropout: f64,
    pub gc_dropout: f64,
    pub het_snp_sensitivity: f64,
    pub het_snp_q: f64,
    // MultilevelMetrics
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

struct Coverage {
    depths: Vec<i32>,
    read_count: i64,
}

impl Coverage {
    fn new(length: i32) -> Self {
        Coverage {
            depths: vec![0; length.max(0) as usize],
            read_count: 0,
        }
    }

    fn add_base(&mut self, offset: i32) {
        if offset >= 0
            && (offset as usize) < self.depths.len()
            && self.depths[offset as usize] < i32::MAX - 1
        {
            self.depths[offset as usize] += 1;
        }
    }

    fn has_coverage(&self) -> bool {
        self.depths.iter().any(|&d| d > 0)
    }

    fn total(&self) -> i64 {
        let mut total: i64 = 0;
        for &d in &self.depths {
            let d = d as i64;
            total += if total < i64::MAX - d {
                d
            } else {
                i64::MAX - total
            };
        }
        total
    }
}

struct Unit {
    metrics: TargetMetrics,
    base_q: [i64; 127],
    uncapped_base_q: [i64; 127],
    high_quality: Vec<Coverage>,
    unfiltered: Vec<Coverage>,
    hq_max_depth: i64,
    adapter_bases: i64,
    mapq_bases: i64,
}

/// The histograms the units share.
struct Shared {
    high_quality_depth: Histogram,
    unfiltered_depth: Histogram,
    unfiltered_base_q: Histogram,
    uncapped_base_q: Histogram,
}

fn aligned_bases(cigar: &Cigar, start: i32) -> i64 {
    alignment_blocks(cigar, start)
        .iter()
        .map(|b| b.length as i64)
        .sum()
}

fn reference_length(cigar: &Cigar) -> i32 {
    cigar.reference_length() as i32
}

/// `SAMUtils.getNumOverlappingAlignedBasesToClip`.
fn overlapping_bases_to_clip(record: &BamRecord) -> i32 {
    if record.flags & PAIRED == 0
        || record.flags & UNMAPPED != 0
        || record.flags & MATE_UNMAPPED != 0
    {
        return 0;
    }
    // Only the left-most end of an overlapping pair is clipped; on a shared start, the first end.
    if record.mate_alignment_start < record.alignment_start
        || (record.mate_alignment_start == record.alignment_start
            && record.flags & FIRST_OF_PAIR != 0)
    {
        return 0;
    }
    let mut clip = 0;
    let ref_start = record.mate_alignment_start;
    let mut ref_pos = record.alignment_start;
    for e in &record.cigar.elements {
        let ref_len = if e.op.consumes_reference_bases() {
            e.length as i32
        } else {
            0
        };
        // `refStartPos <= refPos + refBasesLength - 1`
        if ref_start < ref_pos + ref_len {
            match e.op {
                Op::M => {
                    if ref_start < ref_pos {
                        clip += ref_len;
                    } else {
                        clip += (ref_pos + ref_len) - ref_start;
                    }
                }
                Op::S | Op::H | Op::P | Op::N => {}
                _ => {
                    if e.op.consumes_read_bases() {
                        clip += e.length as i32;
                    }
                }
            }
        }
        ref_pos += ref_len;
    }
    clip.max(0)
}

/// The record a unit walks after `clipOverlappingAlignedBases`, or `None` when the clip made it
/// unmapped (the mate starts at or before it).
fn clipped_cigar(record: &BamRecord, to_clip: i32) -> Option<Cigar> {
    if to_clip <= 0 || record.flags & UNMAPPED != 0 || record.flags & MATE_UNMAPPED != 0 {
        return Some(record.cigar.clone());
    }
    if record.mate_alignment_start <= record.alignment_start {
        return None;
    }
    let mut clip_from = record.read_length() as i32 - to_clip + 1;
    if let Some(last) = record.cigar.elements.last() {
        if last.op == Op::S {
            clip_from -= last.length as i32;
        }
    }
    Some(Cigar::new(soft_clip_end_of_read(
        clip_from,
        &record.cigar.elements,
    )))
}

impl Unit {
    fn new(
        panel: &Panel,
        options: &Options,
        sample: Option<String>,
        library: Option<String>,
        read_group: Option<String>,
    ) -> Self {
        let territory = |spans: &[Span]| spans.iter().map(|s| (s.2 - s.1 + 1) as i64).sum();
        let metrics = TargetMetrics {
            probe_set: Some(options.probe_set_name.clone()),
            probe_territory: territory(&panel.probes),
            target_territory: territory(&panel.targets),
            genome_size: panel.genome_size,
            sample,
            library,
            read_group,
            ..TargetMetrics::default()
        };
        let coverage = || {
            panel
                .targets
                .iter()
                .map(|t| Coverage::new(t.2 - t.1 + 1))
                .collect::<Vec<_>>()
        };
        Unit {
            metrics,
            base_q: [0; 127],
            uncapped_base_q: [0; 127],
            high_quality: coverage(),
            unfiltered: coverage(),
            hq_max_depth: 0,
            adapter_bases: 0,
            mapq_bases: 0,
        }
    }

    fn accept(
        &mut self,
        record: &BamRecord,
        panel: &Panel,
        options: &Options,
        adapter: &AdapterUtility,
    ) {
        let flags = record.flags;
        if flags & SECONDARY != 0 {
            return;
        }
        let unmapped = flags & UNMAPPED != 0;
        let supplementary = flags & SUPPLEMENTARY != 0;
        let duplicate = flags & DUPLICATE != 0;
        let mapped_in_pair =
            flags & PAIRED != 0 && !unmapped && flags & MATE_UNMAPPED == 0 && !supplementary;
        let qualities = &record.base_qualities;
        let bases_aligned = if unmapped {
            0
        } else {
            aligned_bases(&record.cigar, record.alignment_start)
        };
        let m = &mut self.metrics;
        if !supplementary {
            m.total_reads += 1;
            if flags & VENDOR_FAILED == 0 {
                m.pf_reads += 1;
                if !duplicate {
                    m.pf_unique_reads += 1;
                    if !unmapped {
                        m.pf_uq_reads_aligned += 1;
                    }
                }
            }
        }
        if flags & VENDOR_FAILED != 0 {
            return;
        }
        if !supplementary {
            m.pf_bases += record.read_length() as i64;
        }
        if !unmapped {
            m.pf_bases_aligned += bases_aligned;
            if !duplicate {
                m.pf_uq_bases_aligned += bases_aligned;
            }
        }
        if unmapped {
            return;
        }

        let read_start = record.alignment_start;
        let read_end = record.alignment_start + reference_length(&record.cigar) - 1;
        let sequence = record.reference_index;
        let targets: Vec<usize> = panel
            .targets
            .iter()
            .enumerate()
            .filter(|(_, t)| t.0 == sequence && t.1 <= read_end && t.2 >= read_start)
            .map(|(i, _)| i)
            .collect();
        let near = options.near_distance;
        let probes: Vec<&Span> = panel
            .probes
            .iter()
            .filter(|p| p.0 == sequence && p.1 - near <= read_end && p.2 + near >= read_start)
            .collect();

        if !supplementary
            && flags & PAIRED != 0
            && flags & FIRST_OF_PAIR != 0
            && flags & MATE_UNMAPPED == 0
            && !probes.is_empty()
        {
            m.pf_selected_pairs += 1;
            if !duplicate {
                m.pf_selected_unique_pairs += 1;
            }
        }

        if !probes.is_empty() {
            let mut on_bait: i64 = 0;
            for bait in &probes {
                for block in alignment_blocks(&record.cigar, record.alignment_start) {
                    let end = block.reference_start + block.length - 1;
                    for pos in block.reference_start..=end {
                        if pos >= bait.1 && pos <= bait.2 {
                            on_bait += 1;
                        }
                    }
                }
            }
            m.on_probe_bases += on_bait;
            m.near_probe_bases += bases_aligned - on_bait;
        } else {
            m.off_probe_bases += bases_aligned;
        }

        if duplicate {
            m.pct_exc_dupe += bases_aligned as f64;
            return;
        }
        if adapter.is_adapter(record) {
            self.adapter_bases += bases_aligned_any(record);
            return;
        }
        if (record.mapping_quality as i32) < options.minimum_mapping_quality {
            self.mapq_bases += bases_aligned_any(record);
            return;
        }

        let cigar = if options.clip_overlapping_reads {
            let to_clip = overlapping_bases_to_clip(record);
            let clipped = clipped_cigar(record, to_clip);
            self.metrics.pct_exc_overlap += to_clip as f64;
            match clipped {
                Some(c) => c,
                None => return,
            }
        } else {
            record.cigar.clone()
        };

        let m = &mut self.metrics;
        let mut covered_targets: Vec<usize> = Vec::new();
        let mut read_offset: usize = 0;
        let mut ref_pos = record.alignment_start;
        for e in &cigar.elements {
            let op = e.op;
            let alignment = matches!(op, Op::M | Op::Eq | Op::X);
            let indel = matches!(op, Op::I | Op::D);
            for _ in 0..e.length {
                if alignment || (options.include_indels && indel) {
                    let qual = qualities[read_offset] as i8 as i32;
                    let high_qual = qual >= options.minimum_base_quality;
                    let on_target = targets.iter().any(|&t| {
                        let t = panel.targets[t];
                        ref_pos >= t.1 && ref_pos <= t.2
                    });
                    let per_target = op != Op::I;
                    if !high_qual {
                        m.pct_exc_baseq += 1.0;
                    } else if !on_target {
                        m.pct_exc_off_target += 1.0;
                    } else {
                        m.on_target_bases += 1;
                        if mapped_in_pair {
                            m.on_target_from_pair_bases += 1;
                        }
                    }
                    if qual > 2 && per_target && on_target {
                        for &t in &targets {
                            let span = panel.targets[t];
                            if ref_pos >= span.1 && ref_pos <= span.2 {
                                let offset = ref_pos - span.1;
                                let uf = &mut self.unfiltered[t];
                                uf.add_base(offset);
                                if uf.depths[offset as usize] <= options.coverage_cap {
                                    self.base_q[qual as usize] += 1;
                                }
                                self.uncapped_base_q[qual as usize] += 1;
                                if high_qual {
                                    let hq = &mut self.high_quality[t];
                                    hq.add_base(offset);
                                    self.hq_max_depth =
                                        self.hq_max_depth.max(hq.depths[offset as usize] as i64);
                                    if !covered_targets.contains(&t) {
                                        covered_targets.push(t);
                                        hq.read_count += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                if op.consumes_read_bases() {
                    read_offset += 1;
                }
                if op.consumes_reference_bases() {
                    ref_pos += 1;
                }
            }
        }
    }

    fn finish(
        &mut self,
        panel: &Panel,
        options: &Options,
        shared: &mut Shared,
    ) -> Result<(), Thrown> {
        let m = &mut self.metrics;
        m.pct_pf_reads = m.pf_reads as f64 / m.total_reads as f64;
        m.pct_pf_uq_reads = m.pf_unique_reads as f64 / m.total_reads as f64;
        m.pct_pf_uq_reads_aligned = m.pf_uq_reads_aligned as f64 / m.pf_unique_reads as f64;
        let denominator = (m.on_probe_bases + m.near_probe_bases + m.off_probe_bases) as f64;
        m.pct_selected_bases = (m.on_probe_bases + m.near_probe_bases) as f64 / denominator;
        m.pct_off_probe = m.off_probe_bases as f64 / denominator;
        m.on_probe_vs_selected =
            m.on_probe_bases as f64 / (m.on_probe_bases + m.near_probe_bases) as f64;
        m.mean_probe_coverage = m.on_probe_bases as f64 / m.probe_territory as f64;
        m.fold_enrichment = (m.on_probe_bases as f64 / denominator)
            / (m.probe_territory as f64 / m.genome_size as f64);
        let aligned = m.pf_bases_aligned as f64;
        m.pct_exc_dupe /= aligned;
        m.pct_exc_adapter = self.adapter_bases as f64 / aligned;
        m.pct_exc_mapq = self.mapq_bases as f64 / aligned;
        m.pct_exc_baseq /= aligned;
        m.pct_exc_overlap /= aligned;
        m.pct_exc_off_target /= aligned;

        // calculateTargetCoverageMetrics
        let hq = &mut shared.high_quality_depth;
        for i in 0..self.hq_max_depth {
            hq.increment_by(i as f64, 0.0);
        }
        let mut zero_targets = 0;
        let mut total_coverage: i64 = 0;
        let mut min_depth = i64::MAX;
        const DEPTHS: [i64; 18] = [
            0, 1, 2, 10, 20, 30, 40, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000,
            100000,
        ];
        let mut target_bases = [0i64; 18];
        for c in &self.high_quality {
            if !c.has_coverage() {
                zero_targets += 1;
                hq.increment_by(0.0, c.depths.len() as f64);
                target_bases[0] += c.depths.len() as i64;
                min_depth = 0;
                continue;
            }
            for &depth in &c.depths {
                let depth = depth as i64;
                total_coverage += depth;
                hq.increment_by(depth as f64, 1.0);
                min_depth = min_depth.min(depth);
                for (i, &threshold) in DEPTHS.iter().enumerate() {
                    if depth >= threshold {
                        target_bases[i] += 1;
                    } else {
                        break;
                    }
                }
            }
        }
        m.mean_target_coverage = total_coverage as f64 / m.target_territory as f64;
        m.median_target_coverage = hq.median();
        m.max_target_coverage = self.hq_max_depth;
        m.min_target_coverage = min_depth.min(self.hq_max_depth);
        let p20 = hq.percentile(0.2).map_err(|_| {
            "java.lang.IllegalStateException: Cannot calculate percentiles when total is zero."
                .to_string()
        })?;
        m.fold_80_base_penalty = m.mean_target_coverage / p20;
        m.zero_cvg_targets_pct = zero_targets as f64 / panel.raw_target_count as f64;
        m.pct_target_bases = (1..18)
            .map(|i| target_bases[i] as f64 / target_bases[0] as f64)
            .collect();

        // calculateTheoreticalHetSensitivity
        let cap = options.coverage_cap;
        let uf = &mut shared.unfiltered_depth;
        for i in 0..cap.max(0) {
            uf.increment_by(i as f64, 0.0);
        }
        for c in &self.unfiltered {
            if !c.has_coverage() {
                uf.increment_by(0.0, c.depths.len() as f64);
                continue;
            }
            for &depth in &c.depths {
                uf.increment_by(depth.min(cap) as f64, 1.0);
            }
        }
        let base_q_total: i64 = self.base_q.iter().sum();
        if base_q_total as f64 != uf.sum() {
            return Err("picard.PicardException: numbers of bases in the base quality histogram and the coverage histogram are not equal".to_string());
        }
        for (i, &v) in self.base_q.iter().enumerate() {
            shared.unfiltered_base_q.increment_by(i as f64, v as f64);
        }
        for (i, &v) in self.uncapped_base_q.iter().enumerate() {
            shared.uncapped_base_q.increment_by(i as f64, v as f64);
        }
        let depth_values = histogram_as_array(uf);
        let base_q_values = histogram_as_array(&shared.unfiltered_base_q);
        let sensitivity = het_snp_sensitivity(
            &normalize(&depth_values),
            &normalize(&base_q_values),
            options.sample_size,
            LOG_ODDS_THRESHOLD,
        )?;
        m.het_snp_sensitivity = sensitivity;
        m.het_snp_q =
            crate::theoretical_sensitivity::phred_from_error_probability(1.0 - sensitivity) as f64;

        // calculateGcMetrics
        if let Some(gc) = &panel.target_gc {
            let mut target_by_gc = [0i64; 101];
            let mut aligned_by_gc = [0i64; 101];
            for (t, c) in self.high_quality.iter().enumerate() {
                let length = c.depths.len() as i64;
                if length <= 0 {
                    continue;
                }
                let bin = java_round(gc[t] * 100.0) as i32;
                target_by_gc[bin as usize] += length;
                aligned_by_gc[bin as usize] += c.total();
            }
            let total_target: i64 = target_by_gc.iter().sum();
            let total_bases: i64 = aligned_by_gc.iter().sum();
            for i in 0..101 {
                let target_pct = target_by_gc[i] as f64 / total_target as f64;
                let aligned_pct = aligned_by_gc[i] as f64 / total_bases as f64;
                let mut dropout = (aligned_pct - target_pct) * 100.0;
                if dropout < 0.0 {
                    dropout = dropout.abs();
                    if i <= 50 {
                        m.at_dropout += dropout;
                    }
                    if i >= 50 {
                        m.gc_dropout += dropout;
                    }
                }
            }
        }
        Ok(())
    }
}

/// `CountingFilter`'s base count: the record's alignment blocks, read from its cigar whatever
/// its flags say.
fn bases_aligned_any(record: &BamRecord) -> i64 {
    aligned_bases(&record.cigar, record.alignment_start)
}

/// `normalizeHistogram`'s view of a histogram: `get(i)` for `i` in `0..size()`, a missing key
/// reading as zero.
fn histogram_as_array(h: &Histogram) -> Vec<f64> {
    (0..h.size())
        .map(|i| h.get(i as f64).unwrap_or(0.0))
        .collect()
}

/// What a run produces: one converted row per unit, in output order, and the two shared
/// histograms as they stand at the end (each written once per unit).
pub struct Outcome {
    pub rows: Vec<TargetMetrics>,
    pub high_quality_depth: Vec<(i64, f64)>,
    pub uncapped_base_q: Vec<(i64, f64)>,
}

struct Distributor {
    level: Level,
    keys: Vec<Option<String>>,
    units: Vec<Unit>,
}

fn key_of(level: Level, group: &ReadGroup) -> Option<String> {
    match level {
        Level::AllReads => None,
        Level::Sample => group.sample.clone(),
        Level::Library => group.library.clone(),
        Level::ReadGroup => group.platform_unit.clone(),
    }
}

/// Runs the collector over `records`. `group_of` gives each record's `@RG`, if it has one the
/// header knows.
pub fn collect(
    records: &[BamRecord],
    group_of: impl Fn(&BamRecord) -> Option<usize>,
    groups: &[ReadGroup],
    panel: &Panel,
    options: &Options,
) -> Result<Outcome, Thrown> {
    let adapter = AdapterUtility::with_defaults();
    let mut distributors: Vec<Distributor> = Vec::new();
    for level in [
        Level::AllReads,
        Level::Sample,
        Level::Library,
        Level::ReadGroup,
    ] {
        if !options.levels.contains(&level) {
            continue;
        }
        let mut d = Distributor {
            level,
            keys: Vec::new(),
            units: Vec::new(),
        };
        if level == Level::AllReads {
            d.keys.push(None);
            d.units.push(Unit::new(panel, options, None, None, None));
        } else {
            for g in groups {
                let key = key_of(level, g);
                if d.keys.contains(&key) {
                    continue;
                }
                let (sample, library, read_group) = match level {
                    Level::Sample => (g.sample.clone(), None, None),
                    Level::Library => (g.sample.clone(), g.library.clone(), None),
                    _ => (g.sample.clone(), g.library.clone(), g.platform_unit.clone()),
                };
                d.keys.push(key);
                d.units
                    .push(Unit::new(panel, options, sample, library, read_group));
            }
        }
        distributors.push(d);
    }

    for record in records {
        let group = group_of(record).map(|i| &groups[i]);
        for d in distributors.iter_mut() {
            if d.level == Level::AllReads {
                d.units[0].accept(record, panel, options, &adapter);
                continue;
            }
            let key = group
                .and_then(|g| key_of(d.level, g))
                .unwrap_or_else(|| "unknown".to_string());
            let at = d
                .keys
                .iter()
                .position(|k| k.as_deref() == Some(key.as_str()));
            let at = match at {
                Some(at) => at,
                None => {
                    if key != "unknown" {
                        return Err(format!(
                            "picard.PicardException: Could not find collector for {key}"
                        ));
                    }
                    let unknown = Some("unknown".to_string());
                    let unit = match d.level {
                        Level::Sample => Unit::new(panel, options, unknown, None, None),
                        Level::Library => Unit::new(panel, options, unknown.clone(), unknown, None),
                        _ => Unit::new(panel, options, unknown.clone(), unknown.clone(), unknown),
                    };
                    d.keys.push(Some(key));
                    d.units.push(unit);
                    d.units.len() - 1
                }
            };
            d.units[at].accept(record, panel, options, &adapter);
        }
    }

    let mut shared = Shared {
        high_quality_depth: Histogram::new(
            "coverage_or_base_quality",
            "high_quality_coverage_count",
        ),
        unfiltered_depth: Histogram::new("coverage_or_base_quality", "unfiltered_coverage_count"),
        unfiltered_base_q: Histogram::new("baseq", "unfiltered_baseq_count"),
        uncapped_base_q: Histogram::new("baseq", "unfiltered_baseq_count"),
    };
    for d in distributors.iter_mut() {
        for unit in d.units.iter_mut() {
            unit.finish(panel, options, &mut shared)?;
        }
    }
    let rows = distributors
        .into_iter()
        .flat_map(|d| d.units.into_iter().map(|u| u.metrics))
        .collect();
    let bins = |h: &Histogram| h.bins().map(|(k, v)| (k as i64, v)).collect::<Vec<_>>();
    Ok(Outcome {
        rows,
        high_quality_depth: bins(&shared.high_quality_depth),
        uncapped_base_q: bins(&shared.uncapped_base_q),
    })
}

/// `HsMetricCollector.calculateHsPenalty`.
pub fn hs_penalty(library_size: Option<i64>, m: &TargetMetrics, goal: i32) -> f64 {
    let Some(size) = library_size else {
        return 0.0;
    };
    let mean_coverage = m.on_target_from_pair_bases as f64 / m.target_territory as f64;
    let fold80 = m.fold_80_base_penalty;
    let pairs = m.pf_selected_pairs;
    let unique_pairs = m.pf_selected_unique_pairs;
    let on_target_pct = m.on_target_bases as f64 / m.pf_uq_bases_aligned as f64;
    let goal_multiplier = (goal as f64 / mean_coverage) * fold80;
    let mut pair_multiplier = goal_multiplier;
    let mut increment = 1.0;
    let mut going_up = goal_multiplier >= 1.0;
    let mut final_multiplier = -1.0;
    for _ in 0..10000 {
        let unique_multiplier = estimate_roi(size, pair_multiplier, pairs, unique_pairs);
        if (unique_multiplier - goal_multiplier).abs() / goal_multiplier <= 0.001 {
            final_multiplier = pair_multiplier;
            break;
        } else if (unique_multiplier > goal_multiplier && going_up)
            || (unique_multiplier < goal_multiplier && !going_up)
        {
            increment /= 2.0;
            going_up = !going_up;
        }
        pair_multiplier += if going_up { increment } else { -increment };
    }
    if final_multiplier == -1.0 {
        -1.0
    } else {
        let unique_fraction =
            (unique_pairs as f64 * goal_multiplier) / (pairs as f64 * final_multiplier);
        (1.0 / unique_fraction) * fold80 * (1.0 / on_target_pct)
    }
}

/// `DuplicationMetrics.estimateLibrarySize`, re-exported for the HS conversion.
pub fn library_size(m: &TargetMetrics) -> Option<i64> {
    estimate_library_size(m.pf_selected_pairs, m.pf_selected_unique_pairs)
}

/// `SequenceUtil.calculateGc`.
pub fn calculate_gc(bases: &[u8]) -> f64 {
    let gcs = bases
        .iter()
        .filter(|&&b| matches!(b, b'C' | b'G' | b'c' | b'g'))
        .count();
    gcs as f64 / bases.len() as f64
}
