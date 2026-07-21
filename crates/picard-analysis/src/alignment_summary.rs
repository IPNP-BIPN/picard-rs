//! `CollectAlignmentSummaryMetrics`.
//!
//! Ported from `picard.analysis.CollectAlignmentSummaryMetrics`, `AlignmentSummaryMetrics` and
//! `AlignmentSummaryMetricsCollector`, tag 3.4.0, together with `picard.analysis.ChimeraUtil`
//! and `picard.util.MathUtil.divide`.
//!
//! Second member of the `histogram` + `multi_level` + `single_pass` stratum, after
//! `CollectInsertSizeMetrics`, and the tool the archetype delta is measured against at the large
//! end of the size distribution.
//!
//! ## The divergence this port reproduces on purpose
//!
//! `collectQualityData` walks the read's alignment blocks and, on a mismatch, does
//!
//! ```java
//! badCycleHistogram.increment(CoordMath.getCycle(negativeStrand, readBases.length, i));
//! ```
//!
//! where `i` is the offset **within the current alignment block**, while `getCycle`'s third
//! parameter is declared `final int readBaseIndex`. Every other call in the method passes a read
//! index. Two mismatched bases in different blocks that share a block offset therefore land in
//! the same cycle bin, and the read is charged one bad cycle instead of two.
//!
//! Measured in the pinned oracle by `tools/asm-conformance/badcycle_probe.sh`, on one 20-base
//! read with CIGAR `10M5D10M` and two mismatches:
//!
//! ```text
//! collide   read indices 3 and 13, block offsets 3 and 3   BAD_CYCLES=1
//! distinct  read indices 3 and 15, block offsets 3 and 5   BAD_CYCLES=2
//! ```
//!
//! `PF_HQ_MEDIAN_MISMATCHES` is 2 in both, which is the control: both mismatches are counted in
//! both runs, so what differs is the binning.
//!
//! Any read carrying an indel or a splice has more than one block, so this reaches most aligned
//! BAMs. Passing the read index instead - which is what the parameter name asks for, and what a
//! careful reimplementation writes - gives a different `BAD_CYCLES` for real data. It is
//! reproduced here, and the conformance test would fail if it were not.
//!
//! The *unaligned* path a few lines above, which counts no-calls over the whole read, passes a
//! genuine read index. So the same histogram is fed in two different coordinate systems
//! depending on whether the read aligned.
//!
//! ## A comment in the source that is not true of the source
//!
//! `collectReadData` and `collectQualityData` both open with
//!
//! ```text
//! // NB: for read count metrics, do not include supplementary records, but for base count
//! // metrics, do include supplementary records.
//! ```
//!
//! and `collectReadData` returns early on the supplementary flag to implement the first half.
//! The second half never happens. `AlignmentSummaryMetricsCollector.acceptRecord` overrides its
//! parent with
//!
//! ```java
//! if (!rec.isSecondaryOrSupplementary()) {
//!     super.acceptRecord(rec, ref);
//! }
//! ```
//!
//! and `isSecondaryOrSupplementary()` is secondary **or** supplementary, so a supplementary
//! record is dropped before any per-unit collector sees it. The guard inside `collectReadData`
//! and the `!getSupplementaryAlignmentFlag()` test inside `collectQualityData` are unreachable
//! through this tool.
//!
//! Confirmed against the oracle rather than inferred: the `supplementary` corpus case is one
//! ordinary 20-base read plus one supplementary 20-base read, and Picard reports
//! `PF_ALIGNED_BASES = 20`. If the comment were true it would be 40.
//!
//! Both guards are reproduced anyway, because the collector is public and `IndividualCollector`
//! can be driven directly, but the comment is not repeated as though it described behaviour.

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{bases_equal, is_no_call};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, Value};
use htsjdk_metrics::histogram::Histogram;

use crate::adapter::AdapterUtility;
use crate::insert_size::{pair_orientation, PairOrientation};

const READ_PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const FIRST_OF_PAIR: u16 = 0x40;
const SECONDARY: u16 = 0x100;
const VENDOR_FAILED: u16 = 0x200;
const SUPPLEMENTARY: u16 = 0x800;

/// `MAPPING_QUALITY_THRESHOLD`.
const MAPPING_QUALITY_THRESHOLD: u8 = 20;
/// `BASE_QUALITY_THRESHOLD`.
const BASE_QUALITY_THRESHOLD: u8 = 20;
/// `ChimeraUtil.DEFAULT_INSERT_SIZE_LIMIT`.
pub const DEFAULT_INSERT_SIZE_LIMIT: i32 = 100_000;

/// `MathUtil.divide`.
///
/// Not a plain division and not a zero guard either: it returns 0 whenever the denominator's
/// **magnitude** is at most 1e-6. So a genuine ratio with a denominator of 5e-7 comes out as 0
/// rather than as a large number. The denominators here are counts, so the threshold only ever
/// fires at exactly zero - but the rule is the rule, and writing `if d == 0.0` would be a
/// different function.
pub fn divide(numerator: f64, denominator: f64) -> f64 {
    if (0.0 - denominator).abs() > 0.000_001 {
        numerator / denominator
    } else {
        0.0
    }
}

/// `AlignmentSummaryMetrics.Category`, in declaration order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    Unpaired,
    FirstOfPair,
    SecondOfPair,
    Pair,
}

impl Category {
    pub fn name(self) -> &'static str {
        match self {
            Category::Unpaired => "UNPAIRED",
            Category::FirstOfPair => "FIRST_OF_PAIR",
            Category::SecondOfPair => "SECOND_OF_PAIR",
            Category::Pair => "PAIR",
        }
    }
}

/// `AlignmentSummaryMetrics`, in declared field order, with the inherited `MultilevelMetrics`
/// fields last. Same rule as `InsertSizeMetrics`: the column order is what HotSpot's
/// `Class.getFields()` happens to return, which is declaration order then inherited.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentSummaryMetrics {
    pub category: Category,
    pub total_reads: i64,
    pub pf_reads: i64,
    pub pct_pf_reads: f64,
    pub pf_noise_reads: i64,
    pub pf_reads_aligned: i64,
    pub pct_pf_reads_aligned: f64,
    pub pf_aligned_bases: i64,
    pub pf_hq_aligned_reads: i64,
    pub pf_hq_aligned_bases: i64,
    pub pf_hq_aligned_q20_bases: i64,
    pub pf_hq_median_mismatches: f64,
    pub pf_mismatch_rate: f64,
    pub pf_hq_error_rate: f64,
    pub pf_indel_rate: f64,
    pub mean_read_length: f64,
    pub sd_read_length: f64,
    pub median_read_length: f64,
    pub mad_read_length: f64,
    pub min_read_length: f64,
    pub max_read_length: f64,
    pub mean_aligned_read_length: f64,
    pub reads_aligned_in_pairs: i64,
    pub pct_reads_aligned_in_pairs: f64,
    pub pf_reads_improper_pairs: i64,
    pub pct_pf_reads_improper_pairs: f64,
    pub bad_cycles: i64,
    pub strand_balance: f64,
    pub pct_chimeras: f64,
    pub pct_adapter: f64,
    pub pct_softclip: f64,
    pub pct_hardclip: f64,
    pub avg_pos_3prime_softclip_length: f64,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

impl AlignmentSummaryMetrics {
    fn new(category: Category) -> Self {
        AlignmentSummaryMetrics {
            category,
            total_reads: 0,
            pf_reads: 0,
            pct_pf_reads: 0.0,
            pf_noise_reads: 0,
            pf_reads_aligned: 0,
            pct_pf_reads_aligned: 0.0,
            pf_aligned_bases: 0,
            pf_hq_aligned_reads: 0,
            pf_hq_aligned_bases: 0,
            pf_hq_aligned_q20_bases: 0,
            pf_hq_median_mismatches: 0.0,
            pf_mismatch_rate: 0.0,
            pf_hq_error_rate: 0.0,
            pf_indel_rate: 0.0,
            mean_read_length: 0.0,
            sd_read_length: 0.0,
            median_read_length: 0.0,
            mad_read_length: 0.0,
            min_read_length: 0.0,
            max_read_length: 0.0,
            mean_aligned_read_length: 0.0,
            reads_aligned_in_pairs: 0,
            pct_reads_aligned_in_pairs: 0.0,
            pf_reads_improper_pairs: 0,
            pct_pf_reads_improper_pairs: 0.0,
            bad_cycles: 0,
            strand_balance: 0.0,
            pct_chimeras: 0.0,
            pct_adapter: 0.0,
            pct_softclip: 0.0,
            pct_hardclip: 0.0,
            avg_pos_3prime_softclip_length: 0.0,
            sample: None,
            library: None,
            read_group: None,
        }
    }
}

const COLUMNS: &[&str] = &[
    "CATEGORY",
    "TOTAL_READS",
    "PF_READS",
    "PCT_PF_READS",
    "PF_NOISE_READS",
    "PF_READS_ALIGNED",
    "PCT_PF_READS_ALIGNED",
    "PF_ALIGNED_BASES",
    "PF_HQ_ALIGNED_READS",
    "PF_HQ_ALIGNED_BASES",
    "PF_HQ_ALIGNED_Q20_BASES",
    "PF_HQ_MEDIAN_MISMATCHES",
    "PF_MISMATCH_RATE",
    "PF_HQ_ERROR_RATE",
    "PF_INDEL_RATE",
    "MEAN_READ_LENGTH",
    "SD_READ_LENGTH",
    "MEDIAN_READ_LENGTH",
    "MAD_READ_LENGTH",
    "MIN_READ_LENGTH",
    "MAX_READ_LENGTH",
    "MEAN_ALIGNED_READ_LENGTH",
    "READS_ALIGNED_IN_PAIRS",
    "PCT_READS_ALIGNED_IN_PAIRS",
    "PF_READS_IMPROPER_PAIRS",
    "PCT_PF_READS_IMPROPER_PAIRS",
    "BAD_CYCLES",
    "STRAND_BALANCE",
    "PCT_CHIMERAS",
    "PCT_ADAPTER",
    "PCT_SOFTCLIP",
    "PCT_HARDCLIP",
    "AVG_POS_3PRIME_SOFTCLIP_LENGTH",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

impl MetricBean for AlignmentSummaryMetrics {
    fn class_name(&self) -> &str {
        "picard.analysis.AlignmentSummaryMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }

    fn values(&self) -> Vec<Value> {
        let text = |o: &Option<String>| match o {
            Some(s) => Value::Str(s.clone()),
            None => Value::Null,
        };
        vec![
            Value::Str(self.category.name().to_string()),
            Value::Long(self.total_reads),
            Value::Long(self.pf_reads),
            Value::Double(self.pct_pf_reads),
            Value::Long(self.pf_noise_reads),
            Value::Long(self.pf_reads_aligned),
            Value::Double(self.pct_pf_reads_aligned),
            Value::Long(self.pf_aligned_bases),
            Value::Long(self.pf_hq_aligned_reads),
            Value::Long(self.pf_hq_aligned_bases),
            Value::Long(self.pf_hq_aligned_q20_bases),
            Value::Double(self.pf_hq_median_mismatches),
            Value::Double(self.pf_mismatch_rate),
            Value::Double(self.pf_hq_error_rate),
            Value::Double(self.pf_indel_rate),
            Value::Double(self.mean_read_length),
            Value::Double(self.sd_read_length),
            Value::Double(self.median_read_length),
            Value::Double(self.mad_read_length),
            Value::Double(self.min_read_length),
            Value::Double(self.max_read_length),
            Value::Double(self.mean_aligned_read_length),
            Value::Long(self.reads_aligned_in_pairs),
            Value::Double(self.pct_reads_aligned_in_pairs),
            Value::Long(self.pf_reads_improper_pairs),
            Value::Double(self.pct_pf_reads_improper_pairs),
            Value::Long(self.bad_cycles),
            Value::Double(self.strand_balance),
            Value::Double(self.pct_chimeras),
            Value::Double(self.pct_adapter),
            Value::Double(self.pct_softclip),
            Value::Double(self.pct_hardclip),
            Value::Double(self.avg_pos_3prime_softclip_length),
            text(&self.sample),
            text(&self.library),
            text(&self.read_group),
        ]
    }
}

/// The tool's arguments, with Picard's defaults.
#[derive(Debug, Clone)]
pub struct Options {
    /// `COLLECT_ALIGNMENT_INFORMATION`, which the collector receives as `doRefMetrics`.
    pub collect_alignment_information: bool,
    pub max_insert_size: i32,
    /// `EXPECTED_PAIR_ORIENTATIONS`, default `{FR}`.
    pub expected_orientations: Vec<PairOrientation>,
    pub is_bisulfite_sequenced: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            collect_alignment_information: true,
            max_insert_size: DEFAULT_INSERT_SIZE_LIMIT,
            expected_orientations: vec![PairOrientation::Fr],
            is_bisulfite_sequenced: false,
        }
    }
}

/// `ChimeraUtil.isChimeric(rec, maxInsertSize, expectedOrientations)`.
///
/// The `SA` tag is consulted inside `matchesExpectedOrientations`, so a read whose orientation is
/// expected but which carries an `SA` tag is still chimeric. That makes the tag a second, hidden
/// term of the orientation test rather than a separate condition.
fn is_chimeric(rec: &BamRecord, max_insert_size: i32, expected: &[PairOrientation]) -> bool {
    let mapped_pair = rec.flags & READ_PAIRED != 0
        && rec.flags & READ_UNMAPPED == 0
        && rec.flags & MATE_UNMAPPED == 0;
    if !mapped_pair {
        return false;
    }
    let has_sa = rec.tags.get(Tag::new(b"SA")).is_some();
    let orientation_ok = expected.contains(&pair_orientation(rec)) && !has_sa;
    rec.inferred_insert_size.abs() > max_insert_size
        || rec.reference_index != rec.mate_reference_index
        || !orientation_ok
}

/// `getTotalCigarOperatorCount`.
fn total_cigar_operator_count(rec: &BamRecord, op: Op) -> i64 {
    rec.cigar
        .elements
        .iter()
        .filter(|e| e.op == op)
        .map(|e| e.length as i64)
        .sum()
}

/// `get3PrimeSoftClippedBases`.
///
/// The comment in Picard says it returns 0 when there are no non-clipping operators, "as it is
/// unclear which clips should be considered on the 3' end". The mechanism is the
/// `foundNonSoftClipOperator` latch, which also means a leading soft clip is never counted, only
/// clips that come *after* something aligned.
fn three_prime_soft_clipped_bases(rec: &BamRecord) -> i64 {
    let negative_strand = rec.flags & READ_REVERSE != 0;
    let mut elements: Vec<_> = rec.cigar.elements.iter().collect();
    if negative_strand {
        elements.reverse();
    }
    let mut found_non_clip = false;
    let mut soft = 0i64;
    for e in elements {
        // `isClipping()` is true for both S and H, so a hard clip also satisfies the latch's
        // "non-clipping" test in the negative - it does not set the latch, and it does not add.
        if !matches!(e.op, Op::S | Op::H) {
            found_non_clip = true;
            continue;
        }
        if found_non_clip && e.op == Op::S {
            soft += e.length as i64;
        }
    }
    soft
}

/// `getUnclippedBaseCount`: read bases consumed by non-clipping operators.
fn unclipped_base_count(rec: &BamRecord) -> i64 {
    rec.cigar
        .elements
        .iter()
        .filter(|e| matches!(e.op, Op::M | Op::I | Op::Eq | Op::X))
        .map(|e| e.length as i64)
        .sum()
}

/// `CoordMath.getCycle`.
fn get_cycle(is_negative_strand: bool, read_length: usize, read_base_index: usize) -> i32 {
    if is_negative_strand {
        read_length as i32 - read_base_index as i32
    } else {
        read_base_index as i32 + 1
    }
}

/// `IndividualAlignmentSummaryMetricsCollector`.
pub struct IndividualCollector {
    metrics: AlignmentSummaryMetrics,
    num_positive_strand: i64,
    read_length_histogram: Histogram,
    aligned_read_length_histogram: Histogram,
    chimeras: i64,
    chimeras_denominator: i64,
    adapter_reads: i64,
    indels: i64,
    num_soft_clipped: i64,
    num_3prime_soft_clipped_bases: i64,
    num_reads_with_3prime_soft_clips: i64,
    num_hard_clipped: i64,
    non_bisulfite_aligned_bases: i64,
    hq_non_bisulfite_aligned_bases: i64,
    mismatch_histogram: Histogram,
    hq_mismatch_histogram: Histogram,
    bad_cycle_histogram: Histogram,
}

impl IndividualCollector {
    pub fn new(category: Category) -> Self {
        IndividualCollector {
            metrics: AlignmentSummaryMetrics::new(category),
            num_positive_strand: 0,
            read_length_histogram: Histogram::new("count", "readLength"),
            aligned_read_length_histogram: Histogram::new("count", "alignedReadLength"),
            chimeras: 0,
            chimeras_denominator: 0,
            adapter_reads: 0,
            indels: 0,
            num_soft_clipped: 0,
            num_3prime_soft_clipped_bases: 0,
            num_reads_with_3prime_soft_clips: 0,
            num_hard_clipped: 0,
            non_bisulfite_aligned_bases: 0,
            hq_non_bisulfite_aligned_bases: 0,
            mismatch_histogram: Histogram::default(),
            hq_mismatch_histogram: Histogram::default(),
            bad_cycle_histogram: Histogram::default(),
        }
    }

    pub fn metrics(&self) -> &AlignmentSummaryMetrics {
        &self.metrics
    }

    /// `getReadHistogram`.
    pub fn read_histogram(&self) -> &Histogram {
        &self.read_length_histogram
    }

    /// `getAlignedReadHistogram`.
    pub fn aligned_read_histogram(&self) -> &Histogram {
        &self.aligned_read_length_histogram
    }

    /// `acceptRecord`.
    pub fn accept(
        &mut self,
        rec: &BamRecord,
        reference: Option<&[u8]>,
        opts: &Options,
        adapters: &AdapterUtility,
    ) {
        // Secondary alignments are skipped here as well as at the outer collector, so a
        // secondary record is filtered twice and a supplementary one only once.
        if rec.flags & SECONDARY != 0 {
            return;
        }
        self.collect_read_data(rec, opts, adapters);
        self.collect_quality_data(rec, reference, opts);
    }

    /// `collectReadData`.
    ///
    /// The supplementary guard is unreachable through `GroupCollector`, which drops those
    /// records first. See the module note: Picard's comment here claims an asymmetry between
    /// read counts and base counts that the class hierarchy prevents.
    fn collect_read_data(&mut self, rec: &BamRecord, opts: &Options, adapters: &AdapterUtility) {
        if rec.flags & SUPPLEMENTARY != 0 {
            return;
        }
        self.metrics.total_reads += 1;
        if rec.flags & VENDOR_FAILED != 0 {
            return;
        }

        self.metrics.pf_reads += 1;
        if is_noise_read(rec) {
            self.metrics.pf_noise_reads += 1;
        }

        self.read_length_histogram
            .increment(rec.read_bases.len() as f64);
        self.aligned_read_length_histogram
            .increment(unclipped_base_count(rec) as f64);

        if adapters.is_adapter(rec) {
            self.adapter_reads += 1;
        }
        self.num_hard_clipped += total_cigar_operator_count(rec, Op::H);

        if rec.flags & READ_UNMAPPED != 0 {
            return;
        }
        self.num_soft_clipped += total_cigar_operator_count(rec, Op::S);

        let three_prime = three_prime_soft_clipped_bases(rec);
        if three_prime > 0 {
            self.num_3prime_soft_clipped_bases += three_prime;
            self.num_reads_with_3prime_soft_clips += 1;
        }

        if !opts.collect_alignment_information {
            return;
        }

        self.metrics.pf_reads_aligned += 1;
        if rec.flags & READ_PAIRED != 0 && rec.flags & PROPER_PAIR == 0 {
            self.metrics.pf_reads_improper_pairs += 1;
        }
        if rec.flags & READ_REVERSE == 0 {
            self.num_positive_strand += 1;
        }

        if rec.flags & READ_PAIRED != 0 && rec.flags & MATE_UNMAPPED == 0 {
            self.metrics.reads_aligned_in_pairs += 1;

            // Picard's own comment says "check that both ends have mapq > minimum", but the
            // condition is `mateMq == null || mateMq >= T && mapq >= T`, and Java binds && more
            // tightly than ||. So a record with no MQ tag passes regardless of its own mapping
            // quality, and only a record that *has* the tag is tested on both ends.
            let mate_mq = rec.tags.get(Tag::new(b"MQ")).and_then(as_int);
            let passes = match mate_mq {
                None => true,
                Some(mq) => {
                    mq >= MAPPING_QUALITY_THRESHOLD as i64
                        && rec.mapping_quality >= MAPPING_QUALITY_THRESHOLD
                }
            };
            if passes {
                self.chimeras_denominator += 1;
                if is_chimeric(rec, opts.max_insert_size, &opts.expected_orientations) {
                    self.chimeras += 1;
                }
            }
        } else if rec.mapping_quality >= MAPPING_QUALITY_THRESHOLD {
            // Fragments and half-mapped pairs: chimerism is read off the SA tag alone.
            self.chimeras_denominator += 1;
            if rec.tags.get(Tag::new(b"SA")).is_some() {
                self.chimeras += 1;
            }
        }
    }

    /// `collectQualityData`.
    fn collect_quality_data(&mut self, rec: &BamRecord, reference: Option<&[u8]>, opts: &Options) {
        let unmapped = rec.flags & READ_UNMAPPED != 0;
        let vendor_failed = rec.flags & VENDOR_FAILED != 0;
        let negative_strand = rec.flags & READ_REVERSE != 0;
        let read_bases = &rec.read_bases;

        if unmapped || vendor_failed || !opts.collect_alignment_information {
            // The unaligned path passes a genuine read index to getCycle.
            for (i, &base) in read_bases.iter().enumerate() {
                if is_no_call(base) {
                    self.bad_cycle_histogram.increment(get_cycle(
                        negative_strand,
                        read_bases.len(),
                        i,
                    ) as f64);
                }
            }
            return;
        }

        let high_quality_mapping = rec.mapping_quality >= MAPPING_QUALITY_THRESHOLD;
        if high_quality_mapping && rec.flags & SUPPLEMENTARY == 0 {
            self.metrics.pf_hq_aligned_reads += 1;
        }

        let qualities = &rec.base_qualities;
        let ref_length = reference.map_or(i32::MAX as usize, |r| r.len());
        let mut mismatch_count = 0i64;
        let mut hq_mismatch_count = 0i64;

        for block in alignment_blocks(&rec.cigar, rec.alignment_start) {
            let read_index = (block.read_start - 1) as usize;
            let ref_index = (block.reference_start - 1) as usize;

            for i in 0..block.length as usize {
                if ref_index + i >= ref_length {
                    break;
                }
                let read_base_index = read_index + i;
                let mut mismatch = reference
                    .is_some_and(|r| !bases_equal(read_bases[read_base_index], r[ref_index + i]));

                // Picard indexes the *reference* with the **read** index here. With a deletion
                // in the CIGAR the two differ, so the bisulfite test compares a read base
                // against the wrong reference base. Reproduced; it is reachable only with
                // IS_BISULFITE_SEQUENCED, which is not the default.
                let bisulfite_match = reference.is_some_and(|r| {
                    opts.is_bisulfite_sequenced
                        && read_base_index < r.len()
                        && bisulfite_bases_equal(
                            negative_strand,
                            read_bases[read_base_index],
                            r[read_base_index],
                        )
                });

                let bisulfite_base = mismatch && bisulfite_match;
                mismatch = mismatch && !bisulfite_match;

                if mismatch {
                    mismatch_count += 1;
                }

                self.metrics.pf_aligned_bases += 1;
                if !bisulfite_base {
                    self.non_bisulfite_aligned_bases += 1;
                }

                if high_quality_mapping {
                    self.metrics.pf_hq_aligned_bases += 1;
                    if !bisulfite_base {
                        self.hq_non_bisulfite_aligned_bases += 1;
                    }
                    if qualities[read_base_index] >= BASE_QUALITY_THRESHOLD {
                        self.metrics.pf_hq_aligned_q20_bases += 1;
                    }
                    if mismatch {
                        hq_mismatch_count += 1;
                    }
                }

                if mismatch || is_no_call(read_bases[read_base_index]) {
                    // `i` is the offset within the block, not the read. See the module note:
                    // this is the divergence, measured in the oracle and reproduced here.
                    self.bad_cycle_histogram.increment(get_cycle(
                        negative_strand,
                        read_bases.len(),
                        i,
                    ) as f64);
                }
            }
        }

        self.mismatch_histogram.increment(mismatch_count as f64);
        self.hq_mismatch_histogram
            .increment(hq_mismatch_count as f64);

        for e in &rec.cigar.elements {
            if matches!(e.op, Op::I | Op::D) {
                self.indels += 1;
            }
        }
    }

    /// `finish`.
    ///
    /// Returns an error in the one case Picard throws: reads present but none of them PF.
    pub fn finish(&mut self, opts: &Options) -> Result<(), &'static str> {
        if self.metrics.pf_reads == 0 {
            return if self.metrics.total_reads > 0 {
                Err("Input file contains no PF_READS.")
            } else {
                Ok(())
            };
        }

        // Plain division, not MathUtil.divide: Picard writes `/` here and `divide` below, and
        // the two differ when the denominator is zero. TOTAL_READS cannot be zero at this point
        // because PF_READS is positive, so the choice is invisible - but it is still a choice
        // the source made, and inverting it would matter if the guard above ever changed.
        self.metrics.pct_pf_reads = self.metrics.pf_reads as f64 / self.metrics.total_reads as f64;
        self.metrics.pct_adapter = self.adapter_reads as f64 / self.metrics.pf_reads as f64;
        self.metrics.mean_read_length = self.read_length_histogram.mean();
        self.metrics.sd_read_length = self.read_length_histogram.standard_deviation();
        self.metrics.median_read_length = self.read_length_histogram.median();
        self.metrics.mad_read_length = self.read_length_histogram.median_absolute_deviation();
        self.metrics.min_read_length = self.read_length_histogram.min().unwrap_or(0.0);
        self.metrics.max_read_length = self.read_length_histogram.max().unwrap_or(0.0);

        // A cycle is bad when at least 80% of reads have a bad base there. The division is
        // `cycleBin.getValue() / metrics.TOTAL_READS`, i.e. the count over *all* reads including
        // the non-PF ones that never contributed to the histogram.
        self.metrics.bad_cycles = 0;
        for (_id, value) in self.bad_cycle_histogram.bins() {
            if value / self.metrics.total_reads as f64 >= 0.8 {
                self.metrics.bad_cycles += 1;
            }
        }

        if opts.collect_alignment_information {
            let total_bases = self.read_length_histogram.sum();
            self.metrics.pct_pf_reads_aligned = divide(
                self.metrics.pf_reads_aligned as f64,
                self.metrics.pf_reads as f64,
            );
            self.metrics.pct_reads_aligned_in_pairs = divide(
                self.metrics.reads_aligned_in_pairs as f64,
                self.metrics.pf_reads_aligned as f64,
            );
            self.metrics.pct_pf_reads_improper_pairs = divide(
                self.metrics.pf_reads_improper_pairs as f64,
                self.metrics.pf_reads_aligned as f64,
            );
            self.metrics.mean_aligned_read_length = self.aligned_read_length_histogram.mean();
            self.metrics.strand_balance = divide(
                self.num_positive_strand as f64,
                self.metrics.pf_reads_aligned as f64,
            );
            self.metrics.pct_chimeras =
                divide(self.chimeras as f64, self.chimeras_denominator as f64);
            self.metrics.pf_indel_rate =
                divide(self.indels as f64, self.metrics.pf_aligned_bases as f64);
            self.metrics.pf_mismatch_rate = divide(
                self.mismatch_histogram.sum(),
                self.non_bisulfite_aligned_bases as f64,
            );
            self.metrics.pf_hq_error_rate = divide(
                self.hq_mismatch_histogram.sum(),
                self.hq_non_bisulfite_aligned_bases as f64,
            );
            self.metrics.pct_hardclip = divide(self.num_hard_clipped as f64, total_bases);
            self.metrics.pct_softclip = divide(self.num_soft_clipped as f64, total_bases);
            self.metrics.avg_pos_3prime_softclip_length = divide(
                self.num_3prime_soft_clipped_bases as f64,
                self.num_reads_with_3prime_soft_clips as f64,
            );
            self.metrics.pf_hq_median_mismatches = self.hq_mismatch_histogram.median();
        }
        Ok(())
    }
}

/// `SequenceUtil.bisulfiteBasesEqual`.
///
/// On the positive strand a reference `C` matches a read `T`; on the negative strand a reference
/// `G` matches a read `A`. Anything else falls through to the ordinary comparison.
fn bisulfite_bases_equal(negative_strand: bool, read: u8, reference: u8) -> bool {
    if negative_strand {
        if bases_equal(reference, b'G') && bases_equal(read, b'A') {
            return true;
        }
    } else if bases_equal(reference, b'C') && bases_equal(read, b'T') {
        return true;
    }
    bases_equal(read, reference)
}

/// `SAMRecord.getIntegerAttribute`, narrowed to what these two call sites need.
fn as_int(v: &TagValue) -> Option<i64> {
    match v {
        TagValue::Int(i) => Some(*i),
        _ => None,
    }
}

/// `isNoiseRead`: the reserved `XN` tag equal to 1.
fn is_noise_read(rec: &BamRecord) -> bool {
    rec.tags.get(Tag::new(b"XN")).and_then(as_int) == Some(1)
}

/// `GroupAlignmentSummaryMetricsPerUnitMetricCollector`: the four pairing categories.
pub struct GroupCollector {
    pub unpaired: IndividualCollector,
    pub first_of_pair: IndividualCollector,
    pub second_of_pair: IndividualCollector,
    pub pair: IndividualCollector,
    opts: Options,
    adapters: AdapterUtility,
}

impl GroupCollector {
    pub fn new(opts: Options) -> Self {
        GroupCollector {
            unpaired: IndividualCollector::new(Category::Unpaired),
            first_of_pair: IndividualCollector::new(Category::FirstOfPair),
            second_of_pair: IndividualCollector::new(Category::SecondOfPair),
            pair: IndividualCollector::new(Category::Pair),
            adapters: AdapterUtility::with_defaults(),
            opts,
        }
    }

    /// `AlignmentSummaryMetricsCollector.acceptRecord`, which filters secondary **and**
    /// supplementary records before the per-unit collectors see them, and
    /// `GroupAlignmentSummaryMetricsPerUnitMetricCollector.acceptRecord`, which routes by pairing.
    ///
    /// A paired record goes to **two** collectors: its own end and the PAIR aggregate. So
    /// `TOTAL_READS` summed over the categories double-counts every paired read, by design.
    pub fn accept(&mut self, rec: &BamRecord, reference: Option<&[u8]>) {
        if rec.flags & SECONDARY != 0 || rec.flags & SUPPLEMENTARY != 0 {
            return;
        }
        if rec.flags & READ_PAIRED != 0 {
            if rec.flags & FIRST_OF_PAIR != 0 {
                self.first_of_pair
                    .accept(rec, reference, &self.opts, &self.adapters);
            } else {
                self.second_of_pair
                    .accept(rec, reference, &self.opts, &self.adapters);
            }
            self.pair.accept(rec, reference, &self.opts, &self.adapters);
        } else {
            self.unpaired
                .accept(rec, reference, &self.opts, &self.adapters);
        }
    }

    pub fn finish(&mut self) -> Result<(), &'static str> {
        self.unpaired.finish(&self.opts)?;
        self.first_of_pair.finish(&self.opts)?;
        self.second_of_pair.finish(&self.opts)?;
        self.pair.finish(&self.opts)?;
        Ok(())
    }

    /// The four read-length histograms `CollectAlignmentSummaryMetrics.finish` appends, in the
    /// order it appends them, with the labels it gives them.
    ///
    /// `addHistogramToMetrics` **mutates** the collector's histogram to set the labels rather
    /// than copying it, so the labels are a property of the tool and not of the collector. And
    /// the PAIRED pair comes first even for an unpaired file, where both are empty and
    /// `printHistogram` drops them - which is why an unpaired file's histogram table has two
    /// columns rather than four, and why the column order is not evidence of the append order.
    pub fn read_length_histograms(&self) -> Vec<(&'static str, &Histogram)> {
        vec![
            ("PAIRED_TOTAL_LENGTH_COUNT", self.pair.read_histogram()),
            (
                "PAIRED_ALIGNED_LENGTH_COUNT",
                self.pair.aligned_read_histogram(),
            ),
            (
                "UNPAIRED_TOTAL_LENGTH_COUNT",
                self.unpaired.read_histogram(),
            ),
            (
                "UNPAIRED_ALIGNED_LENGTH_COUNT",
                self.unpaired.aligned_read_histogram(),
            ),
        ]
    }

    /// `addMetricsToFile`.
    ///
    /// Two rules that decide which rows exist at all:
    ///
    /// * `PAIR`'s `BAD_CYCLES` is **overwritten** with the sum of the two ends', discarding what
    ///   the PAIR collector computed from its own histogram. So the one metric is a sum where
    ///   every other metric on the row is an independent accumulation.
    /// * `UNPAIRED` is emitted when it has reads *or* when `FIRST_OF_PAIR` has none. An
    ///   all-paired file therefore has three rows, and an empty file has one row of zeros.
    pub fn rows(&mut self) -> Vec<AlignmentSummaryMetrics> {
        let mut out = Vec::new();
        if self.first_of_pair.metrics.total_reads > 0 {
            self.pair.metrics.bad_cycles =
                self.first_of_pair.metrics.bad_cycles + self.second_of_pair.metrics.bad_cycles;
            out.push(self.first_of_pair.metrics.clone());
            out.push(self.second_of_pair.metrics.clone());
            out.push(self.pair.metrics.clone());
        }
        if self.unpaired.metrics.total_reads > 0 || self.first_of_pair.metrics.total_reads == 0 {
            out.push(self.unpaired.metrics.clone());
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::{Cigar, CigarElement};

    fn rec(cigar: &[(u32, Op)], bases: &[u8], flags: u16) -> BamRecord {
        BamRecord {
            read_name: "r".to_string(),
            flags,
            reference_index: 0,
            alignment_start: 1,
            mapping_quality: 60,
            cigar: Cigar::new(
                cigar
                    .iter()
                    .map(|&(length, op)| CigarElement { length, op })
                    .collect(),
            ),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: bases.to_vec(),
            base_qualities: vec![40; bases.len()],
            tags: Default::default(),
        }
    }

    /// The probe from `tools/asm-conformance/badcycle_probe.sh`, in Rust. Two mismatches in
    /// different alignment blocks sharing a block offset collapse into one bad cycle.
    #[test]
    fn bad_cycles_bins_by_block_offset_not_read_index() {
        let opts = Options::default();
        let adapters = AdapterUtility::with_defaults();
        let read = rec(&[(10, Op::M), (5, Op::D), (10, Op::M)], &[b'A'; 20], 0);

        // Reference position 4 aligns to read index 3 (block offset 3); reference position 19
        // aligns to read index 13 (also block offset 3).
        let mut reference = vec![b'A'; 100];
        reference[3] = b'C';
        reference[18] = b'C';
        let mut c = IndividualCollector::new(Category::Unpaired);
        c.accept(&read, Some(&reference), &opts, &adapters);
        c.finish(&opts).unwrap();
        assert_eq!(c.metrics.bad_cycles, 1, "two mismatches, one bad cycle");
        assert_eq!(c.metrics.pf_hq_median_mismatches, 2.0, "both were counted");

        // Reference position 21 aligns to read index 15, block offset 5: no collision.
        let mut reference = vec![b'A'; 100];
        reference[3] = b'C';
        reference[20] = b'C';
        let mut c = IndividualCollector::new(Category::Unpaired);
        c.accept(&read, Some(&reference), &opts, &adapters);
        c.finish(&opts).unwrap();
        assert_eq!(c.metrics.bad_cycles, 2);
        assert_eq!(c.metrics.pf_hq_median_mismatches, 2.0);
    }

    /// The 3' soft clip latch: a leading clip never counts, a trailing one does, and on the
    /// negative strand the roles swap because the CIGAR is walked backwards.
    #[test]
    fn only_soft_clips_after_an_aligned_operator_count() {
        let leading = rec(&[(5, Op::S), (10, Op::M)], &[b'A'; 15], 0);
        assert_eq!(three_prime_soft_clipped_bases(&leading), 0);

        let trailing = rec(&[(10, Op::M), (5, Op::S)], &[b'A'; 15], 0);
        assert_eq!(three_prime_soft_clipped_bases(&trailing), 5);

        let leading_rev = rec(&[(5, Op::S), (10, Op::M)], &[b'A'; 15], READ_REVERSE);
        assert_eq!(three_prime_soft_clipped_bases(&leading_rev), 5);
    }

    #[test]
    fn a_cigar_with_no_aligned_operator_gives_no_three_prime_clip() {
        let all_clip = rec(&[(5, Op::S), (5, Op::S)], &[b'A'; 10], 0);
        assert_eq!(three_prime_soft_clipped_bases(&all_clip), 0);
    }

    /// A paired read is accepted by two collectors, so summing TOTAL_READS across rows
    /// double-counts. That is the tool's design and not a bug in the port.
    #[test]
    fn a_paired_read_is_counted_in_two_categories() {
        let mut g = GroupCollector::new(Options::default());
        let read = rec(&[(10, Op::M)], &[b'A'; 10], READ_PAIRED | FIRST_OF_PAIR);
        g.accept(&read, None);
        g.finish().unwrap();
        assert_eq!(g.first_of_pair.metrics.total_reads, 1);
        assert_eq!(g.pair.metrics.total_reads, 1);
        assert_eq!(g.second_of_pair.metrics.total_reads, 0);
        assert_eq!(g.unpaired.metrics.total_reads, 0);
    }

    /// An all-paired file still gets an UNPAIRED row only if FIRST_OF_PAIR is empty, so here it
    /// does not; and PAIR's BAD_CYCLES is the sum of the two ends.
    #[test]
    fn the_pair_row_takes_its_bad_cycles_from_the_two_ends() {
        let mut g = GroupCollector::new(Options::default());
        g.accept(
            &rec(&[(10, Op::M)], &[b'A'; 10], READ_PAIRED | FIRST_OF_PAIR),
            None,
        );
        g.accept(&rec(&[(10, Op::M)], &[b'A'; 10], READ_PAIRED), None);
        g.finish().unwrap();
        g.first_of_pair.metrics.bad_cycles = 3;
        g.second_of_pair.metrics.bad_cycles = 4;
        let rows = g.rows();
        assert_eq!(rows.len(), 3, "no UNPAIRED row when there are paired reads");
        assert_eq!(rows[2].category, Category::Pair);
        assert_eq!(rows[2].bad_cycles, 7);
    }

    /// An empty input still produces one row: the UNPAIRED zeros.
    #[test]
    fn an_empty_input_gives_one_row_of_zeros() {
        let mut g = GroupCollector::new(Options::default());
        g.finish().unwrap();
        let rows = g.rows();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].category, Category::Unpaired);
        assert_eq!(rows[0].total_reads, 0);
    }

    /// Reads present but none PF is the one condition Picard throws on.
    #[test]
    fn reads_with_none_passing_filter_is_an_error() {
        let opts = Options::default();
        let adapters = AdapterUtility::with_defaults();
        let mut c = IndividualCollector::new(Category::Unpaired);
        c.accept(
            &rec(&[(10, Op::M)], &[b'A'; 10], VENDOR_FAILED),
            None,
            &opts,
            &adapters,
        );
        assert_eq!(c.finish(&opts), Err("Input file contains no PF_READS."));
    }

    /// `MathUtil.divide` is a magnitude threshold, not a zero test.
    #[test]
    fn divide_returns_zero_below_a_millionth() {
        assert_eq!(divide(1.0, 0.0), 0.0);
        assert_eq!(divide(1.0, 5e-7), 0.0, "not 2000000.0");
        assert_eq!(divide(1.0, 2e-6), 500000.0);
        assert_eq!(divide(1.0, -5e-7), 0.0, "the magnitude is what is tested");
    }

    /// A supplementary record reaches nothing: the group collector drops it before any
    /// per-unit collector runs, so neither the read counts nor the base counts move. This is
    /// what makes Picard's "do include supplementary records for base counts" comment false of
    /// its own tool; the oracle agrees, at `PF_ALIGNED_BASES = 20` rather than 40 on the
    /// `supplementary` corpus case.
    #[test]
    fn a_supplementary_record_reaches_neither_counter_through_the_group_collector() {
        let mut g = GroupCollector::new(Options::default());
        let reference = vec![b'A'; 100];
        g.accept(&rec(&[(10, Op::M)], &[b'A'; 10], 0), Some(&reference));
        g.accept(
            &rec(&[(10, Op::M)], &[b'A'; 10], SUPPLEMENTARY),
            Some(&reference),
        );
        assert_eq!(g.unpaired.metrics.total_reads, 1);
        assert_eq!(g.unpaired.metrics.pf_aligned_bases, 10, "not 20");
    }

    /// Driven directly, bypassing the group filter, the inner guard does what its author meant.
    /// The two tests together are the finding: the behaviour exists and the tool cannot reach it.
    #[test]
    fn the_unreachable_guard_still_behaves_as_written() {
        let opts = Options::default();
        let adapters = AdapterUtility::with_defaults();
        let mut c = IndividualCollector::new(Category::Unpaired);
        let reference = vec![b'A'; 100];
        c.accept(
            &rec(&[(10, Op::M)], &[b'A'; 10], SUPPLEMENTARY),
            Some(&reference),
            &opts,
            &adapters,
        );
        assert_eq!(c.metrics.total_reads, 0, "not counted as a read");
        assert_eq!(c.metrics.pf_aligned_bases, 10, "counted as bases");
    }
}
