//! GC content over a sliding window.
//!
//! Ported from `picard.analysis.GcBiasUtils` and its inner `CalculateGcState`, tag 3.4.0. This
//! is the core of `CollectGcBiasMetrics`, and it is a dense little function: five behaviours in
//! thirty lines, four of which a reimplementation gets wrong by writing the obvious thing.
//!
//! ## The window range excludes both ends
//!
//! Both callers loop `for (int i = 1; i < lastWindowStart; ++i)`. So the window starting at
//! reference position **0 is never computed**, and neither is the one starting at
//! `lastWindowStart`. A port that looped `0..=lastWindowStart` would produce two more windows
//! than Picard and a different GC histogram, on every reference.
//!
//! ## The comment and the code disagree about the no-call threshold
//!
//! ```java
//! // If the window includes more than five no-calls then -1 is returned.
//! ...
//! if (state.nCount > 4) return -1;
//! ```
//!
//! More than *four*. The code is what runs, so four is what is ported.
//!
//! ## GC is truncating integer division
//!
//! `(gcCount * 100) / (endIndex - startIndex)`, in `int`. A window that is 50.9% GC bins as 50.
//!
//! ## The no-call test is case-insensitive on three of its four uses
//!
//! The initialising branch counts no-calls with `SequenceUtil.basesEqual(base, 'N')`, which folds
//! case. The incremental branch counts the **incoming** base with a raw byte comparison:
//!
//! ```java
//! else if (newBase == 'N') ++state.nCount;
//! ```
//!
//! while still *decrementing* on the outgoing base with the case-insensitive test. So over a
//! sequence containing lowercase `n`, no-calls are removed from the count that were never added,
//! and `nCount` can go negative — after which a window that should be rejected is accepted and
//! its GC value enters the histogram.
//!
//! **Measured in the oracle**, `tools/gcbias-conformance/GcAsymmetryProbe.java`, over a sequence
//! with seven no-calls entering the window:
//!
//! ```text
//! upper: [.., 60, -1, -1, -1, -1, -1, -1, -1, -1, 60, ..]   8 windows rejected
//! lower: [.., 60, 50, 40, 30, 30, 30, 30, 40, 50, 60, ..]   0 windows rejected
//! ```
//!
//! **And it is unreachable through the tool.** Both callers in Picard run
//! `StringUtil.toUpperCase(refBases)` before calling in, so no lowercase base ever arrives. The
//! asymmetry is a latent bug in a public static method, not a divergence in
//! `CollectGcBiasMetrics`. It is reproduced here because `GcBiasUtils` is public and a future
//! caller that skipped the uppercase would depend on it, and it is documented as unreachable
//! rather than presented as a live divergence.
//!
//! ## Where the uncomputed windows go
//!
//! The zeros left at index 0 and at or past `lastWindowStart` are not inert.
//! `GcBiasMetricsCollector.addRead` reads `gc[pos]` and bins on the result without checking
//! whether that entry was ever written, so **a read whose window start falls outside the
//! computed range is charged to GC bin 0**, not skipped. Measured in the oracle
//! (`tools/gcbias-conformance/WindowBinningProbe.java`), a read at position 320 of a 400-base
//! contig with a 100-base window reports `gc=0`.
//!
//! Two more from the same probe, both about which window a read is charged to:
//!
//! ```java
//! final int pos = rec.getReadNegativeStrandFlag()
//!     ? rec.getAlignmentEnd() - scanWindowSize
//!     : rec.getAlignmentStart();
//! ```
//!
//! A forward read is charged to the window at its alignment *start*, a reverse read to the one
//! at its alignment *end* minus the window size. Two reads covering **the same reference bases**
//! therefore land in different bins: over a pure-GC stretch the forward read reports `gc=99` and
//! the reverse `gc=100`. riker's errata describes "a forward-strand window-binning fix", so this
//! is a place an independent reimplementation chose to differ.
//!
//! And `if (pos > 0)` silently drops a reverse read near the contig start from the GC bins while
//! still counting it in `totalAlignedReads`, so the per-bin read counts do not sum to the
//! aligned-read total.
//!
//! Getting to the case asymmetry took two probes. The first put the no-calls at the start of the sequence,
//! where they only ever *left* the window, so the buggy branch never ran and the two cases
//! agreed. That probe proved nothing; the second one, with the no-calls entering, is the
//! evidence.

use htsjdk_bam::sequence::bases_equal;
use htsjdk_metrics::file::{MetricBean, Value};

/// `GcBiasUtils.CalculateGcState`.
///
/// The state is what makes the window incremental, and the whole asymmetry above lives in the
/// difference between its first update and its later ones.
#[derive(Debug, Clone)]
pub struct CalculateGcState {
    pub init: bool,
    pub n_count: i32,
    pub gc_count: i32,
    pub prior_base: u8,
}

impl Default for CalculateGcState {
    fn default() -> Self {
        CalculateGcState {
            init: true,
            n_count: 0,
            gc_count: 0,
            prior_base: 0,
        }
    }
}

/// `GcBiasUtils.calculateGc(bases, startIndex, endIndex, state)`.
///
/// Returns the GC percentage of `bases[start..end]`, or `-1` when the window holds more than
/// four no-calls.
///
/// `n_count` is `i32` rather than an unsigned type on purpose: htsjdk's is an `int` and the
/// asymmetry above can drive it negative. A `u32` here would panic in debug and wrap in release,
/// and either would be a different function.
pub fn calculate_gc(bases: &[u8], start: usize, end: usize, state: &mut CalculateGcState) -> i32 {
    if state.init {
        state.init = false;
        state.gc_count = 0;
        state.n_count = 0;
        for &base in &bases[start..end] {
            if bases_equal(base, b'G') || bases_equal(base, b'C') {
                state.gc_count += 1;
            } else if bases_equal(base, b'N') {
                state.n_count += 1;
            }
        }
    } else {
        let new_base = bases[end - 1];
        if bases_equal(new_base, b'G') || bases_equal(new_base, b'C') {
            state.gc_count += 1;
        } else if new_base == b'N' {
            // Raw byte comparison, not `bases_equal`. This is the asymmetry; see the module note.
            state.n_count += 1;
        }

        if bases_equal(state.prior_base, b'G') || bases_equal(state.prior_base, b'C') {
            state.gc_count -= 1;
        } else if bases_equal(state.prior_base, b'N') {
            state.n_count -= 1;
        }
    }
    state.prior_base = bases[start];

    if state.n_count > 4 {
        -1
    } else {
        // Integer division, truncating, as the Java's int arithmetic does.
        (state.gc_count * 100) / (end - start) as i32
    }
}

/// `GcBiasUtils.calculateAllGcs(refBases, lastWindowStart, windowSize)`.
///
/// The returned array is `refBases.len() + 1` long and is indexed by window start, so index 0 and
/// every index at or past `last_window_start` are left at zero rather than computed. Those zeros
/// are indistinguishable from a genuine 0% GC window, which is why the range matters: a caller
/// cannot tell the two apart afterwards.
pub fn calculate_all_gcs(
    ref_bases: &[u8],
    last_window_start: usize,
    window_size: usize,
) -> Vec<i8> {
    let mut state = CalculateGcState::default();
    let mut gc = vec![0i8; ref_bases.len() + 1];
    let mut i = 1;
    while i < last_window_start {
        let window_end = i + window_size;
        gc[i] = calculate_gc(ref_bases, i, window_end, &mut state) as i8;
        i += 1;
    }
    gc
}

/// `GcBiasUtils.calculateRefWindowsByGc(windows, referenceSequence, windowSize)`, taking the
/// contigs already read rather than a path.
///
/// **This one uppercases and `calculate_all_gcs` does not**, which is the whole of the
/// difference between them. htsjdk's FASTA reader returns bases with their case intact, so a
/// soft-masked reference arrives lowercase; this function folds it and the other does not.
pub fn calculate_ref_windows_by_gc(
    windows: usize,
    contigs: &[Vec<u8>],
    window_size: usize,
) -> Vec<i32> {
    let mut windows_by_gc = vec![0i32; windows];
    for contig in contigs {
        let ref_bases: Vec<u8> = contig.to_ascii_uppercase();
        if ref_bases.len() < window_size {
            continue;
        }
        let last_window_start = ref_bases.len() - window_size;
        let mut state = CalculateGcState::default();
        let mut i = 1;
        while i < last_window_start {
            let gc_bin = calculate_gc(&ref_bases, i, i + window_size, &mut state);
            if gc_bin != -1 {
                windows_by_gc[gc_bin as usize] += 1;
            }
            i += 1;
        }
    }
    windows_by_gc
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gcs(seq: &str, window: usize) -> Vec<i8> {
        let bases = seq.as_bytes();
        calculate_all_gcs(bases, bases.len() - window, window)
    }

    /// The oracle's numbers for the probe sequence, uppercase, pinned directly.
    #[test]
    fn the_probe_sequence_matches_the_oracle_uppercase() {
        let seq = "GCGCGCGCGCGNNNNNNNGCGCGCGCGCGCATATATATAT";
        let got = gcs(seq, 10);
        assert_eq!(
            &got[..22],
            &[
                0, 100, 90, 80, 70, 60, -1, -1, -1, -1, -1, -1, -1, -1, 60, 70, 80, 90, 100, 100,
                100, 90
            ]
        );
    }

    /// ...and lowercase, where the incoming-base test stops recognising the no-calls and every
    /// window is accepted. Measured in the oracle before being written here.
    #[test]
    fn the_probe_sequence_matches_the_oracle_lowercase() {
        let seq = "GCGCGCGCGCGnnnnnnnGCGCGCGCGCGCATATATATAT";
        let got = gcs(seq, 10);
        assert_eq!(
            &got[..22],
            &[
                0, 100, 90, 80, 70, 60, 50, 40, 30, 30, 30, 30, 40, 50, 60, 70, 80, 90, 100, 100,
                100, 90
            ]
        );
        assert!(
            !got[1..30].contains(&-1),
            "no window is rejected once the no-calls are lowercase"
        );
    }

    /// Window 0 is never computed, and neither is the last one. Both stay at the initial zero,
    /// which is indistinguishable from a genuine 0% GC window.
    #[test]
    fn the_first_and_last_windows_are_never_computed() {
        let seq = "GCGCGCGCGCGCGCGCGCGC";
        let got = gcs(seq, 10);
        assert_eq!(got[0], 0, "window 0 is skipped, not computed as 100");
        assert_eq!(got[1], 100, "window 1 is computed");
        assert_eq!(got[10], 0, "lastWindowStart is excluded by the `<` bound");
    }

    /// The threshold is four, not the five the comment claims.
    #[test]
    fn the_threshold_is_four_no_calls_not_five() {
        let mut state = CalculateGcState::default();
        // Exactly four no-calls: accepted.
        assert_ne!(calculate_gc(b"NNNNGCGCGC", 0, 10, &mut state), -1);
        let mut state = CalculateGcState::default();
        // Five: rejected.
        assert_eq!(calculate_gc(b"NNNNNGCGCG", 0, 10, &mut state), -1);
    }

    /// Truncating integer division: 5 of 9 bases is 55.5%, reported as 55.
    #[test]
    fn gc_is_truncated_not_rounded() {
        let mut state = CalculateGcState::default();
        assert_eq!(calculate_gc(b"GCGCGATAT", 0, 9, &mut state), 55);
    }

    /// The uppercasing is the difference between the two entry points, and it is
    /// `calculate_ref_windows_by_gc` that has it.
    #[test]
    fn only_the_windows_by_gc_entry_point_folds_case() {
        let seq = b"GCGCGCGCGCGnnnnnnnGCGCGCGCGCGCATATATATAT".to_vec();
        let folded = calculate_ref_windows_by_gc(101, std::slice::from_ref(&seq), 10);
        let unfolded = calculate_all_gcs(&seq, seq.len() - 10, 10);
        // Folded: the no-call windows are rejected, so they contribute to no bin.
        let total: i32 = folded.iter().sum();
        assert!(
            total < (seq.len() - 10) as i32,
            "some windows were rejected"
        );
        // Unfolded: none were.
        assert!(!unfolded[1..30].contains(&-1));
    }

    /// A contig shorter than the window produces no windows rather than panicking on the
    /// underflow that `refLength - windowSize` would give.
    #[test]
    fn a_contig_shorter_than_the_window_is_skipped() {
        assert_eq!(
            calculate_ref_windows_by_gc(101, &[b"GCGC".to_vec()], 10)
                .iter()
                .sum::<i32>(),
            0
        );
    }
}

// ---------------------------------------------------------------------------------------------
// The collector
// ---------------------------------------------------------------------------------------------

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{count_deleted_bases, count_inserted_bases, count_mismatches};

/// `GcBiasMetricsCollector.BINS`: 101 bins, one per whole GC percentage from 0 to 100.
pub const BINS: usize = 101;
/// `CollectGcBiasMetrics.SCAN_WINDOW_SIZE` default.
pub const DEFAULT_SCAN_WINDOW_SIZE: i32 = 100;

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const READ_REVERSE: u16 = 0x10;
const FIRST_OF_PAIR: u16 = 0x40;

/// `QualityUtil.getPhredScoreFromObsAndErrors(observations, errors)`.
///
/// `(int) Math.round(-10 * Math.log10(errors / observations))`. Three things to keep:
///
///   * `Math.log10`, which decision 0006 in htsjdk-rs establishes is targeted at correct
///     rounding, so `jmath::math::log10` and not the system libm;
///   * `Math.round(double)` is `floor(x + 0.5)`, **not** Rust's `f64::round`, which rounds half
///     away from zero. The two differ on negative halves: Java gives 0 for -0.5, Rust gives -1.
///     The argument here is non-negative whenever `errors <= observations`, so the difference is
///     unreachable through this caller — but writing `.round()` would be a different function;
///   * the cast to `int` truncates toward zero *after* the rounding, which for a non-negative
///     value is a no-op and is kept only so the shape matches.
pub fn phred_score_from_obs_and_errors(observations: f64, errors: f64) -> i32 {
    htsjdk_bam::quality_util::phred_score_from_obs_and_errors(observations, errors)
}

/// `GcBiasMetricsCollector.GcObject`: one accumulation level's counters.
#[derive(Debug, Clone)]
pub struct GcObject {
    pub total_clusters: i64,
    pub total_aligned_reads: i64,
    pub reads_by_gc: Vec<i32>,
    pub bases_by_gc: Vec<i64>,
    pub errors_by_gc: Vec<i64>,
}

impl Default for GcObject {
    fn default() -> Self {
        GcObject {
            total_clusters: 0,
            total_aligned_reads: 0,
            reads_by_gc: vec![0; BINS],
            bases_by_gc: vec![0; BINS],
            errors_by_gc: vec![0; BINS],
        }
    }
}

/// `GcBiasDetailMetrics`, in declared field order with the inherited `MultilevelMetrics` fields
/// last, as `Class.getFields()` returns them.
#[derive(Debug, Clone, PartialEq)]
pub struct GcBiasDetailMetrics {
    pub accumulation_level: String,
    pub reads_used: String,
    pub gc: i32,
    pub windows: i32,
    pub read_starts: i64,
    pub mean_base_quality: i32,
    pub normalized_coverage: f64,
    pub error_bar_width: f64,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

const DETAIL_COLUMNS: &[&str] = &[
    "ACCUMULATION_LEVEL",
    "READS_USED",
    "GC",
    "WINDOWS",
    "READ_STARTS",
    "MEAN_BASE_QUALITY",
    "NORMALIZED_COVERAGE",
    "ERROR_BAR_WIDTH",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

impl MetricBean for GcBiasDetailMetrics {
    fn class_name(&self) -> &str {
        "picard.analysis.GcBiasDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        DETAIL_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let text = |o: &Option<String>| match o {
            Some(s) => Value::Str(s.clone()),
            None => Value::Null,
        };
        vec![
            Value::Str(self.accumulation_level.clone()),
            Value::Str(self.reads_used.clone()),
            Value::Long(self.gc as i64),
            Value::Long(self.windows as i64),
            Value::Long(self.read_starts),
            Value::Long(self.mean_base_quality as i64),
            Value::Double(self.normalized_coverage),
            Value::Double(self.error_bar_width),
            text(&self.sample),
            text(&self.library),
            text(&self.read_group),
        ]
    }
}

/// `GcBiasSummaryMetrics`.
#[derive(Debug, Clone, PartialEq)]
pub struct GcBiasSummaryMetrics {
    pub accumulation_level: String,
    pub reads_used: String,
    pub window_size: i32,
    pub total_clusters: i64,
    pub aligned_reads: i64,
    pub at_dropout: f64,
    pub gc_dropout: f64,
    pub gc_nc_0_19: f64,
    pub gc_nc_20_39: f64,
    pub gc_nc_40_59: f64,
    pub gc_nc_60_79: f64,
    pub gc_nc_80_100: f64,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

const SUMMARY_COLUMNS: &[&str] = &[
    "ACCUMULATION_LEVEL",
    "READS_USED",
    "WINDOW_SIZE",
    "TOTAL_CLUSTERS",
    "ALIGNED_READS",
    "AT_DROPOUT",
    "GC_DROPOUT",
    "GC_NC_0_19",
    "GC_NC_20_39",
    "GC_NC_40_59",
    "GC_NC_60_79",
    "GC_NC_80_100",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

impl MetricBean for GcBiasSummaryMetrics {
    fn class_name(&self) -> &str {
        "picard.analysis.GcBiasSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        SUMMARY_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let text = |o: &Option<String>| match o {
            Some(s) => Value::Str(s.clone()),
            None => Value::Null,
        };
        vec![
            Value::Str(self.accumulation_level.clone()),
            Value::Str(self.reads_used.clone()),
            Value::Long(self.window_size as i64),
            Value::Long(self.total_clusters),
            Value::Long(self.aligned_reads),
            Value::Double(self.at_dropout),
            Value::Double(self.gc_dropout),
            Value::Double(self.gc_nc_0_19),
            Value::Double(self.gc_nc_20_39),
            Value::Double(self.gc_nc_40_59),
            Value::Double(self.gc_nc_60_79),
            Value::Double(self.gc_nc_80_100),
            text(&self.sample),
            text(&self.library),
            text(&self.read_group),
        ]
    }
}

/// `GcBiasMetricsCollector` restricted to the ALL_READS level, which is the default, as
/// `InsertSizeMetricsCollector` is.
pub struct GcBiasMetricsCollector {
    scan_window_size: i32,
    bisulfite: bool,
    all_reads: GcObject,
    /// The reference windows by GC bin, computed once from the whole reference before any read
    /// is seen. This is `calculate_ref_windows_by_gc`, which uppercases; the per-contig `gc`
    /// array used for binning reads does not. See the module note.
    windows_by_gc: Vec<i32>,
    /// The per-contig GC array from `calculate_all_gcs`, recomputed when the contig changes.
    gc: Vec<i8>,
    reference_bases: Vec<u8>,
}

impl GcBiasMetricsCollector {
    pub fn new(contigs: &[Vec<u8>], scan_window_size: i32, bisulfite: bool) -> Self {
        GcBiasMetricsCollector {
            scan_window_size,
            bisulfite,
            all_reads: GcObject::default(),
            windows_by_gc: calculate_ref_windows_by_gc(BINS, contigs, scan_window_size as usize),
            gc: Vec::new(),
            reference_bases: Vec::new(),
        }
    }

    /// `acceptRecord`, for a record whose contig's bases are supplied.
    ///
    /// An unmapped read still bumps `totalClusters`, through `updateTotalClusters`, so the
    /// cluster count includes reads that never reach a GC bin.
    pub fn accept(&mut self, rec: &BamRecord, reference: &[u8]) {
        // A read with `*` in SEQ is dropped with a warning and never counted at all, not even as
        // a cluster.
        if rec.read_bases.is_empty() {
            return;
        }
        if rec.flags & READ_UNMAPPED != 0 {
            self.update_total_clusters(rec);
            return;
        }
        if self.reference_bases != reference {
            // The collector uppercases here, which is what makes the lowercase asymmetry in
            // `calculate_gc` unreachable through this tool.
            self.reference_bases = reference.to_ascii_uppercase();
            let last_window_start = self
                .reference_bases
                .len()
                .saturating_sub(self.scan_window_size as usize);
            self.gc = calculate_all_gcs(
                &self.reference_bases,
                last_window_start,
                self.scan_window_size as usize,
            );
        }
        self.add_read(rec);
    }

    fn update_total_clusters(&mut self, rec: &BamRecord) {
        if rec.flags & READ_PAIRED == 0 || rec.flags & FIRST_OF_PAIR != 0 {
            self.all_reads.total_clusters += 1;
        }
    }

    /// `addRead`.
    ///
    /// The window a read is charged to depends on its strand: forward reads use their alignment
    /// start, reverse reads their alignment end minus the window size. Two reads over the same
    /// bases therefore land in different bins. See the module note; this is measured, and it is
    /// a place riker chose to differ.
    fn add_read(&mut self, rec: &BamRecord) {
        self.update_total_clusters(rec);
        let pos = if rec.flags & READ_REVERSE != 0 {
            alignment_end(rec) - self.scan_window_size
        } else {
            rec.alignment_start
        };
        self.all_reads.total_aligned_reads += 1;
        if pos <= 0 {
            // A reverse read near the contig start is dropped from the GC bins while still
            // counting as an aligned read, so the per-bin counts do not sum to the total.
            return;
        }
        // `gc[pos]` is read without checking whether that entry was ever computed. Index 0 and
        // everything at or past lastWindowStart are still zero, so a read whose window falls
        // there is charged to GC bin 0.
        let window_gc = *self.gc.get(pos as usize).unwrap_or(&0);
        if window_gc >= 0 {
            let bin = window_gc as usize;
            self.all_reads.reads_by_gc[bin] += 1;
            self.all_reads.bases_by_gc[bin] += rec.read_bases.len() as i64;
            let blocks = alignment_blocks(&rec.cigar, rec.alignment_start);
            self.all_reads.errors_by_gc[bin] += count_mismatches(
                &rec.read_bases,
                &blocks,
                &self.reference_bases,
                0,
                rec.flags & READ_REVERSE != 0,
                self.bisulfite,
            ) as i64
                + count_inserted_bases(&rec.cigar) as i64
                + count_deleted_bases(&rec.cigar) as i64;
        }
    }

    /// `addGcDataToFile`: the 101 detail rows and the summary, or nothing at all when no read
    /// aligned.
    pub fn rows(&self) -> Option<(Vec<GcBiasDetailMetrics>, GcBiasSummaryMetrics)> {
        if self.all_reads.total_aligned_reads == 0 {
            return None;
        }
        let total_windows: f64 = self.windows_by_gc.iter().map(|&w| w as f64).sum();
        let total_reads: f64 = self.all_reads.reads_by_gc.iter().map(|&r| r as f64).sum();
        let mean_reads_per_window = total_reads / total_windows;

        let mut details = Vec::with_capacity(BINS);
        for i in 0..BINS {
            let windows = self.windows_by_gc[i];
            let read_starts = self.all_reads.reads_by_gc[i] as i64;
            let mut detail = GcBiasDetailMetrics {
                accumulation_level: "All Reads".to_string(),
                reads_used: "ALL".to_string(),
                gc: i as i32,
                windows,
                read_starts,
                mean_base_quality: 0,
                normalized_coverage: 0.0,
                error_bar_width: 0.0,
                sample: None,
                library: None,
                read_group: None,
            };
            // Only set when there were errors, so a bin with none keeps the default 0 rather
            // than the infinite phred score the formula would give.
            if self.all_reads.errors_by_gc[i] > 0 {
                detail.mean_base_quality = phred_score_from_obs_and_errors(
                    self.all_reads.bases_by_gc[i] as f64,
                    self.all_reads.errors_by_gc[i] as f64,
                );
            }
            if windows != 0 {
                detail.normalized_coverage =
                    (read_starts as f64 / windows as f64) / mean_reads_per_window;
                // Math.sqrt is the one Math function IEEE-754 specifies exactly, so Rust's
                // f64::sqrt is bit-identical to Java's and jmath is not needed here.
                detail.error_bar_width =
                    ((read_starts as f64).sqrt() / windows as f64) / mean_reads_per_window;
            }
            details.push(detail);
        }

        let (at_dropout, gc_dropout) = dropout_metrics(&details);
        let summary = GcBiasSummaryMetrics {
            accumulation_level: "All Reads".to_string(),
            reads_used: "ALL".to_string(),
            window_size: self.scan_window_size,
            total_clusters: self.all_reads.total_clusters,
            aligned_reads: self.all_reads.total_aligned_reads,
            at_dropout,
            gc_dropout,
            gc_nc_0_19: self.gc_norm_coverage(mean_reads_per_window, 0, 19),
            gc_nc_20_39: self.gc_norm_coverage(mean_reads_per_window, 20, 39),
            gc_nc_40_59: self.gc_norm_coverage(mean_reads_per_window, 40, 59),
            gc_nc_60_79: self.gc_norm_coverage(mean_reads_per_window, 60, 79),
            gc_nc_80_100: self.gc_norm_coverage(mean_reads_per_window, 80, 100),
            sample: None,
            library: None,
            read_group: None,
        };
        Some((details, summary))
    }

    /// `calculateGcNormCoverage`.
    ///
    /// Bins with no reference windows are skipped from **both** the numerator and the window
    /// total, so reads in a GC bin the reference never reaches contribute nothing at all.
    fn gc_norm_coverage(&self, mean_reads_per_window: f64, start: usize, end: usize) -> f64 {
        let mut windows_total = 0i64;
        let mut sum = 0.0;
        for i in start..=end {
            if self.windows_by_gc[i] != 0 {
                sum += self.all_reads.reads_by_gc[i] as f64;
                windows_total += self.windows_by_gc[i] as i64;
            }
        }
        if windows_total == 0 {
            0.0
        } else {
            sum / (windows_total as f64 * mean_reads_per_window)
        }
    }
}

/// `SAMRecord.getAlignmentEnd`: the last reference position the CIGAR covers.
fn alignment_end(rec: &BamRecord) -> i32 {
    let blocks = alignment_blocks(&rec.cigar, rec.alignment_start);
    match blocks.last() {
        Some(b) => b.reference_start + b.length - 1,
        None => rec.alignment_start,
    }
}

/// `calculateDropoutMetrics`.
///
/// Only **positive** dropouts are summed; a bin with more reads than its window share
/// contributes nothing rather than offsetting. And the split is `GC <= 50` for AT, so bin 50
/// itself counts as AT.
fn dropout_metrics(details: &[GcBiasDetailMetrics]) -> (f64, f64) {
    let total_reads: f64 = details.iter().map(|d| d.read_starts as f64).sum();
    let total_windows: f64 = details.iter().map(|d| d.windows as f64).sum();
    let mut at_dropout = 0.0;
    let mut gc_dropout = 0.0;
    for d in details {
        let relative_reads = d.read_starts as f64 / total_reads;
        let relative_windows = d.windows as f64 / total_windows;
        let dropout = (relative_windows - relative_reads) * 100.0;
        if dropout > 0.0 {
            if d.gc <= 50 {
                at_dropout += dropout;
            } else {
                gc_dropout += dropout;
            }
        }
    }
    (at_dropout, gc_dropout)
}
