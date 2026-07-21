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
//! Getting to that took two probes. The first put the no-calls at the start of the sequence,
//! where they only ever *left* the window, so the buggy branch never ran and the two cases
//! agreed. That probe proved nothing; the second one, with the no-calls entering, is the
//! evidence.

use htsjdk_bam::sequence::bases_equal;

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
