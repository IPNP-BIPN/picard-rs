//! `ExtractIlluminaBarcodes`: which declared barcode a cluster is, and when it is none of them.
//!
//! The match is by EDIT DISTANCE rather than equality, and it is two tests rather than one. A
//! cluster is a barcode when it is within `MAX_MISMATCHES` of it AND the next best barcode is
//! further away by `MIN_MISMATCH_DELTA`. Two barcodes equidistant from a cluster therefore match
//! neither, however close they are.
//!
//! Ported from `picard.illumina.ExtractIlluminaBarcodes` and
//! `picard.illumina.ExtractBarcodesProgram` in Picard 3.4.0.

/// The arguments the match reads.
#[derive(Debug, Clone)]
pub struct Options {
    pub max_mismatches: i32,
    pub min_mismatch_delta: i32,
    pub max_no_calls: i32,
    pub minimum_base_quality: u8,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_mismatches: 1,
            min_mismatch_delta: 1,
            max_no_calls: 2,
            minimum_base_quality: 0,
        }
    }
}

/// What one cluster's barcode cycles carry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Observed {
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
}

/// What the tool decided about one cluster.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Match {
    /// Whether a barcode was matched at all, which is the `Y` or `N` of the per-tile file.
    pub matched: bool,
    /// The barcode's index in the declared list, where one was matched or nearly was.
    pub barcode: Option<usize>,
    pub mismatches: i32,
    pub mismatch_delta: i32,
}

/// `countMismatches`: a no-call counts as a mismatch, and a base below the floor is a no-call.
pub fn mismatches(observed: &Observed, barcode: &[u8], minimum_quality: u8) -> i32 {
    let mut count = 0;
    for (index, expected) in barcode.iter().enumerate() {
        let base = observed.bases.get(index).copied().unwrap_or(b'.');
        let quality = observed.qualities.get(index).copied().unwrap_or(0);
        if base == b'.' || quality < minimum_quality || base != *expected {
            count += 1;
        }
    }
    count
}

/// How many of a cluster's barcode bases are no-calls, which `MAX_NO_CALLS` caps.
pub fn no_calls(observed: &Observed) -> i32 {
    observed.bases.iter().filter(|base| **base == b'.').count() as i32
}

/// `findBestBarcode`: the nearest barcode, and whether it wins by enough.
///
/// The comparison keeps the FIRST of an equal-distance pair as the best and the second as the
/// runner-up, so the delta between them is zero and neither is matched.
pub fn best_match(observed: &Observed, barcodes: &[Vec<u8>], options: &Options) -> Match {
    let mut best: Option<usize> = None;
    let mut best_distance = i32::MAX;
    let mut next_distance = i32::MAX;
    for (index, barcode) in barcodes.iter().enumerate() {
        let distance = mismatches(observed, barcode, options.minimum_base_quality);
        if distance < best_distance {
            next_distance = best_distance;
            best_distance = distance;
            best = Some(index);
        } else if distance < next_distance {
            next_distance = distance;
        }
    }
    let delta = next_distance.saturating_sub(best_distance);
    let matched = best.is_some()
        && no_calls(observed) <= options.max_no_calls
        && best_distance <= options.max_mismatches
        && delta >= options.min_mismatch_delta;
    Match {
        matched,
        barcode: best,
        mismatches: best_distance,
        mismatch_delta: delta,
    }
}

/// One barcode's row of the metrics file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BarcodeMetrics {
    pub barcode: String,
    pub reads: i64,
    pub pf_reads: i64,
    pub perfect_matches: i64,
    pub pf_perfect_matches: i64,
    pub one_mismatch_matches: i64,
    pub pf_one_mismatch_matches: i64,
}

/// The metrics over a tile: one row per declared barcode, and one for everything else.
///
/// The unmatched row's barcode is a string of `N`s as long as the declared ones, which is what the
/// golden's `NN` row is, and the PF columns count only the clusters the filter file passed.
pub fn metrics(
    observations: &[(Observed, bool)],
    barcodes: &[Vec<u8>],
    options: &Options,
) -> Vec<BarcodeMetrics> {
    let length = barcodes.first().map(Vec::len).unwrap_or(0);
    let mut rows: Vec<BarcodeMetrics> = barcodes
        .iter()
        .map(|barcode| BarcodeMetrics {
            barcode: String::from_utf8_lossy(barcode).to_string(),
            ..BarcodeMetrics::default()
        })
        .collect();
    rows.push(BarcodeMetrics {
        barcode: "N".repeat(length),
        ..BarcodeMetrics::default()
    });
    let unmatched = rows.len() - 1;
    for (observed, passed) in observations {
        let decision = best_match(observed, barcodes, options);
        let row = if decision.matched {
            decision.barcode.expect("a matched barcode has an index")
        } else {
            unmatched
        };
        rows[row].reads += 1;
        if *passed {
            rows[row].pf_reads += 1;
        }
        if decision.matched && decision.mismatches == 0 {
            rows[row].perfect_matches += 1;
            if *passed {
                rows[row].pf_perfect_matches += 1;
            }
        }
        if decision.matched && decision.mismatches == 1 {
            rows[row].one_mismatch_matches += 1;
            if *passed {
                rows[row].pf_one_mismatch_matches += 1;
            }
        }
    }
    rows
}
