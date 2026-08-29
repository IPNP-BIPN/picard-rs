//! `MarkIlluminaAdapters`: where an adapter starts, as the `XT` tag counts it.
//!
//! The tool marks rather than clips. The tag is the ONE-BASED position an adapter begins at, and
//! the histogram counts the bases a read would lose if it were clipped there.
//!
//! Ported from `picard.illumina.MarkIlluminaAdapters`, `picard.util.ClippingUtility` and
//! `picard.util.IlluminaUtil.IlluminaAdapterPair`.

/// `ClippingUtility.MIN_MATCH_BASES`, for a single-end run.
pub const MIN_MATCH_BASES: usize = 12;
/// `ClippingUtility.MIN_MATCH_PE_BASES`, which is HALF of it: a pair is matched twice, so each
/// read is allowed to carry less of the adapter than a single read must.
pub const MIN_MATCH_PE_BASES: usize = 6;
/// `ClippingUtility.MAX_ERROR_RATE` and its paired twin, which are the same number.
pub const MAX_ERROR_RATE: f64 = 0.10;
pub const MAX_PE_ERROR_RATE: f64 = 0.10;

/// `ClippingUtility.NO_MATCH`.
pub const NO_MATCH: i32 = -1;

/// The three adapter pairs the tool tries by default, in the order it tries them.
///
/// The enum declares nine and the default list is three of them: the first that matches wins, so
/// the order is part of the answer.
pub const DEFAULT_ADAPTERS: [&str; 3] = ["INDEXED", "DUAL_INDEXED", "PAIRED_END"];

/// The three-prime sequence of each pair the default list names.
pub fn three_prime(name: &str) -> Option<&'static str> {
    match name {
        "PAIRED_END" => Some("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"),
        "INDEXED" => Some("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),
        "SINGLE_END" => Some("AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"),
        "DUAL_INDEXED" => {
            Some("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
        }
        _ => None,
    }
}

/// `SequenceUtil.isNoCall`, which is what makes an `N` in an adapter match anything.
fn is_no_call(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'N' | b'.')
}

fn bases_equal(left: u8, right: u8) -> bool {
    left.eq_ignore_ascii_case(&right)
}

/// `findIndexOfClipSequence`: the ZERO-based start of the adapter, or [`NO_MATCH`].
///
/// Three things decide the answer and none of them is obvious from the name.
///
/// The loop runs from `read.len() - min_match` DOWN to zero and returns the first start it can,
/// which is the LAST position in the read: a repeated adapter prefix is found at its last
/// occurrence and not at its first.
///
/// The comparison length is the OVERLAP, `min(read.len() - start, adapter.len())`, so a match near
/// the end of the read compares fewer bases than one further in, and the error allowance is
/// computed from that shorter length.
///
/// The allowance itself is `(int)(length * rate)`, truncated, so twelve bases at a tenth allow one
/// mismatch and nine bases allow none.
pub fn find_index_of_clip_sequence(
    read: &[u8],
    adapter: &[u8],
    min_match: usize,
    max_error_rate: f64,
) -> i32 {
    if read.len() < min_match {
        return NO_MATCH;
    }
    let mut start = read.len() - min_match;
    loop {
        let length = (read.len() - start).min(adapter.len());
        let allowed = (length as f64 * max_error_rate) as usize;
        let mut mismatches = 0;
        let mut matched = true;
        for index in 0..length {
            if !is_no_call(adapter[index]) && !bases_equal(adapter[index], read[start + index]) {
                mismatches += 1;
                if mismatches > allowed {
                    matched = false;
                    break;
                }
            }
        }
        if matched {
            return start as i32;
        }
        if start == 0 {
            return NO_MATCH;
        }
        start -= 1;
    }
}

/// The `XT` tag's value, which is the index PLUS ONE.
pub fn xt_tag(index: i32) -> Option<i32> {
    if index == NO_MATCH {
        None
    } else {
        Some(index + 1)
    }
}

/// The first adapter of a list that matches, which is what the default list's order decides.
pub fn first_matching_adapter(
    read: &[u8],
    adapters: &[&str],
    min_match: usize,
    max_error_rate: f64,
) -> Option<(usize, i32)> {
    for (position, name) in adapters.iter().enumerate() {
        let Some(sequence) = three_prime(name) else {
            continue;
        };
        let index =
            find_index_of_clip_sequence(read, sequence.as_bytes(), min_match, max_error_rate);
        if index != NO_MATCH {
            return Some((position, index));
        }
    }
    None
}

/// How many bases a read marked at `tag` would lose, which is what the histogram counts.
pub fn clipped_bases(read_length: usize, tag: i32) -> usize {
    read_length - (tag as usize - 1)
}

/// Whether a pair whose reads disagree is marked read by read, which it is not.
///
/// The paired path finds ONE index for the pair and writes it onto both reads, so a pair where one
/// read carries an adapter and the other does not comes back with two tags at the same position.
pub const A_PAIR_IS_MARKED_AS_A_PAIR: bool = true;

/// Whether a tag the input already carried survives a run that finds no adapter, which it does
/// not.
///
/// The tag is SET from the search's answer rather than merged with what was there, so a file
/// marked twice with two adapter lists carries the second run's answer and not the union.
pub const AN_EXISTING_TAG_SURVIVES: bool = false;
