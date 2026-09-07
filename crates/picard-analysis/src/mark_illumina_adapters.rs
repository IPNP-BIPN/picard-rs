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

/// The nine adapter pairs the enum declares, as `(name, five prime, three prime)`, in the order
/// the enum declares them.
///
/// `ALTERNATIVE_SINGLE_END` is last on purpose: its three prime sequence is a suffix of several of
/// the others, so a list that tried it earlier would answer it for reads that belong to another
/// pair. The comment in `IlluminaUtil` says so.
pub const ADAPTER_PAIRS: [(&str, &str, &str); 9] = [
    (
        "PAIRED_END",
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "INDEXED",
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "SINGLE_END",
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "NEXTERA_V1",
        "AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG",
        "CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCGNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "NEXTERA_V2",
        "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        "CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "DUAL_INDEXED",
        "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "FLUIDIGM",
        "AATGATACGGCGACCACCGAGATCTACACTGACGACATGGTTCTACA",
        "AGACCAAGTCTCTGCTACCGTANNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TRUSEQ_SMALLRNA",
        "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC",
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "ALTERNATIVE_SINGLE_END",
        "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACGATC",
        "TCGTATGCCGTCTTCTGCTTG",
    ),
];

/// `SequenceUtil.reverseComplement`, which is what makes a five prime adapter searchable.
///
/// The five prime sequence is declared as the oligo, and a read that runs into it reads its
/// complement backwards, so the sequence looked for in read two is the reverse complement of the
/// declared one. The enum computes it in its constructor; nothing searches the declared form.
pub fn reverse_complement(sequence: &str) -> String {
    let mut bases = sequence.as_bytes().to_vec();
    htsjdk_bam::sequence::reverse_complement(&mut bases);
    String::from_utf8(bases).expect("a complemented base is still ASCII")
}

/// One adapter pair as the search uses it: both sequences in READ order.
///
/// `AdapterMarker` never searches with the pairs it was given. It searches with truncated copies,
/// and this is one of those: the name is the truncated pair's name, and the sequences are already
/// cut and already in read order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    pub name: String,
    pub three_prime_read_order: String,
    pub five_prime_read_order: String,
}

/// The pair a name stands for, with the five prime side turned into read order.
pub fn adapter_pair(name: &str) -> Option<Adapter> {
    ADAPTER_PAIRS
        .iter()
        .find(|(pair_name, _, _)| *pair_name == name)
        .map(|(pair_name, five, three)| Adapter {
            name: (*pair_name).to_string(),
            three_prime_read_order: (*three).to_string(),
            five_prime_read_order: reverse_complement(five),
        })
}

/// `AdapterMarker.DEFAULT_ADAPTER_LENGTH`: every adapter is cut to this before anything is
/// searched.
pub const DEFAULT_ADAPTER_LENGTH: usize = 30;

/// `AdapterMarker.substringAndRemoveTrailingNs`: cut to `length`, then keep cutting while the last
/// base is a no-call.
///
/// The second half is what makes the cut position depend on the sequence: `NEXTERA_V2`'s five
/// prime in read order has its Ns near the front, so cutting at thirty leaves a real base last and
/// nothing more is removed, while a pair whose Ns straddle the cut comes back shorter than thirty.
pub fn substring_and_remove_trailing_ns(sequence: &str, length: usize) -> String {
    let bases = sequence.as_bytes();
    let mut length = length.min(bases.len());
    while length > 0 && is_no_call(bases[length - 1]) {
        length -= 1;
    }
    sequence[..length].to_string()
}

/// `new AdapterMarker(length, adapters)`: truncate every pair, then COLLAPSE the ones that became
/// the same.
///
/// The collapse is the part with teeth. `INDEXED` and `DUAL_INDEXED` differ only past base thirty
/// on both sides, so the default list of three becomes a list of TWO, and a read that would have
/// been attributed to the second pair is attributed to the first. The name records the merge
/// (`truncated INDEXED|DUAL_INDEXED`) and is the only place the difference survives.
pub fn truncated_adapters(names: &[&str], length: usize) -> Vec<Adapter> {
    let mut truncated: Vec<Adapter> = Vec::new();
    for name in names {
        let Some(pair) = adapter_pair(name) else {
            continue;
        };
        let candidate = Adapter {
            name: format!("truncated {}", pair.name),
            three_prime_read_order: substring_and_remove_trailing_ns(
                &pair.three_prime_read_order,
                length,
            ),
            five_prime_read_order: substring_and_remove_trailing_ns(
                &pair.five_prime_read_order,
                length,
            ),
        };
        // `equals` on the truncated pair ignores the name, so this compares the two sequences and
        // nothing else.
        match truncated.iter_mut().find(|existing| {
            existing.three_prime_read_order == candidate.three_prime_read_order
                && existing.five_prime_read_order == candidate.five_prime_read_order
        }) {
            Some(existing) => existing.name = format!("{}|{}", existing.name, name),
            None => truncated.push(candidate),
        }
    }
    truncated
}

/// The bases as the search sees them: reverse complemented when the record is on the negative
/// strand.
///
/// `ClippingUtility.getReadBases` does this, so an aligned input is searched in read order however
/// it was mapped, and a copy is complemented rather than the record: what is written out is still
/// the record's own bases.
pub fn read_bases_in_read_order(bases: &[u8], negative_strand: bool) -> Vec<u8> {
    let mut copy = bases.to_vec();
    if negative_strand {
        htsjdk_bam::sequence::reverse_complement(&mut copy);
    }
    copy
}

/// `adapterTrimIlluminaSingleRead`: the first pair whose THREE prime side matches, and where.
pub fn trim_single_read(
    bases: &[u8],
    adapters: &[Adapter],
    min_match: usize,
    max_error_rate: f64,
) -> Option<(usize, i32)> {
    for (position, adapter) in adapters.iter().enumerate() {
        let index = find_index_of_clip_sequence(
            bases,
            adapter.three_prime_read_order.as_bytes(),
            min_match,
            max_error_rate,
        );
        if index != NO_MATCH {
            return Some((position, index));
        }
    }
    None
}

/// What a pair's search decided: a tag for each read, and which pair was credited.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PairedTrim {
    pub first_tag: Option<i32>,
    pub second_tag: Option<i32>,
    pub matched: Option<usize>,
}

/// `attemptOneSidedMatch`: re-check the read that DID match against a stricter minimum, and mark
/// both reads at that position if it survives.
///
/// The stricter minimum is not applied by searching again. It is applied to the match already
/// found, as a length: the matched read must carry at least `stricter_min_match` bases from the
/// match onwards. The other read is then marked at the same index whether or not anything was
/// found in it, which is how a pair comes back with two tags when only one read carried an
/// adapter.
fn attempt_one_sided_match(
    first: &[u8],
    second: &[u8],
    index1: i32,
    index2: i32,
    stricter_min_match: usize,
    trim: &mut PairedTrim,
) -> bool {
    let matched_index = if index1 == NO_MATCH { index2 } else { index1 };
    let matched_length = if index1 == NO_MATCH {
        second.len()
    } else {
        first.len()
    };
    if matched_length as i32 - matched_index >= stricter_min_match as i32 {
        if first.len() as i32 > matched_index {
            trim.first_tag = Some(matched_index + 1);
        }
        if second.len() as i32 > matched_index {
            trim.second_tag = Some(matched_index + 1);
        }
        return true;
    }
    false
}

/// `adapterTrimIlluminaPairedReads`: one index for the pair, from two searches that have to agree.
///
/// Read one is searched for the THREE prime sequence and read two for the FIVE prime one in read
/// order, and the two indices are compared. Three outcomes, and only the first ends the loop:
///
///  * **equal and found**: both reads are marked there and the pair is the answer, immediately;
///  * **exactly one found**: the one-sided path may mark both reads, but the loop KEEPS GOING,
///    so a later pair that matches exactly overwrites what it wrote;
///  * **both found, at different places**: nothing at all, and the tags an earlier pair wrote
///    stay.
///
/// The tags are therefore cumulative across the list rather than decided once, and the returned
/// pair is the LAST one-sided match when no pair matched exactly.
pub fn trim_paired_reads(
    first: &[u8],
    second: &[u8],
    adapters: &[Adapter],
    min_match: usize,
    max_error_rate: f64,
) -> PairedTrim {
    let mut trim = PairedTrim::default();
    for (position, adapter) in adapters.iter().enumerate() {
        let index1 = find_index_of_clip_sequence(
            first,
            adapter.three_prime_read_order.as_bytes(),
            min_match,
            max_error_rate,
        );
        let index2 = find_index_of_clip_sequence(
            second,
            adapter.five_prime_read_order.as_bytes(),
            min_match,
            max_error_rate,
        );
        if index1 == index2 {
            if index1 != NO_MATCH {
                trim.first_tag = Some(index1 + 1);
                trim.second_tag = Some(index2 + 1);
                trim.matched = Some(position);
                return trim;
            }
        } else if (index1 == NO_MATCH || index2 == NO_MATCH)
            && attempt_one_sided_match(first, second, index1, index2, 2 * min_match, &mut trim)
        {
            trim.matched = Some(position);
        }
    }
    trim
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The default list of three is a list of TWO by the time anything is searched.
    ///
    /// `INDEXED` and `DUAL_INDEXED` differ only past base thirty on both sides -- the dual pair's
    /// eight Ns sit at bases 29 to 36 of its five prime sequence, which the cut removes -- so the
    /// truncation makes them the same adapter and the second is folded into the first.
    #[test]
    fn the_default_list_collapses_to_two_adapters() {
        let adapters = truncated_adapters(&DEFAULT_ADAPTERS, DEFAULT_ADAPTER_LENGTH);
        assert_eq!(adapters.len(), 2);
        assert_eq!(adapters[0].name, "truncated INDEXED|DUAL_INDEXED");
        assert_eq!(adapters[1].name, "truncated PAIRED_END");
        assert_eq!(
            adapters[0].three_prime_read_order,
            "AGATCGGAAGAGCACACGTCTGAACTCCAG"
        );
        assert_eq!(
            adapters[1].three_prime_read_order,
            "AGATCGGAAGAGCGGTTCAGCAGGAATGCC"
        );
    }

    /// The five prime sequence is searched reverse complemented, and every default pair's is the
    /// same one: the three pairs are told apart by their three prime side alone.
    #[test]
    fn the_five_prime_side_is_searched_in_read_order() {
        let indexed = adapter_pair("INDEXED").expect("a known pair");
        assert!(indexed
            .five_prime_read_order
            .starts_with("AGATCGGAAGAGCGTCGTGTAGGGAAAGAG"));
        assert_eq!(
            substring_and_remove_trailing_ns(&indexed.five_prime_read_order, 30),
            substring_and_remove_trailing_ns(
                &adapter_pair("PAIRED_END")
                    .expect("a known pair")
                    .five_prime_read_order,
                30
            )
        );
    }

    /// The cut is not thirty bases: it is thirty bases and then back past the Ns.
    ///
    /// `FLUIDIGM`'s three prime sequence carries ten Ns from base 22, so cutting at thirty leaves
    /// eight of them last and all eight go: the sequence searched is TWENTY-TWO bases long, which
    /// is what the error allowance is computed from.
    #[test]
    fn trailing_no_calls_shorten_the_cut() {
        let fluidigm = adapter_pair("FLUIDIGM").expect("a known pair");
        assert_eq!(
            substring_and_remove_trailing_ns(&fluidigm.three_prime_read_order, 30),
            "AGACCAAGTCTCTGCTACCGTA"
        );
    }

    /// A pair where one read carries the adapter comes back with TWO tags at the same position.
    #[test]
    fn one_sided_match_marks_both_reads() {
        let adapters = truncated_adapters(&["PAIRED_END"], DEFAULT_ADAPTER_LENGTH);
        let first = b"TTTTTTTTTTTTTTTTTTTTTTTTTTAGATCGGAAGAGCGGTTCAGCAGG".to_vec();
        let second = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec();
        let trim = trim_paired_reads(
            &first,
            &second,
            &adapters,
            MIN_MATCH_PE_BASES,
            MAX_PE_ERROR_RATE,
        );
        assert_eq!(trim.first_tag, Some(27));
        assert_eq!(trim.second_tag, Some(27));
    }
}
