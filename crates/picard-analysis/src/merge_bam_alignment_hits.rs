//! Which alignment `MergeBamAlignment` calls primary.
//!
//! An aligner may report a read several times. One of those alignments has to be made primary, the
//! rest marked secondary, and four strategies make that choice differently: three of them take the
//! mapping quality and one takes the alignment that maps the read's earliest base.
//!
//! The strategies also differ over the ALIGNER's own choice. `BestMapq` looks first at whether the
//! aligner named a primary and leaves it alone if it named exactly one; `EarliestFragment` never
//! looks, and chooses again from scratch.
//!
//! Where two alignments tie, the reference picks one of them with a seeded generator it does not
//! expose. That draw is not ported: the functions here return every alignment that tied, and the
//! measurement is built so that nothing ties.
//!
//! Ported from `picard.sam.BestMapqPrimaryAlignmentSelectionStrategy`,
//! `picard.sam.EarliestFragmentPrimaryAlignmentSelectionStrategy`,
//! `picard.sam.MostDistantPrimaryAlignmentSelectionStrategy`,
//! `picard.sam.BestEndMapqPrimaryAlignmentStrategy` and `picard.sam.HitsForInsert` in Picard 3.4.0,
//! and from `htsjdk.samtools.SAMUtils` in htsjdk 4.2.0.

/// One alignment of one end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub reference: String,
    pub start: i32,
    pub mapping_quality: i32,
    pub cigar: Vec<(usize, char)>,
    pub negative_strand: bool,
    /// Whether the aligner marked this one as secondary.
    pub secondary: bool,
}

impl Alignment {
    pub fn end(&self) -> i32 {
        self.start
            + self
                .cigar
                .iter()
                .filter(|(_, operator)| matches!(operator, 'M' | 'D' | 'N' | '=' | 'X'))
                .map(|(length, _)| *length as i32)
                .sum::<i32>()
            - 1
    }
}

/// One hit: the alignments of one end, or of both.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hit {
    pub first: Option<Alignment>,
    pub second: Option<Alignment>,
}

/// How many primaries the aligner named on one end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumPrimaryAlignmentState {
    None,
    One,
    MoreThanOne,
}

/// `HitsForInsert.tallyPrimaryAlignments`.
pub fn tally_primary_alignments(hits: &[Hit], first_end: bool) -> NumPrimaryAlignmentState {
    let mut primaries = 0;
    for hit in hits {
        let alignment = if first_end { &hit.first } else { &hit.second };
        if let Some(alignment) = alignment {
            if !alignment.secondary {
                primaries += 1;
            }
        }
    }
    match primaries {
        0 => NumPrimaryAlignmentState::None,
        1 => NumPrimaryAlignmentState::One,
        _ => NumPrimaryAlignmentState::MoreThanOne,
    }
}

/// `SAMUtils.combineMapqs`: a hundred times each, and the unknown quality of 255 counted as one.
pub fn combine_mapqs(first: i32, second: i32) -> i32 {
    let scale = |quality: i32| if quality == 255 { 1 } else { quality * 100 };
    scale(first) + scale(second)
}

/// The mapping quality of one hit, both ends combined; a missing end counts as nought.
pub fn hit_mapq(hit: &Hit) -> i32 {
    combine_mapqs(
        hit.first.as_ref().map_or(0, |end| end.mapping_quality),
        hit.second.as_ref().map_or(0, |end| end.mapping_quality),
    )
}

/// The four strategies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    BestMapq,
    EarliestFragment,
    BestEndMapq,
    MostDistant,
}

/// Whether the aligner's own choice stands, which is what `BestMapq` asks before choosing.
///
/// It stands only when exactly one primary was named on each end that has alignments: none on
/// both, or more than one on either, and the strategy chooses again.
pub fn aligners_choice_stands(hits: &[Hit]) -> bool {
    let first = tally_primary_alignments(hits, true);
    let second = tally_primary_alignments(hits, false);
    let neither =
        first == NumPrimaryAlignmentState::None && second == NumPrimaryAlignmentState::None;
    let too_many = first == NumPrimaryAlignmentState::MoreThanOne
        || second == NumPrimaryAlignmentState::MoreThanOne;
    !(neither || too_many)
}

/// The hits with the best combined mapping quality.
pub fn best_mapq_candidates(hits: &[Hit]) -> Vec<usize> {
    let mut best = -1;
    let mut candidates = Vec::new();
    for (index, hit) in hits.iter().enumerate() {
        let quality = hit_mapq(hit);
        if quality > best {
            best = quality;
            candidates.clear();
        }
        if quality == best {
            candidates.push(index);
        }
    }
    candidates
}

/// The read base an alignment starts at, counted from the five-prime end of the read.
///
/// A leading soft or hard clip is what pushes it up, and on a negative-strand record the clip that
/// counts is the trailing one, the read having been stored the other way round.
pub fn index_of_first_aligned_base(alignment: &Alignment) -> usize {
    let elements: Vec<(usize, char)> = if alignment.negative_strand {
        alignment.cigar.iter().rev().copied().collect()
    } else {
        alignment.cigar.clone()
    };
    let mut index = 0;
    for (length, operator) in elements {
        if matches!(operator, 'S' | 'H') {
            index += length;
        } else {
            break;
        }
    }
    index
}

/// The hits whose fragment maps the earliest base of the read, the mapping quality breaking ties.
///
/// This strategy reads a FRAGMENT, so a paired read is refused rather than chosen from.
pub fn earliest_fragment_candidates(hits: &[Hit]) -> Vec<usize> {
    let mut earliest = usize::MAX;
    let mut best_mapq = -1;
    let mut candidates = Vec::new();
    for (index, hit) in hits.iter().enumerate() {
        let Some(fragment) = hit.first.as_ref() else {
            continue;
        };
        let first_base = index_of_first_aligned_base(fragment);
        let quality = fragment.mapping_quality;
        if first_base < earliest || (first_base == earliest && quality > best_mapq) {
            candidates.clear();
            candidates.push(index);
            earliest = first_base;
            best_mapq = quality;
        } else if first_base == earliest && quality == best_mapq {
            candidates.push(index);
        }
    }
    candidates
}

/// The refusal `EarliestFragment` answers a paired read with.
pub fn earliest_fragment_refusal(read_name: &str, description: &str) -> String {
    format!("getFragment called for paired read: {read_name} {description}")
}

/// How far apart a pairing reaches: the first base of one end to the last of the other.
pub fn pair_distance(first: &Alignment, second: &Alignment) -> i32 {
    let start = first.start.min(second.start);
    let end = first.end().max(second.end());
    end - start + 1
}

/// The pairings that reach furthest, the combined mapping quality breaking ties.
///
/// Only pairings on one contig are considered: a chimeric pairing has no distance to compare.
pub fn most_distant_candidates(hits: &[Hit]) -> Vec<usize> {
    let mut best_distance = -1;
    let mut best_mapq = -1;
    let mut candidates = Vec::new();
    for (index, hit) in hits.iter().enumerate() {
        let (Some(first), Some(second)) = (hit.first.as_ref(), hit.second.as_ref()) else {
            continue;
        };
        if first.reference != second.reference {
            continue;
        }
        let distance = pair_distance(first, second);
        let quality = combine_mapqs(first.mapping_quality, second.mapping_quality);
        if distance > best_distance || (distance == best_distance && quality > best_mapq) {
            best_distance = distance;
            best_mapq = quality;
            candidates.clear();
            candidates.push(index);
        } else if distance == best_distance && quality == best_mapq {
            candidates.push(index);
        }
    }
    // Where every pairing would be chimeric there is nothing to compare, and the strategy falls
    // back to each end's own best mapping quality.
    if candidates.is_empty() {
        return best_mapq_candidates(hits);
    }
    candidates
}

/// The hits any one strategy would accept as primary.
///
/// A single answer is the usual case; more than one means the reference would have drawn between
/// them.
pub fn primary_candidates(hits: &[Hit], strategy: Strategy) -> Vec<usize> {
    match strategy {
        // Three of the four look at what the aligner said first; only this one is written to.
        Strategy::BestMapq => {
            if aligners_choice_stands(hits) {
                return hits
                    .iter()
                    .enumerate()
                    .filter(|(_, hit)| {
                        hit.first.as_ref().is_some_and(|end| !end.secondary)
                            || hit.second.as_ref().is_some_and(|end| !end.secondary)
                    })
                    .map(|(index, _)| index)
                    .collect();
            }
            best_mapq_candidates(hits)
        }
        Strategy::EarliestFragment => earliest_fragment_candidates(hits),
        Strategy::BestEndMapq => best_mapq_candidates(hits),
        Strategy::MostDistant => most_distant_candidates(hits),
    }
}

/// What is written out of a read's alignments, given which one was made primary.
///
/// The losers are written as secondary rather than dropped, unless the run says otherwise.
pub fn written(hits: &[Hit], primary: usize, include_secondary: bool) -> Vec<usize> {
    if include_secondary {
        (0..hits.len()).collect()
    } else {
        vec![primary]
    }
}
