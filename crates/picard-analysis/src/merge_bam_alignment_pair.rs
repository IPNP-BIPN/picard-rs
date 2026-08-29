//! What `MergeBamAlignment` does to a PAIR, which is not what it does to two reads.
//!
//! The mate fields of each end are written from the other end's alignment, whatever the aligner
//! put there; the insert size is computed from both ends and signed by which starts first; the
//! proper-pair flag is the merger's own decision unless the caller asks for the aligner's; and
//! where the two ends overlap, one of them is clipped so the template's bases are not counted
//! twice.
//!
//! Ported from `picard.sam.AbstractAlignmentMerger` in Picard 3.4.0 and
//! `htsjdk.samtools.SamPairUtil` in htsjdk 4.2.0.

/// One end of a pair, as much of it as the fix-up touches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct End {
    pub reference: Option<String>,
    /// One-based; `0` where there is no alignment.
    pub start: i32,
    pub mapping_quality: i32,
    pub negative_strand: bool,
    pub unmapped: bool,
    /// The cigar as `(length, operator)`.
    pub cigar: Vec<(usize, char)>,
    pub bases: Vec<u8>,
    /// Phred, not ASCII.
    pub qualities: Vec<u8>,
    pub proper_pair: bool,
    pub mate_reference: Option<String>,
    pub mate_start: i32,
    pub mate_negative_strand: bool,
    pub mate_unmapped: bool,
    pub insert_size: i32,
    /// `MQ`, the mate's mapping quality.
    pub mate_mapping_quality: Option<i32>,
    /// `MC`, the mate's cigar.
    pub mate_cigar: Option<String>,
    /// `XB` and `XQ`, where hard clipping put what it removed.
    pub clipped_bases: Option<String>,
    pub clipped_qualities: Option<String>,
}

impl End {
    /// The cigar as a string.
    pub fn cigar_string(&self) -> String {
        if self.cigar.is_empty() {
            return "*".to_string();
        }
        self.cigar
            .iter()
            .map(|(length, operator)| format!("{length}{operator}"))
            .collect()
    }

    /// The last reference position the alignment covers.
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

    /// The first position the read would cover if its soft clips were aligned.
    pub fn unclipped_start(&self) -> i32 {
        self.start
            - self
                .cigar
                .iter()
                .take_while(|(_, operator)| matches!(operator, 'S' | 'H'))
                .filter(|(_, operator)| *operator == 'S')
                .map(|(length, _)| *length as i32)
                .sum::<i32>()
    }

    /// The last, likewise.
    pub fn unclipped_end(&self) -> i32 {
        self.end()
            + self
                .cigar
                .iter()
                .rev()
                .take_while(|(_, operator)| matches!(operator, 'S' | 'H'))
                .filter(|(_, operator)| *operator == 'S')
                .map(|(length, _)| *length as i32)
                .sum::<i32>()
    }

    /// How many bases the read holds.
    pub fn read_length(&self) -> usize {
        self.bases.len()
    }

    /// Whether the two alignments cover a common position.
    pub fn overlaps(&self, other: &End) -> bool {
        self.reference == other.reference && self.start <= other.end() && other.start <= self.end()
    }
}

/// `SamPairUtil.computeInsertSize`: the distance between the two five-prime ends, and one more,
/// signed by which end starts first.
pub fn insert_size(first: &End, second: &End) -> i32 {
    if first.unmapped || second.unmapped || first.reference != second.reference {
        return 0;
    }
    let five_prime = |end: &End| {
        if end.negative_strand {
            end.end()
        } else {
            end.start
        }
    };
    let first_position = five_prime(first);
    let second_position = five_prime(second);
    let adjustment = if second_position >= first_position {
        1
    } else {
        -1
    };
    second_position - first_position + adjustment
}

/// The orientations a pair may be read in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairOrientation {
    /// The two ends point at each other, which is what a library normally gives.
    FR,
    /// They point away from each other.
    RF,
    /// They point the same way.
    Tandem,
}

/// The orientations a proper pair may be in, unless the command line names others.
pub fn default_expected_orientations() -> Vec<PairOrientation> {
    vec![PairOrientation::FR]
}

/// `SamPairUtil.getPairOrientation`, read off one end and its mate fields.
pub fn pair_orientation(record: &End) -> PairOrientation {
    if record.negative_strand == record.mate_negative_strand {
        return PairOrientation::Tandem;
    }
    let positive_five_prime = if record.negative_strand {
        i64::from(record.mate_start)
    } else {
        i64::from(record.start)
    };
    let negative_five_prime = if record.negative_strand {
        i64::from(record.end())
    } else {
        i64::from(record.start) + i64::from(record.insert_size)
    };
    if positive_five_prime < negative_five_prime {
        PairOrientation::FR
    } else {
        PairOrientation::RF
    }
}

/// `SamPairUtil.isProperPair`: both ends placed, on one contig, in an expected orientation.
pub fn is_proper_pair(first: &End, second: &End, expected: &[PairOrientation]) -> bool {
    if first.unmapped || second.unmapped {
        return false;
    }
    if first.reference.is_none() || first.reference != second.reference {
        return false;
    }
    expected.contains(&pair_orientation(first))
}

/// `SamPairUtil.setMateInfo`: each end's mate fields from the other end, and the insert size from
/// both.
///
/// An end with no alignment is given its mate's coordinates rather than left at nought, and it is
/// the one that keeps `MC` and `MQ`: the mapped end has no mate cigar to carry, its mate having no
/// cigar.
pub fn set_mate_info(first: &mut End, second: &mut End, add_mate_cigar: bool) {
    if !first.unmapped && !second.unmapped {
        first.mate_reference = second.reference.clone();
        first.mate_start = second.start;
        first.mate_negative_strand = second.negative_strand;
        first.mate_unmapped = false;
        first.mate_mapping_quality = Some(second.mapping_quality);

        second.mate_reference = first.reference.clone();
        second.mate_start = first.start;
        second.mate_negative_strand = first.negative_strand;
        second.mate_unmapped = false;
        second.mate_mapping_quality = Some(first.mapping_quality);

        if add_mate_cigar {
            first.mate_cigar = Some(second.cigar_string());
            second.mate_cigar = Some(first.cigar_string());
        } else {
            first.mate_cigar = None;
            second.mate_cigar = None;
        }
    } else if first.unmapped && second.unmapped {
        let (first_negative, second_negative) = (first.negative_strand, second.negative_strand);
        for (end, other_negative) in [
            (&mut *first, second_negative),
            (&mut *second, first_negative),
        ] {
            end.reference = None;
            end.start = 0;
            end.mate_reference = None;
            end.mate_start = 0;
            end.mate_negative_strand = other_negative;
            end.mate_unmapped = true;
            end.mate_mapping_quality = None;
            end.mate_cigar = None;
            end.insert_size = 0;
        }
    } else {
        let (mapped, unmapped) = if first.unmapped {
            (&mut *second, &mut *first)
        } else {
            (&mut *first, &mut *second)
        };
        unmapped.reference = mapped.reference.clone();
        unmapped.start = mapped.start;

        mapped.mate_reference = unmapped.reference.clone();
        mapped.mate_start = unmapped.start;
        mapped.mate_negative_strand = unmapped.negative_strand;
        mapped.mate_unmapped = true;
        mapped.mate_mapping_quality = None;
        mapped.mate_cigar = None;
        mapped.insert_size = 0;

        unmapped.mate_reference = mapped.reference.clone();
        unmapped.mate_start = mapped.start;
        unmapped.mate_negative_strand = mapped.negative_strand;
        unmapped.mate_unmapped = false;
        unmapped.mate_mapping_quality = Some(mapped.mapping_quality);
        unmapped.mate_cigar = if add_mate_cigar {
            Some(mapped.cigar_string())
        } else {
            None
        };
        unmapped.insert_size = 0;
    }

    let size = insert_size(first, second);
    first.insert_size = size;
    second.insert_size = -size;
}

/// `SamPairUtil.setProperPairFlags`: the flag both ends carry, decided once.
pub fn set_proper_pair_flags(first: &mut End, second: &mut End, expected: &[PairOrientation]) {
    let proper = !first.unmapped && !second.unmapped && is_proper_pair(first, second, expected);
    first.proper_pair = proper;
    second.proper_pair = proper;
}

/// The read position a reference position falls at, with soft clips counted as aligned bases.
///
/// One-based, and nought where the position is outside the read.
pub fn read_position_at_reference_position_ignoring_soft_clips(
    record: &End,
    position: i32,
) -> usize {
    if position <= 0 {
        return 0;
    }
    // The soft clips become matches, and the alignment starts that many bases earlier.
    let mut leading = 0;
    for (length, operator) in &record.cigar {
        match operator {
            'S' => leading += *length as i32,
            'H' => {}
            _ => break,
        }
    }
    let mut reference = record.start - leading;
    let mut read = 1usize;
    for (length, operator) in &record.cigar {
        match operator {
            'M' | 'S' | '=' | 'X' => {
                for _ in 0..*length {
                    if reference == position {
                        return read;
                    }
                    reference += 1;
                    read += 1;
                }
            }
            'I' => read += *length,
            'D' | 'N' => {
                for _ in 0..*length {
                    if reference == position {
                        // A position inside a deletion answers with the base before it.
                        return read - 1;
                    }
                    reference += 1;
                }
            }
            _ => {}
        }
    }
    0
}

/// Clip an end from a read position to its three-prime end, softly or hard.
///
/// `clip_from` is one-based and counted from the five-prime end of the read as the machine read
/// it, so on a negative-strand record it is counted from the far end of the stored bases.
pub fn clip_3_prime_end(record: &mut End, clip_from: usize, hard: bool) {
    if hard {
        move_clipped_bases_to_tag(record, clip_from);
    }
    let operator = if hard { 'H' } else { 'S' };
    let length = record.read_length() - clip_from + 1;
    if record.negative_strand {
        clip_start(record, length, operator);
    } else {
        clip_end(record, length, operator);
    }
}

/// What hard clipping keeps: the bases and qualities it is about to remove, in the order the
/// sequencer produced them.
fn move_clipped_bases_to_tag(record: &mut End, clip_from: usize) {
    let read_length = record.read_length();
    let (from, to) = if record.negative_strand {
        (0, read_length - clip_from + 1)
    } else {
        (clip_from - 1, read_length)
    };
    let bases = String::from_utf8(record.bases[from..to].to_vec()).expect("bases are ASCII");
    let qualities: String = record.qualities[from..to]
        .iter()
        .map(|quality| (quality + 33) as char)
        .collect();
    if record.negative_strand {
        record.clipped_bases = Some(reverse_complement(&bases));
        record.clipped_qualities = Some(qualities.chars().rev().collect());
    } else {
        record.clipped_bases = Some(bases);
        record.clipped_qualities = Some(qualities);
    }
}

fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

/// Clip the end of the alignment, moving nothing.
fn clip_end(record: &mut End, length: usize, operator: char) {
    let mut remaining = length;
    let mut cigar: Vec<(usize, char)> = Vec::new();
    let mut kept: Vec<(usize, char)> = Vec::new();
    for element in record.cigar.iter().rev() {
        if remaining == 0 {
            kept.push(*element);
            continue;
        }
        let (element_length, element_operator) = *element;
        if element_operator == operator || element_operator == 'S' {
            remaining = remaining.saturating_sub(element_length);
            cigar.push((element_length, operator));
            continue;
        }
        if element_length <= remaining {
            remaining -= element_length;
            cigar.push((element_length, operator));
        } else {
            cigar.push((remaining, operator));
            kept.push((element_length - remaining, element_operator));
            remaining = 0;
        }
    }
    kept.reverse();
    let clipped: usize = cigar.iter().map(|(length, _)| length).sum();
    kept.push((clipped, operator));
    record.cigar = kept;
    if operator == 'H' {
        record.bases.truncate(record.bases.len() - length);
        record.qualities.truncate(record.qualities.len() - length);
    }
}

/// Clip the start of the alignment, which moves it forward.
fn clip_start(record: &mut End, length: usize, operator: char) {
    let mut remaining = length;
    let mut cigar: Vec<(usize, char)> = Vec::new();
    let mut consumed_reference = 0i32;
    let mut clipped = 0usize;
    for element in record.cigar.iter() {
        let (element_length, element_operator) = *element;
        if remaining == 0 {
            cigar.push(*element);
            continue;
        }
        if element_operator == operator || element_operator == 'S' {
            remaining = remaining.saturating_sub(element_length);
            clipped += element_length;
            continue;
        }
        if element_length <= remaining {
            remaining -= element_length;
            clipped += element_length;
            if matches!(element_operator, 'M' | 'D' | 'N' | '=' | 'X') {
                consumed_reference += element_length as i32;
            }
        } else {
            clipped += remaining;
            if matches!(element_operator, 'M' | 'D' | 'N' | '=' | 'X') {
                consumed_reference += remaining as i32;
            }
            cigar.push((element_length - remaining, element_operator));
            remaining = 0;
        }
    }
    let mut result = vec![(clipped, operator)];
    result.extend(cigar);
    record.cigar = result;
    record.start += consumed_reference;
    if operator == 'H' {
        record.bases.drain(0..length);
        record.qualities.drain(0..length);
    }
}

/// Clip whichever end reaches past its mate, so the template's bases are counted once.
///
/// The two ends have to be mapped, on opposite strands, and overlapping. Soft clipping equalises
/// the ALIGNED ends; hard clipping is applied on top of it, and equalises the unclipped ones.
pub fn clip_for_overlapping_reads(first: &mut End, second: &mut End, hard: bool) {
    if first.unmapped
        || second.unmapped
        || first.negative_strand == second.negative_strand
        || !first.overlaps(second)
    {
        return;
    }
    let first_is_negative = first.negative_strand;
    {
        let (positive, negative) = if first_is_negative {
            (&mut *second, &mut *first)
        } else {
            (&mut *first, &mut *second)
        };
        clip_3_prime_ends_to_5_prime_ends(positive, negative, false, false);
        if hard {
            clip_3_prime_ends_to_5_prime_ends(positive, negative, true, true);
        }
    }
}

fn clip_3_prime_ends_to_5_prime_ends(
    positive: &mut End,
    negative: &mut End,
    hard: bool,
    use_unclipped_ends: bool,
) {
    let negative_end = if use_unclipped_ends {
        negative.unclipped_end()
    } else {
        negative.end()
    };
    let positive_start = if use_unclipped_ends {
        positive.unclipped_start()
    } else {
        positive.start
    };

    let last_kept = read_position_at_reference_position_ignoring_soft_clips(positive, negative_end);
    if last_kept > 0 && last_kept < positive.read_length() {
        clip_3_prime_end(positive, last_kept + 1, hard);
    }

    let from_start =
        read_position_at_reference_position_ignoring_soft_clips(negative, positive_start - 1);
    let from_five_prime_end = if from_start > 0 {
        negative.read_length() + 1 - from_start
    } else {
        0
    };
    if from_five_prime_end > 0 {
        clip_3_prime_end(negative, from_five_prime_end, hard);
    }
}

/// What a run was asked for, as far as a pair is concerned.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    pub clip_overlapping_reads: bool,
    pub hard_clip_overlapping_reads: bool,
    pub add_mate_cigar: bool,
    /// Keep the aligner's proper-pair flags rather than deciding them again.
    pub aligner_proper_pair_flags: bool,
    pub expected_orientations: Vec<PairOrientation>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            clip_overlapping_reads: true,
            hard_clip_overlapping_reads: false,
            add_mate_cigar: true,
            aligner_proper_pair_flags: false,
            expected_orientations: default_expected_orientations(),
        }
    }
}

/// The whole fix-up, in the order the reference applies it.
///
/// The clipping happens FIRST, so the mate fields and the insert size are written from the clipped
/// alignments rather than from the aligner's.
pub fn fix_up_pair(first: &mut End, second: &mut End, options: &Options) {
    if options.clip_overlapping_reads {
        clip_for_overlapping_reads(first, second, options.hard_clip_overlapping_reads);
    }
    set_mate_info(first, second, options.add_mate_cigar);
    if !options.aligner_proper_pair_flags {
        set_proper_pair_flags(first, second, &options.expected_orientations);
    }
}
