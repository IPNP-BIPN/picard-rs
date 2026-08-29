//! What `MergeBamAlignment` clips, and what it unmaps.
//!
//! Three things shorten or remove an alignment on the way through: the contig it runs off the end
//! of, the adapter the unmapped bam marked in it, and the contamination filter that decides an
//! alignment is too short to believe.
//!
//! The filter counts two things rather than one. An alignment is a contaminant when it has fewer
//! aligned bases than `--MIN_UNCLIPPED_BASES` AND is soft-clipped at both ends, so the same ten
//! aligned bases clipped at one end only are never a contaminant however short they are.
//!
//! Ported from `picard.sam.AbstractAlignmentMerger`, `picard.filter.OverclippedReadFilter` and
//! `htsjdk.samtools.util.CigarUtil` in Picard 3.4.0 and htsjdk 4.2.0.

/// One cigar element.
pub type Element = (usize, char);

/// How many bases of the read a cigar consumes.
pub fn read_length(cigar: &[Element]) -> usize {
    cigar
        .iter()
        .filter(|(_, operator)| matches!(operator, 'M' | 'I' | 'S' | '=' | 'X'))
        .map(|(length, _)| length)
        .sum()
}

/// How many bases of the reference it consumes.
pub fn reference_length(cigar: &[Element]) -> usize {
    cigar
        .iter()
        .filter(|(_, operator)| matches!(operator, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|(length, _)| length)
        .sum()
}

/// `CigarUtil.softClipEndOfRead`: everything from a read position on becomes a soft clip.
///
/// The position is one-based and counted along the read.
pub fn soft_clip_end_of_read(clip_from: usize, cigar: &[Element]) -> Vec<Element> {
    let mut kept: Vec<Element> = Vec::new();
    let mut position = 1usize;
    let mut clipped = 0usize;
    for (length, operator) in cigar {
        let consumes_read = matches!(operator, 'M' | 'I' | 'S' | '=' | 'X');
        if !consumes_read {
            if position < clip_from {
                kept.push((*length, *operator));
            }
            continue;
        }
        let element_end = position + length - 1;
        if element_end < clip_from {
            kept.push((*length, *operator));
        } else if position >= clip_from {
            clipped += length;
        } else {
            let keep = clip_from - position;
            kept.push((keep, *operator));
            clipped += length - keep;
        }
        position += length;
    }
    // A deletion or a skip left at the end of the kept cigar makes no sense once what followed it
    // is clipped away.
    while matches!(kept.last(), Some((_, 'D')) | Some((_, 'N'))) {
        kept.pop();
    }
    if clipped > 0 {
        kept.push((clipped, 'S'));
    }
    kept
}

/// The cigar an alignment that runs off the end of its contig is given.
///
/// The clip is counted from the READ: the overhang is the number of reference bases past the end,
/// and the clip starts that many bases from the end of the read, LESS whatever soft clip the read
/// already carried there, so a read clipped once is not clipped twice.
pub fn clip_off_the_end(start: i32, cigar: &[Element], contig_length: i32) -> Option<Vec<Element>> {
    let alignment_end = start + reference_length(cigar) as i32 - 1;
    let overhang = alignment_end - contig_length;
    if overhang <= 0 {
        return None;
    }
    let mut clip_from = read_length(cigar) as i32 - overhang + 1;
    if let Some((length, 'S')) = cigar.last() {
        clip_from -= *length as i32;
    }
    Some(soft_clip_end_of_read(clip_from as usize, cigar))
}

/// `CigarUtil.softClip3PrimeEndOfRead`: the adapter's clip, from the base `XT` names.
///
/// On a negative-strand record the three-prime end of the read is the START of the stored bases,
/// so the clip is applied from the other end.
pub fn clip_adapter(cigar: &[Element], negative_strand: bool, xt: usize) -> Vec<Element> {
    if negative_strand {
        let reversed: Vec<Element> = cigar.iter().rev().copied().collect();
        let clipped = soft_clip_end_of_read(read_length(cigar) - xt + 2, &reversed);
        return clipped.iter().rev().copied().collect();
    }
    soft_clip_end_of_read(xt, cigar)
}

/// `OverclippedReadFilter.filterOut`: whether an alignment is too short to believe.
///
/// Consecutive soft clips count as one block, and it takes two blocks: a read clipped at one end
/// only is never a contaminant, whatever is left of it.
pub fn is_contaminant(cigar: &[Element], minimum_unclipped_bases: usize) -> bool {
    let mut aligned = 0usize;
    let mut soft_clip_blocks = 0usize;
    let mut last: Option<char> = None;
    for (length, operator) in cigar {
        if *operator == 'S' {
            if last != Some('S') {
                soft_clip_blocks += 1;
            }
        } else if matches!(operator, 'M' | 'I' | '=' | 'X') {
            aligned += length;
        }
        last = Some(*operator);
    }
    aligned < minimum_unclipped_bases && soft_clip_blocks >= 2
}

/// What an unmapped contaminant is left holding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnmappingReadStrategy {
    /// Keep the coordinates, drop the cigar and the mapping quality, write no tag.
    DoNotChange,
    /// Keep the coordinates and the tag, drop the cigar and the mapping quality.
    CopyToTag,
    /// Drop the coordinates too.
    MoveToTag,
    /// Keep everything, on a record flagged unmapped, which is not a valid record.
    DoNotChangeInvalid,
}

impl UnmappingReadStrategy {
    pub fn reset_mapping_information(&self) -> bool {
        matches!(self, UnmappingReadStrategy::MoveToTag)
    }
    pub fn populate_oa_tag(&self) -> bool {
        matches!(
            self,
            UnmappingReadStrategy::CopyToTag | UnmappingReadStrategy::MoveToTag
        )
    }
    /// Whether the record is left valid: an unmapped record may carry neither a cigar nor a
    /// mapping quality, and one strategy declines to enforce that.
    pub fn keep_valid(&self) -> bool {
        !matches!(self, UnmappingReadStrategy::DoNotChangeInvalid)
    }
}

/// The comment every unmapped contaminant carries.
pub const CONTAMINATION_COMMENT: &str = "Cross-species contamination";

/// One record, as much of it as the unmapping touches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub reference: Option<String>,
    pub start: i32,
    pub mapping_quality: i32,
    pub cigar: Vec<Element>,
    pub unmapped: bool,
    pub edit_distance: Option<i32>,
    pub original_alignment: Option<String>,
    pub comment: Option<String>,
}

/// `encodeMappingInformation`: the `OA` tag's contents, in the SA tag's format.
pub fn encode_mapping_information(record: &Record) -> String {
    format!(
        "{},{},{},{},{};",
        record.reference.clone().unwrap_or_default(),
        record.start,
        cigar_string(&record.cigar),
        record.mapping_quality,
        record
            .edit_distance
            .map(|distance| distance.to_string())
            .unwrap_or_default()
    )
}

/// A cigar as a string.
pub fn cigar_string(cigar: &[Element]) -> String {
    if cigar.is_empty() {
        return "*".to_string();
    }
    cigar
        .iter()
        .map(|(length, operator)| format!("{length}{operator}"))
        .collect()
}

/// Unmap a contaminant, keeping whatever the strategy says it keeps.
///
/// An existing comment is not replaced: the new one is appended after a pipe.
pub fn unmap_contaminant(record: &mut Record, strategy: UnmappingReadStrategy) {
    if strategy.populate_oa_tag() {
        record.original_alignment = Some(encode_mapping_information(record));
    }
    if strategy.reset_mapping_information() {
        record.reference = None;
        record.start = 0;
        record.edit_distance = None;
    }
    record.unmapped = true;
    if strategy.keep_valid() {
        record.mapping_quality = 0;
        record.cigar = Vec::new();
    }
    record.comment = Some(match &record.comment {
        Some(existing) => format!("{existing} | {CONTAMINATION_COMMENT}"),
        None => CONTAMINATION_COMMENT.to_string(),
    });
}
