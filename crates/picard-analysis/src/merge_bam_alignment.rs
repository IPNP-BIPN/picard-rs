//! `MergeBamAlignment`, the core alignment-info transfer.
//!
//! `MergeBamAlignment` merges an unmapped BAM (the original reads, bases/qualities/tags) with an
//! aligned BAM (the aligner's output for the same reads), producing reads that carry the alignment
//! while keeping their original sequence and user tags. The engine is `AbstractAlignmentMerger`; this
//! module opens the port with its innermost named symbol, `setValuesFromAlignment`
//! (`AbstractAlignmentMerger.java:946`), which copies one aligned record's alignment onto its unmapped
//! mate, plus the `isReservedTag` predicate that decides which aligner tags carry over.
//!
//! `setValuesFromAlignment`, for default options (`ATTRIBUTES_TO_RETAIN` / `ATTRIBUTES_TO_REMOVE`
//! empty), does:
//!   - copy every **non-reserved** aligner attribute onto the read (a reserved tag starts with a
//!     lower-case letter or one of `X`/`Y`/`Z`, so aligner-private tags like `X0` do not leak, while
//!     standard tags like `NM` do);
//!   - set the unmapped, reference, start, negative-strand, secondary and supplementary state from the
//!     alignment (the reference is set by **name**, so the two files' dictionaries may differ in
//!     order);
//!   - set the cigar and mapping quality when the alignment is mapped;
//!   - set the proper-pair flag when the read is paired;
//!   - reverse-complement the bases and reverse the qualities (and reverse-complement `E2`/`SQ`,
//!     reverse `OQ`/`U2`) when the alignment is on the negative strand, so the stored read is in
//!     reference orientation.
//!
//! Scope of this slice: the single-record transfer only. The surrounding engine, all subsequent
//! slices, is the read matching (`MultiHitAlignedReadIterator`), the paired mate-info and clipped-pair
//! fixing, the off-end-of-reference cigar clipping, the reference-based `NM`/`MD`/`UQ` recomputation
//! (which htsjdk-rs already provides via [`htsjdk_bam::md_nm`]), the program-record (`PG`) linkage, and
//! the merged-header construction. This function reproduces the transfer those stages build on.

use htsjdk_bam::header::SequenceRecord;
use htsjdk_bam::pair::set_mate_info;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{reverse, reverse_complement, reverse_qualities};
use htsjdk_bam::tag::{Tag, TagValue};

use crate::set_nm_md_and_uq_tags::fix_nm_md_and_uq;

const READ_PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const READ_UNMAPPED: u16 = 0x4;
const READ_REVERSE_STRAND: u16 = 0x10;
const MATE_NEGATIVE_STRAND: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const NOT_PRIMARY_ALIGNMENT: u16 = 0x100;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;

/// Tags whose values are reverse-complemented on a negative-strand read (`SAMRecord.TAGS_TO_REVERSE_COMPLEMENT`).
const TAGS_TO_REVERSE_COMPLEMENT: [&[u8; 2]; 2] = [b"E2", b"SQ"];
/// Tags whose values are merely reversed on a negative-strand read (`SAMRecord.TAGS_TO_REVERSE`).
const TAGS_TO_REVERSE: [&[u8; 2]; 2] = [b"OQ", b"U2"];

/// Why the transfer could not run.
#[derive(Debug)]
pub enum MergeAlignmentError {
    /// The read taken from the unmapped BAM is itself mapped (`setValuesFromAlignment` throws
    /// `UNMAPPED_BAM contains mapped reads`).
    UnmappedBamContainsMappedRead(String),
}

/// `AbstractAlignmentMerger.isReservedTag`: a tag is reserved (and so is not copied from the aligner
/// unless explicitly retained) when its first character is lower-case or one of `X`/`Y`/`Z`.
pub fn is_reserved_tag(code: &[u8; 2]) -> bool {
    let first = code[0];
    first.is_ascii_lowercase() || matches!(first, b'X' | b'Y' | b'Z')
}

/// Set or clear a flag bit to match `on`.
fn set_flag(flags: &mut u16, mask: u16, on: bool) {
    if on {
        *flags |= mask;
    } else {
        *flags &= !mask;
    }
}

/// Reverse-complement or reverse the string-valued reverse tags in place, as
/// `SAMRecord.reverseComplement` does alongside the bases and qualities.
fn reverse_complement_tags(record: &mut BamRecord) {
    for code in TAGS_TO_REVERSE_COMPLEMENT {
        if let Some(TagValue::Str(v)) = record.tags.get(Tag::new(code)) {
            let mut bytes = v.clone().into_bytes();
            reverse_complement(&mut bytes);
            let value = String::from_utf8(bytes).expect("reverse-complement of ASCII stays ASCII");
            record.tags.insert(Tag::new(code), TagValue::Str(value));
        }
    }
    for code in TAGS_TO_REVERSE {
        if let Some(TagValue::Str(v)) = record.tags.get(Tag::new(code)) {
            let mut bytes = v.clone().into_bytes();
            let len = bytes.len();
            reverse(&mut bytes, 0, len);
            let value = String::from_utf8(bytes).expect("reversal keeps ASCII");
            record.tags.insert(Tag::new(code), TagValue::Str(value));
        }
    }
}

/// `AbstractAlignmentMerger.setValuesFromAlignment` for default options: copy `aligned`'s alignment
/// onto `unmapped` (which keeps its own bases, qualities and user tags). `aligned_reference_name` is
/// the aligned record's `RNAME` (as resolved against the aligned file's dictionary, `"*"` when
/// unmapped), and `out_sequences` is the output header's dictionary, against which the reference is
/// re-resolved by name so the two files' dictionaries may differ in order.
pub fn transfer_alignment_info_to_fragment(
    unmapped: &mut BamRecord,
    aligned: &BamRecord,
    aligned_reference_name: &str,
    out_sequences: &[SequenceRecord],
) -> Result<(), MergeAlignmentError> {
    if unmapped.flags & READ_UNMAPPED == 0 {
        return Err(MergeAlignmentError::UnmappedBamContainsMappedRead(
            unmapped.read_name.clone(),
        ));
    }

    // Copy over any non-reserved aligner attributes (attributesToRetain/Remove default empty).
    for (tag, value) in aligned.tags.iter() {
        if !is_reserved_tag(&tag.name()) {
            unmapped.tags.insert(*tag, value.clone());
        }
    }

    let aligned_unmapped = aligned.flags & READ_UNMAPPED != 0;
    set_flag(&mut unmapped.flags, READ_UNMAPPED, aligned_unmapped);

    // Reference is set by name, not index, in case the dictionaries differ in order.
    unmapped.reference_index = if aligned_reference_name == "*" {
        NO_ALIGNMENT_REFERENCE_INDEX
    } else {
        out_sequences
            .iter()
            .position(|s| s.name == aligned_reference_name)
            .map(|i| i as i32)
            .unwrap_or(NO_ALIGNMENT_REFERENCE_INDEX)
    };

    unmapped.alignment_start = aligned.alignment_start;
    set_flag(
        &mut unmapped.flags,
        READ_REVERSE_STRAND,
        aligned.flags & READ_REVERSE_STRAND != 0,
    );
    set_flag(
        &mut unmapped.flags,
        NOT_PRIMARY_ALIGNMENT,
        aligned.flags & NOT_PRIMARY_ALIGNMENT != 0,
    );
    set_flag(
        &mut unmapped.flags,
        SUPPLEMENTARY_ALIGNMENT,
        aligned.flags & SUPPLEMENTARY_ALIGNMENT != 0,
    );

    if !aligned_unmapped {
        // Only aligned reads carry a cigar and mapping quality.
        unmapped.cigar = aligned.cigar.clone();
        unmapped.mapping_quality = aligned.mapping_quality;
    }

    if unmapped.flags & READ_PAIRED != 0 {
        set_flag(
            &mut unmapped.flags,
            PROPER_PAIR,
            aligned.flags & PROPER_PAIR != 0,
        );
    }

    if unmapped.flags & READ_REVERSE_STRAND != 0 {
        reverse_complement(&mut unmapped.read_bases);
        reverse_qualities(&mut unmapped.base_qualities);
        reverse_complement_tags(unmapped);
    }

    Ok(())
}

/// `AbstractAlignmentMerger.maybeSetPgTag`: link the read to the merge's program record by setting its
/// `PG` tag (default `ADD_PG_TAG_TO_READS=true`). `program_id` is `None` when no program record is in
/// play.
pub fn maybe_set_pg_tag(record: &mut BamRecord, program_id: Option<&str>) {
    if let Some(id) = program_id {
        record
            .tags
            .insert(Tag::new(b"PG"), TagValue::Str(id.to_string()));
    }
}

/// One merged fragment for the default coordinate-sorted output: transfer the alignment, link the
/// program record, and (for a mapped read) recompute `NM`/`MD`/`UQ` against the reference, exactly as
/// `AbstractAlignmentMerger` does across `setValuesFromAlignment`, `maybeSetPgTag`, and the
/// coordinate-sorted `fixNmMdAndUq` pass. `reference_bases` is the bases of the contig the read maps
/// to; the caller resolves it from the merged record's reference name.
///
/// Scope: a single unpaired primary hit. Paired mate-fixing, off-end clipping and multi-hit selection
/// are later slices.
pub fn merge_aligned_fragment(
    unmapped: &BamRecord,
    aligned: &BamRecord,
    aligned_reference_name: &str,
    out_sequences: &[SequenceRecord],
    reference_bases: &[u8],
    program_id: Option<&str>,
) -> Result<BamRecord, MergeAlignmentError> {
    let mut merged = unmapped.clone();
    transfer_alignment_info_to_fragment(
        &mut merged,
        aligned,
        aligned_reference_name,
        out_sequences,
    )?;
    maybe_set_pg_tag(&mut merged, program_id);
    if merged.flags & READ_UNMAPPED == 0 {
        fix_nm_md_and_uq(&mut merged, reference_bases);
    }
    Ok(merged)
}

/// `SamPairUtil.PairOrientation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PairOrientation {
    Fr,
    Rf,
    Tandem,
}

/// `SamPairUtil.getPairOrientation`: the orientation of a paired, both-ends-mapped read from its own
/// and its mate's strand and 5' positions.
fn get_pair_orientation(rec: &BamRecord) -> PairOrientation {
    let read_reverse = rec.flags & READ_REVERSE_STRAND != 0;
    let mate_reverse = rec.flags & MATE_NEGATIVE_STRAND != 0;
    if read_reverse == mate_reverse {
        return PairOrientation::Tandem;
    }
    let positive_5prime = if read_reverse {
        rec.mate_alignment_start
    } else {
        rec.alignment_start
    };
    let negative_5prime = if read_reverse {
        rec.alignment_end()
    } else {
        rec.alignment_start + rec.inferred_insert_size
    };
    if positive_5prime < negative_5prime {
        PairOrientation::Fr
    } else {
        PairOrientation::Rf
    }
}

/// `SamPairUtil.isProperPair`: both ends mapped to the same reference with an expected orientation.
fn is_proper_pair(first: &BamRecord, second: &BamRecord, expected: &[PairOrientation]) -> bool {
    if first.flags & READ_UNMAPPED != 0 || second.flags & READ_UNMAPPED != 0 {
        return false;
    }
    if first.reference_index < 0 || first.reference_index != second.reference_index {
        return false;
    }
    expected.contains(&get_pair_orientation(first))
}

/// `SamPairUtil.setProperPairFlags`: set the proper-pair flag on both ends to whether they form a
/// proper pair (both mapped, expected orientation), or clear it otherwise.
fn set_proper_pair_flags(rec1: &mut BamRecord, rec2: &mut BamRecord, expected: &[PairOrientation]) {
    let proper = rec1.flags & READ_UNMAPPED == 0
        && rec2.flags & READ_UNMAPPED == 0
        && is_proper_pair(rec1, rec2, expected);
    set_flag(&mut rec1.flags, PROPER_PAIR, proper);
    set_flag(&mut rec2.flags, PROPER_PAIR, proper);
}

/// The merged pair for a paired template, for the default coordinate-sorted output:
/// `transferAlignmentInfoToPairedRead` transfers each end, `SamPairUtil.setMateInfo` cross-sets the
/// mate coordinates / strand / `MQ` / `MC` and insert size, `setProperPairFlags` recomputes the
/// proper-pair flag against the default `FR` orientation, and each end then gets its `PG` and
/// reference-based `NM`/`MD`/`UQ`. `*_reference_bases` is the bases of the contig each end maps to.
///
/// Scope: both ends primary and mapped, non-overlapping (so `CLIP_OVERLAPPING_READS`, though on by
/// default, does nothing). Overlap clipping and off-end clipping are later slices.
#[allow(clippy::too_many_arguments)]
pub fn merge_aligned_pair(
    first_unmapped: &BamRecord,
    second_unmapped: &BamRecord,
    first_aligned: &BamRecord,
    second_aligned: &BamRecord,
    first_reference_name: &str,
    second_reference_name: &str,
    out_sequences: &[SequenceRecord],
    first_reference_bases: &[u8],
    second_reference_bases: &[u8],
    program_id: Option<&str>,
) -> Result<(BamRecord, BamRecord), MergeAlignmentError> {
    let mut first = first_unmapped.clone();
    let mut second = second_unmapped.clone();
    transfer_alignment_info_to_fragment(
        &mut first,
        first_aligned,
        first_reference_name,
        out_sequences,
    )?;
    transfer_alignment_info_to_fragment(
        &mut second,
        second_aligned,
        second_reference_name,
        out_sequences,
    )?;

    // htsjdk calls setMateInfo(second, first, addMateCigar) then setProperPairFlags(second, first).
    set_mate_info(&mut second, &mut first, true);
    set_proper_pair_flags(&mut second, &mut first, &[PairOrientation::Fr]);

    for (rec, bases) in [
        (&mut first, first_reference_bases),
        (&mut second, second_reference_bases),
    ] {
        maybe_set_pg_tag(rec, program_id);
        if rec.flags & READ_UNMAPPED == 0 {
            fix_nm_md_and_uq(rec, bases);
        }
    }
    Ok((first, second))
}

/// The merged, coordinate-sorted records for the default unpaired single-hit path: match each aligned
/// record to its unmapped read by name, [`merge_aligned_fragment`] the pair, then sort by coordinate
/// (`MergingSamRecordIterator` walks the two queryname-sorted inputs and the output writer re-sorts to
/// `SORT_ORDER=coordinate`).
///
/// `aligned_sequences` is the aligned file's dictionary (to resolve each aligned record's reference
/// name), `out_sequences` the output dictionary, and `reference_bases` maps each output contig name to
/// its bases (for the `NM`/`MD`/`UQ` recomputation). `program_id` is linked into each read's `PG`.
///
/// Scope: one primary hit per read, every unmapped read having exactly one aligned record of the same
/// name. The merged-header construction, paired mate-fixing, clipping and multi-hit selection are
/// later slices.
pub fn merge_bam_alignment_records(
    unmapped: &[BamRecord],
    aligned: &[BamRecord],
    aligned_sequences: &[SequenceRecord],
    out_sequences: &[SequenceRecord],
    reference_bases: &std::collections::HashMap<String, Vec<u8>>,
    program_id: Option<&str>,
) -> Result<Vec<BamRecord>, MergeAlignmentError> {
    let unmapped_by_name: std::collections::HashMap<&str, &BamRecord> =
        unmapped.iter().map(|r| (r.read_name.as_str(), r)).collect();

    let mut merged: Vec<BamRecord> = Vec::with_capacity(aligned.len());
    for a in aligned {
        let reference_name = if a.reference_index < 0 {
            "*".to_string()
        } else {
            aligned_sequences[a.reference_index as usize].name.clone()
        };
        let contig_bases: &[u8] = reference_bases
            .get(&reference_name)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let u = unmapped_by_name
            .get(a.read_name.as_str())
            .expect("every aligned read has a same-named unmapped read");
        merged.push(merge_aligned_fragment(
            u,
            a,
            &reference_name,
            out_sequences,
            contig_bases,
            program_id,
        )?);
    }

    merged.sort_by(htsjdk_bam::coordinate::compare);
    Ok(merged)
}

/// The paired analogue of [`merge_bam_alignment_records`]: group each template's first- and
/// second-of-pair ends in both files by name, [`merge_aligned_pair`] them, and coordinate-sort.
///
/// Scope: every template has both ends primary and mapped, non-overlapping. Singletons, overlap and
/// off-end clipping, and multi-hit selection are later slices.
pub fn merge_bam_alignment_paired(
    unmapped: &[BamRecord],
    aligned: &[BamRecord],
    aligned_sequences: &[SequenceRecord],
    out_sequences: &[SequenceRecord],
    reference_bases: &std::collections::HashMap<String, Vec<u8>>,
    program_id: Option<&str>,
) -> Result<Vec<BamRecord>, MergeAlignmentError> {
    use std::collections::HashMap;

    // Index each file's ends by (name, is-first-of-pair).
    let key = |r: &BamRecord| (r.read_name.clone(), r.flags & FIRST_OF_PAIR != 0);
    let unmapped_by: HashMap<(String, bool), &BamRecord> =
        unmapped.iter().map(|r| (key(r), r)).collect();
    let aligned_by: HashMap<(String, bool), &BamRecord> =
        aligned.iter().map(|r| (key(r), r)).collect();

    let reference_name = |a: &BamRecord| {
        if a.reference_index < 0 {
            "*".to_string()
        } else {
            aligned_sequences[a.reference_index as usize].name.clone()
        }
    };
    let bases_for = |name: &str| {
        reference_bases
            .get(name)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    };

    // Each template appears once as a first-of-pair aligned record.
    let mut merged: Vec<BamRecord> = Vec::with_capacity(aligned.len());
    for a_first in aligned.iter().filter(|r| r.flags & FIRST_OF_PAIR != 0) {
        let name = a_first.read_name.as_str();
        let a_second = aligned_by[&(name.to_string(), false)];
        let u_first = unmapped_by[&(name.to_string(), true)];
        let u_second = unmapped_by[&(name.to_string(), false)];
        let first_ref = reference_name(a_first);
        let second_ref = reference_name(a_second);
        let (first, second) = merge_aligned_pair(
            u_first,
            u_second,
            a_first,
            a_second,
            &first_ref,
            &second_ref,
            out_sequences,
            bases_for(&first_ref),
            bases_for(&second_ref),
            program_id,
        )?;
        merged.push(first);
        merged.push(second);
    }

    merged.sort_by(htsjdk_bam::coordinate::compare);
    Ok(merged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::SequenceRecord;
    use htsjdk_bam::sam_file::{read_sam, write_sam};

    fn one(sam: &str) -> BamRecord {
        read_sam(sam).unwrap().1.into_iter().next().unwrap()
    }

    fn chr1() -> Vec<SequenceRecord> {
        vec![SequenceRecord::new("chr1", 40)]
    }

    // The transferred record rendered as its SAM data line, for readable assertions.
    fn line(rec: &BamRecord) -> String {
        let header = {
            let mut h = htsjdk_bam::header::SamHeader::new();
            h.sequences = chr1();
            h
        };
        write_sam(&header, std::slice::from_ref(rec))
            .unwrap()
            .lines()
            .find(|l| !l.starts_with('@'))
            .unwrap()
            .to_string()
    }

    const UNMAPPED: &str = "@HD\tVN:1.6\n@RG\tID:rg1\tSM:s\n\
        r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tab:Z:keepme\n";
    const ALIGNED: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:s\n\
        r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\tX0:i:1\n";

    #[test]
    fn a_forward_hit_transfers_alignment_and_keeps_the_original_read() {
        let mut unmapped = one(UNMAPPED);
        let aligned = one(ALIGNED);
        transfer_alignment_info_to_fragment(&mut unmapped, &aligned, "chr1", &chr1()).unwrap();
        // Flag 0, mapped to chr1:1, 8M, MAPQ 60; bases/quals kept; NM copied, X0 (reserved) dropped;
        // tags in htsjdk's packed-code order (RG, NM, ab).
        assert_eq!(
            line(&unmapped),
            "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\tab:Z:keepme"
        );
    }

    #[test]
    fn a_reserved_aligner_tag_is_not_copied() {
        let mut unmapped = one(UNMAPPED);
        let aligned = one(ALIGNED);
        transfer_alignment_info_to_fragment(&mut unmapped, &aligned, "chr1", &chr1()).unwrap();
        assert!(unmapped.tags.get(Tag::new(b"X0")).is_none());
        assert!(unmapped.tags.get(Tag::new(b"NM")).is_some());
    }

    #[test]
    fn a_negative_strand_hit_reverse_complements_bases_and_reverses_quals() {
        // A non-palindromic read so the reverse-complement is visible.
        let mut unmapped = one("@HD\tVN:1.6\n\
            r\t4\t*\t0\t0\t*\t*\t0\t0\tAAACGGGG\t01234567\tOQ:Z:HGFEDCBA\n");
        let aligned = one("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n\
            r\t16\tchr1\t1\t60\t8M\t*\t0\t0\tCCCCGTTT\t76543210\n");
        transfer_alignment_info_to_fragment(&mut unmapped, &aligned, "chr1", &chr1()).unwrap();
        // AAACGGGG -> reverse-complement CCCCGTTT; quals 01234567 -> reversed 76543210; OQ reversed.
        assert_eq!(unmapped.read_bases, b"CCCCGTTT");
        assert_eq!(
            unmapped.base_qualities,
            one("@HD\tVN:1.6\n\
            x\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCGTTT\t76543210\n")
            .base_qualities
        );
        assert_eq!(
            unmapped.tags.get(Tag::new(b"OQ")),
            Some(&TagValue::Str("ABCDEFGH".to_string()))
        );
    }

    #[test]
    fn a_mapped_unmapped_bam_read_is_rejected() {
        let mut mapped = one("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n\
            r1\t0\tchr1\t5\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n");
        let aligned = one(ALIGNED);
        assert!(matches!(
            transfer_alignment_info_to_fragment(&mut mapped, &aligned, "chr1", &chr1()),
            Err(MergeAlignmentError::UnmappedBamContainsMappedRead(_))
        ));
    }

    #[test]
    fn is_reserved_tag_matches_the_rule() {
        assert!(is_reserved_tag(b"X0")); // aligner-private
        assert!(is_reserved_tag(b"ab")); // lower-case user tag
        assert!(!is_reserved_tag(b"NM")); // standard tag carries over
        assert!(!is_reserved_tag(b"AS"));
    }
}
