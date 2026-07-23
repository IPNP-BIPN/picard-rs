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
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{reverse, reverse_complement, reverse_qualities};
use htsjdk_bam::tag::{Tag, TagValue};

const READ_PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const READ_UNMAPPED: u16 = 0x4;
const READ_REVERSE_STRAND: u16 = 0x10;
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
