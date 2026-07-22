//! `RevertSam`.
//!
//! Ports `picard.sam.RevertSam.revertSamRecord` + `createOutHeader` + `doWork` at tag 3.4.0, for the
//! **default option path**: undo the alignment of a SAM/BAM so it can be re-aligned. With the
//! defaults (`REMOVE_DUPLICATE_INFORMATION`, `REMOVE_ALIGNMENT_INFORMATION`,
//! `RESTORE_ORIGINAL_QUALITIES` all true, `SANITIZE`/`OUTPUT_BY_READGROUP`/`RESTORE_HARDCLIPS` false,
//! `SORT_ORDER=queryname`) each record is stripped back to an unmapped read and the file is
//! queryname-sorted.
//!
//! The output header, from `createOutHeader` with `removeAlignmentInformation=true`, is a **fresh**
//! `SAMFileHeader`: `@HD VN:1.6 SO:queryname` plus the input's `@RG` lines (added verbatim,
//! `doWork` l.289), and **no `@SQ`, no `@PG`, no `@CO`**. RevertSam adds no `@PG` and no timestamp, so
//! the whole file is comparable raw. Every record ends unmapped, so a missing sequence dictionary is
//! fine: every RNAME/RNEXT resolves to `*`.
//!
//! The per-record revert (`revertSamRecord`) is independent, so it runs on all cores and stays
//! byte-identical (decision 0006); the queryname sort must be a **stable** in-memory sort for
//! byte-identity (decision 0021).
//!
//! Only the default path is claimed. `SANITIZE` (a two-pass discard of unpaired/duplicate-name
//! reads), `OUTPUT_BY_READGROUP` (a file per read group), `RESTORE_HARDCLIPS`, `SAMPLE_ALIAS`,
//! `LIBRARY_NAME`, a non-default `SORT_ORDER`, and a customized `ATTRIBUTE_TO_CLEAR` are separate
//! surfaces.

use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::fastq::fastq_to_phred;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::query_name;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::sequence::{reverse_complement, reverse_qualities};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use rayon::prelude::*;

const READ_PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const MATE_NEGATIVE_STRAND: u16 = 0x20;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const DUPLICATE_READ: u16 = 0x400;

const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;
const NO_ALIGNMENT_START: i32 = 0;

/// `RevertSam.ATTRIBUTE_TO_CLEAR` default: tags calculated from the alignment.
const ATTRIBUTE_TO_CLEAR: [&[u8; 2]; 8] = [b"NM", b"UQ", b"PG", b"MD", b"MQ", b"SA", b"MC", b"AS"];
/// `SAMRecord.TAGS_TO_REVERSE_COMPLEMENT`.
const TAGS_TO_REVERSE_COMPLEMENT: [&[u8; 2]; 2] = [b"E2", b"SQ"];
/// `SAMRecord.TAGS_TO_REVERSE`.
const TAGS_TO_REVERSE: [&[u8; 2]; 2] = [b"OQ", b"U2"];

/// `SAMRecord.reverseComplement(TAGS_TO_REVERSE_COMPLEMENT, TAGS_TO_REVERSE, inplace=true)` for the
/// fields that survive the revert: the bases (reverse-complemented), the qualities (reversed), and
/// the string tags in the two default lists. The alignment is dropped right after, so the CIGAR the
/// full htsjdk method would also reverse is not reproduced.
fn reverse_complement_record(rec: &mut BamRecord) {
    reverse_complement(&mut rec.read_bases);
    reverse_qualities(&mut rec.base_qualities);

    for name in TAGS_TO_REVERSE_COMPLEMENT {
        if let Some(TagValue::Str(s)) = rec.tags.get(Tag::new(name)) {
            let mut bytes = s.clone().into_bytes();
            reverse_complement(&mut bytes);
            let value =
                String::from_utf8(bytes).expect("a reverse-complemented base string is ASCII");
            rec.tags.insert(Tag::new(name), TagValue::Str(value));
        }
    }
    for name in TAGS_TO_REVERSE {
        if let Some(TagValue::Str(s)) = rec.tags.get(Tag::new(name)) {
            let value: String = s.chars().rev().collect();
            rec.tags.insert(Tag::new(name), TagValue::Str(value));
        }
    }
}

/// `RevertSam.revertSamRecord` with the default options.
fn revert_record(rec: &mut BamRecord) {
    // RESTORE_ORIGINAL_QUALITIES: move OQ back into QUAL and drop OQ.
    if let Some(TagValue::Str(oq)) = rec.tags.get(Tag::new(b"OQ")) {
        rec.base_qualities = fastq_to_phred(oq);
        rec.tags.remove(Tag::new(b"OQ"));
    }

    // REMOVE_DUPLICATE_INFORMATION.
    rec.flags &= !DUPLICATE_READ;

    // REMOVE_ALIGNMENT_INFORMATION.
    if rec.flags & READ_NEGATIVE_STRAND != 0 {
        reverse_complement_record(rec);
        rec.flags &= !READ_NEGATIVE_STRAND;
    }

    // Remove all alignment-based information about the read itself.
    rec.reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
    rec.alignment_start = NO_ALIGNMENT_START;
    rec.cigar = Cigar::default(); // NO_ALIGNMENT_CIGAR, "*"
    rec.mapping_quality = 0; // NO_MAPPING_QUALITY
    rec.inferred_insert_size = 0;
    rec.flags &= !SECONDARY_ALIGNMENT; // setNotPrimaryAlignmentFlag(false)
    rec.flags &= !PROPER_PAIR; // setProperPairFlag(false)
    rec.flags |= READ_UNMAPPED;

    // Then remove any mate flags and info related to alignment.
    rec.mate_alignment_start = NO_ALIGNMENT_START;
    rec.flags &= !MATE_NEGATIVE_STRAND;
    rec.mate_reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
    // setMateUnmappedFlag(getReadPairedFlag()): a paired read's mate becomes unmapped too.
    if rec.flags & READ_PAIRED != 0 {
        rec.flags |= MATE_UNMAPPED;
    } else {
        rec.flags &= !MATE_UNMAPPED;
    }

    // And then remove any tags that are calculated from the alignment.
    for name in ATTRIBUTE_TO_CLEAR {
        rec.tags.remove(Tag::new(name));
    }
}

/// `RevertSam.doWork` for SAM input and output, default options.
pub fn revert_sam(input_sam: &str) -> Result<String, ParseError> {
    // RevertSam opens the input EAGERLY_DECODE at VALIDATION_STRINGENCY.SILENT; stringency does not
    // reach the bytes.
    let (input_header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    // createOutHeader(removeAlignmentInformation=true): a fresh header, sort order queryname, plus
    // the input read groups verbatim; no @SQ, no @PG, no @CO.
    let mut out_header = SamHeader::new();
    out_header.set_sort_order("queryname");
    out_header.read_groups = input_header.read_groups.clone();

    // The per-record revert is independent and order-preserving (decision 0006).
    records.par_iter_mut().for_each(revert_record);

    // The presorted=false writer sorts by SORT_ORDER (queryname); a stable sort keeps records that
    // compare equal in input order (decision 0021).
    records.sort_by(query_name::compare);

    Ok(write_sam(&out_header, &records).expect("unmapped records always encode as SAM text"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::sam_file::read_sam;

    // Coordinate-sorted input (zeb@100, amy@200, mid@300); queryname order is amy, mid, zeb, so the
    // sort is load-bearing.
    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @RG\tID:rg1\tSM:s\n\
        zeb\t1024\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tOQ:Z:5555\tNM:i:1\tMD:Z:4\tAS:i:40\tRG:Z:rg1\n\
        amy\t16\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tABCD\tRG:Z:rg1\n\
        mid\t99\tchr1\t300\t60\t4M\t=\t350\t54\tACGT\tIIII\tMC:Z:4M\tRG:Z:rg1\n";

    fn row<'a>(sam: &'a str, name: &str) -> Vec<&'a str> {
        sam.lines()
            .find(|l| l.starts_with(name))
            .unwrap()
            .split('\t')
            .collect()
    }

    #[test]
    fn the_header_is_bare_with_read_groups_and_queryname_order() {
        let out = revert_sam(INPUT).unwrap();
        assert!(out.starts_with("@HD\tVN:1.6\tSO:queryname\n"), "got {out}");
        assert!(out.contains("@RG\tID:rg1\tSM:s"), "read groups kept: {out}");
        assert!(!out.contains("@SQ"), "no sequence dictionary: {out}");
        assert!(!out.contains("@PG"), "no program record: {out}");
    }

    #[test]
    fn records_come_out_in_queryname_order() {
        let out = revert_sam(INPUT).unwrap();
        let names: Vec<&str> = out
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(names, ["amy", "mid", "zeb"]);
    }

    #[test]
    fn a_duplicate_read_has_its_alignment_and_calculated_tags_removed_and_oq_restored() {
        let out = revert_sam(INPUT).unwrap();
        let f = row(&out, "zeb");
        assert_eq!(f[1], "4"); // dup(0x400) cleared, unmapped(0x4) set
        assert_eq!(f[2], "*"); // RNAME
        assert_eq!(f[3], "0"); // POS
        assert_eq!(f[4], "0"); // MAPQ
        assert_eq!(f[5], "*"); // CIGAR
        assert_eq!(f[10], "5555"); // QUAL restored from OQ:Z:5555
        let tags = &f[11..];
        assert!(tags.contains(&"RG:Z:rg1"), "RG kept: {tags:?}");
        assert!(
            !tags.iter().any(|t| t.starts_with("OQ")
                || t.starts_with("NM")
                || t.starts_with("MD")
                || t.starts_with("AS")),
            "calculated tags and OQ removed: {tags:?}"
        );
    }

    #[test]
    fn a_negative_strand_read_is_reverse_complemented() {
        let out = revert_sam(INPUT).unwrap();
        let f = row(&out, "amy");
        assert_eq!(f[1], "4"); // negative-strand(0x10) cleared, unmapped set
                               // revcomp(ACGT) = ACGT (palindrome); quals ABCD reversed to DCBA.
        assert_eq!(f[9], "ACGT");
        assert_eq!(f[10], "DCBA");
    }

    #[test]
    fn a_proper_pair_loses_its_mate_alignment_and_mc() {
        let out = revert_sam(INPUT).unwrap();
        let f = row(&out, "mid");
        // 99 = paired|proper|mate-neg|first; -> paired|unmapped|mate-unmapped|first = 77.
        assert_eq!(f[1], "77");
        assert_eq!(f[6], "*"); // RNEXT
        assert_eq!(f[7], "0"); // PNEXT
        assert_eq!(f[8], "0"); // TLEN
        assert!(!out
            .lines()
            .any(|l| l.starts_with("mid") && l.contains("MC:Z")));
    }

    /// The parallel revert must produce the same bytes as a serial one (decision 0006).
    #[test]
    fn parallel_and_serial_reverts_agree() {
        let parallel = revert_sam(INPUT).unwrap();

        let (input_header, mut records) = read_sam(INPUT).unwrap();
        let mut out_header = SamHeader::new();
        out_header.set_sort_order("queryname");
        out_header.read_groups = input_header.read_groups.clone();
        for rec in &mut records {
            revert_record(rec);
        }
        records.sort_by(query_name::compare);
        let serial = write_sam(&out_header, &records).unwrap();

        assert_eq!(parallel, serial);
    }
}
