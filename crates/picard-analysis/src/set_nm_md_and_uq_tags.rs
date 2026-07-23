//! `SetNmMdAndUqTags`.
//!
//! Ports `picard.sam.SetNmMdAndUqTags.doWork` + `AbstractAlignmentMerger.fixNmMdAndUq` at tag 3.4.0,
//! for the **default (non-bisulfite) path**: for every mapped read, recompute the `MD`, `NM`, and `UQ`
//! tags against the reference, and leave unmapped reads untouched. The first reference tool that walks
//! a read's alignment against the reference bases.
//!
//! - `MD`/`NM` come from [`calculate_md_and_nm`](htsjdk_bam::md_nm::calculate_md_and_nm).
//! - `UQ` is `SequenceUtil.sumQualitiesOfMismatches`: the sum of a read's base qualities at the
//!   positions where it mismatches the reference (only when the read carries qualities).
//!
//! The input must be **coordinate-sorted** (the tool throws otherwise), the header is written
//! unchanged, and no `@PG` and no timestamp are added, so the whole output is comparable raw. The
//! records keep their input order (`makeWriter(header, presorted=true, ...)`). Bisulfite handling
//! (`IS_BISULFITE_SEQUENCE`) and `SET_ONLY_UQ` are separate surfaces.

use std::collections::HashMap;

use htsjdk_bam::cigar::Op;
use htsjdk_bam::fasta::{read_fasta, FastaError};
use htsjdk_bam::md_nm::calculate_md_and_nm;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::sequence::bases_equal;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

const READ_UNMAPPED: u16 = 0x4;

/// Why `SetNmMdAndUqTags` could not run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SetTagsError {
    Parse(ParseError),
    Fasta(String),
    /// The input was not coordinate-sorted.
    NotCoordinateSorted(String),
    /// A mapped read references a contig with no reference bases.
    MissingReferenceContig(String),
}

impl From<ParseError> for SetTagsError {
    fn from(e: ParseError) -> Self {
        SetTagsError::Parse(e)
    }
}

impl From<FastaError> for SetTagsError {
    fn from(e: FastaError) -> Self {
        SetTagsError::Fasta(format!("{e:?}"))
    }
}

/// `SequenceUtil.sumQualitiesOfMismatches` (non-bisulfite, referenceOffset 0): the sum of the read's
/// base qualities at positions where it mismatches the reference, walking the alignment (`M`/`=`/`X`)
/// blocks.
fn sum_qualities_of_mismatches(rec: &BamRecord, ref_bases: &[u8]) -> i32 {
    let mut qualities: i32 = 0;
    let mut block_ref = (rec.alignment_start - 1) as isize;
    let mut block_read: usize = 0;
    for ce in &rec.cigar.elements {
        let len = ce.length as usize;
        match ce.op {
            Op::M | Op::Eq | Op::X => {
                for i in 0..len {
                    let ref_idx = block_ref + i as isize;
                    if ref_idx < 0 || ref_idx as usize >= ref_bases.len() {
                        continue;
                    }
                    let read_base = rec.read_bases[block_read + i];
                    if !bases_equal(read_base, ref_bases[ref_idx as usize]) {
                        qualities += rec.base_qualities[block_read + i] as i32;
                    }
                }
                block_ref += len as isize;
                block_read += len;
            }
            Op::D | Op::N => block_ref += len as isize,
            Op::I | Op::S => block_read += len,
            Op::H | Op::P => {}
        }
    }
    qualities
}

/// `SetNmMdAndUqTags.doWork` up to the write: the coordinate-sorted header and the records with their
/// recomputed `MD`/`NM`/`UQ` tags.
fn set_tags(
    input_sam: &str,
    fasta: &[u8],
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), SetTagsError> {
    let (header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    if header.attributes.get("SO") != Some("coordinate") {
        return Err(SetTagsError::NotCoordinateSorted(
            header
                .attributes
                .get("SO")
                .unwrap_or("unsorted")
                .to_string(),
        ));
    }

    // Map each contig name to its bases, and each reference index to that name via the @SQ order.
    let contigs = read_fasta(fasta)?;
    let by_name: HashMap<&str, &[u8]> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.bases.as_slice()))
        .collect();

    for rec in &mut records {
        if rec.flags & READ_UNMAPPED != 0 {
            continue; // fixRecord leaves unmapped reads alone
        }
        let name = &header.sequences[rec.reference_index as usize].name;
        let ref_bases = *by_name
            .get(name.as_str())
            .ok_or_else(|| SetTagsError::MissingReferenceContig(name.clone()))?;
        fix_nm_md_and_uq(rec, ref_bases);
    }

    Ok((header, records))
}

/// `AbstractAlignmentMerger.fixNmMdAndUq` (non-bisulfite): recompute `MD`/`NM` from the reference and,
/// when the read carries qualities, `UQ`. Shared with the alignment merger, which recomputes the same
/// tags in its coordinate-sorted final pass.
pub(crate) fn fix_nm_md_and_uq(rec: &mut BamRecord, ref_bases: &[u8]) {
    let (md, nm) = calculate_md_and_nm(rec.alignment_start, &rec.cigar, &rec.read_bases, ref_bases);
    rec.tags.insert(Tag::new(b"MD"), TagValue::Str(md));
    rec.tags.insert(Tag::new(b"NM"), TagValue::Int(nm as i64));

    // fixUq: UQ only when the read carries qualities.
    if !rec.base_qualities.is_empty() {
        let uq = sum_qualities_of_mismatches(rec, ref_bases);
        rec.tags.insert(Tag::new(b"UQ"), TagValue::Int(uq as i64));
    }
}

/// `SetNmMdAndUqTags.doWork` for SAM input and output, default (non-bisulfite) options. `fasta` is the
/// `REFERENCE_SEQUENCE` bytes.
pub fn set_nm_md_and_uq_tags(input_sam: &str, fasta: &[u8]) -> Result<String, SetTagsError> {
    let (header, records) = set_tags(input_sam, fasta)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// `SetNmMdAndUqTags.doWork` for SAM input and **BAM** output, its default output format. Same
/// recomputation, written through htsjdk-rs's byte-identical `BamWriter`; the tool adds no `@PG`.
/// Byte-identity to Picard's `USE_JDK_DEFLATER=true` output follows transitively: the tagged records
/// are the ones `set_nm_md_and_uq_tags` already reproduces (its oracle), and the `BamWriter` is proven
/// byte-identical over arbitrary records (the SamFormatConverter oracle and htsjdk-rs's
/// `every_file_is_byte_identical_to_htsjdks`).
pub fn set_nm_md_and_uq_tags_to_bam(
    input_sam: &str,
    fasta: &[u8],
) -> Result<Vec<u8>, SetTagsError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, records) = set_tags(input_sam, fasta)?;
    let mut writer = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        writer
            .write(rec)
            .expect("records that parsed re-encode as BAM");
    }
    Ok(writer
        .finish()
        .expect("in-memory BAM writer never fails to finish"))
}

/// `SetNmAndUqTags.doWork`. `SetNmAndUqTags` is a `@Deprecated` subclass of `SetNmMdAndUqTags` that
/// overrides nothing (`SetNmAndUqTags.java` at Picard 3.4.0 is only annotations and usage strings), so
/// it inherits `doWork` unchanged and, despite its name, still writes `MD` as well as `NM`/`UQ`. Its
/// output is therefore byte-for-byte the same as [`set_nm_md_and_uq_tags`]; this alias exists so the
/// deprecated tool name is covered by the port.
pub fn set_nm_and_uq_tags(input_sam: &str, fasta: &[u8]) -> Result<String, SetTagsError> {
    set_nm_md_and_uq_tags(input_sam, fasta)
}

/// `SetNmAndUqTags.doWork` for **BAM** output, identical to [`set_nm_md_and_uq_tags_to_bam`] for the
/// same reason.
pub fn set_nm_and_uq_tags_to_bam(input_sam: &str, fasta: &[u8]) -> Result<Vec<u8>, SetTagsError> {
    set_nm_md_and_uq_tags_to_bam(input_sam, fasta)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA: &[u8] = b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:40\n\
        r\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACCTACGT\t##II##II\n";

    /// The BAM output decodes back to exactly the SAM output; the writer's byte-identity to htsjdk is
    /// proven elsewhere.
    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        let sam = set_nm_md_and_uq_tags(INPUT, FASTA).unwrap();
        let bam = set_nm_md_and_uq_tags_to_bam(INPUT, FASTA).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bam).expect("bam decompresses");
        let reader = htsjdk_bam::reader::BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(
            htsjdk_bam::sam_file::write_sam(&header, &records).unwrap(),
            sam
        );
    }

    /// A sanity check that the recomputation actually fires: the one mismatch yields NM:i:1 and a UQ
    /// summing the two low qualities at that position (here a single mismatch at offset 2).
    #[test]
    fn a_mismatch_gets_md_nm_and_uq() {
        let sam = set_nm_md_and_uq_tags(INPUT, FASTA).unwrap();
        let row = sam.lines().find(|l| l.starts_with('r')).unwrap();
        assert!(row.contains("MD:Z:2G5"), "got {row}");
        assert!(row.contains("NM:i:1"), "got {row}");
        assert!(row.contains("UQ:i:"), "got {row}");
    }

    /// The deprecated `SetNmAndUqTags` inherits `doWork`, so its output is identical to
    /// `SetNmMdAndUqTags` in both formats.
    #[test]
    fn the_deprecated_alias_is_identical_to_the_base_tool() {
        assert_eq!(
            set_nm_and_uq_tags(INPUT, FASTA).unwrap(),
            set_nm_md_and_uq_tags(INPUT, FASTA).unwrap()
        );
        assert_eq!(
            set_nm_and_uq_tags_to_bam(INPUT, FASTA).unwrap(),
            set_nm_md_and_uq_tags_to_bam(INPUT, FASTA).unwrap()
        );
    }
}
