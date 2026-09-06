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
//! records keep their input order (`makeWriter(header, presorted=true, ...)`).
//!
//! [`Options`] carries the tool's other two arguments, and neither is a variation on the default
//! path:
//!
//! * `SET_ONLY_UQ` writes `UQ` and leaves `MD` and `NM` as the input had them, so a read whose
//!   input tags are wrong stays wrong in two of the three;
//! * `IS_BISULFITE_SEQUENCE` changes what counts as a mismatch, but not everywhere. `MD` is built
//!   by `calculateMdAndNmTags`, whose comparison is `bases[read] == bases[ref]` with no bisulfite
//!   branch at all, so `MD` is IDENTICAL under both settings. What changes is `NM`, recomputed
//!   through `calculateSamNmTag`, and `UQ`, both of which forgive C to T on the positive strand
//!   and G to A on the negative one.

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
const READ_REVERSE_STRAND: u16 = 0x10;

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

/// `SequenceUtil.isBisulfiteConverted`: C to T on the positive strand, G to A on the negative one.
///
/// The direction is the reference's, not the read's: a `C` in the reference read as `T` is the
/// conversion, and the reverse is a mismatch like any other.
fn is_bisulfite_converted(read: u8, reference: u8, negative_strand: bool) -> bool {
    if negative_strand {
        bases_equal(reference, b'G') && bases_equal(read, b'A')
    } else {
        bases_equal(reference, b'C') && bases_equal(read, b'T')
    }
}

/// `SequenceUtil.bisulfiteBasesEqual`.
fn bisulfite_bases_equal(negative_strand: bool, read: u8, reference: u8) -> bool {
    bases_equal(read, reference) || is_bisulfite_converted(read, reference, negative_strand)
}

/// `SequenceUtil.countMismatches(read, ref, 0, bisulfiteSequence, matchAmbiguousRef = false)`.
///
/// Walks the same alignment blocks as [`sum_qualities_of_mismatches`] and counts, rather than sums.
fn count_mismatches(rec: &BamRecord, ref_bases: &[u8], bisulfite: bool) -> i32 {
    let negative_strand = rec.flags & READ_REVERSE_STRAND != 0;
    let mut mismatches: i32 = 0;
    walk_alignment(rec, |read_offset, ref_index| {
        if ref_index >= ref_bases.len() {
            return;
        }
        let read_base = rec.read_bases[read_offset];
        let ref_base = ref_bases[ref_index];
        let matched = if bisulfite {
            bisulfite_bases_equal(negative_strand, read_base, ref_base)
        } else {
            bases_equal(read_base, ref_base)
        };
        if !matched {
            mismatches += 1;
        }
    });
    mismatches
}

/// `SequenceUtil.calculateSamNmTag(read, ref, 0, bisulfiteSequence)`: the mismatches, plus every
/// inserted and deleted base.
fn calculate_sam_nm_tag(rec: &BamRecord, ref_bases: &[u8], bisulfite: bool) -> i32 {
    let mut nm = count_mismatches(rec, ref_bases, bisulfite);
    for ce in &rec.cigar.elements {
        if matches!(ce.op, Op::I | Op::D) {
            nm += ce.length as i32;
        }
    }
    nm
}

/// The read's aligned positions, as `(read offset, reference index)` pairs.
///
/// `SAMRecord.getAlignmentBlocks` builds these from the `M`, `=` and `X` elements; `I` and `S`
/// advance the read, `D` and `N` advance the reference, and `H` and `P` advance neither.
fn walk_alignment(rec: &BamRecord, mut visit: impl FnMut(usize, usize)) {
    let mut block_ref = (rec.alignment_start - 1) as isize;
    let mut block_read: usize = 0;
    for ce in &rec.cigar.elements {
        let len = ce.length as usize;
        match ce.op {
            Op::M | Op::Eq | Op::X => {
                for i in 0..len {
                    let ref_index = block_ref + i as isize;
                    if ref_index >= 0 {
                        visit(block_read + i, ref_index as usize);
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
}

/// `SequenceUtil.sumQualitiesOfMismatches` (non-bisulfite, referenceOffset 0): the sum of the read's
/// base qualities at positions where it mismatches the reference, walking the alignment (`M`/`=`/`X`)
/// blocks.
fn sum_qualities_of_mismatches_with(rec: &BamRecord, ref_bases: &[u8], bisulfite: bool) -> i32 {
    if !bisulfite {
        return sum_qualities_of_mismatches(rec, ref_bases);
    }
    let negative_strand = rec.flags & READ_REVERSE_STRAND != 0;
    let mut qualities: i32 = 0;
    walk_alignment(rec, |read_offset, ref_index| {
        if ref_index >= ref_bases.len() {
            return;
        }
        if !bisulfite_bases_equal(
            negative_strand,
            rec.read_bases[read_offset],
            ref_bases[ref_index],
        ) {
            qualities += rec.base_qualities[read_offset] as i32;
        }
    });
    qualities
}

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
    fix_record(rec, ref_bases, Options::default());
}

/// The tool's own two arguments, with Picard's defaults.
#[derive(Debug, Clone, Copy, Default)]
pub struct Options {
    /// `IS_BISULFITE_SEQUENCE`: C to T on the positive strand and G to A on the negative one stop
    /// counting as mismatches, for `NM` and `UQ` but not for `MD`.
    pub is_bisulfite_sequence: bool,
    /// `SET_ONLY_UQ`: write `UQ` and leave `MD` and `NM` as the input had them.
    pub set_only_uq: bool,
}

/// `SetNmMdAndUqTags.fixRecord`, which is `fixUq` or `fixNmMdAndUq` depending on `SET_ONLY_UQ`.
///
/// The unmapped guard is the caller's in Picard (`if (!record.getReadUnmappedFlag())`) and is kept
/// here so a caller cannot forget it: an unmapped read has no reference bases to compare against.
pub fn fix_record(rec: &mut BamRecord, ref_bases: &[u8], options: Options) {
    if rec.flags & READ_UNMAPPED != 0 {
        return;
    }
    if !options.set_only_uq {
        // `calculateMdAndNmTags(record, ref, true, !isBisulfiteSequence)`: MD always, NM only when
        // this is not bisulfite data, because bisulfite NM is recomputed just below by a function
        // that forgives the conversion.
        let (md, nm) =
            calculate_md_and_nm(rec.alignment_start, &rec.cigar, &rec.read_bases, ref_bases);
        rec.tags.insert(Tag::new(b"MD"), TagValue::Str(md));
        let nm = if options.is_bisulfite_sequence {
            calculate_sam_nm_tag(rec, ref_bases, true)
        } else {
            nm
        };
        rec.tags.insert(Tag::new(b"NM"), TagValue::Int(nm as i64));
    }

    // fixUq: UQ only when the read carries qualities.
    if !rec.base_qualities.is_empty() {
        let uq = sum_qualities_of_mismatches_with(rec, ref_bases, options.is_bisulfite_sequence);
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
