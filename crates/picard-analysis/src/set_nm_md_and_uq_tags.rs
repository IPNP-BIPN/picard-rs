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

/// `SetNmMdAndUqTags.doWork` for SAM input and output, default (non-bisulfite) options. `fasta` is the
/// `REFERENCE_SEQUENCE` bytes.
pub fn set_nm_md_and_uq_tags(input_sam: &str, fasta: &[u8]) -> Result<String, SetTagsError> {
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

        let (md, nm) =
            calculate_md_and_nm(rec.alignment_start, &rec.cigar, &rec.read_bases, ref_bases);
        rec.tags.insert(Tag::new(b"MD"), TagValue::Str(md));
        rec.tags.insert(Tag::new(b"NM"), TagValue::Int(nm as i64));

        // fixUq: UQ only when the read carries qualities.
        if !rec.base_qualities.is_empty() {
            let uq = sum_qualities_of_mismatches(rec, ref_bases);
            rec.tags.insert(Tag::new(b"UQ"), TagValue::Int(uq as i64));
        }
    }

    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}
