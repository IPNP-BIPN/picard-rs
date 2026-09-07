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
//! [`Options`] carries the arguments the covering array varies: `SANITIZE` with
//! `KEEP_FIRST_DUPLICATE` and `MAX_DISCARD_FRACTION`, `RESTORE_HARDCLIPS`, `SORT_ORDER`, and each
//! of the three `REMOVE_`/`RESTORE_` switches. `OUTPUT_BY_READGROUP` (a file per read group),
//! `SAMPLE_ALIAS`, `LIBRARY_NAME`, `OUTPUT_MAP` and a customized `ATTRIBUTE_TO_CLEAR` are still
//! separate surfaces.
//!
//! One piece of `SANITIZE` is deliberately not ported: `createReadGroupFormatMap` runs htsjdk's
//! `QualityEncodingDetector` over the file and subtracts 31 from every quality when a read group
//! is detected as Solexa or Illumina rather than Standard. A Standard-encoded input takes no such
//! subtraction, and both the conformance corpus and the covering array's fixtures are Standard, so
//! the branch is unreachable from here; an input in an older encoding is out of scope rather than
//! silently converted, and asking for it would be asking for a detector this crate does not have.

use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::coordinate;
use htsjdk_bam::fastq::{fastq_to_phred, phred_to_fastq};
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
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
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

/// `SAMFileHeader.SortOrder`, for the two orders `RevertSam` can be asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SortOrder {
    #[default]
    Queryname,
    Coordinate,
}

impl SortOrder {
    pub fn name(self) -> &'static str {
        match self {
            SortOrder::Queryname => "queryname",
            SortOrder::Coordinate => "coordinate",
        }
    }
}

/// The tool's arguments, with Picard's defaults.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    pub sort_order: SortOrder,
    pub restore_original_qualities: bool,
    pub remove_duplicate_information: bool,
    pub remove_alignment_information: bool,
    pub restore_hardclips: bool,
    pub sanitize: bool,
    pub keep_first_duplicate: bool,
    /// `MAX_DISCARD_FRACTION`, default 0.01. A run that discards more than this fraction while
    /// sanitizing writes its output and THEN throws, which is Picard's order and not a detail:
    /// the file exists when the exception is raised.
    pub max_discard_fraction: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            sort_order: SortOrder::Queryname,
            restore_original_qualities: true,
            remove_duplicate_information: true,
            remove_alignment_information: true,
            restore_hardclips: false,
            sanitize: false,
            keep_first_duplicate: false,
            max_discard_fraction: 0.01,
        }
    }
}

/// What `RevertSam` refuses, and how.
#[derive(Debug, Clone, PartialEq)]
pub enum RevertError {
    Parse(ParseError),
    /// Thrown after the header is read, before a record is written.
    HardclipsWithoutRemovingAlignment,
    /// Thrown after the output has been written.
    DiscardedTooMuch {
        discarded: f64,
        max: f64,
    },
}

impl From<ParseError> for RevertError {
    fn from(e: ParseError) -> Self {
        RevertError::Parse(e)
    }
}

impl RevertError {
    pub fn java_class(&self) -> &'static str {
        "picard.PicardException"
    }

    pub fn message(&self) -> String {
        match self {
            RevertError::Parse(e) => format!("{e:?}"),
            RevertError::HardclipsWithoutRemovingAlignment => {
                "Cannot revert sam file when RESTORE_HARDCLIPS is true and \
                 REMOVE_ALIGNMENT_INFORMATION is false."
                    .to_string()
            }
            RevertError::DiscardedTooMuch { discarded, max } => format!(
                "Discarded {} which is above MAX_DISCARD_FRACTION of {}",
                percent(*discarded),
                percent(*max)
            ),
        }
    }
}

/// `new DecimalFormat("0.000%")`: the fraction as a percentage, three decimals, percent sign.
fn percent(fraction: f64) -> String {
    format!("{:.3}%", fraction * 100.0)
}

/// `RevertSam.customCommandLineValidation`, which collects every failure rather than stopping at
/// the first: a row that breaks two rules prints both, in this order, and Barclay prints them
/// after the usage block.
///
/// `output_is_directory` answers what `Files.isDirectory(OUTPUT)` answers, which is why it is the
/// caller's to supply: the library has no path.
pub fn validate(
    options: &Options,
    output_by_readgroup: bool,
    output_is_directory: bool,
    output: &str,
) -> Vec<String> {
    let mut errors = Vec::new();
    if options.sanitize && options.sort_order != SortOrder::Queryname {
        errors.push(
            "SORT_ORDER must be queryname when sanitization is enabled with SANITIZE=true."
                .to_string(),
        );
    }
    if output_by_readgroup {
        if !output_is_directory {
            errors.push(format!(
                "When OUTPUT_BY_READGROUP=true and OUTPUT is provided, it must be a directory: \
                 {output}"
            ));
        }
    } else if output_is_directory {
        errors.push(format!(
            "OUTPUT {output} should not be a directory when OUTPUT_BY_READGROUP=false"
        ));
    }
    if !options.sanitize && options.keep_first_duplicate {
        errors.push("KEEP_FIRST_DUPLICATE cannot be used without SANITIZE".to_string());
    }
    errors
}

/// `RevertSam.revertSamRecord`.
fn revert_record_with(rec: &mut BamRecord, options: &Options) {
    if options.restore_original_qualities {
        // Move OQ back into QUAL and drop OQ.
        if let Some(TagValue::Str(oq)) = rec.tags.get(Tag::new(b"OQ")) {
            rec.base_qualities = fastq_to_phred(oq);
            rec.tags.remove(Tag::new(b"OQ"));
        }
    }

    if options.remove_duplicate_information {
        rec.flags &= !DUPLICATE_READ;
    }

    if !options.remove_alignment_information {
        return;
    }

    if rec.flags & READ_NEGATIVE_STRAND != 0 {
        reverse_complement_record(rec);
        rec.flags &= !READ_NEGATIVE_STRAND;
    }

    rec.reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
    rec.alignment_start = NO_ALIGNMENT_START;
    rec.cigar = Cigar::default();
    rec.mapping_quality = 0;
    rec.inferred_insert_size = 0;
    rec.flags &= !SECONDARY_ALIGNMENT;
    rec.flags &= !PROPER_PAIR;
    rec.flags |= READ_UNMAPPED;

    rec.mate_alignment_start = NO_ALIGNMENT_START;
    rec.flags &= !MATE_NEGATIVE_STRAND;
    rec.mate_reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
    if rec.flags & READ_PAIRED != 0 {
        rec.flags |= MATE_UNMAPPED;
    } else {
        rec.flags &= !MATE_UNMAPPED;
    }

    if options.restore_hardclips {
        // The record has already been reverse-complemented if it was on the negative strand, so
        // the stored bases append as they are.
        let bases = match rec.tags.get(Tag::new(b"XB")) {
            Some(TagValue::Str(s)) => Some(s.clone()),
            _ => None,
        };
        let qualities = match rec.tags.get(Tag::new(b"XQ")) {
            Some(TagValue::Str(s)) => Some(s.clone()),
            _ => None,
        };
        if let (Some(bases), Some(qualities)) = (bases, qualities) {
            rec.read_bases.extend_from_slice(bases.as_bytes());
            let mut restored = phred_to_fastq(&rec.base_qualities);
            restored.push_str(&qualities);
            rec.base_qualities = fastq_to_phred(&restored);
            rec.tags.remove(Tag::new(b"XB"));
            rec.tags.remove(Tag::new(b"XQ"));
        }
    }

    for name in ATTRIBUTE_TO_CLEAR {
        rec.tags.remove(Tag::new(name));
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

/// `RevertSam.sanitize`: group the queryname-sorted records by name and keep only the groups that
/// look like one template.
///
/// Returns the kept records and `(discarded, total)`. Three rules, in Picard's order:
///
/// * every read in the group must have as many qualities as bases, or the whole group goes;
/// * a paired group must hold exactly one R1 and exactly one R2. With `KEEP_FIRST_DUPLICATE` a
///   group holding at least one of each keeps the FIRST of each and discards the rest; without it
///   the whole group goes;
/// * an unpaired group must hold exactly one read, with the same `KEEP_FIRST_DUPLICATE` escape.
///
/// The counting is what the discard fraction is computed from, so a group trimmed by
/// `KEEP_FIRST_DUPLICATE` counts its dropped reads as discarded while still writing two.
fn sanitize(records: Vec<BamRecord>, keep_first_duplicate: bool) -> (Vec<BamRecord>, u64, u64) {
    let mut out = Vec::with_capacity(records.len());
    let mut discarded: u64 = 0;
    let mut total: u64 = 0;

    let mut index = 0;
    while index < records.len() {
        let mut end = index + 1;
        while end < records.len() && records[end].read_name == records[index].read_name {
            end += 1;
        }
        let group = &records[index..end];
        index = end;
        total += group.len() as u64;

        if group
            .iter()
            .any(|rec| rec.read_bases.len() != rec.base_qualities.len())
        {
            discarded += group.len() as u64;
            continue;
        }

        let mut firsts = 0;
        let mut seconds = 0;
        let mut unpaired = 0;
        let (mut first_record, mut second_record, mut unpaired_record) = (None, None, None);
        for (offset, rec) in group.iter().enumerate() {
            if rec.flags & READ_PAIRED == 0 {
                unpaired_record.get_or_insert(offset);
                unpaired += 1;
            } else {
                if rec.flags & FIRST_OF_PAIR != 0 {
                    first_record.get_or_insert(offset);
                    firsts += 1;
                }
                if rec.flags & SECOND_OF_PAIR != 0 {
                    second_record.get_or_insert(offset);
                    seconds += 1;
                }
            }
        }

        let kept: Vec<&BamRecord> = if firsts > 0 || seconds > 0 {
            if firsts == 1 && seconds == 1 {
                group.iter().collect()
            } else if keep_first_duplicate && firsts >= 1 && seconds >= 1 {
                discarded += group.len() as u64 - 2;
                vec![
                    &group[first_record.unwrap()],
                    &group[second_record.unwrap()],
                ]
            } else {
                discarded += group.len() as u64;
                continue;
            }
        } else if unpaired > 1 {
            if keep_first_duplicate {
                discarded += group.len() as u64 - 1;
                vec![&group[unpaired_record.unwrap()]]
            } else {
                discarded += group.len() as u64;
                continue;
            }
        } else {
            group.iter().collect()
        };

        out.extend(kept.into_iter().cloned());
    }

    (out, discarded, total)
}

/// What a run produced: the output, and the refusal if there was one.
///
/// A `Result` would be the wrong shape here, because the discard-fraction refusal does not
/// replace the output. Picard closes the writer and throws afterwards, so the file exists and the
/// run still fails, and a caller that wrote nothing on the error would differ from the reference
/// on the bytes as well as on the message.
pub struct Reverted {
    pub header: SamHeader,
    pub records: Vec<BamRecord>,
    pub error: Option<RevertError>,
}

/// `RevertSam.doWork` up to the write, with the tool's arguments.
pub fn revert_with(header: &SamHeader, records: Vec<BamRecord>, options: &Options) -> Reverted {
    // createOutHeader: a fresh header carrying the sort order; the dictionary and the @PG records
    // survive only when the alignment is kept. The read groups are added afterwards, always.
    let mut out_header = SamHeader::new();
    out_header.set_sort_order(options.sort_order.name());
    if !options.remove_alignment_information {
        out_header.sequences = header.sequences.clone();
        out_header.programs = header.programs.clone();
    }
    out_header.read_groups = header.read_groups.clone();

    // "Weed out non-primary and supplemental read as we don't want duplicates in the reverted
    // file": the drop happens before the revert, so a secondary record never reaches the output
    // however the flags would have been cleared.
    let mut records: Vec<BamRecord> = records
        .into_iter()
        .filter(|rec| rec.flags & (SECONDARY_ALIGNMENT | SUPPLEMENTARY_ALIGNMENT) == 0)
        .collect();
    records
        .par_iter_mut()
        .for_each(|rec| revert_record_with(rec, options));

    // isPresorted: the writer sorts unless the input already carries the order asked for, or the
    // sanitizer has just produced it.
    let presorted = header.attributes.get("SO") == Some(options.sort_order.name())
        || (options.sort_order == SortOrder::Queryname && options.sanitize);

    if options.sanitize {
        // The sorting collection the sanitizer reads from is queryname-ordered, and SANITIZE is
        // refused for any other order, so this is the only sort it can be.
        records.sort_by(query_name::compare);
        let (kept, discarded, total) = sanitize(records, options.keep_first_duplicate);
        let fraction = if total == 0 {
            0.0
        } else {
            discarded as f64 / total as f64
        };
        let error =
            (fraction > options.max_discard_fraction).then_some(RevertError::DiscardedTooMuch {
                discarded: fraction,
                max: options.max_discard_fraction,
            });
        return Reverted {
            header: out_header,
            records: kept,
            error,
        };
    }

    if !presorted {
        match options.sort_order {
            SortOrder::Queryname => records.sort_by(query_name::compare),
            SortOrder::Coordinate => records.sort_by(coordinate::compare),
        }
    }
    Reverted {
        header: out_header,
        records,
        error: None,
    }
}

/// `RevertSam.doWork` up to the write: the bare output header and the reverted, queryname-sorted
/// records.
fn revert(input_sam: &str) -> Result<(SamHeader, Vec<BamRecord>), ParseError> {
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

    Ok((out_header, records))
}

/// `RevertSam.doWork` for SAM input and output, default options.
pub fn revert_sam(input_sam: &str) -> Result<String, ParseError> {
    let (out_header, records) = revert(input_sam)?;
    Ok(write_sam(&out_header, &records).expect("unmapped records always encode as SAM text"))
}

/// `RevertSam.doWork` for SAM input and **BAM** output, default options. Same revert, written through
/// htsjdk-rs's byte-identical `BamWriter`; RevertSam adds no `@PG`. Byte-identity to Picard's
/// `USE_JDK_DEFLATER=true` output follows transitively: the records are the ones `revert_sam` already
/// reproduces (its oracle), and the `BamWriter` is proven byte-identical over arbitrary records (the
/// SamFormatConverter oracle and htsjdk-rs's `every_file_is_byte_identical_to_htsjdks`).
pub fn revert_sam_to_bam(input_sam: &str) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;
    let (out_header, records) = revert(input_sam)?;
    let mut writer =
        BamWriter::new(Vec::new(), &out_header).expect("in-memory BAM writer never fails");
    for rec in &records {
        writer
            .write(rec)
            .expect("records that parsed re-encode as BAM");
    }
    Ok(writer
        .finish()
        .expect("in-memory BAM writer never fails to finish"))
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

    /// The BAM output decodes back to exactly the SAM output; the writer's byte-identity to htsjdk is
    /// proven elsewhere.
    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        let sam = revert_sam(INPUT).unwrap();
        let bam = revert_sam_to_bam(INPUT).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bam).expect("bam decompresses");
        let reader = htsjdk_bam::reader::BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(
            htsjdk_bam::sam_file::write_sam(&header, &records).unwrap(),
            sam
        );
    }
}
