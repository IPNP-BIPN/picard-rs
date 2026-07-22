//! `ReorderSam`.
//!
//! Ports `picard.sam.ReorderSam.doWork` at tag 3.4.0, for **SAM input** (the unindexed path). Reorders
//! a file's reads to a second reference's contig ordering, matched by exact contig name: each read's
//! reference index and mate reference index are remapped to the new dictionary, a read (or mate) on a
//! contig absent from the new reference is turned unmapped, and the output header carries the new
//! dictionary in place of the old.
//!
//! Two things fix the byte layout of the output. First, `doWork` **clones the input header** and only
//! swaps in the new sequence dictionary, adding no `@PG` and no timestamp, so the header is the input
//! header with a different `@SQ` block and the whole file is comparable raw. Second, the writer is
//! `makeWriter(outHeader, presorted=false, ...)`, so for a `SO:coordinate` or `SO:queryname` header it
//! **re-sorts** the remapped records through a `SortingCollection` (`SAMFileWriterImpl.init`,
//! l.172-174); for `SO:unsorted` it writes them in input order. The remapped indices feed that sort,
//! so an unmapped read lands last under the coordinate order.
//!
//! The remap is a per-record transform reading only the shared index map, so it runs on all cores and
//! stays byte-identical (decision 0006); the re-sort must be a **stable** in-memory sort for
//! byte-identity (decision 0021), which it is.
//!
//! Only the unindexed (SAM) path is claimed here. The indexed-BAM path (`in.hasIndex()`) writes reads
//! grouped by the *input* contig order via per-contig queries, a different output order, and is a
//! separate surface.

use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, read_sam_with, write_sam};
use htsjdk_bam::tag::Tag;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_bam::{coordinate, query_name};
use rayon::prelude::*;

/// `SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX`.
const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;
/// `SAMRecord.NO_ALIGNMENT_START`.
const NO_ALIGNMENT_START: i32 = 0;
/// `SAMFlag.READ_UNMAPPED`.
const READ_UNMAPPED: u16 = 0x4;
/// `SAMFlag.MATE_UNMAPPED`.
const MATE_UNMAPPED: u16 = 0x8;

/// `ReorderSam`'s two switches, both defaulting to Picard's defaults (false).
#[derive(Debug, Clone, Copy, Default)]
pub struct Options {
    /// `ALLOW_INCOMPLETE_DICT_CONCORDANCE`: if true, a read contig with no match in the new reference
    /// maps to unmapped instead of aborting.
    pub allow_incomplete_dict_concordance: bool,
    /// `ALLOW_CONTIG_LENGTH_DISCORDANCE`: if true, a name match with a different length warns instead
    /// of aborting.
    pub allow_contig_length_discordance: bool,
}

/// Why `ReorderSam` could not produce an output, mirroring the `PicardException`s it throws plus the
/// sort orders this port does not cover.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReorderError {
    /// The input SAM did not parse.
    Parse(ParseError),
    /// `buildSequenceDictionaryMap`: a name match with a different length and no
    /// `ALLOW_CONTIG_LENGTH_DISCORDANCE`.
    DiscordantContigLength {
        read_name: String,
        read_length: i32,
        ref_name: String,
        ref_length: i32,
    },
    /// `buildSequenceDictionaryMap`: a read contig with no match in the new reference and no
    /// `ALLOW_INCOMPLETE_DICT_CONCORDANCE`.
    MissingContig(String),
    /// The output header's sort order is one this port does not write (only `coordinate`,
    /// `queryname`, and `unsorted` are handled).
    UnsupportedSortOrder(String),
}

impl From<ParseError> for ReorderError {
    fn from(e: ParseError) -> Self {
        ReorderError::Parse(e)
    }
}

/// `buildSequenceDictionaryMap(refDict, readsDict)`: for each read-dictionary index, the new
/// reference index, or `-1` (unmapped) when the contig is absent and incompleteness is allowed.
///
/// The read-dictionary indices are `0..reads.len()`, so the map is a dense `Vec` rather than a
/// `HashMap`; every entry is filled exactly as Java fills its map, so the lookup never misses.
fn build_new_index(
    reads: &[SequenceRecord],
    reference: &[SequenceRecord],
    opts: &Options,
) -> Result<Vec<i32>, ReorderError> {
    let mut new_index = vec![i32::MIN; reads.len()];

    // First pass, over the reference dictionary: fill in every read contig that has a name match,
    // checking lengths as htsjdk does.
    for (ref_idx, ref_rec) in reference.iter().enumerate() {
        if let Some((read_idx, read_rec)) = reads
            .iter()
            .enumerate()
            .find(|(_, r)| r.name == ref_rec.name)
        {
            if ref_rec.length != read_rec.length && !opts.allow_contig_length_discordance {
                return Err(ReorderError::DiscordantContigLength {
                    read_name: read_rec.name.clone(),
                    read_length: read_rec.length,
                    ref_name: ref_rec.name.clone(),
                    ref_length: ref_rec.length,
                });
            }
            new_index[read_idx] = ref_idx as i32;
        }
    }

    // Second pass, over the read dictionary: any contig still unmapped is either sent to unmapped
    // (incompleteness allowed) or aborts the run.
    for (read_idx, read_rec) in reads.iter().enumerate() {
        if new_index[read_idx] == i32::MIN {
            if opts.allow_incomplete_dict_concordance {
                new_index[read_idx] = NO_ALIGNMENT_REFERENCE_INDEX;
            } else {
                return Err(ReorderError::MissingContig(read_rec.name.clone()));
            }
        }
    }

    Ok(new_index)
}

/// `ReorderSam.newOrderIndex`: `-1` stays `-1` (unmapped); otherwise look up the new index.
fn new_order_index(old_index: i32, new_index: &[i32]) -> i32 {
    if old_index == NO_ALIGNMENT_REFERENCE_INDEX {
        NO_ALIGNMENT_REFERENCE_INDEX
    } else {
        new_index[old_index as usize]
    }
}

/// `ReorderSam.writeReads`' per-record update: remap the reference and mate reference indices, and
/// turn a read (or its mate) unmapped when its contig has left the reference.
fn reorder_record(rec: &mut BamRecord, new_index: &[i32]) {
    let old_ref = rec.reference_index;
    let old_mate = rec.mate_reference_index;

    let new_ref = new_order_index(old_ref, new_index);
    rec.reference_index = new_ref;

    // Read becoming unmapped.
    if old_ref != NO_ALIGNMENT_REFERENCE_INDEX && new_ref == NO_ALIGNMENT_REFERENCE_INDEX {
        rec.alignment_start = NO_ALIGNMENT_START;
        rec.flags |= READ_UNMAPPED;
        rec.cigar = Cigar::default(); // SAMRecord.NO_ALIGNMENT_CIGAR, "*"
        rec.mapping_quality = 0; // SAMRecord.NO_MAPPING_QUALITY
    }

    let new_mate = new_order_index(old_mate, new_index);
    // Mate becoming unmapped.
    if old_mate != NO_ALIGNMENT_REFERENCE_INDEX && new_mate == NO_ALIGNMENT_REFERENCE_INDEX {
        rec.mate_alignment_start = NO_ALIGNMENT_START;
        rec.flags |= MATE_UNMAPPED;
        rec.tags.remove(Tag::new(b"MC")); // Set the Mate Cigar String to null
    }
    rec.mate_reference_index = new_mate;
}

/// A record comparator, as the sort orders expose one.
type RecordComparator = fn(&BamRecord, &BamRecord) -> std::cmp::Ordering;

/// The comparator the `presorted=false` writer would apply for the output header's sort order, or
/// `None` for `unsorted` (write in input order). Errors on an order this port does not write.
fn writer_sort(header: &SamHeader) -> Result<Option<RecordComparator>, ReorderError> {
    // `SAMFileHeader.getSortOrder()` defaults an absent `SO` to `unsorted`.
    match header.attributes.get("SO").unwrap_or("unsorted") {
        "coordinate" => Ok(Some(coordinate::compare)),
        "queryname" => Ok(Some(query_name::compare)),
        "unsorted" => Ok(None),
        other => Err(ReorderError::UnsupportedSortOrder(other.to_string())),
    }
}

/// `ReorderSam.doWork` for SAM input and output. `dict_text` is the `SEQUENCE_DICTIONARY` file (a
/// `.dict`, or any header-bearing text); its `@SQ` records are the new reference dictionary.
pub fn reorder_sam(
    input_sam: &str,
    dict_text: &str,
    opts: &Options,
) -> Result<String, ReorderError> {
    // The input is read leniently, as the tool reads it (VALIDATION_STRINGENCY does not reach the
    // bytes; it only decides which records are admitted).
    let (input_header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    let (dict_header, _) = read_sam(dict_text)?;
    let reference = dict_header.sequences;

    let new_index = build_new_index(&input_header.sequences, &reference, opts)?;

    // outHeader = in.getFileHeader().clone(); outHeader.setSequenceDictionary(outputDictionary).
    let mut out_header = input_header.clone();
    out_header.sequences = reference;

    // The remap is independent per record and order-preserving (decision 0006).
    records
        .par_iter_mut()
        .for_each(|rec| reorder_record(rec, &new_index));

    // The presorted=false writer re-sorts by the header's sort order; a stable sort keeps records
    // that compare equal in input order (decision 0021).
    if let Some(cmp) = writer_sort(&out_header)? {
        records.sort_by(cmp);
    }

    Ok(write_sam(&out_header, &records).expect("records that parsed re-encode as SAM text"))
}

#[cfg(test)]
mod tests {
    use super::*;

    // A read dictionary chr1, chr2, chr3; the new reference swaps chr1/chr2 and drops chr3.
    const DICT: &str = "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n";

    // Coordinate-sorted in the input order: chr1, chr2, chr3.
    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @SQ\tSN:chr2\tLN:1000\n\
        @SQ\tSN:chr3\tLN:1000\n\
        r1\t1\tchr1\t100\t60\t4M\tchr3\t500\t0\tACGT\tIIII\tMC:Z:4M\n\
        r2\t0\tchr2\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        r3\t0\tchr3\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

    fn incomplete() -> Options {
        Options {
            allow_incomplete_dict_concordance: true,
            allow_contig_length_discordance: false,
        }
    }

    fn field(line: &str, n: usize) -> &str {
        line.split('\t').nth(n).unwrap()
    }

    fn record<'a>(sam: &'a str, name: &str) -> &'a str {
        sam.lines().find(|l| l.starts_with(name)).unwrap()
    }

    #[test]
    fn the_output_header_carries_the_new_dictionary_in_the_new_order() {
        let out = reorder_sam(INPUT, DICT, &incomplete()).unwrap();
        // The @SQ block is the new reference: chr2 then chr1, and chr3 is gone.
        assert!(
            out.contains("@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n"),
            "got {out}"
        );
        assert!(!out.contains("chr3\tLN"), "chr3 must be dropped: {out}");
        // The @HD is the input's, unchanged.
        assert!(out.starts_with("@HD\tVN:1.6\tSO:coordinate\n"));
    }

    #[test]
    fn reads_are_re_sorted_into_the_new_coordinate_order() {
        let out = reorder_sam(INPUT, DICT, &incomplete()).unwrap();
        let names: Vec<&str> = out
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| field(l, 0))
            .collect();
        // New indices: chr2=0 (r2), chr1=1 (r1); r3 and u1 are unmapped and sort last, r3 before u1
        // by name.
        assert_eq!(names, ["r2", "r1", "r3", "u1"]);
    }

    #[test]
    fn a_read_on_a_dropped_contig_becomes_unmapped() {
        let out = reorder_sam(INPUT, DICT, &incomplete()).unwrap();
        let r3 = record(&out, "r3");
        assert_eq!(field(r3, 1), "4", "unmapped flag set");
        assert_eq!(field(r3, 2), "*", "RNAME cleared");
        assert_eq!(field(r3, 3), "0", "POS cleared");
        assert_eq!(field(r3, 4), "0", "MAPQ zeroed");
        assert_eq!(field(r3, 5), "*", "CIGAR cleared");
    }

    #[test]
    fn a_mate_on_a_dropped_contig_becomes_unmapped_and_loses_its_mate_cigar() {
        let out = reorder_sam(INPUT, DICT, &incomplete()).unwrap();
        let r1 = record(&out, "r1");
        // r1 was paired (0x1); its mate on chr3 is dropped, so MATE_UNMAPPED (0x8) is set => 0x9.
        assert_eq!(field(r1, 1), "9");
        assert_eq!(field(r1, 6), "*", "RNEXT cleared");
        assert_eq!(field(r1, 7), "0", "PNEXT cleared");
        assert!(!r1.contains("MC:Z"), "mate cigar removed: {r1}");
        // r1 itself stays mapped on chr1 (new index 1).
        assert_eq!(field(r1, 2), "chr1");
        assert_eq!(field(r1, 3), "100");
    }

    #[test]
    fn a_missing_contig_without_allow_incomplete_is_an_error() {
        let err = reorder_sam(INPUT, DICT, &Options::default()).unwrap_err();
        assert_eq!(err, ReorderError::MissingContig("chr3".to_string()));
    }

    #[test]
    fn a_length_mismatch_without_allow_discordance_is_an_error() {
        let dict =
            "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:999\n@SQ\tSN:chr3\tLN:1000\n";
        let err = reorder_sam(INPUT, dict, &Options::default()).unwrap_err();
        assert_eq!(
            err,
            ReorderError::DiscordantContigLength {
                read_name: "chr1".to_string(),
                read_length: 1000,
                ref_name: "chr1".to_string(),
                ref_length: 999,
            }
        );
    }
}
