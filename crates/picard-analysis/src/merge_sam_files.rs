//! `MergeSamFiles`.
//!
//! Ports `picard.sam.MergeSamFiles.doWork` at tag 3.4.0 for **inputs that already share a header**:
//! merge several SAM files into one, sorted by `SORT_ORDER`. The tool builds the output header with
//! `SamFileHeaderMerger` and merges the records with a `MergingSamRecordIterator` (a k-way merge on
//! the `SORT_ORDER` comparator), and it adds **no** `@PG` (only optional `@CO`, default none), so the
//! output compares raw.
//!
//! Scope of this slice: the inputs share a **sequence dictionary** and their `@RG`/`@PG` IDs do not
//! collide with differing content. The merged header always gains one change `SamFileHeaderMerger`
//! makes (`SamFileHeaderMerger.java:185`): the group order is set to `none`. htsjdk builds the merged
//! `@HD` fresh, so its attributes come out in insertion order `VN`, `GO`, `SO` (not the input's `VN`,
//! `SO`). The `@RG` and `@PG` records are the **union** across inputs in file order (an ID reused with
//! identical content is deduped; distinct IDs are all kept), and `@CO` is every input's comments
//! concatenated (the merger does not dedupe comments). Records from all inputs are concatenated and
//! **stably** sorted by the `SORT_ORDER` comparator; the coordinate/queryname comparators fully order
//! any two distinct records, so this reproduces the k-way merge except for wholly-identical records,
//! where the stable sort keeps input order (first file first), as the merge does.
//!
//! Deferred to a further slice: `SamFileHeaderMerger`'s collision **renaming** (a reused `@RG`/`@PG`
//! ID with different content becomes `ID.1` and its records are remapped), and the sequence-dictionary
//! **union** (`MERGE_SEQUENCE_DICTIONARIES`). Both are reported as errors here rather than merged.

use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::text_parse::ParseError;

use crate::sort_sam::SortOrder;

/// Why a merge could not run.
#[derive(Debug)]
pub enum MergeError {
    Parse(ParseError),
    /// Two inputs declare the same `@RG` ID with different content. `SamFileHeaderMerger` renames the
    /// later one (`ID` -> `ID.1`) and remaps its records; that renaming is a separate slice.
    ReadGroupCollision(String),
    /// The same for a `@PG` ID.
    ProgramCollision(String),
    /// The inputs' sequence dictionaries differ. `MERGE_SEQUENCE_DICTIONARIES` and the general
    /// dictionary union are a separate slice.
    SequenceDictionaryMismatch,
}

impl From<ParseError> for MergeError {
    fn from(e: ParseError) -> Self {
        MergeError::Parse(e)
    }
}

/// Union the `@RG` records across the input headers in file order: an ID seen again with identical
/// content is deduped; a new ID is appended; an ID reused with **different** content is a collision
/// (renaming deferred). This covers both the identical-header case and distinct read groups.
fn union_read_groups(headers: &[SamHeader]) -> Result<Vec<ReadGroup>, MergeError> {
    let mut out: Vec<ReadGroup> = Vec::new();
    for h in headers {
        for rg in &h.read_groups {
            match out.iter().find(|e| e.id == rg.id) {
                Some(existing) if existing == rg => {}
                Some(_) => return Err(MergeError::ReadGroupCollision(rg.id.clone())),
                None => out.push(rg.clone()),
            }
        }
    }
    Ok(out)
}

/// The same union for `@PG` records.
fn union_programs(headers: &[SamHeader]) -> Result<Vec<ProgramRecord>, MergeError> {
    let mut out: Vec<ProgramRecord> = Vec::new();
    for h in headers {
        for pg in &h.programs {
            match out.iter().find(|e| e.id == pg.id) {
                Some(existing) if existing == pg => {}
                Some(_) => return Err(MergeError::ProgramCollision(pg.id.clone())),
                None => out.push(pg.clone()),
            }
        }
    }
    Ok(out)
}

/// The merged output header: a fresh `@HD` carrying `VN`, then `GO:none`, then the output `SO`; the
/// `@SQ` dictionary of the (identical) inputs; the union of their `@RG` and `@PG`; and every input's
/// `@CO` concatenated in order (`SamFileHeaderMerger` does not dedupe comments). Ports
/// `SamFileHeaderMerger`, which always sets `GroupOrder.none`, for inputs whose sequence dictionaries
/// match and whose read-group/program IDs do not collide with differing content.
fn merged_header(
    headers: &[SamHeader],
    comments: Vec<String>,
    order: SortOrder,
) -> Result<SamHeader, MergeError> {
    let sequences: Vec<SequenceRecord> = headers
        .first()
        .map(|h| h.sequences.clone())
        .unwrap_or_default();
    if headers.iter().any(|h| h.sequences != sequences) {
        return Err(MergeError::SequenceDictionaryMismatch);
    }

    let mut h = SamHeader::new(); // @HD VN:<current>
    h.set_group_order("none"); // VN, GO
    h.set_sort_order(order.name()); // VN, GO, SO
    h.sequences = sequences;
    h.read_groups = union_read_groups(headers)?;
    h.programs = union_programs(headers)?;
    h.comments = comments;
    Ok(h)
}

/// The merge itself: the merged header and the merged, sorted records. Shared by the SAM and BAM
/// renderers so the header construction and record ordering cannot drift.
fn merge_records(
    inputs: &[&str],
    order: SortOrder,
) -> Result<(SamHeader, Vec<BamRecord>), MergeError> {
    let mut headers: Vec<SamHeader> = Vec::with_capacity(inputs.len());
    let mut comments: Vec<String> = Vec::new();
    let mut records: Vec<BamRecord> = Vec::new();
    for input in inputs {
        let (header, recs) = read_sam(input)?;
        comments.extend(header.comments.iter().cloned());
        headers.push(header);
        records.extend(recs);
    }
    let header = merged_header(&headers, comments, order)?;

    // A stable sort so wholly-identical records keep input order (decision 0021).
    records.sort_by(order.comparator());
    Ok((header, records))
}

/// `MergeSamFiles.doWork` for SAM inputs whose dictionaries match and read groups do not collide: the
/// merged, sorted SAM.
pub fn merge_sam_files(inputs: &[&str], order: SortOrder) -> Result<String, MergeError> {
    let (header, records) = merge_records(inputs, order)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// The same merge with **BAM** output, written through htsjdk-rs's byte-identical `BamWriter`.
/// Byte-identical to Picard with `USE_JDK_DEFLATER=true` (the merged records are those the SAM path
/// already reproduces).
pub fn merge_sam_files_to_bam(inputs: &[&str], order: SortOrder) -> Result<Vec<u8>, MergeError> {
    use htsjdk_bam::writer::BamWriter;

    let (header, records) = merge_records(inputs, order)?;
    let mut w = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        w.write(rec).expect("record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: &str =
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s\tLB:lib1\n";

    fn rec(name: &str, start: i32) -> String {
        format!("{name}\t0\tchr1\t{start}\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n")
    }

    #[test]
    fn merges_and_sorts_under_a_group_ordered_header() {
        let a = format!("{H}{}{}", rec("a", 10), rec("c", 30));
        let b = format!("{H}{}{}", rec("b", 20), rec("d", 40));
        let out = merge_sam_files(&[&a, &b], SortOrder::Coordinate).unwrap();
        assert_eq!(
            out,
            "@HD\tVN:1.6\tGO:none\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:1000\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\n\
             a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
             b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
             c\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
             d\t0\tchr1\t40\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
        );
    }

    #[test]
    fn identical_read_groups_are_not_duplicated() {
        let a = format!("{H}{}", rec("a", 10));
        let b = format!("{H}{}", rec("b", 20));
        let out = merge_sam_files(&[&a, &b], SortOrder::Coordinate).unwrap();
        assert_eq!(out.matches("@RG\t").count(), 1);
    }

    #[test]
    fn bam_output_decodes_to_the_same_as_the_sam_output() {
        use htsjdk_bam::reader::BamReader;
        let a = format!("{H}{}{}", rec("a", 10), rec("c", 30));
        let b = format!("{H}{}{}", rec("b", 20), rec("d", 40));
        let sam = merge_sam_files(&[&a, &b], SortOrder::Coordinate).unwrap();
        let bam = merge_sam_files_to_bam(&[&a, &b], SortOrder::Coordinate).unwrap();
        let decoded = htsjdk_bgzf::decompress_all(&bam).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        let header = reader.header.text.clone();
        let recs: Vec<_> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(write_sam(&header, &recs).unwrap(), sam);
    }

    #[test]
    fn distinct_read_groups_are_unioned() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s1\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg2\tSM:s2\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n";
        let out = merge_sam_files(&[a, b], SortOrder::Coordinate).unwrap();
        assert!(out.contains("@RG\tID:rg1\tSM:s1"));
        assert!(out.contains("@RG\tID:rg2\tSM:s2"));
    }

    #[test]
    fn a_colliding_read_group_id_with_different_content_is_an_error() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s1\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s2\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        assert!(matches!(
            merge_sam_files(&[a, b], SortOrder::Coordinate),
            Err(MergeError::ReadGroupCollision(_))
        ));
    }

    #[test]
    fn a_mismatched_sequence_dictionary_is_an_error() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:2000\n\
            b\t0\tchr2\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        assert!(matches!(
            merge_sam_files(&[a, b], SortOrder::Coordinate),
            Err(MergeError::SequenceDictionaryMismatch)
        ));
    }
}
