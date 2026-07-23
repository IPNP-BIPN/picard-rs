//! `MergeSamFiles`.
//!
//! Ports `picard.sam.MergeSamFiles.doWork` at tag 3.4.0 for **inputs that already share a header**:
//! merge several SAM files into one, sorted by `SORT_ORDER`. The tool builds the output header with
//! `SamFileHeaderMerger` and merges the records with a `MergingSamRecordIterator` (a k-way merge on
//! the `SORT_ORDER` comparator), and it adds **no** `@PG` (only optional `@CO`, default none), so the
//! output compares raw.
//!
//! Scope of this slice: the inputs share a header (same `@SQ`/`@RG`/`@PG`/`@CO`). Then the merged
//! header is that header with one change, which `SamFileHeaderMerger` always makes
//! (`SamFileHeaderMerger.java:185`): the group order is set to `none`. htsjdk builds the merged `@HD`
//! fresh, so its attributes come out in insertion order `VN`, `GO`, `SO` (not the input's `VN`, `SO`).
//! Records from all inputs are concatenated and **stably** sorted by the `SORT_ORDER` comparator; the
//! coordinate/queryname comparators fully order any two distinct records, so this reproduces the k-way
//! merge except for wholly-identical records, where the stable sort keeps input order (first file
//! first), as the merge does. `SamFileHeaderMerger`'s read-group/program collision renaming and
//! sequence-dictionary union for **differing** headers are a separate slice, deferred.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::text_parse::ParseError;

use crate::sort_sam::SortOrder;

/// The merged output header: a fresh `@HD` carrying `VN`, then `GO:none`, then the output `SO`, with
/// the shared inputs' `@SQ`/`@RG`/`@PG` copied verbatim from the first (identical inputs dedupe to
/// one) and every input's `@CO` concatenated in order (`SamFileHeaderMerger` does not dedupe
/// comments). Ports the identical-header result of `SamFileHeaderMerger`, which always sets
/// `GroupOrder.none`.
fn merged_header(first: &SamHeader, comments: Vec<String>, order: SortOrder) -> SamHeader {
    let mut h = SamHeader::new(); // @HD VN:<current>
    h.set_group_order("none"); // VN, GO
    h.set_sort_order(order.name()); // VN, GO, SO
    h.sequences = first.sequences.clone();
    h.read_groups = first.read_groups.clone();
    h.programs = first.programs.clone();
    h.comments = comments;
    h
}

/// The merge itself: the merged header and the merged, sorted records. Shared by the SAM and BAM
/// renderers so the header construction and record ordering cannot drift.
fn merge_records(
    inputs: &[&str],
    order: SortOrder,
) -> Result<(SamHeader, Vec<BamRecord>), ParseError> {
    let mut first_header: Option<SamHeader> = None;
    let mut comments: Vec<String> = Vec::new();
    let mut records: Vec<BamRecord> = Vec::new();
    for input in inputs {
        let (header, recs) = read_sam(input)?;
        comments.extend(header.comments.iter().cloned());
        if first_header.is_none() {
            first_header = Some(header);
        }
        records.extend(recs);
    }
    let header = merged_header(&first_header.unwrap_or_default(), comments, order);

    // A stable sort so wholly-identical records keep input order (decision 0021).
    records.sort_by(order.comparator());
    Ok((header, records))
}

/// `MergeSamFiles.doWork` for shared-header SAM inputs: the merged, sorted SAM.
pub fn merge_sam_files(inputs: &[&str], order: SortOrder) -> Result<String, ParseError> {
    let (header, records) = merge_records(inputs, order)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// `MergeSamFiles.doWork` for shared-header inputs with **BAM** output: the same merge, written as a
/// BAM through htsjdk-rs's byte-identical `BamWriter`. Byte-identical to Picard with
/// `USE_JDK_DEFLATER=true` (the merged records are those the SAM path already reproduces).
pub fn merge_sam_files_to_bam(inputs: &[&str], order: SortOrder) -> Result<Vec<u8>, ParseError> {
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
}
