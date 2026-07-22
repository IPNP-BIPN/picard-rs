//! `ReplaceSamHeader`.
//!
//! Ports `picard.sam.ReplaceSamHeader.standardReheader` at tag 3.4.0, the **SAM (non-BAM) path**:
//! write the records of `INPUT` under the header of `HEADER`, verbatim. The tool replaces a file's
//! header wholesale, which is how a caller fixes `@RG`/`@CO`/`@PG` lines while keeping the reads.
//!
//! Two invariants make the output comparable raw. `standardReheader` writes the replacement header
//! **unchanged** (adds no `@PG` and no timestamp) and then writes each record with
//! `makeWriter(replacementHeader, presorted=true, ...)`, so the records are emitted **in input order
//! with no re-sort**. A record carries an integer reference index, not a name, and `rec.setHeader`
//! keeps that index: the RNAME/RNEXT are re-resolved against the *new* dictionary at the same index.
//! When the two dictionaries share their `@SQ` block (the usual case) this is a no-op on the reads,
//! and the whole output is the new header followed by the input records byte-for-byte.
//!
//! There is no per-record transform here (the record bytes are unchanged), so there is nothing to
//! parallelize; this port stays serial by nature, not by omission.
//!
//! Before writing, `standardReheader` requires the two headers to declare the **same sort order**,
//! throwing a `PicardException` otherwise. The BAM path (`BamFileIoUtils.reheaderBamFile`, a BGZF
//! block copy) is a separate surface and is not ported here.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::sam_file::{read_sam, read_sam_with, write_sam};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

/// `SAMFileHeader.SortOrder`, in the normalization `getSortOrder()` applies: an absent `SO` is
/// `unsorted`, and an unrecognized value is `unknown`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SortOrder {
    Coordinate,
    Queryname,
    Unsorted,
    Duplicate,
    Unknown,
}

impl SortOrder {
    /// `SAMFileHeader.getSortOrder()`.
    fn of(header: &SamHeader) -> SortOrder {
        match header.attributes.get("SO") {
            None => SortOrder::Unsorted,
            Some("coordinate") => SortOrder::Coordinate,
            Some("queryname") => SortOrder::Queryname,
            Some("unsorted") => SortOrder::Unsorted,
            Some("duplicate") => SortOrder::Duplicate,
            Some(_) => SortOrder::Unknown,
        }
    }

    /// `SortOrder.name()`, as it appears in the mismatch message.
    fn name(self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
            SortOrder::Unsorted => "unsorted",
            SortOrder::Duplicate => "duplicate",
            SortOrder::Unknown => "unknown",
        }
    }
}

/// Why `ReplaceSamHeader` could not produce an output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReplaceHeaderError {
    /// The `INPUT` SAM did not parse.
    Parse(ParseError),
    /// `standardReheader`: the sort orders of `INPUT` and `HEADER` do not agree.
    SortOrderMismatch {
        input: &'static str,
        header: &'static str,
    },
}

impl From<ParseError> for ReplaceHeaderError {
    fn from(e: ParseError) -> Self {
        ReplaceHeaderError::Parse(e)
    }
}

/// `ReplaceSamHeader.standardReheader` for SAM input and output. `header_text` is the `HEADER` file
/// (a stub SAM or any header-bearing text); its header replaces `INPUT`'s.
pub fn replace_sam_header(
    input_sam: &str,
    header_text: &str,
) -> Result<String, ReplaceHeaderError> {
    // The tool opens INPUT with VALIDATION_STRINGENCY.SILENT; stringency does not reach the bytes.
    let (input_header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    let (replacement_header, _) = read_sam(header_text)?;

    let input_order = SortOrder::of(&input_header);
    let header_order = SortOrder::of(&replacement_header);
    if input_order != header_order {
        return Err(ReplaceHeaderError::SortOrderMismatch {
            input: input_order.name(),
            header: header_order.name(),
        });
    }

    // makeWriter(replacementHeader, presorted=true): the records are written in input order, and
    // each keeps its integer reference index, re-resolved against the replacement dictionary.
    Ok(
        write_sam(&replacement_header, &records)
            .expect("records that parsed re-encode as SAM text"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @SQ\tSN:chr2\tLN:1000\n\
        @RG\tID:old\tSM:sampleA\n\
        r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n\
        r2\t0\tchr2\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n";

    // Same @SQ block and sort order, but a different @RG and an added @CO.
    const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @SQ\tSN:chr2\tLN:1000\n\
        @RG\tID:old\tSM:sampleB\tLB:lib1\n\
        @CO\tedited by hand\n";

    #[test]
    fn the_output_carries_the_replacement_header() {
        let out = replace_sam_header(INPUT, HEADER).unwrap();
        // The new @RG (SM:sampleB, LB:lib1) and the @CO are present.
        assert!(
            out.contains("@RG\tID:old\tSM:sampleB\tLB:lib1"),
            "got {out}"
        );
        assert!(out.contains("@CO\tedited by hand"), "got {out}");
        // The old @RG (SM:sampleA) is gone.
        assert!(!out.contains("SM:sampleA"), "old header leaked: {out}");
    }

    #[test]
    fn the_records_are_written_verbatim_in_input_order() {
        let out = replace_sam_header(INPUT, HEADER).unwrap();
        let body: Vec<&str> = out.lines().filter(|l| !l.starts_with('@')).collect();
        assert_eq!(
            body,
            [
                "r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old",
                "r2\t0\tchr2\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old",
            ]
        );
    }

    #[test]
    fn a_sort_order_mismatch_is_an_error() {
        let header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n";
        let err = replace_sam_header(INPUT, header).unwrap_err();
        assert_eq!(
            err,
            ReplaceHeaderError::SortOrderMismatch {
                input: "coordinate",
                header: "queryname",
            }
        );
    }

    #[test]
    fn an_absent_sort_order_reads_as_unsorted_on_both_sides() {
        let input = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@CO\tnote\n";
        // Both default to unsorted, so they agree and the reheader succeeds.
        let out = replace_sam_header(input, header).unwrap();
        assert!(out.contains("@CO\tnote"));
    }
}
