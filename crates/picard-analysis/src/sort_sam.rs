//! `SortSam`.
//!
//! Ports `picard.sam.SortSam.doWork` at tag 3.4.0: read a SAM/BAM, set the header's sort order, and
//! rewrite the records in that order. A record transform whose byte-identity is the whole output
//! file.
//!
//! `doWork` keeps the **input header unchanged except for the `SO` field**, and adds no `@PG` and no
//! timestamp, so the output is comparable raw. The records are then sorted by the chosen order's
//! comparator: `SAMRecordCoordinateComparator` for `coordinate`, `SAMRecordQueryNameComparator` for
//! `queryname`. Both comparators treat records that agree on every field as equal, so the sort must
//! be **stable** to be byte-identical (decision 0021); Rust's `sort_by` is.
//!
//! This ports the `coordinate` and `queryname` orders. The `unsorted` and `duplicate` orders are
//! separate surfaces.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::text_parse::ParseError;
use htsjdk_bam::{coordinate, query_name};

/// The `SORT_ORDER` values this port handles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Queryname,
}

impl SortOrder {
    /// `SAMFileHeader.SortOrder.getSortOrder()`: the string written into the `@HD SO` field.
    fn name(self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
        }
    }

    fn comparator(self) -> fn(&BamRecord, &BamRecord) -> std::cmp::Ordering {
        match self {
            SortOrder::Coordinate => coordinate::compare,
            SortOrder::Queryname => query_name::compare,
        }
    }
}

/// `SortSam.doWork` for SAM input and output: set the sort order on the input header, sort the
/// records, and write. The header is otherwise the input header, byte for byte.
pub fn sort_sam(input_sam: &str, order: SortOrder) -> Result<String, ParseError> {
    let (mut header, mut records) = read_sam(input_sam)?;
    header.set_sort_order(order.name());
    // A stable sort so equal records keep input order (decision 0021).
    records.sort_by(order.comparator());
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "@HD\tVN:1.6\tSO:unsorted\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @SQ\tSN:chr2\tLN:1000\n\
        r3\t0\tchr1\t500\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        r1\t0\tchr2\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        r2\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    fn body(sam: &str) -> Vec<String> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    }

    #[test]
    fn coordinate_sorts_by_reference_then_position() {
        let out = sort_sam(INPUT, SortOrder::Coordinate).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:coordinate"));
        // chr1:100 (r2), chr1:500 (r3), then chr2:100 (r1).
        assert_eq!(body(&out), ["r2", "r3", "r1"]);
    }

    #[test]
    fn queryname_sorts_by_name() {
        let out = sort_sam(INPUT, SortOrder::Queryname).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:queryname"));
        assert_eq!(body(&out), ["r1", "r2", "r3"]);
    }

    #[test]
    fn only_the_sort_order_field_of_the_header_changes() {
        let out = sort_sam(INPUT, SortOrder::Coordinate).unwrap();
        // The sequence dictionary is preserved verbatim and in order.
        assert!(out.contains("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000"));
    }
}
