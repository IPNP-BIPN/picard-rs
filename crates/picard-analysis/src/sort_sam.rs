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
//! This ports the `coordinate` and `queryname` orders, and the precondition the `duplicate` order
//! is refused on. A duplicate sort compares records by, among other things, their MATE's unclipped
//! coordinate, which it reads out of the `MC` tag; a file whose records do not carry one cannot be
//! sorted that way at all, and the reference says so before it writes a single record. The
//! ordering itself is a separate surface, and unreachable from a corpus without mate cigars.

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

/// What `--SORT_ORDER` may name.
///
/// Two of the three are sorts this port performs. The third is the duplicate order, which is
/// refused on any file whose paired records carry no mate cigar, and that is every file this
/// corpus has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequestedOrder {
    Coordinate,
    Queryname,
    Duplicate,
}

impl RequestedOrder {
    pub fn parse(name: &str) -> Option<RequestedOrder> {
        match name {
            "coordinate" => Some(RequestedOrder::Coordinate),
            "queryname" => Some(RequestedOrder::Queryname),
            "duplicate" => Some(RequestedOrder::Duplicate),
            _ => None,
        }
    }

    /// The `@HD SO` the output header carries.
    pub fn name(self) -> &'static str {
        match self {
            RequestedOrder::Coordinate => "coordinate",
            RequestedOrder::Queryname => "queryname",
            RequestedOrder::Duplicate => "duplicate",
        }
    }

    /// The sort this port performs for it, where it performs one.
    pub fn sort_order(self) -> Option<SortOrder> {
        match self {
            RequestedOrder::Coordinate => Some(SortOrder::Coordinate),
            RequestedOrder::Queryname => Some(SortOrder::Queryname),
            RequestedOrder::Duplicate => None,
        }
    }
}

/// One record, as much of it as the mate-cigar precondition looks at.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MateCigarCheck {
    pub read_name: String,
    pub paired: bool,
    pub first_of_pair: bool,
    pub unmapped: bool,
    pub mate_unmapped: bool,
    pub read_length: usize,
    pub reference: Option<String>,
    pub start: i32,
    pub end: i32,
    pub mate_cigar: Option<String>,
}

impl MateCigarCheck {
    /// `SAMRecord.toString()`, which is what the refusal names the record with.
    pub fn description(&self) -> String {
        let mut text = self.read_name.clone();
        if self.paired {
            text.push_str(if self.first_of_pair { " 1/2" } else { " 2/2" });
        }
        text.push(' ');
        text.push_str(&format!("{}b", self.read_length));
        if self.unmapped {
            text.push_str(" unmapped read.");
        } else {
            text.push_str(&format!(
                " aligned to {}:{}-{}.",
                self.reference.clone().unwrap_or_default(),
                self.start,
                self.end
            ));
        }
        text
    }

    /// Whether sorting by duplicate order would ask this record for its mate's coordinate.
    pub fn needs_mate_cigar(&self) -> bool {
        self.paired && !self.unmapped && !self.mate_unmapped
    }
}

/// The order the comparator touches the records in, which is what decides WHICH record a refusal
/// names.
///
/// The sort compares the second record against the first before anything else, so the second is
/// the first one asked for its mate's coordinate; after that each record is touched as it is
/// inserted, in file order.
pub fn touch_order(count: usize) -> Vec<usize> {
    match count {
        0 => Vec::new(),
        1 => vec![0],
        _ => {
            let mut order = vec![1, 0];
            order.extend(2..count);
            order
        }
    }
}

/// The refusal a duplicate sort answers a file with no mate cigars with.
///
/// `SAMUtils.getMateUnclippedStart` throws on the first record that needs a mate cigar and has
/// none, and the message carries the record as `SAMRecord.toString()` writes it. The space after
/// the colon is the unclipped-coordinate methods'; the mate-alignment-end method next to them
/// writes the same sentence WITHOUT it, which is a difference a reader would not invent.
pub fn mate_cigar_refusal(records: &[MateCigarCheck]) -> Option<String> {
    for index in touch_order(records.len()) {
        let record = &records[index];
        if record.needs_mate_cigar() && record.mate_cigar.is_none() {
            return Some(format!(
                "Mate CIGAR (Tag MC) not found: {}",
                record.description()
            ));
        }
    }
    None
}

impl SortOrder {
    /// `SAMFileHeader.SortOrder.getSortOrder()`: the string written into the `@HD SO` field.
    pub fn name(self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
        }
    }

    pub fn comparator(self) -> fn(&BamRecord, &BamRecord) -> std::cmp::Ordering {
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

/// `SortSam.doWork` for SAM input and **BAM** output: the same sort, written as a BAM through
/// htsjdk-rs's byte-identical `BamWriter`. The returned bytes are the whole BAM file. Byte-identical
/// to Picard run with `USE_JDK_DEFLATER=true` (the BGZF blocks come from the port's zlib writer);
/// Picard's default GKL deflater is a separate surface.
pub fn sort_sam_to_bam(input_sam: &str, order: SortOrder) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;

    let (mut header, mut records) = read_sam(input_sam)?;
    header.set_sort_order(order.name());
    records.sort_by(order.comparator());

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
