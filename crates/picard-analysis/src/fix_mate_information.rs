//! `FixMateInformation`.
//!
//! Ports `picard.sam.FixMateInformation.doWork` and the `SamPairUtil.SetMateInfoIterator` it drives,
//! at tag 3.4.0, for a single SAM input with default options. Makes every pair self-consistent:
//! group the reads of each template, call [`set_mate_info`](htsjdk_bam::pair::set_mate_info) on the
//! two primary ends, and write the result in the output sort order.
//!
//! The pipeline is: sort into queryname order (to group templates), fix each template, then re-sort
//! to the output order, which defaults to the input's sort order (`SORT_ORDER` unset). `doWork` adds
//! no `@PG` and no timestamp, so the whole output is comparable raw. `ADD_MATE_CIGAR` defaults true.
//!
//! This ports templates of one or two **primary** reads. A template carrying secondary or
//! supplementary records needs `setMateInformationOnSupplementalAlignment` too, a separate surface;
//! the port asserts there are none rather than emitting a half-fixed template.

use htsjdk_bam::pair::set_mate_info;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_bam::{coordinate, query_name};

const READ_PAIRED: u16 = 0x1;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;

/// The output sort order options this port handles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Queryname,
}

impl SortOrder {
    fn name(self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
        }
    }

    fn from_str(s: &str) -> Option<SortOrder> {
        match s {
            "coordinate" => Some(SortOrder::Coordinate),
            "queryname" => Some(SortOrder::Queryname),
            _ => None,
        }
    }
}

fn is_primary(rec: &BamRecord) -> bool {
    rec.flags & (SECONDARY_ALIGNMENT | SUPPLEMENTARY_ALIGNMENT) == 0
}

/// `SamPairUtil.SetMateInfoIterator`: over queryname-sorted records, fix each template's mate info.
fn fix_templates(records: &mut [BamRecord], set_mate_cigar: bool) {
    let mut start = 0;
    while start < records.len() {
        // The run of records sharing this read name.
        let mut end = start + 1;
        while end < records.len() && records[end].read_name == records[start].read_name {
            end += 1;
        }

        let mut first_primary: Option<usize> = None;
        let mut second_primary: Option<usize> = None;
        for (offset, r) in records[start..end].iter().enumerate() {
            let i = start + offset;
            assert!(
                is_primary(r),
                "FixMateInformation: secondary/supplementary records are not ported"
            );
            if r.flags & READ_PAIRED != 0 {
                if r.flags & FIRST_OF_PAIR != 0 {
                    assert!(first_primary.is_none(), "two first-of-pair primaries");
                    first_primary = Some(i);
                } else if r.flags & SECOND_OF_PAIR != 0 {
                    assert!(second_primary.is_none(), "two second-of-pair primaries");
                    second_primary = Some(i);
                }
            }
        }

        if let (Some(f), Some(s)) = (first_primary, second_primary) {
            // Borrow both ends mutably: split at the higher index.
            let (lo, hi) = (f.min(s), f.max(s));
            let (left, right) = records.split_at_mut(hi);
            let (a, b) = (&mut left[lo], &mut right[0]);
            if f < s {
                set_mate_info(a, b, set_mate_cigar);
            } else {
                set_mate_info(b, a, set_mate_cigar);
            }
        }
        // A lone paired read (missing mate) would throw unless IGNORE_MISSING_MATES; the corpora do
        // not exercise that, and a fragment (unpaired) read simply passes through.

        start = end;
    }
}

/// `FixMateInformation.doWork` for a single SAM input.
pub fn fix_mate_information(
    input_sam: &str,
    sort_order: Option<SortOrder>,
) -> Result<String, ParseError> {
    let (mut header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    // The output order defaults to the input header's sort order.
    let output_order = sort_order.unwrap_or_else(|| {
        header
            .attributes
            .get("SO")
            .and_then(SortOrder::from_str)
            .unwrap_or(SortOrder::Coordinate)
    });

    // Sort into queryname order to group templates, fix them, then re-sort to the output order.
    records.sort_by(query_name::compare);
    fix_templates(&mut records, true);
    match output_order {
        SortOrder::Coordinate => records.sort_by(coordinate::compare),
        SortOrder::Queryname => {} // already queryname-sorted
    }

    header.set_sort_order(output_order.name());
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

#[cfg(test)]
mod tests {
    use super::*;

    // A coordinate-sorted pair with the mate fields absent (RNEXT *, PNEXT 0, TLEN 0, no MC/MQ).
    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        p\t65\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
        p\t145\tchr1\t200\t50\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";

    #[test]
    fn a_pair_gets_consistent_mate_fields() {
        let out = fix_mate_information(INPUT, None).unwrap();
        let rows: Vec<Vec<&str>> = out
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').collect())
            .collect();
        // First end: RNEXT '=', PNEXT 200, TLEN 110, MC/MQ set.
        let first = rows.iter().find(|r| r[3] == "100").unwrap();
        assert_eq!(first[6], "="); // RNEXT
        assert_eq!(first[7], "200"); // PNEXT
        assert_eq!(first[8], "110"); // TLEN
        assert!(first.contains(&"MC:Z:10M"));
        assert!(first.contains(&"MQ:i:50"));
        // Second end: TLEN -110.
        let second = rows.iter().find(|r| r[3] == "200").unwrap();
        assert_eq!(second[8], "-110");
    }

    #[test]
    fn the_output_keeps_the_input_sort_order() {
        let out = fix_mate_information(INPUT, None).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:coordinate"));
    }
}
