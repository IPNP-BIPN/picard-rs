//! `RevertOriginalBaseQualitiesAndAddMateCigar`.
//!
//! Ports `picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar.doWork` at tag 3.4.0, for a single
//! SAM input with default options. Two passes over the reads: restore each read's original base
//! qualities from its `OQ` tag, then group the reads of each template and stamp the mate cigar (`MC`)
//! and mate info onto the pair.
//!
//! The pipeline is `doWork`'s: revert `OQ` while pushing every record into a **queryname** sorting
//! collection, then walk that collection through `SamPairUtil.SetMateInfoIterator(setMateCigar=true)`
//! into an output writer whose sort order is `SORT_ORDER` (unset, so the input's). `doWork` clones
//! the input header, sets only its sort order, and adds no `@PG` and no timestamp, so the whole output
//! is comparable raw.
//!
//! The `OQ` revert is an independent per-record transform, so it runs on all cores and stays
//! byte-identical (decision 0006); the two sorts are stable in-memory sorts (decision 0021) and the
//! mate-info pass is a sequential grouped walk.
//!
//! Scope: on-reference reads with one or two **primary** ends. `createNewCigarsIfMapsOffEndOfReference`
//! (a no-op unless a read hangs off the contig end), the `setMateInformationOnSupplementalAlignment`
//! path for secondary/supplementary records, and the `canSkipSAMFile` shortcut (which suppresses the
//! output entirely when there is nothing to do) are separate surfaces; the port asserts there are no
//! secondary/supplementary records rather than emitting a half-fixed template.

use htsjdk_bam::fastq::fastq_to_phred;
use htsjdk_bam::pair::set_mate_info;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_bam::{coordinate, query_name};
use rayon::prelude::*;

const READ_PAIRED: u16 = 0x1;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;

/// The output sort orders this port writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SortOrder {
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

/// `RESTORE_ORIGINAL_QUALITIES`: move the `OQ` tag back into `QUAL` and drop it.
fn restore_original_qualities(rec: &mut BamRecord) {
    if let Some(TagValue::Str(oq)) = rec.tags.get(Tag::new(b"OQ")) {
        rec.base_qualities = fastq_to_phred(oq);
        rec.tags.remove(Tag::new(b"OQ"));
    }
}

/// `SamPairUtil.SetMateInfoIterator(setMateCigar=true)`: over queryname-sorted records, set mate info
/// and the mate cigar on each template's primary pair.
fn add_mate_info(records: &mut [BamRecord]) {
    let mut start = 0;
    while start < records.len() {
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
                "RevertOriginalBaseQualitiesAndAddMateCigar: secondary/supplementary records are not ported"
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
            let (lo, hi) = (f.min(s), f.max(s));
            let (left, right) = records.split_at_mut(hi);
            let (a, b) = (&mut left[lo], &mut right[0]);
            if f < s {
                set_mate_info(a, b, true);
            } else {
                set_mate_info(b, a, true);
            }
        }

        start = end;
    }
}

/// `RevertOriginalBaseQualitiesAndAddMateCigar.doWork` for a single SAM input, default options.
pub fn revert_original_base_qualities_and_add_mate_cigar(
    input_sam: &str,
) -> Result<String, ParseError> {
    // The tool opens the input EAGERLY_DECODE at whatever VALIDATION_STRINGENCY; stringency does not
    // reach the bytes.
    let (mut header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    // SORT_ORDER defaults to the input header's sort order.
    let output_order = header
        .attributes
        .get("SO")
        .and_then(SortOrder::from_str)
        .unwrap_or(SortOrder::Coordinate);

    // Restore original qualities: independent per record, so parallel (decision 0006).
    records.par_iter_mut().for_each(restore_original_qualities);

    // Queryname sort to group templates, add mate info + mate cigar, then re-sort to the output order.
    records.sort_by(query_name::compare);
    add_mate_info(&mut records);
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

    // A coordinate-sorted proper pair, both mapped on-reference, each with an OQ to restore and no MC.
    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        p1\t99\tchr1\t100\t60\t4M\t=\t300\t204\tACGT\tIIII\tOQ:Z:5555\n\
        p1\t147\tchr1\t300\t50\t4M\t=\t100\t-204\tACGT\tJJJJ\tOQ:Z:AAAA\n";

    fn rows(sam: &str) -> Vec<Vec<&str>> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').collect())
            .collect()
    }

    #[test]
    fn original_qualities_are_restored_and_oq_dropped() {
        let out = revert_original_base_qualities_and_add_mate_cigar(INPUT).unwrap();
        let r = rows(&out);
        let first = r.iter().find(|x| x[3] == "100").unwrap();
        assert_eq!(first[10], "5555"); // QUAL restored from OQ:Z:5555
        assert!(!first.iter().any(|t| t.starts_with("OQ")), "OQ dropped");
        let second = r.iter().find(|x| x[3] == "300").unwrap();
        assert_eq!(second[10], "AAAA");
    }

    #[test]
    fn the_mate_cigar_and_mate_mapping_quality_are_added() {
        let out = revert_original_base_qualities_and_add_mate_cigar(INPUT).unwrap();
        let r = rows(&out);
        let first = r.iter().find(|x| x[3] == "100").unwrap();
        // MC is the mate's cigar; MQ the mate's mapping quality (50). MC sorts before MQ by tag code.
        assert!(first.contains(&"MC:Z:4M"), "got {first:?}");
        assert!(first.contains(&"MQ:i:50"), "got {first:?}");
        let second = r.iter().find(|x| x[3] == "300").unwrap();
        assert!(second.contains(&"MQ:i:60"), "got {second:?}");
    }

    #[test]
    fn the_output_keeps_the_input_coordinate_order() {
        let out = revert_original_base_qualities_and_add_mate_cigar(INPUT).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:coordinate"));
        let names: Vec<&str> = rows(&out).iter().map(|r| r[3]).collect();
        assert_eq!(names, ["100", "300"]);
    }
}
