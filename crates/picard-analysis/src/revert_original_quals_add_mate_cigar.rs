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
const MATE_UNMAPPED: u16 = 0x8;
const READ_UNMAPPED: u16 = 0x4;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const MATE_NEGATIVE_STRAND: u16 = 0x20;

/// Why a run stops: the input did not parse, or the reference threw while it walked the templates.
#[derive(Debug)]
pub enum RevertError {
    Parse(ParseError),
    /// The exception, class and message, as `Exception in thread "main"` prints it.
    Thrown(String),
}

impl From<ParseError> for RevertError {
    fn from(e: ParseError) -> Self {
        RevertError::Parse(e)
    }
}

/// The output sort orders, which are `SAMFileHeader.SortOrder`'s five.
///
/// The writer is built with `presorted = false`, so it SORTS into the order the header names --
/// and two of the five name no comparator at all, which is not "leave it alone by accident" but
/// the documented behaviour of `SortOrder.getComparatorInstance` answering null.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
    Unknown,
}

impl SortOrder {
    pub fn name(self) -> &'static str {
        match self {
            SortOrder::Unsorted => "unsorted",
            SortOrder::Queryname => "queryname",
            SortOrder::Coordinate => "coordinate",
            SortOrder::Duplicate => "duplicate",
            SortOrder::Unknown => "unknown",
        }
    }

    pub fn parse(s: &str) -> Option<SortOrder> {
        match s {
            "unsorted" => Some(SortOrder::Unsorted),
            "queryname" => Some(SortOrder::Queryname),
            "coordinate" => Some(SortOrder::Coordinate),
            "duplicate" => Some(SortOrder::Duplicate),
            "unknown" => Some(SortOrder::Unknown),
            _ => None,
        }
    }
}

/// `RevertOriginalBaseQualitiesAndAddMateCigar`'s own arguments.
#[derive(Debug, Clone)]
pub struct Options {
    pub restore_original_qualities: bool,
    /// `null` means the input header's order, which is what the tool substitutes before it writes.
    pub sort_order: Option<SortOrder>,
    pub max_records_to_examine: i32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            restore_original_qualities: true,
            sort_order: None,
            max_records_to_examine: 10_000,
        }
    }
}

/// `CanSkipSamFile`: what the first records of the input say about whether there is work to do.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanSkip {
    /// A record with no `OQ` whose mate is mapped and whose `MC` is already there: nothing to do,
    /// and the tool returns having written NO output file at all.
    CanSkip,
    FoundOq,
    FoundNoMateCigar,
    FoundNoEvidence,
}

impl CanSkip {
    pub fn skips(self) -> bool {
        self == CanSkip::CanSkip
    }
}

/// `canSkipSAMFile`, which decides on the FIRST record that answers either way.
///
/// The loop stops at the first record carrying an `OQ` (cannot skip), and otherwise at the first
/// paired record whose mate is mapped: that record's `MC` decides it for the whole file. Only
/// records that answer neither count against `MAX_RECORDS_TO_EXAMINE`, so a file whose first
/// record is unpaired and has no `OQ` is examined further while one whose first record is a mapped
/// pair is decided immediately.
pub fn can_skip(records: &[BamRecord], options: &Options) -> CanSkip {
    // Only records that answer neither question count against the limit, and the loop stops at
    // the first that answers either.
    for record in records
        .iter()
        .take(options.max_records_to_examine.max(0) as usize)
    {
        if options.restore_original_qualities && record.tags.get(Tag::new(b"OQ")).is_some() {
            return CanSkip::FoundOq;
        }
        if record.flags & READ_PAIRED != 0 && record.flags & MATE_UNMAPPED == 0 {
            return match record.tags.get(Tag::new(b"MC")) {
                Some(_) => CanSkip::CanSkip,
                None => CanSkip::FoundNoMateCigar,
            };
        }
    }
    CanSkip::FoundNoEvidence
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

/// `SamPairUtil.SetMateInfoIterator(setMateCigar=true, ignoreMissingMates=false)`: over
/// queryname-sorted records, set mate info and the mate cigar on each template's primary pair, and
/// on its supplementary records from the other end's primary.
///
/// Secondary records are carried through untouched. A template whose paired primaries are not
/// exactly one first and one second end is the `SAMException` the iterator throws: a second primary
/// on the same end as soon as it is met, a missing end once the template is read.
fn add_mate_info(records: &mut [BamRecord]) -> Result<(), String> {
    let mut start = 0;
    while start < records.len() {
        let mut end = start + 1;
        while end < records.len() && records[end].read_name == records[start].read_name {
            end += 1;
        }

        let mut first_primary: Option<usize> = None;
        let mut second_primary: Option<usize> = None;
        let mut contains_supplementary = false;
        for (offset, r) in records[start..end].iter().enumerate() {
            let i = start + offset;
            if r.flags & READ_PAIRED == 0 {
                continue;
            }
            if is_primary(r) {
                if r.flags & FIRST_OF_PAIR != 0 {
                    if first_primary.is_some() {
                        return Err(format!(
                            "htsjdk.samtools.SAMException: Found two records that are paired, not \
                             supplementary, and first of the pair: {}",
                            r.read_name
                        ));
                    }
                    first_primary = Some(i);
                } else if r.flags & SECOND_OF_PAIR != 0 {
                    if second_primary.is_some() {
                        return Err(format!(
                            "htsjdk.samtools.SAMException: Found two records that are paired, not \
                             supplementary, and second of the pair: {}",
                            r.read_name
                        ));
                    }
                    second_primary = Some(i);
                }
            }
            if r.flags & SUPPLEMENTARY_ALIGNMENT != 0 {
                contains_supplementary = true;
            }
        }

        match (first_primary, second_primary) {
            (Some(f), Some(s)) => {
                let (lo, hi) = (f.min(s), f.max(s));
                let (left, right) = records.split_at_mut(hi);
                let (a, b) = (&mut left[lo], &mut right[0]);
                if f < s {
                    set_mate_info(a, b, true);
                } else {
                    set_mate_info(b, a, true);
                }
                if contains_supplementary {
                    let first = records[f].clone();
                    let second = records[s].clone();
                    for r in &mut records[start..end] {
                        if r.flags & READ_PAIRED != 0 && r.flags & SUPPLEMENTARY_ALIGNMENT != 0 {
                            let mate = if r.flags & FIRST_OF_PAIR != 0 {
                                &second
                            } else {
                                &first
                            };
                            set_mate_info_on_supplementary(r, mate);
                        }
                    }
                }
            }
            (Some(f), None) => {
                return Err(format!(
                    "htsjdk.samtools.SAMException: Missing second read of pair: {}",
                    records[f].read_name
                ));
            }
            (None, Some(s)) => {
                return Err(format!(
                    "htsjdk.samtools.SAMException: Missing first read of pair: {}",
                    records[s].read_name
                ));
            }
            (None, None) => {}
        }

        start = end;
    }
    Ok(())
}

/// `SamPairUtil.setMateInformationOnSupplementalAlignment(supplemental, matePrimary, true)`, in its
/// setter order: the tags are `MC` then `MQ`, which is the order a record that had neither gains them.
fn set_mate_info_on_supplementary(supplemental: &mut BamRecord, mate: &BamRecord) {
    supplemental.mate_reference_index = mate.reference_index;
    supplemental.mate_alignment_start = mate.alignment_start;
    set_flag(
        supplemental,
        MATE_NEGATIVE_STRAND,
        mate.flags & READ_NEGATIVE_STRAND != 0,
    );
    let mate_unmapped = mate.flags & READ_UNMAPPED != 0;
    set_flag(supplemental, MATE_UNMAPPED, mate_unmapped);
    supplemental.inferred_insert_size = mate.inferred_insert_size.wrapping_neg();
    if mate_unmapped {
        supplemental.tags.remove(Tag::new(b"MC"));
    } else {
        supplemental
            .tags
            .insert(Tag::new(b"MC"), TagValue::Str(mate.cigar.to_text()));
    }
    supplemental.tags.insert(
        Tag::new(b"MQ"),
        TagValue::Int(i64::from(mate.mapping_quality)),
    );
}

fn set_flag(rec: &mut BamRecord, bit: u16, value: bool) {
    if value {
        rec.flags |= bit;
    } else {
        rec.flags &= !bit;
    }
}

/// `RevertOriginalBaseQualitiesAndAddMateCigar.doWork` for a single SAM input, default options.
pub fn revert_original_base_qualities_and_add_mate_cigar(
    input_sam: &str,
) -> Result<String, RevertError> {
    let (header, records) = revert_original_records(input_sam)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// The transform up to the write: the header (with the output sort order) and the reverted,
/// mate-info'd, sorted records. Shared by the SAM and BAM renderers so they cannot drift.
fn revert_original_records(
    input_sam: &str,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), RevertError> {
    revert_original_records_with(input_sam, &Options::default())
}

/// The transform with the tool's own arguments.
///
/// `SORT_ORDER` is the header's when the caller names none, and it decides TWO things: the `SO`
/// the output header carries, and how the records are ordered -- the writer is built with
/// `presorted = false`, so it sorts. `unsorted` and `unknown` name no comparator, and the records
/// then come out in the order the mate-info pass left them, which is query name.
pub fn revert_original_records_with(
    input_sam: &str,
    options: &Options,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), RevertError> {
    // The tool opens the input EAGERLY_DECODE at whatever VALIDATION_STRINGENCY; stringency does not
    // reach the bytes.
    let (mut header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    let output_order = options.sort_order.unwrap_or_else(|| {
        header
            .attributes
            .get("SO")
            .and_then(SortOrder::parse)
            .unwrap_or(SortOrder::Unsorted)
    });

    // Restore original qualities: independent per record, so parallel (decision 0006).
    if options.restore_original_qualities {
        records.par_iter_mut().for_each(restore_original_qualities);
    }

    // Queryname sort to group templates, add mate info + mate cigar, then re-sort to the output order.
    records.sort_by(query_name::compare);
    add_mate_info(&mut records).map_err(RevertError::Thrown)?;
    match output_order {
        SortOrder::Coordinate => records.sort_by(coordinate::compare),
        // `SO:duplicate`'s comparator is htsjdk's `SAMRecordDuplicateComparator`, the same one the
        // duplicate-set iterator orders a set by.
        SortOrder::Duplicate => sort_by_duplicate_order(&header, &mut records),
        // Already query-name sorted, and the two orders that name no comparator leave it there.
        SortOrder::Queryname | SortOrder::Unsorted | SortOrder::Unknown => {}
    }

    header.set_sort_order(output_order.name());
    Ok((header, records))
}

/// `SO:duplicate`: htsjdk's duplicate comparator over the records, which needs each record's
/// library and read group as well as its own fields.
fn sort_by_duplicate_order(header: &htsjdk_bam::header::SamHeader, records: &mut [BamRecord]) {
    use crate::mark_duplicates::{Record, ScoringStrategy};
    let library_of = |record: &BamRecord| -> (String, i32) {
        let id = match record.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(value)) => Some(value.clone()),
            _ => None,
        };
        match id.and_then(|id| {
            header
                .read_groups
                .iter()
                .position(|group| group.id == id)
                .map(|position| (position, &header.read_groups[position]))
        }) {
            Some((position, group)) => (
                group
                    .attributes
                    .get("LB")
                    .unwrap_or("Unknown Library")
                    .to_string(),
                position as i32,
            ),
            None => ("Unknown Library".to_string(), -1),
        }
    };
    let as_record = |record: &BamRecord| -> Record {
        let (library, read_group) = library_of(record);
        Record {
            name: record.read_name.clone(),
            flags: record.flags,
            reference_index: record.reference_index,
            alignment_start: record.alignment_start,
            cigar: record.cigar.clone(),
            qualities: record.base_qualities.clone(),
            mate_reference_index: record.mate_reference_index,
            library,
            read_group,
            barcode: None,
            existing_dt: None,
            mate_cigar: match record.tags.get(Tag::new(b"MC")) {
                Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
                _ => None,
            },
            mate_alignment_start: record.mate_alignment_start,
        }
    };
    let mut libraries: Vec<String> = Vec::new();
    let decorated: Vec<(Record, i32)> = records
        .iter()
        .map(|record| {
            let converted = as_record(record);
            let library = match libraries.iter().position(|k| *k == converted.library) {
                Some(at) => at as i32,
                None => {
                    libraries.push(converted.library.clone());
                    (libraries.len() - 1) as i32
                }
            };
            (converted, library)
        })
        .collect();
    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_by(|a, b| {
        crate::duplicate_set::compare(
            &decorated[*a].0,
            &decorated[*b].0,
            decorated[*a].1,
            decorated[*b].1,
            ScoringStrategy::SumOfBaseQualities,
        )
    });
    let sorted: Vec<BamRecord> = order.iter().map(|index| records[*index].clone()).collect();
    records.clone_from_slice(&sorted);
}

/// The same transform for **BAM** output, byte-identical to Picard with `USE_JDK_DEFLATER=true` via
/// `BamWriter`. The tool adds no `@PG`, so byte-identity follows transitively (the records are those
/// the SAM path already reproduces, and `BamWriter` is oracle-gated over arbitrary records).
pub fn revert_original_base_qualities_and_add_mate_cigar_to_bam(
    input_sam: &str,
) -> Result<Vec<u8>, RevertError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, records) = revert_original_records(input_sam)?;
    let mut w = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        w.write(rec).expect("record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
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
    fn the_bam_output_round_trips_to_the_sam_output() {
        use htsjdk_bam::reader::BamReader;
        let sam = revert_original_base_qualities_and_add_mate_cigar(INPUT).unwrap();
        let bam = revert_original_base_qualities_and_add_mate_cigar_to_bam(INPUT).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bam).unwrap();
        let reader = BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(write_sam(&header, &records).unwrap(), sam);
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
