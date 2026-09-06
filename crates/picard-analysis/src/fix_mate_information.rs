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
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_bam::{coordinate, query_name};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const READ_REVERSE_STRAND: u16 = 0x10;
const MATE_UNMAPPED: u16 = 0x8;
const MATE_REVERSE_STRAND: u16 = 0x20;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;

/// `FixMateInformation`'s arguments, defaulting to Picard's defaults.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `SORT_ORDER`. `None` is Picard's default and means "the input header's own order", which
    /// is not the same as `coordinate`.
    pub sort_order: Option<SortOrder>,
    /// `ASSUME_SORTED`: take the input as queryname-sorted whatever its header says, which skips
    /// the re-sort. On a coordinate-sorted file that means templates are NOT grouped, so only
    /// records that happen to be adjacent are fixed together -- and that is the reference's
    /// behaviour, not a degradation of it.
    pub assume_sorted: bool,
    /// `ADD_MATE_CIGAR`: write `MC`, or clear it.
    pub add_mate_cigar: bool,
    /// `IGNORE_MISSING_MATES`: pass a lone end through, or refuse.
    pub ignore_missing_mates: bool,
    /// `CREATE_INDEX`, which is refused unless the output order is coordinate.
    pub create_index: bool,
}

impl Default for Options {
    /// Picard's defaults, which are not `bool::default()` twice: `ADD_MATE_CIGAR` and
    /// `IGNORE_MISSING_MATES` are true.
    fn default() -> Options {
        Options {
            sort_order: None,
            assume_sorted: false,
            add_mate_cigar: true,
            ignore_missing_mates: true,
            create_index: false,
        }
    }
}

/// What `FixMateInformation` refuses, with the reference's own messages.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FixMateError {
    /// The input did not parse.
    Parse(String),
    /// `SAMException("Missing second read of pair: " + name)`, and its first-read twin. The
    /// reference names the end it FOUND, so a template with only a second read reports "Missing
    /// first read of pair".
    MissingMate { first: bool, name: String },
    /// `PicardException("Can't CREATE_INDEX unless sort order is coordinate")`.
    IndexWithoutCoordinate,
}

impl FixMateError {
    /// The message the reference prints, without the `Exception in thread "main"` prefix.
    pub fn message(&self) -> String {
        match self {
            FixMateError::Parse(detail) => detail.clone(),
            FixMateError::MissingMate { first, name } => {
                let which = if *first { "first" } else { "second" };
                format!("Missing {which} read of pair: {name}")
            }
            FixMateError::IndexWithoutCoordinate => {
                "Can't CREATE_INDEX unless sort order is coordinate".to_string()
            }
        }
    }

    /// The exception class, which decides how the reference's own handler prints it.
    pub fn java_class(&self) -> &'static str {
        match self {
            FixMateError::Parse(_) => "htsjdk.samtools.SAMFormatException",
            FixMateError::MissingMate { .. } => "htsjdk.samtools.SAMException",
            FixMateError::IndexWithoutCoordinate => "picard.PicardException",
        }
    }
}

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
fn fix_templates(
    records: &mut [BamRecord],
    set_mate_cigar: bool,
    ignore_missing_mates: bool,
) -> Result<(), FixMateError> {
    let mut start = 0;
    while start < records.len() {
        // The run of records sharing this read name.
        let mut end = start + 1;
        while end < records.len() && records[end].read_name == records[start].read_name {
            end += 1;
        }

        let mut first_primary: Option<usize> = None;
        let mut second_primary: Option<usize> = None;
        let mut supplementals: Vec<usize> = Vec::new();
        for (offset, r) in records[start..end].iter().enumerate() {
            let i = start + offset;
            if r.flags & READ_PAIRED != 0 {
                // `isSecondaryOrSupplementary()`: only a PRIMARY record can be the end whose mate
                // info the others are set from. A secondary alignment is neither that nor a
                // record the reference touches afterwards; a supplementary one is touched, from
                // the OTHER end's primary, which is what `supplementals` collects.
                if r.flags & SUPPLEMENTARY_ALIGNMENT != 0 {
                    supplementals.push(i);
                }
                if is_primary(r) {
                    if r.flags & FIRST_OF_PAIR != 0 {
                        assert!(first_primary.is_none(), "two first-of-pair primaries");
                        first_primary = Some(i);
                    } else if r.flags & SECOND_OF_PAIR != 0 {
                        assert!(second_primary.is_none(), "two second-of-pair primaries");
                        second_primary = Some(i);
                    }
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

            // `setMateInformationOnSupplementalAlignment`, from the OTHER end's primary: a
            // first-of-pair supplementary takes the second primary's, and the reverse. The
            // reference does this only when the template has both primaries, which is why it sits
            // inside this branch.
            for index in &supplementals {
                let from = if records[*index].flags & FIRST_OF_PAIR != 0 {
                    s
                } else {
                    f
                };
                let mate = records[from].clone();
                set_mate_info_on_supplemental(&mut records[*index], &mate, set_mate_cigar);
            }
        } else if !ignore_missing_mates {
            // The reference names the end it FOUND: a template with only a second read reports
            // "Missing first read of pair". A fragment read is not a missing mate at all, which is
            // why the paired flag is tested again here.
            if let Some(f) = first_primary {
                if records[f].flags & READ_PAIRED != 0 {
                    return Err(FixMateError::MissingMate {
                        first: false,
                        name: records[f].read_name.clone(),
                    });
                }
            } else if let Some(s) = second_primary {
                if records[s].flags & READ_PAIRED != 0 {
                    return Err(FixMateError::MissingMate {
                        first: true,
                        name: records[s].read_name.clone(),
                    });
                }
            }
        }
        // A lone paired read (missing mate) would throw unless IGNORE_MISSING_MATES; the corpora do
        // not exercise that, and a fragment (unpaired) read simply passes through.

        start = end;
    }
    Ok(())
}

/// `SAMRecord.setMateNegativeStrandFlag` and friends: one bit, set or cleared.
fn set_flag(record: &mut BamRecord, bit: u16, value: bool) {
    if value {
        record.flags |= bit;
    } else {
        record.flags &= !bit;
    }
}

/// `SamPairUtil.setMateInformationOnSupplementalAlignment`.
///
/// Six fields and two tags, and the insert size is the mate's NEGATED: a supplementary alignment
/// describes the same fragment from the other side. `MC` is set when the caller asked for a mate
/// cigar AND the mate is mapped, and CLEARED otherwise -- the else branch writes null rather than
/// leaving whatever was there, which is the difference between a record that carries a stale mate
/// cigar and one that carries none.
fn set_mate_info_on_supplemental(
    supplemental: &mut BamRecord,
    mate: &BamRecord,
    set_mate_cigar: bool,
) {
    supplemental.mate_reference_index = mate.reference_index;
    supplemental.mate_alignment_start = mate.alignment_start;
    set_flag(
        supplemental,
        MATE_REVERSE_STRAND,
        mate.flags & READ_REVERSE_STRAND != 0,
    );
    set_flag(supplemental, MATE_UNMAPPED, mate.flags & READ_UNMAPPED != 0);
    supplemental.inferred_insert_size = -mate.inferred_insert_size;
    if set_mate_cigar && mate.flags & READ_UNMAPPED == 0 {
        supplemental
            .tags
            .insert(Tag::new(b"MC"), TagValue::Str(mate.cigar.to_text()));
    } else {
        supplemental.tags.remove(Tag::new(b"MC"));
    }
    supplemental
        .tags
        .insert(Tag::new(b"MQ"), TagValue::Int(mate.mapping_quality as i64));
}

/// `FixMateInformation.doWork` for a single SAM input, with Picard's defaults.
pub fn fix_mate_information(
    input_sam: &str,
    sort_order: Option<SortOrder>,
) -> Result<String, ParseError> {
    let options = Options {
        sort_order,
        ..Options::default()
    };
    // The old signature answered with a `ParseError`, and its callers only ever met the parse
    // one: with Picard's defaults the two refusals below cannot be raised at all.
    fix_mate_information_with(input_sam, &options).map_err(|error| match error {
        FixMateError::Parse(_) => read_sam_with(input_sam, ValidationStringency::Lenient)
            .expect_err("the parse failed once already"),
        other => panic!("unreachable with Picard's defaults: {}", other.message()),
    })
}

/// `doWork` with the arguments the tool declares.
pub fn fix_mate_information_with(
    input_sam: &str,
    options: &Options,
) -> Result<String, FixMateError> {
    let (header, records) = fix_mates_with(input_sam, options)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// `FixMateInformation.doWork` up to the write: the header (with its final sort order) and the
/// mate-fixed, re-sorted records.
fn fix_mates_with(
    input_sam: &str,
    options: &Options,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), FixMateError> {
    let (mut header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)
        .map_err(|error| FixMateError::Parse(format!("{error:?}")))?;

    // The output order defaults to the input header's sort order.
    let output_order = options.sort_order.unwrap_or_else(|| {
        header
            .attributes
            .get("SO")
            .and_then(SortOrder::from_str)
            .unwrap_or(SortOrder::Coordinate)
    });

    // `if (CREATE_INDEX && header.getSortOrder() != coordinate) throw`, and the header's order at
    // that point is the OUTPUT order rather than the input's.
    if options.create_index && output_order != SortOrder::Coordinate {
        return Err(FixMateError::IndexWithoutCoordinate);
    }

    // `if (ASSUME_SORTED || allQueryNameSorted)` the input is taken as already grouped and the
    // re-sort is skipped. On a coordinate-sorted file that changes the answer rather than only the
    // speed: templates are not adjacent, so a pair whose ends sit apart is never fixed together.
    let already_queryname = header.attributes.get("SO") == Some("queryname");
    if !(options.assume_sorted || already_queryname) {
        records.sort_by(query_name::compare);
    }
    fix_templates(
        &mut records,
        options.add_mate_cigar,
        options.ignore_missing_mates,
    )?;
    match output_order {
        SortOrder::Coordinate => records.sort_by(coordinate::compare),
        SortOrder::Queryname => records.sort_by(query_name::compare),
    }

    header.set_sort_order(output_order.name());
    Ok((header, records))
}

/// `FixMateInformation.doWork` for SAM input and **BAM** output. The same fix, written through
/// htsjdk-rs's byte-identical `BamWriter`; FixMateInformation adds no `@PG`. Byte-identity to Picard's
/// `USE_JDK_DEFLATER=true` output follows transitively: the records are the ones
/// `fix_mate_information` already reproduces (its oracle), and the `BamWriter` is proven byte-identical
/// over arbitrary records (the SamFormatConverter oracle and htsjdk-rs's
/// `every_file_is_byte_identical_to_htsjdks`).
pub fn fix_mate_information_to_bam(
    input_sam: &str,
    sort_order: Option<SortOrder>,
) -> Result<Vec<u8>, ParseError> {
    let options = Options {
        sort_order,
        ..Options::default()
    };
    fix_mate_information_to_bam_with(input_sam, &options).map_err(|error| match error {
        FixMateError::Parse(_) => read_sam_with(input_sam, ValidationStringency::Lenient)
            .expect_err("the parse failed once already"),
        other => panic!("unreachable with Picard's defaults: {}", other.message()),
    })
}

/// The same, with the arguments the tool declares.
pub fn fix_mate_information_to_bam_with(
    input_sam: &str,
    options: &Options,
) -> Result<Vec<u8>, FixMateError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, records) = fix_mates_with(input_sam, options)?;
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

    /// The BAM output decodes back to exactly the SAM output; the writer's byte-identity to htsjdk is
    /// proven elsewhere.
    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        let sam = fix_mate_information(INPUT, None).unwrap();
        let bam = fix_mate_information_to_bam(INPUT, None).unwrap();
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
