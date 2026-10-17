//! `FilterSamReads`, the read-list filters.
//!
//! Ports `picard.sam.FilterSamReads` driving `htsjdk.samtools.filter.ReadNameFilter` at Picard 3.4.0,
//! for `FILTER=includeReadList` and `FILTER=excludeReadList`: keep (or drop) the reads whose name
//! appears in a `READ_LIST_FILE`. `ReadNameFilter.filterOut` is exact set membership on the whole
//! read name (`readNameFilterSet.contains(name) != includeReads`), so a read is kept when
//! `contains(name) == include`.
//!
//! With `SORT_ORDER` unset the output keeps the input's sort order, so the writer is `presorted` and
//! there is **no re-sort**; `doWork` adds no `@PG` and no timestamp, so the whole output is comparable
//! raw. The membership test is a pure per-record predicate reading the shared, immutable name set, so
//! the filter runs on all cores with rayon's ordered `collect` preserving input order (decision 0006).
//!
//! The aligned, tag-value and paired-interval filters and `WRITE_READS_FILES` live here too; the
//! binary drives them and the JavaScript filter.

use std::collections::HashSet;

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use rayon::prelude::*;

/// Whether the `READ_LIST_FILE` names are the reads to keep or the reads to drop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    /// `includeReadList`: keep the reads whose name is in the list.
    IncludeReadList,
    /// `excludeReadList`: keep the reads whose name is **not** in the list.
    ExcludeReadList,
}

impl Filter {
    /// `ReadNameFilter`'s `includeReads` flag.
    fn include(self) -> bool {
        matches!(self, Filter::IncludeReadList)
    }
}

/// Parses a `READ_LIST_FILE`: one read name per line, blank lines ignored.
pub fn read_name_set(read_list_text: &str) -> HashSet<&str> {
    read_list_text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect()
}

/// `ReadNameFilter.filterOut` over records already in memory: a read is kept when its membership
/// in the set matches the filter's `includeReads` flag.
///
/// The text entry points below parse and re-encode SAM; a caller holding records (a BAM reader,
/// say) needs the predicate and not the parsing.
pub fn keep_by_read_list(
    records: &[BamRecord],
    names: &HashSet<&str>,
    filter: Filter,
) -> Vec<BamRecord> {
    let include = filter.include();
    records
        .par_iter()
        .filter(|rec| names.contains(rec.read_name.as_str()) == include)
        .cloned()
        .collect()
}

/// `FilteringSamIterator` with `filterReadPairs` and an `AlignedFilter`, record by record as
/// htsjdk runs it: a first-of-pair is judged with the record after it, which must be its mate
/// (`SamPairUtil.assertMate`), and a second-of-pair is passed through on the assumption that its
/// first was. The error is the assertion's, for a file whose mates do not follow each other.
pub fn keep_by_alignment_pairwise(
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    filter: AlignedFilter,
) -> Result<Vec<BamRecord>, String> {
    let (kept, error) = alignment_pairwise_indices(header, records, filter);
    match error {
        Some((_, message)) => Err(message),
        None => Ok(kept.into_iter().map(|i| records[i].clone()).collect()),
    }
}

/// [`keep_by_alignment_pairwise`] as indices into `records`, with the error, if any, beside the
/// index of the record being judged when it was thrown and the records kept before it.
pub fn alignment_pairwise_indices(
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    filter: AlignedFilter,
) -> (Vec<usize>, Option<(usize, String)>) {
    const PAIRED: u16 = 0x1;
    const FIRST: u16 = 0x40;
    const SECOND: u16 = 0x80;
    let include = matches!(filter, AlignedFilter::IncludeAligned);
    let unmapped = |rec: &BamRecord| !is_aligned(rec);
    let mut kept = Vec::new();
    let mut index = 0;
    while index < records.len() {
        let at = index;
        let record = &records[index];
        index += 1;
        let paired = record.flags & PAIRED != 0;
        if paired && record.flags & FIRST != 0 && index < records.len() {
            let mate = &records[index];
            if mate.flags & PAIRED == 0 {
                return (
                    kept,
                    Some((
                        at,
                        format!(
                    "htsjdk.samtools.SAMException: Second record is not marked as paired: {}",
                    sam_line(header, mate)
                ),
                    )),
                );
            }
            if mate.flags & SECOND == 0 {
                return (
                    kept,
                    Some((
                        at,
                        format!(
                    "htsjdk.samtools.SAMException: Second record is not marked as second of \
                     pair: {}",
                    sam_line(header, mate)
                ),
                    )),
                );
            }
            if mate.read_name != record.read_name {
                return (
                    kept,
                    Some((
                        at,
                        format!(
                    "htsjdk.samtools.SAMException: First [{}] and Second [{}] readnames do not \
                     match!",
                    record.read_name, mate.read_name
                ),
                    )),
                );
            }
            let filtered_out = if include {
                unmapped(record) || unmapped(mate)
            } else {
                !(unmapped(record) || unmapped(mate))
            };
            if filtered_out {
                index += 1;
            } else {
                kept.push(at);
            }
        } else if (paired && record.flags & SECOND != 0) || unmapped(record) != include {
            // A second-of-pair passes on the assumption that its first did; anything else is
            // judged alone.
            kept.push(at);
        }
    }
    (kept, None)
}

/// `IntervalList.fromFile(...).getIntervals()` as `IntervalKeepPairFilter` holds them: the
/// intervals per contig, as an `OverlapDetector` with no buffer keeps them. An interval with no
/// bases (`start == end + 1`) is never stored, as `addLhs` drops it.
pub struct PairedIntervals {
    by_contig: std::collections::HashMap<String, Vec<(i32, i32)>>,
}

impl PairedIntervals {
    /// `IntervalList.fromReader` and `IntervalListCodec.decode`: the header is required, a line
    /// has exactly five tab-separated fields (Java's `split` drops trailing empty ones), and an
    /// interval on a contig the header does not name is ignored. The error is the exception's
    /// message, which `doWork` catches.
    pub fn parse(text: &str) -> Result<PairedIntervals, String> {
        let mut lengths: std::collections::HashMap<&str, i32> = std::collections::HashMap::new();
        let mut lines = text.lines().peekable();
        let mut header = false;
        while let Some(line) = lines.peek() {
            if !line.starts_with('@') {
                break;
            }
            header = true;
            if line.starts_with("@SQ") {
                let field = |key: &str| line.split('\t').find_map(|f| f.strip_prefix(key));
                if let Some(name) = field("SN:") {
                    let length = field("LN:").and_then(|l| l.parse().ok()).unwrap_or(0);
                    lengths.insert(name, length);
                }
            }
            lines.next();
        }
        if !header {
            return Err(
                "java.lang.IllegalStateException: Interval list file must contain header. "
                    .to_string(),
            );
        }
        let mut by_contig: std::collections::HashMap<String, Vec<(i32, i32)>> =
            std::collections::HashMap::new();
        for line in lines {
            if line.starts_with('@') || line.trim().is_empty() {
                continue;
            }
            let mut fields: Vec<&str> = line.split('\t').collect();
            while fields.last() == Some(&"") {
                fields.pop();
            }
            if fields.len() != 5 {
                return Err(format!(
                    "htsjdk.tribble.TribbleException: Invalid interval record contains {} fields: {line}",
                    fields.len()
                ));
            }
            let number = |s: &str| {
                s.parse::<i32>().map_err(|_| {
                    format!("java.lang.NumberFormatException: For input string: \"{s}\"")
                })
            };
            let (start, end) = (number(fields[1])?, number(fields[2])?);
            if start < 1 {
                return Err(format!(
                    "java.lang.IllegalArgumentException: Coordinate less than 1: start value of {start} is less than 1 and thus illegal"
                ));
            }
            if start as i64 > end as i64 + 1 {
                return Err(format!(
                    "java.lang.IllegalArgumentException: Start value of {start} is greater than end + 1 for end of value: {end}. \
                     I'm afraid I cannot let you do that."
                ));
            }
            if !matches!(fields[3], "+" | "-") {
                return Err(format!(
                    "java.lang.IllegalArgumentException: Invalid strand field: {}",
                    fields[3]
                ));
            }
            let Some(&length) = lengths.get(fields[0]) else {
                continue;
            };
            if length > 0 && length < end {
                return Err(format!(
                    "java.lang.IllegalArgumentException: interval with end: {end} extends beyond end of sequence with length: {length}"
                ));
            }
            if start <= end {
                by_contig
                    .entry(fields[0].to_string())
                    .or_default()
                    .push((start, end));
            }
        }
        Ok(PairedIntervals { by_contig })
    }

    /// `hasOverlaps`: whether any interval shares a base with `contig:start-end`.
    fn overlaps(&self, contig: &str, start: i32, end: i32) -> bool {
        start <= end
            && self
                .by_contig
                .get(contig)
                .is_some_and(|spans| spans.iter().any(|&(s, e)| s <= end && start <= e))
    }
}

/// `IntervalKeepPairFilter` driven record by record (`FilterSamReads` builds its
/// `FilteringSamIterator` without `filterReadPairs`): a paired, primary read is kept when it or its
/// mate overlaps an interval, the mate's extent read from its `MC` tag. A read whose mate is
/// mapped but carries no `MC` is `SAMUtils.getMateAlignmentEnd`'s exception.
pub fn keep_by_paired_intervals(
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    intervals: &PairedIntervals,
) -> Result<Vec<BamRecord>, String> {
    let (kept, error) = paired_interval_indices(header, records, intervals);
    match error {
        Some((_, message)) => Err(message),
        None => Ok(kept.into_iter().map(|i| records[i].clone()).collect()),
    }
}

/// [`keep_by_paired_intervals`] as indices into `records`, with the error, if any, beside the
/// index of the record that threw it.
pub fn paired_interval_indices(
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    intervals: &PairedIntervals,
) -> (Vec<usize>, Option<(usize, String)>) {
    let sequence_names: Vec<&str> = header.sequences.iter().map(|s| s.name.as_str()).collect();
    const PAIRED: u16 = 0x1;
    const MATE_UNMAPPED: u16 = 0x8;
    const SECONDARY_OR_SUPPLEMENTARY: u16 = 0x100 | 0x800;
    let name = |index: i32| {
        usize::try_from(index)
            .ok()
            .and_then(|i| sequence_names.get(i).copied())
            .unwrap_or("*")
    };
    let mut kept = Vec::new();
    for (at, record) in records.iter().enumerate() {
        if record.flags & SECONDARY_OR_SUPPLEMENTARY != 0 || record.flags & PAIRED == 0 {
            continue;
        }
        let read_overlaps = is_aligned(record)
            && intervals.overlaps(
                name(record.reference_index),
                record.alignment_start,
                record.alignment_end(),
            );
        let keep = read_overlaps
            || (record.flags & MATE_UNMAPPED == 0 && {
                let Some(length) = crate::insert_size::mate_cigar_reference_length(record) else {
                    let message = format!(
                        "htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found:{}",
                        sam_line(header, record)
                    );
                    return (kept, Some((at, message)));
                };
                intervals.overlaps(
                    name(record.mate_reference_index),
                    record.mate_alignment_start,
                    record.mate_alignment_start + length - 1,
                )
            });
        if keep {
            kept.push(at);
        }
    }
    (kept, None)
}

/// `SAMRecord.toString()` in htsjdk 5.0.0: the record's SAM line, without its newline.
pub fn sam_line(header: &htsjdk_bam::header::SamHeader, record: &BamRecord) -> String {
    let mut text = reads_file_text(header, std::slice::from_ref(record));
    text.truncate(text.trim_end_matches('\n').len());
    text
}

/// A `WRITE_READS_FILES` `.reads` file: `rec.toString()` and a newline per record. htsjdk 5.0.0,
/// which Picard 3.5.0 ships, renders a record's `toString()` as its SAM line (4.2.0 printed a
/// one-line summary instead).
pub fn reads_file_text(header: &htsjdk_bam::header::SamHeader, records: &[BamRecord]) -> String {
    let header_only = write_sam(header, &[]).unwrap_or_default();
    let mut text = write_sam(header, records).unwrap_or_default();
    text.drain(..header_only.len().min(text.len()));
    text
}

/// `FilterSamReads` with a read-list filter, for SAM input and output. `read_list_text` is the
/// `READ_LIST_FILE` contents.
pub fn filter_sam_reads(
    input_sam: &str,
    read_list_text: &str,
    filter: Filter,
) -> Result<String, ParseError> {
    let (header, kept) = filter_read_list(input_sam, read_list_text, filter)?;
    Ok(write_sam(&header, &kept).expect("records that parsed re-encode as SAM text"))
}

/// The read-list filter's work up to the write: the header and the kept records.
fn filter_read_list(
    input_sam: &str,
    read_list_text: &str,
    filter: Filter,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    let names = read_name_set(read_list_text);
    let include = filter.include();

    // Keep a read iff its membership in the set matches the include flag. Independent per record, so
    // the filter is parallel and order-preserving.
    let kept: Vec<BamRecord> = records
        .par_iter()
        .filter(|rec| names.contains(rec.read_name.as_str()) == include)
        .cloned()
        .collect();

    Ok((header, kept))
}

/// `FilterSamReads` with a read-list filter, for SAM input and **BAM** output. Same filter, written
/// through htsjdk-rs's byte-identical `BamWriter`; FilterSamReads adds no `@PG`. Byte-identity to
/// Picard's `USE_JDK_DEFLATER=true` output follows transitively: the kept records are the ones
/// `filter_sam_reads` already reproduces (its oracle), and the `BamWriter` is proven byte-identical
/// over arbitrary records (the SamFormatConverter oracle and htsjdk-rs's
/// `every_file_is_byte_identical_to_htsjdks`).
pub fn filter_sam_reads_to_bam(
    input_sam: &str,
    read_list_text: &str,
    filter: Filter,
) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, kept) = filter_read_list(input_sam, read_list_text, filter)?;
    let mut writer = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &kept {
        writer
            .write(rec)
            .expect("records that parsed re-encode as BAM");
    }
    Ok(writer
        .finish()
        .expect("in-memory BAM writer never fails to finish"))
}

/// Whether the reads matching `TAG_VALUE` are the reads to keep or the reads to drop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TagFilter {
    /// `includeTagValues`: keep the reads whose `TAG` value is one of the listed values.
    IncludeTagValues,
    /// `excludeTagValues`: keep the reads whose `TAG` value is **not** one of the listed values.
    ExcludeTagValues,
}

impl TagFilter {
    /// `TagFilter`'s `includeReads` flag.
    fn include(self) -> bool {
        matches!(self, TagFilter::IncludeTagValues)
    }
}

/// `FilterSamReads` with `FILTER=includeTagValues` / `excludeTagValues`, ported from
/// `htsjdk.samtools.filter.TagFilter`: a read matches when its `tag` attribute equals one of the
/// `values`. `TagFilter.filterOut` is `values.contains(record.getAttribute(tag)) != includeReads`,
/// so a read is kept when `matches == include`.
///
/// `getAttribute` returns the typed value, and the `TAG_VALUE` command-line arguments are strings, so
/// a value can only match a **string** (`Z`) tag; an integer tag or an absent tag never matches, just
/// as an `Integer` (or `null`) is never in a `List<String>`. Like the read-list filter this is a pure
/// per-record predicate, so it runs on all cores in input order (decision 0006), adds no `@PG`, and
/// leaves the sort order untouched.
pub fn filter_sam_reads_by_tag(
    input_sam: &str,
    tag: &[u8; 2],
    values: &[&str],
    filter: TagFilter,
) -> Result<String, ParseError> {
    let (header, kept) = filter_by_tag(input_sam, tag, values, filter)?;
    Ok(write_sam(&header, &kept).expect("records that parsed re-encode as SAM text"))
}

/// The tag-value filter's work up to the write: the header and the kept records.
fn filter_by_tag(
    input_sam: &str,
    tag: &[u8; 2],
    values: &[&str],
    filter: TagFilter,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    let value_set: HashSet<&str> = values.iter().copied().collect();
    let tag = Tag::new(tag);
    let include = filter.include();

    let matches = |rec: &BamRecord| match rec.tags.get(tag) {
        Some(TagValue::Str(s)) => value_set.contains(s.as_str()),
        _ => false, // absent, or a non-string tag, never matches a string value
    };

    let kept: Vec<BamRecord> = records
        .par_iter()
        .filter(|rec| matches(rec) == include)
        .cloned()
        .collect();
    Ok((header, kept))
}

/// The tag-value filter for **BAM** output, byte-identical to Picard's `USE_JDK_DEFLATER=true` via
/// `BamWriter` (transitive, as for the read-list `_to_bam`).
pub fn filter_sam_reads_by_tag_to_bam(
    input_sam: &str,
    tag: &[u8; 2],
    values: &[&str],
    filter: TagFilter,
) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, kept) = filter_by_tag(input_sam, tag, values, filter)?;
    let mut w = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &kept {
        w.write(rec).expect("record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

const READ_UNMAPPED: u16 = 0x4;

/// Whether `includeAligned` (keep aligned templates) or `excludeAligned` (keep templates with an
/// unmapped read).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedFilter {
    /// `includeAligned`: keep a template only when **all** its reads are aligned.
    IncludeAligned,
    /// `excludeAligned`: keep a template when **any** of its reads is unmapped.
    ExcludeAligned,
}

fn is_aligned(rec: &BamRecord) -> bool {
    rec.flags & READ_UNMAPPED == 0
}

/// `FilterSamReads` with `FILTER=includeAligned` / `excludeAligned`, ported from
/// `htsjdk.samtools.filter.AlignedFilter` driven pairwise by `FilteringSamIterator` (constructed with
/// `filterByPair = true`).
///
/// The input must be **queryname-sorted**, so the reads of a template are consecutive. A template is
/// kept or dropped as a unit: for `includeAligned`, a two-read template is kept only when both reads
/// are aligned and a lone read only when it is aligned; for `excludeAligned`, a two-read template is
/// kept when either read is unmapped and a lone read only when it is unmapped. The kept reads pass
/// through in input order; no `@PG` is added and the sort order is untouched, so the output is
/// comparable raw.
///
/// This is a **sequential** pass (the pairing depends on encounter order), so unlike the read-list
/// and tag-value filters it is not parallelized. A template of more than two reads goes through
/// [`keep_by_alignment_pairwise`] as htsjdk drives it: a first-of-pair followed by anything but its
/// mate is `SamPairUtil.assertMate`'s error, returned as [`ParseError::Inconsistent`].
pub fn filter_sam_reads_aligned(
    input_sam: &str,
    filter: AlignedFilter,
) -> Result<String, ParseError> {
    let (header, kept) = filter_aligned(input_sam, filter)?;
    Ok(write_sam(&header, &kept).expect("records that parsed re-encode as SAM text"))
}

/// The aligned filter's work up to the write: the header and the kept records.
fn filter_aligned(
    input_sam: &str,
    filter: AlignedFilter,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    let kept =
        keep_by_alignment_pairwise(&header, &records, filter).map_err(ParseError::Inconsistent)?;
    Ok((header, kept))
}

/// The aligned filter for **BAM** output, byte-identical to Picard's `USE_JDK_DEFLATER=true` via
/// `BamWriter` (transitive, as for the read-list `_to_bam`).
pub fn filter_sam_reads_aligned_to_bam(
    input_sam: &str,
    filter: AlignedFilter,
) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, kept) = filter_aligned(input_sam, filter)?;
    let mut w = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &kept {
        w.write(rec).expect("record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:100000\n\
        keepA\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        dropB\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        keepC\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        dropD\t0\tchr1\t400\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        keepE\t0\tchr1\t500\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    const LIST: &str = "keepA\nkeepC\nkeepE\n";

    fn names(sam: &str) -> Vec<&str> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap())
            .collect()
    }

    #[test]
    fn include_keeps_the_listed_reads_in_input_order() {
        let out = filter_sam_reads(INPUT, LIST, Filter::IncludeReadList).unwrap();
        assert_eq!(names(&out), ["keepA", "keepC", "keepE"]);
        // The header is unchanged and no @PG is added.
        assert!(out.starts_with("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"));
        assert!(!out.contains("@PG"));
    }

    #[test]
    fn exclude_keeps_the_unlisted_reads_in_input_order() {
        let out = filter_sam_reads(INPUT, LIST, Filter::ExcludeReadList).unwrap();
        assert_eq!(names(&out), ["dropB", "dropD"]);
    }

    #[test]
    fn blank_lines_in_the_list_are_ignored() {
        let out = filter_sam_reads(INPUT, "keepA\n\n  \nkeepE\n", Filter::IncludeReadList).unwrap();
        assert_eq!(names(&out), ["keepA", "keepE"]);
    }

    const TAGGED: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:100000\n\
        @RG\tID:rg1\tSM:s\n@RG\tID:rg2\tSM:t\n\
        a\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
        b\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n\
        c\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
        d\t0\tchr1\t400\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    #[test]
    fn include_tag_values_keeps_reads_whose_tag_matches() {
        let out =
            filter_sam_reads_by_tag(TAGGED, b"RG", &["rg1"], TagFilter::IncludeTagValues).unwrap();
        assert_eq!(names(&out), ["a", "c"]);
    }

    #[test]
    fn exclude_tag_values_keeps_reads_whose_tag_does_not_match_including_absent() {
        // b has rg2 and d has no RG tag; both are kept, a and c (rg1) are dropped.
        let out =
            filter_sam_reads_by_tag(TAGGED, b"RG", &["rg1"], TagFilter::ExcludeTagValues).unwrap();
        assert_eq!(names(&out), ["b", "d"]);
    }

    #[test]
    fn multiple_tag_values_are_a_set() {
        let out =
            filter_sam_reads_by_tag(TAGGED, b"RG", &["rg1", "rg2"], TagFilter::IncludeTagValues)
                .unwrap();
        assert_eq!(names(&out), ["a", "b", "c"]);
    }

    // Queryname-sorted: a both-aligned pair, a one-unmapped pair, a both-unmapped pair, an aligned
    // singleton, and an unmapped singleton.
    const PAIRS: &str = "@HD\tVN:1.6\tSO:queryname\n\
        @SQ\tSN:chr1\tLN:100000\n\
        pairAA\t99\tchr1\t100\t60\t4M\t=\t200\t104\tACGT\tIIII\n\
        pairAA\t147\tchr1\t200\t60\t4M\t=\t100\t-104\tACGT\tIIII\n\
        pairAU\t97\tchr1\t300\t60\t4M\t=\t300\t0\tACGT\tIIII\n\
        pairAU\t141\t*\t0\t0\t*\t=\t300\t0\tACGT\tIIII\n\
        pairUU\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
        pairUU\t141\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
        singleA\t0\tchr1\t500\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        singleU\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

    #[test]
    fn include_aligned_keeps_fully_aligned_templates() {
        let out = filter_sam_reads_aligned(PAIRS, AlignedFilter::IncludeAligned).unwrap();
        assert_eq!(names(&out), ["pairAA", "pairAA", "singleA"]);
    }

    #[test]
    fn exclude_aligned_keeps_templates_with_an_unmapped_read() {
        let out = filter_sam_reads_aligned(PAIRS, AlignedFilter::ExcludeAligned).unwrap();
        assert_eq!(
            names(&out),
            ["pairAU", "pairAU", "pairUU", "pairUU", "singleU"]
        );
    }

    /// The BAM output decodes back to exactly the SAM output; the writer's byte-identity to htsjdk is
    /// proven elsewhere.
    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        let sam = filter_sam_reads(INPUT, LIST, Filter::IncludeReadList).unwrap();
        let bam = filter_sam_reads_to_bam(INPUT, LIST, Filter::IncludeReadList).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bam).expect("bam decompresses");
        let reader = htsjdk_bam::reader::BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(
            htsjdk_bam::sam_file::write_sam(&header, &records).unwrap(),
            sam
        );
    }

    fn bam_matches_sam(sam: &str, bam: &[u8]) {
        let plain = htsjdk_bgzf::decompress_all(bam).expect("bam decompresses");
        let reader = htsjdk_bam::reader::BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(
            htsjdk_bam::sam_file::write_sam(&header, &records).unwrap(),
            sam
        );
    }

    #[test]
    fn the_tag_bam_output_round_trips_to_the_sam_output() {
        let sam =
            filter_sam_reads_by_tag(TAGGED, b"RG", &["rg1"], TagFilter::IncludeTagValues).unwrap();
        let bam =
            filter_sam_reads_by_tag_to_bam(TAGGED, b"RG", &["rg1"], TagFilter::IncludeTagValues)
                .unwrap();
        bam_matches_sam(&sam, &bam);
    }

    #[test]
    fn the_aligned_bam_output_round_trips_to_the_sam_output() {
        let sam = filter_sam_reads_aligned(PAIRS, AlignedFilter::IncludeAligned).unwrap();
        let bam = filter_sam_reads_aligned_to_bam(PAIRS, AlignedFilter::IncludeAligned).unwrap();
        bam_matches_sam(&sam, &bam);
    }
}
