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
//! The interval, tag, mapping-quality, and JavaScript filters, `WRITE_READS_FILES`, and a non-default
//! `SORT_ORDER` are separate surfaces.

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
fn read_name_set(read_list_text: &str) -> HashSet<&str> {
    read_list_text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect()
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
/// This is a **sequential** grouped pass (the pairing depends on encounter order), so unlike the
/// read-list and tag-value filters it is not parallelized. Templates of one or two reads are handled;
/// a template with secondary/supplementary records is a separate surface and is asserted against.
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
    let include = matches!(filter, AlignedFilter::IncludeAligned);

    let mut kept: Vec<BamRecord> = Vec::with_capacity(records.len());
    let mut start = 0;
    while start < records.len() {
        let mut end = start + 1;
        while end < records.len() && records[end].read_name == records[start].read_name {
            end += 1;
        }
        let group = &records[start..end];
        assert!(
            group.len() <= 2,
            "filter_sam_reads_aligned: templates with more than two reads are not ported"
        );

        // A template is "aligned" when all its reads are aligned; keep it iff that matches include.
        let all_aligned = group.iter().all(is_aligned);
        if all_aligned == include {
            kept.extend_from_slice(group);
        }
        start = end;
    }
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
