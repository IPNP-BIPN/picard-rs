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

    Ok(write_sam(&header, &kept).expect("records that parsed re-encode as SAM text"))
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

    Ok(write_sam(&header, &kept).expect("records that parsed re-encode as SAM text"))
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
}
