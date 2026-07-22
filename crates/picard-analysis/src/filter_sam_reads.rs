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
}
