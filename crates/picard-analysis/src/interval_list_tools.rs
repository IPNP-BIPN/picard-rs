//! `IntervalListTools` (CONCAT / UNION / INTERSECT / OVERLAPS slice).
//!
//! Ports `picard.util.IntervalListTools.doWork` at tag 3.4.0 for the four actions that do not need
//! `IntervalList.invert` (and therefore no contig lengths):
//! * `CONCAT`: concatenate every `INPUT`'s intervals (file order), then optionally `SORT` / `UNIQUE`.
//! * `UNION`: `CONCAT` with `SORT` and `UNIQUE` forced on.
//! * `INTERSECT`: `reduceEach = IntervalList.intersection`, the sorted-and-merged set of loci
//!   contained in all `INPUT`s. For two inputs, `intersection(a, b)` builds an overlap detector over
//!   `a`, and for each interval `i` of `b` intersects `i` with every overlapping `a` interval
//!   (`i.intersect(j)`, so the merged name reads `"{b-name} intersection {a-name}"`), then
//!   `uniqued()`.
//! * `OVERLAPS`: keep each whole `INPUT` interval that overlaps any interval of `SECOND_INPUT` (an
//!   overlap detector over `SECOND_INPUT.sorted().uniqued()`).
//!
//! `SUBTRACT` and `SYMDIFF` both route through `invert` (so they need contig lengths) and are
//! deferred with `INVERT` / `BREAK_BANDS_AT_MULTIPLES_OF` / `SCATTER` / `PADDING > 0` /
//! multi-dictionary union. This slice assumes one shared sequence dictionary across the inputs and
//! `PADDING = 0`.
//!
//! ## The output header, per action, nailed against the oracle
//!
//! The tool always adds a `@PG` whose `CL` is the command line; it is non-reproducible, so it is
//! canonicalized away exactly as for `DownsampleSam` (the port emits none, the conformance strips
//! `@PG`). The `@HD` sort order, however, differs by action: `CONCAT`/`UNION` (via
//! `IntervalList.addOther`) and `OVERLAPS` (via `overlaps`) both call `setSortOrder(unsorted)`, so
//! their `@HD` is `SO:unsorted`; `INTERSECT` clones the first input's header unchanged, so its `@HD`
//! is emitted verbatim. `sorted()` never rewrites the header's sort order, which is why even a sorted
//! result keeps `SO:unsorted`.

use htsjdk_bam::interval::{Interval, IntervalList, ParseError};
use htsjdk_bam::overlap::OverlapDetector;

/// `IntervalListTools.Action`, restricted to the set-algebra actions this slice covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// `CONCAT`: concatenate all inputs, no implied sort or merge.
    Concat,
    /// `UNION`: `CONCAT` with `SORT` and `UNIQUE` forced on.
    Union,
    /// `INTERSECT`: the sorted, merged set of loci contained in all inputs.
    Intersect,
    /// `OVERLAPS`: whole `INPUT` intervals that overlap any `SECOND_INPUT` interval.
    Overlaps,
}

/// `IntervalListTools`'s options for this slice.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `ACTION`.
    pub action: Action,
    /// `SORT`: coordinate-sort the result (ignored, treated as `true`, when `action` is `UNION`).
    pub sort: bool,
    /// `UNIQUE`: merge overlapping (and, unless `dont_merge_abutting`, abutting) intervals.
    pub unique: bool,
    /// `DONT_MERGE_ABUTTING`: when uniquing, do not merge intervals that only abut.
    pub dont_merge_abutting: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            action: Action::Concat,
            sort: true,
            unique: false,
            dont_merge_abutting: false,
        }
    }
}

/// The `@SQ` `SN:` names in order (for the coordinate sort) and the verbatim `@SQ` lines (for the
/// output header).
fn dictionary(interval_list: &str) -> (Vec<String>, Vec<String>) {
    let mut names = Vec::new();
    let mut sq_lines = Vec::new();
    for line in interval_list.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }
        for field in line.split('\t') {
            if let Some(v) = field.strip_prefix("SN:") {
                names.push(v.to_string());
            }
        }
        sq_lines.push(line.to_string());
    }
    (names, sq_lines)
}

/// The output `@HD` line. `INTERSECT` clones the first input's header verbatim; every other action
/// calls `setSortOrder(unsorted)`, which rewrites (or adds) `SO:unsorted` after `VN`.
fn output_hd(first_input: &str, force_unsorted: bool) -> String {
    let hd = first_input
        .lines()
        .find(|l| l.starts_with("@HD"))
        .unwrap_or("@HD\tVN:1.6");
    if !force_unsorted {
        return hd.to_string();
    }
    let mut fields: Vec<&str> = hd.split('\t').filter(|f| !f.starts_with("SO:")).collect();
    fields.push("SO:unsorted");
    fields.join("\t")
}

fn parse(interval_list: &str) -> Result<Vec<Interval>, ParseError> {
    Ok(IntervalList::parse_body(Vec::new(), interval_list)?.intervals)
}

/// `IntervalListTools.uniqued(list, mergeAbutting)` = `getUniqueIntervals(mergeAbutting,
/// concatenateNames = true, enforceSameStrands = false)`: sort, then fold each overlapping (or, when
/// `merge_abutting`, abutting) run into one interval whose names are joined with `|`.
fn uniqued(list: &IntervalList, merge_abutting: bool) -> Vec<Interval> {
    let sorted = list.sorted();
    let mut out: Vec<Interval> = Vec::new();
    let mut current: Option<Interval> = None;
    let mut names: Vec<String> = Vec::new();

    let finish = |cur: &Interval, names: &[String]| Interval {
        contig: cur.contig.clone(),
        start: cur.start,
        end: cur.end,
        negative_strand: cur.negative_strand,
        name: if names.is_empty() {
            None
        } else {
            Some(names.join("|"))
        },
    };

    for next in &sorted.intervals {
        match &mut current {
            None => {
                current = Some(next.clone());
                names.clear();
                if let Some(n) = &next.name {
                    names.push(n.clone());
                }
            }
            Some(cur)
                if (merge_abutting && cur.within_distance_of(next, 1)) || cur.intersects(next) =>
            {
                cur.end = cur.end.max(next.end);
                if let Some(n) = &next.name {
                    names.push(n.clone());
                }
            }
            Some(cur) => {
                out.push(finish(cur, &names));
                current = Some(next.clone());
                names.clear();
                if let Some(n) = &next.name {
                    names.push(n.clone());
                }
            }
        }
    }
    if let Some(cur) = &current {
        out.push(finish(cur, &names));
    }
    out
}

/// `IntervalList.intersection(list1, list2)`: for each interval of `list2`, intersect it with every
/// overlapping interval of `list1`, then `uniqued()`. The result names read
/// `"{list2-name} intersection {list1-name}"`.
fn intersection(list1: &[Interval], list2: &[Interval], names: &[String]) -> Vec<Interval> {
    let mut detector: OverlapDetector<Interval> = OverlapDetector::create();
    for i in list1 {
        detector.add(&i.contig, i.start, i.end, i.clone());
    }
    let mut hits: Vec<Interval> = Vec::new();
    for i in list2 {
        for j in detector.get_overlaps(&i.contig, i.start, i.end) {
            hits.push(i.intersect(j));
        }
    }
    uniqued(
        &IntervalList {
            dictionary: names.to_vec(),
            intervals: hits,
        },
        true,
    )
}

/// `IntervalList.overlaps(lhs, rhs)`: every whole `lhs` interval that overlaps any `rhs` interval,
/// in `lhs` order. `rhs` is `sorted().uniqued()` before the overlap detector is built.
fn overlaps(lhs: &[Interval], rhs: &[Interval], names: &[String]) -> Vec<Interval> {
    let rhs_unique = uniqued(
        &IntervalList {
            dictionary: names.to_vec(),
            intervals: rhs.to_vec(),
        },
        true,
    );
    let mut detector: OverlapDetector<()> = OverlapDetector::create();
    for i in &rhs_unique {
        detector.add(&i.contig, i.start, i.end, ());
    }
    lhs.iter()
        .filter(|i| detector.overlaps_any(&i.contig, i.start, i.end))
        .cloned()
        .collect()
}

/// `IntervalListTools.doWork` over the `INPUT` and `SECOND_INPUT` interval lists, returning the
/// output interval list.
///
/// All inputs are assumed to share one sequence dictionary; its `@SQ` lines (taken from the first
/// `INPUT`) are emitted verbatim. `second_inputs` is only consulted for `OVERLAPS`.
pub fn interval_list_tools(
    inputs: &[&str],
    second_inputs: &[&str],
    opts: &Options,
) -> Result<String, ParseError> {
    let (names, sq_lines) = inputs.first().map(|s| dictionary(s)).unwrap_or_default();

    // openIntervalLists reduces the INPUT files with the action's reduceEach: concatenate for every
    // action except INTERSECT, which reduces by intersection. PADDING is 0, so no padding.
    let lists: Vec<Vec<Interval>> = inputs.iter().map(|s| parse(s)).collect::<Result<_, _>>()?;
    let combined: Vec<Interval> = match opts.action {
        Action::Intersect => {
            let mut acc = lists.first().cloned().unwrap_or_default();
            for next in lists.iter().skip(1) {
                acc = intersection(&acc, next, &names);
            }
            acc
        }
        _ => lists.into_iter().flatten().collect(),
    };

    // act: OVERLAPS combines with the reduced SECOND_INPUT; the others ignore it.
    let result: Vec<Interval> = match opts.action {
        Action::Overlaps => {
            let second: Vec<Interval> = second_inputs
                .iter()
                .map(|s| parse(s))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .flatten()
                .collect();
            overlaps(&combined, &second, &names)
        }
        _ => combined,
    };

    // UNION forces SORT and UNIQUE on.
    let (sort, unique) = match opts.action {
        Action::Union => (true, true),
        _ => (opts.sort, opts.unique),
    };

    let list = IntervalList {
        dictionary: names.clone(),
        intervals: result,
    };
    let sorted = if sort { list.sorted() } else { list };
    let final_intervals = if unique {
        uniqued(&sorted, !opts.dont_merge_abutting)
    } else {
        sorted.intervals
    };

    let force_unsorted = opts.action != Action::Intersect;
    let mut output = output_hd(inputs.first().copied().unwrap_or(""), force_unsorted);
    output.push('\n');
    for sq in &sq_lines {
        output.push_str(sq);
        output.push('\n');
    }
    for interval in &final_intervals {
        output.push_str(&interval.to_file_line());
        output.push('\n');
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    const IN1: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t10\t+\tA\nchr1\t50\t60\t+\tC\n";
    const IN2: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t5\t15\t+\tB\nchr1\t61\t70\t+\tD\n";
    const UNSORTED: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n";
    const NOSORT: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";

    #[test]
    fn concat_sorted_orders_by_coordinate() {
        let out = interval_list_tools(&[IN1, IN2], &[], &Options::default()).unwrap();
        assert_eq!(
            out,
            format!("{UNSORTED}chr1\t1\t10\t+\tA\nchr1\t5\t15\t+\tB\nchr1\t50\t60\t+\tC\nchr1\t61\t70\t+\tD\n")
        );
    }

    #[test]
    fn union_merges_overlaps_and_abutting_and_joins_names() {
        let opts = Options {
            action: Action::Union,
            ..Options::default()
        };
        let out = interval_list_tools(&[IN1, IN2], &[], &opts).unwrap();
        assert_eq!(
            out,
            format!("{UNSORTED}chr1\t1\t15\t+\tA|B\nchr1\t50\t70\t+\tC|D\n")
        );
    }

    #[test]
    fn intersect_keeps_the_overlap_with_the_second_over_first_name() {
        // in1 A(1-20) C(50-60); in2 B(10-30) D(55-70). B∩A=10-20, D∩C=55-60. Header verbatim (no SO).
        let a = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t20\t+\tA\nchr1\t50\t60\t+\tC\n";
        let b = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t10\t30\t+\tB\nchr1\t55\t70\t+\tD\n";
        let opts = Options {
            action: Action::Intersect,
            ..Options::default()
        };
        let out = interval_list_tools(&[a, b], &[], &opts).unwrap();
        assert_eq!(
            out,
            format!(
                "{NOSORT}chr1\t10\t20\t+\tB intersection A\nchr1\t55\t60\t+\tD intersection C\n"
            )
        );
    }

    #[test]
    fn overlaps_keeps_whole_input_intervals_that_hit_the_second() {
        let a = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t20\t+\tA\nchr1\t50\t60\t+\tC\nchr1\t80\t90\t+\tE\n";
        let second = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t10\t30\t+\tB\n";
        let opts = Options {
            action: Action::Overlaps,
            ..Options::default()
        };
        let out = interval_list_tools(&[a], &[second], &opts).unwrap();
        // Only A(1-20) overlaps B(10-30); C and E do not.
        assert_eq!(out, format!("{UNSORTED}chr1\t1\t20\t+\tA\n"));
    }

    #[test]
    fn dont_merge_abutting_keeps_abutting_intervals_apart() {
        let opts = Options {
            unique: true,
            dont_merge_abutting: true,
            ..Options::default()
        };
        let out = interval_list_tools(&[IN1, IN2], &[], &opts).unwrap();
        assert_eq!(
            out,
            format!("{UNSORTED}chr1\t1\t15\t+\tA|B\nchr1\t50\t60\t+\tC\nchr1\t61\t70\t+\tD\n")
        );
    }
}
