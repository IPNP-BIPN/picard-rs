//! `IntervalListTools` (all actions plus the `INVERT` option; no scatter/padding/break-bands).
//!
//! Ports `picard.util.IntervalListTools.doWork` at tag 3.4.0 for every `ACTION`:
//! * `CONCAT`: concatenate every `INPUT`'s intervals (file order), then optionally `SORT` / `UNIQUE`.
//! * `UNION`: `CONCAT` with `SORT` and `UNIQUE` forced on.
//! * `INTERSECT`: `reduceEach = IntervalList.intersection`, the sorted-and-merged set of loci
//!   contained in all `INPUT`s. For two inputs, `intersection(a, b)` builds an overlap detector over
//!   `a`, and for each interval `i` of `b` intersects `i` with every overlapping `a` interval
//!   (`i.intersect(j)`, so the merged name reads `"{b-name} intersection {a-name}"`), then
//!   `uniqued()`.
//! * `OVERLAPS`: keep each whole `INPUT` interval that overlaps any interval of `SECOND_INPUT` (an
//!   overlap detector over `SECOND_INPUT.sorted().uniqued()`).
//! * `SUBTRACT`: `intersection(INPUT, invert(SECOND_INPUT))`, the loci in `INPUT` not in
//!   `SECOND_INPUT`.
//! * `SYMDIFF`: `union(subtract(a, b), subtract(b, a))`, the loci in exactly one of the two.
//!
//! The `INVERT` option (orthogonal to `ACTION`) complements the result against the dictionary and
//! forces `SORT=false`, `UNIQUE=true`. `SUBTRACT`, `SYMDIFF` and `INVERT` all use
//! [`IntervalList::invert`](htsjdk_bam::interval::IntervalList::invert), which needs the `(name,
//! length)` dictionary read from the `@SQ LN` fields. `PADDING` pads each input (via
//! [`IntervalList::padded`](htsjdk_bam::interval::IntervalList::padded), clamped to the contig)
//! before the action, and `BREAK_BANDS_AT_MULTIPLES_OF` splits the final intervals via
//! [`break_intervals_at_band_multiples`](htsjdk_bam::interval::break_intervals_at_band_multiples).
//! Still deferred: `SCATTER` (a multi-directory output shape) and multi-dictionary union; this slice
//! assumes one shared dictionary.
//!
//! ## The output header, per action, nailed against the oracle
//!
//! The tool always adds a `@PG` whose `CL` is the command line; it is non-reproducible, so it is
//! canonicalized away exactly as for `DownsampleSam` (the port emits none, the conformance strips
//! `@PG`). `doWork` takes the output header from `result.getHeader()`, whose sort order is
//! `SO:unsorted` only when `act` ran `setSortOrder(unsorted)`: a **multi-`INPUT`** `concatenate`
//! (`addOther`), `OVERLAPS`, or `SYMDIFF`'s `union`. A single-`INPUT` `reduce` returns the input
//! header untouched, and `INTERSECT`/`SUBTRACT` clone the left operand's header, so those are emitted
//! verbatim (no forced `SO`). `sorted()` and `invert()` never rewrite the header's sort order.

use htsjdk_bam::interval::{break_intervals_at_band_multiples, Interval, IntervalList, ParseError};
use htsjdk_bam::overlap::OverlapDetector;

/// `IntervalListTools.Action`.
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
    /// `SUBTRACT`: loci in `INPUT` but not `SECOND_INPUT` (`intersection(lhs, invert(rhs))`).
    Subtract,
    /// `SYMDIFF`: loci in `INPUT` or `SECOND_INPUT` but not both
    /// (`union(subtract(a, b), subtract(b, a))`).
    Symdiff,
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
    /// `INVERT`: complement the result against the dictionary (forces `SORT=false`, `UNIQUE=true`).
    pub invert: bool,
    /// `PADDING`: pad every input interval by this many bases on each side (clamped to the contig),
    /// applied before the action. `0` is a no-op.
    pub padding: i32,
    /// `BREAK_BANDS_AT_MULTIPLES_OF`: split the final intervals at integer multiples of this value.
    /// `0` disables.
    pub break_bands_at_multiples_of: i32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            action: Action::Concat,
            sort: true,
            unique: false,
            dont_merge_abutting: false,
            invert: false,
            padding: 0,
            break_bands_at_multiples_of: 0,
        }
    }
}

/// The ordered `(name, length)` dictionary (from `@SQ SN:`/`LN:`, for the coordinate sort and for
/// `invert`) and the verbatim `@SQ` lines (for the output header).
fn dictionary(interval_list: &str) -> (Vec<(String, i32)>, Vec<String>) {
    let mut sequences = Vec::new();
    let mut sq_lines = Vec::new();
    for line in interval_list.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }
        let mut name = None;
        let mut length = None;
        for field in line.split('\t') {
            if let Some(v) = field.strip_prefix("SN:") {
                name = Some(v.to_string());
            } else if let Some(v) = field.strip_prefix("LN:") {
                length = v.parse().ok();
            }
        }
        if let (Some(name), Some(length)) = (name, length) {
            sequences.push((name, length));
        }
        sq_lines.push(line.to_string());
    }
    (sequences, sq_lines)
}

/// The output `@HD` line, which `doWork` takes from `result.getHeader()`. When `unsorted` (the
/// header went through `setSortOrder(unsorted)`), any existing `SO` is rewritten to `SO:unsorted`
/// after `VN`; otherwise the first input's `@HD` is emitted verbatim (a single-`INPUT` reduce, or
/// `INTERSECT`/`SUBTRACT` whose `intersection` clones the left header unchanged).
fn output_hd(first_input: &str, unsorted: bool) -> String {
    let hd = first_input
        .lines()
        .find(|l| l.starts_with("@HD"))
        .unwrap_or("@HD\tVN:1.6");
    if !unsorted {
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

/// `IntervalList.subtract(lhs, rhs)` = `intersection(lhs, invert(rhs))`: the loci in `lhs` not in
/// `rhs`. `invert` needs the `(name, length)` dictionary.
fn subtract(
    lhs: &[Interval],
    rhs: &[Interval],
    names: &[String],
    sequences: &[(String, i32)],
) -> Vec<Interval> {
    let inverted = IntervalList {
        dictionary: names.to_vec(),
        intervals: rhs.to_vec(),
    }
    .invert(sequences);
    intersection(lhs, &inverted.intervals, names)
}

/// `IntervalList.padded(PADDING)`: pad each interval by `padding` on both sides, clamped to
/// `[1, contig length]`. A `padding` of 0 is the identity (and needs no lengths).
fn pad(
    intervals: &[Interval],
    padding: i32,
    sequences: &[(String, i32)],
    names: &[String],
) -> Vec<Interval> {
    if padding == 0 {
        return intervals.to_vec();
    }
    IntervalList {
        dictionary: names.to_vec(),
        intervals: intervals.to_vec(),
    }
    .padded(padding, padding, |contig| {
        sequences
            .iter()
            .find(|(n, _)| n == contig)
            .map(|(_, l)| *l)
            .unwrap_or(0)
    })
    .intervals
}

/// `IntervalList.union(a, b)` = `concatenate(a, b).uniqued()`: merge the two, then unique.
fn union2(a: &[Interval], b: &[Interval], names: &[String]) -> Vec<Interval> {
    let mut both = a.to_vec();
    both.extend_from_slice(b);
    uniqued(
        &IntervalList {
            dictionary: names.to_vec(),
            intervals: both,
        },
        true,
    )
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
    let (sequences, sq_lines) = inputs.first().map(|s| dictionary(s)).unwrap_or_default();
    let names: Vec<String> = sequences.iter().map(|(n, _)| n.clone()).collect();

    // openIntervalLists pads each file (PADDING) then reduces with the action's reduceEach:
    // concatenate for every action except INTERSECT, which reduces by intersection.
    let lists: Vec<Vec<Interval>> = inputs
        .iter()
        .map(|s| parse(s).map(|iv| pad(&iv, opts.padding, &sequences, &names)))
        .collect::<Result<_, _>>()?;
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

    // The reduced SECOND_INPUT (padded then concatenated), used by the actions that take a second
    // input.
    let second: Vec<Interval> = second_inputs
        .iter()
        .map(|s| parse(s).map(|iv| pad(&iv, opts.padding, &sequences, &names)))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
        .collect();

    // act: combine the reduced INPUT with the reduced SECOND_INPUT per the action.
    let result: Vec<Interval> = match opts.action {
        Action::Overlaps => overlaps(&combined, &second, &names),
        Action::Subtract => subtract(&combined, &second, &names, &sequences),
        Action::Symdiff => union2(
            &subtract(&combined, &second, &names, &sequences),
            &subtract(&second, &combined, &names, &sequences),
            &names,
        ),
        _ => combined,
    };

    // UNION forces SORT and UNIQUE on; INVERT forces SORT off and UNIQUE on.
    let mut sort = opts.sort;
    let mut unique = opts.unique;
    if opts.invert {
        sort = false;
        unique = true;
    }
    if opts.action == Action::Union {
        sort = true;
        unique = true;
    }

    let list = IntervalList {
        dictionary: names.clone(),
        intervals: result,
    };
    let sorted = if sort { list.sorted() } else { list };
    let inverted = if opts.invert {
        sorted.invert(&sequences)
    } else {
        sorted
    };
    let mut final_intervals = if unique {
        uniqued(&inverted, !opts.dont_merge_abutting)
    } else {
        inverted.intervals
    };

    // BREAK_BANDS splits the final intervals at integer multiples of the band.
    if opts.break_bands_at_multiples_of > 0 {
        final_intervals =
            break_intervals_at_band_multiples(&final_intervals, opts.break_bands_at_multiples_of);
    }

    // The output header is result.getHeader(): SO:unsorted when act ran setSortOrder(unsorted) - a
    // multi-INPUT concatenate (addOther), OVERLAPS, or SYMDIFF's union - else the input header
    // verbatim (single-INPUT reduce, or INTERSECT/SUBTRACT whose intersection clones the left).
    let unsorted = match opts.action {
        Action::Overlaps | Action::Symdiff => true,
        Action::Intersect => false,
        Action::Concat | Action::Union | Action::Subtract => inputs.len() >= 2,
    };
    let mut output = output_hd(inputs.first().copied().unwrap_or(""), unsorted);
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
    fn invert_complements_against_the_dictionary() {
        // CONCAT (single input) + INVERT: complement of {10-20, 50-60} on chr1[1..100].
        let a = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t10\t20\t+\tA\nchr1\t50\t60\t+\tC\n";
        let opts = Options {
            invert: true,
            ..Options::default()
        };
        let out = interval_list_tools(&[a], &[], &opts).unwrap();
        // Single INPUT -> header verbatim (no SO). Gaps named interval-1..3.
        assert_eq!(
            out,
            format!("{NOSORT}chr1\t1\t9\t+\tinterval-1\nchr1\t21\t49\t+\tinterval-2\nchr1\t61\t100\t+\tinterval-3\n")
        );
    }

    #[test]
    fn subtract_removes_the_second_from_the_first() {
        // {1-50} minus {20-30} = {1-19, 31-50}. Single INPUT -> header verbatim.
        let big = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t50\t+\tA\n";
        let mid = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t20\t30\t+\tB\n";
        let opts = Options {
            action: Action::Subtract,
            ..Options::default()
        };
        let out = interval_list_tools(&[big], &[mid], &opts).unwrap();
        assert_eq!(
            out,
            format!("{NOSORT}chr1\t1\t19\t+\tinterval-1 intersection A\nchr1\t31\t50\t+\tinterval-2 intersection A\n")
        );
    }

    #[test]
    fn symdiff_keeps_loci_in_exactly_one_input() {
        // {1-30} xor {20-50} = {1-19, 31-50}. SYMDIFF's union -> SO:unsorted.
        let left = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t30\t+\tA\n";
        let right = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t20\t50\t+\tB\n";
        let opts = Options {
            action: Action::Symdiff,
            ..Options::default()
        };
        let out = interval_list_tools(&[left], &[right], &opts).unwrap();
        assert_eq!(
            out,
            format!("{UNSORTED}chr1\t1\t19\t+\tinterval-1 intersection A\nchr1\t31\t50\t+\tinterval-1 intersection B\n")
        );
    }

    #[test]
    fn padding_pads_each_interval_clamped_to_the_contig() {
        // PADDING=10: A(5-20) -> 1-30 (start clamped), C(95-100) -> 85-100 (end clamped).
        let a = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t5\t20\t+\tA\nchr1\t95\t100\t+\tC\n";
        let opts = Options {
            padding: 10,
            ..Options::default()
        };
        let out = interval_list_tools(&[a], &[], &opts).unwrap();
        assert_eq!(
            out,
            format!("{NOSORT}chr1\t1\t30\t+\tA\nchr1\t85\t100\t+\tC\n")
        );
    }

    #[test]
    fn break_bands_splits_the_final_intervals() {
        // BREAK_BANDS=10: A(5-25) -> A.1(5-9), A.2(10-19), A.3(20-25).
        let a = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t5\t25\t+\tA\n";
        let opts = Options {
            break_bands_at_multiples_of: 10,
            ..Options::default()
        };
        let out = interval_list_tools(&[a], &[], &opts).unwrap();
        assert_eq!(
            out,
            format!("{NOSORT}chr1\t5\t9\t+\tA.1\nchr1\t10\t19\t+\tA.2\nchr1\t20\t25\t+\tA.3\n")
        );
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
