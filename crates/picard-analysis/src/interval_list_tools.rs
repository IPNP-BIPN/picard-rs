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
//! [`invert`](crate::compact_interval::invert), which needs the `(name,
//! length)` dictionary read from the `@SQ LN` fields. `PADDING` pads each input (via
//! [`padded`](crate::compact_interval::padded), clamped to the contig)
//! before the action, and `BREAK_BANDS_AT_MULTIPLES_OF` splits the final intervals via
//! [`break_at_band_multiples`](crate::compact_interval::break_at_band_multiples).
//! Still deferred: `SCATTER` (a multi-directory output shape) and multi-dictionary union; this slice
//! assumes one shared dictionary.
//!
//! ## The output header, per action, nailed against the oracle
//!
//! The tool always adds a `@PG` whose `CL` is the command line; the binary adds it (see
//! `command_line`), after the input's own `@PG` records, which the output header keeps along with
//! its read groups and comments. `doWork` takes the output header from `result.getHeader()`, whose sort order is
//! `SO:unsorted` only when `act` ran `setSortOrder(unsorted)`: a **multi-`INPUT`** `concatenate`
//! (`addOther`), `OVERLAPS`, or `SYMDIFF`'s `union`. A single-`INPUT` `reduce` returns the input
//! header untouched, and `INTERSECT`/`SUBTRACT` clone the left operand's header, so those are emitted
//! verbatim (no forced `SO`). `sorted()` and `invert()` never rewrite the header's sort order.

use std::io::{self, BufRead};

use htsjdk_bam::interval::ParseError;
use htsjdk_bam::overlap::OverlapDetector;

use crate::compact_interval::{self as ci, Contigs, Header, Iv};

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
/// `IntervalList.intersection(list1, list2)`: for each interval of `list2`, intersect it with every
/// overlapping interval of `list1`, then `uniqued()`. The result names read
/// `"{list2-name} intersection {list1-name}"`.
fn intersection(contigs: &Contigs, list1: &[Iv], list2: &[Iv]) -> Vec<Iv> {
    let mut detector: OverlapDetector<usize> = OverlapDetector::create();
    for (index, i) in list1.iter().enumerate() {
        detector.add(contigs.name(i.contig), i.start, i.end, index);
    }
    let mut hits = Vec::new();
    for i in list2 {
        for &j in detector.get_overlaps(contigs.name(i.contig), i.start, i.end) {
            hits.push(i.intersect(&list1[j]));
        }
    }
    ci::uniqued(contigs, hits, true)
}

/// `IntervalList.overlaps(lhs, rhs)`: every whole `lhs` interval that overlaps any `rhs` interval,
/// in `lhs` order. `rhs` is `sorted().uniqued()` before the overlap detector is built.
fn overlaps(contigs: &Contigs, mut lhs: Vec<Iv>, rhs: Vec<Iv>) -> Vec<Iv> {
    let mut detector: OverlapDetector<()> = OverlapDetector::create();
    for i in &ci::uniqued(contigs, rhs, true) {
        detector.add(contigs.name(i.contig), i.start, i.end, ());
    }
    lhs.retain(|i| detector.overlaps_any(contigs.name(i.contig), i.start, i.end));
    lhs
}

/// `IntervalList.subtract(lhs, rhs)` = `intersection(lhs, invert(rhs))`: the loci in `lhs` not in
/// `rhs`. `invert` needs the `(name, length)` dictionary.
fn subtract(contigs: &mut Contigs, lhs: &[Iv], rhs: &[Iv], sequences: &[(String, i32)]) -> Vec<Iv> {
    let inverted = ci::invert(contigs, sequences, rhs.to_vec());
    intersection(contigs, lhs, &inverted)
}

/// What `IntervalListTools.doWork` writes, before the binary adds its `@PG` and comments.
pub struct Output {
    /// The header lines, each ending in a newline.
    pub header: String,
    pub contigs: Contigs,
    pub intervals: Vec<Iv>,
}

/// `IntervalListTools.doWork` over the `INPUT` and `SECOND_INPUT` interval lists, opened in that
/// order as they are reached.
///
/// All inputs are assumed to share one sequence dictionary; its `@SQ` lines (taken from the first
/// `INPUT`) are emitted verbatim. A line refused in one list is reported only once every later list
/// has read cleanly, the order a whole-file read of each gave.
pub fn run<R: BufRead>(
    inputs: impl IntoIterator<Item = io::Result<R>>,
    second_inputs: impl IntoIterator<Item = io::Result<R>>,
    opts: &Options,
) -> io::Result<Result<Output, ParseError>> {
    let mut contigs = Contigs::default();
    let mut first: Option<Header> = None;
    let mut sequences: Vec<(String, i32)> = Vec::new();
    let mut input_count = 0;
    let mut combined: Option<Vec<Iv>> = None;
    // The reduced SECOND_INPUT (padded then concatenated), used by the actions that take one.
    let mut second: Vec<Iv> = Vec::new();
    let mut refusal = None;
    let all = inputs.into_iter().map(|r| (false, r));
    for (is_second, reader) in all.chain(second_inputs.into_iter().map(|r| (true, r))) {
        let reader = reader?;
        if refusal.is_some() {
            ci::drain(reader)?;
            continue;
        }
        let mut list = Vec::new();
        let header = match ci::read_list(reader, &mut contigs, &mut list)? {
            Ok(header) => header,
            Err(error) => {
                refusal = Some(error);
                continue;
            }
        };
        if !is_second && first.is_none() {
            sequences = dictionary(&header.text()).0;
            let names: Vec<String> = sequences.iter().map(|(n, _)| n.clone()).collect();
            contigs.set_dictionary(&names);
            first = Some(header);
        }
        // openIntervalLists pads each file (PADDING) then reduces with the action's reduceEach:
        // concatenate for every action except INTERSECT, which reduces by intersection.
        let list = ci::padded(&contigs, &sequences, list, opts.padding);
        if is_second {
            second.extend(list);
            continue;
        }
        input_count += 1;
        combined = Some(match combined {
            None => list,
            Some(acc) if opts.action == Action::Intersect => intersection(&contigs, &acc, &list),
            Some(mut acc) => {
                acc.extend(list);
                acc
            }
        });
    }
    if let Some(error) = refusal {
        return Ok(Err(error));
    }
    let combined = combined.unwrap_or_default();

    // act: combine the reduced INPUT with the reduced SECOND_INPUT per the action.
    let mut result: Vec<Iv> = match opts.action {
        Action::Overlaps => overlaps(&contigs, combined, second),
        Action::Subtract => subtract(&mut contigs, &combined, &second, &sequences),
        Action::Symdiff => {
            let mut both = subtract(&mut contigs, &combined, &second, &sequences);
            both.extend(subtract(&mut contigs, &second, &combined, &sequences));
            ci::uniqued(&contigs, both, true)
        }
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
    if sort {
        ci::sort(&contigs, &mut result);
    }
    if opts.invert {
        result = ci::invert(&mut contigs, &sequences, result);
    }
    if unique {
        result = ci::uniqued(&contigs, result, !opts.dont_merge_abutting);
    }
    // BREAK_BANDS splits the final intervals at integer multiples of the band.
    if opts.break_bands_at_multiples_of > 0 {
        result = ci::break_at_band_multiples(result, opts.break_bands_at_multiples_of);
    }

    // The output header is result.getHeader(): SO:unsorted when act ran setSortOrder(unsorted) - a
    // multi-INPUT concatenate (addOther), OVERLAPS, or SYMDIFF's union - else the input header
    // verbatim (single-INPUT reduce, or INTERSECT/SUBTRACT whose intersection clones the left).
    let unsorted = match opts.action {
        Action::Overlaps | Action::Symdiff => true,
        Action::Intersect => false,
        Action::Concat | Action::Union | Action::Subtract => input_count >= 2,
    };
    let first = first.unwrap_or_default();
    let first_text = first.text();
    let mut header = output_hd(&first_text, unsorted);
    header.push('\n');
    for sq in &dictionary(&first_text).1 {
        header.push_str(sq);
        header.push('\n');
    }
    // The rest of that header comes along: its read groups, its program records and its comments,
    // in the order the text codec writes them.
    for prefix in ["@RG\t", "@PG\t", "@CO\t"] {
        for line in first.leading_lines() {
            if line.starts_with(prefix) {
                header.push_str(line);
                header.push('\n');
            }
        }
    }
    Ok(Ok(Output {
        header,
        contigs,
        intervals: result,
    }))
}

/// [`run`] over the interval lists' texts, returning the output interval list.
pub fn interval_list_tools(
    inputs: &[&str],
    second_inputs: &[&str],
    opts: &Options,
) -> Result<String, ParseError> {
    fn open<'a>(texts: &[&'a str]) -> Vec<io::Result<&'a [u8]>> {
        texts.iter().map(|text| Ok(text.as_bytes())).collect()
    }
    let output = run(open(inputs), open(second_inputs), opts).expect("a text reads")?;
    let mut text = output.header.into_bytes();
    ci::write_lines(&output.contigs, &output.intervals, &mut text).expect("a Vec writes");
    Ok(String::from_utf8(text).expect("the inputs were text"))
}

/// `IntervalListScatterMode`: how `SCATTER_COUNT` divides the final list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatterMode {
    /// `INTERVAL_SUBDIVISION`: equal base counts, cutting intervals where it must.
    IntervalSubdivision,
    /// `BALANCING_WITHOUT_INTERVAL_SUBDIVISION`: whole intervals, a list closing before it would
    /// pass the ideal size.
    BalancingWithoutIntervalSubdivision,
    /// `BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW`: the same, but a list may pass the
    /// ideal size while the lists still to come would be larger than it.
    BalancingWithoutIntervalSubdivisionWithOverflow,
    /// `INTERVAL_COUNT`: equal interval counts, the remainder in the last list.
    IntervalCount,
    /// `INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER`: equal interval counts, the remainder spread.
    IntervalCountWithDistributedRemainder,
}

impl ScatterMode {
    pub fn parse(name: &str) -> Option<ScatterMode> {
        Some(match name {
            "INTERVAL_SUBDIVISION" => ScatterMode::IntervalSubdivision,
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION" => {
                ScatterMode::BalancingWithoutIntervalSubdivision
            }
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW" => {
                ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow
            }
            "INTERVAL_COUNT" => ScatterMode::IntervalCount,
            "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER" => {
                ScatterMode::IntervalCountWithDistributedRemainder
            }
            _ => return None,
        })
    }

    fn by_count(self) -> bool {
        matches!(
            self,
            ScatterMode::IntervalCount | ScatterMode::IntervalCountWithDistributedRemainder
        )
    }

    /// `listWeight`: the bases of a list, or its interval count.
    pub fn list_weight(self, intervals: &[Iv]) -> i64 {
        if self.by_count() {
            intervals.len() as i64
        } else {
            intervals.iter().map(|i| i64::from(i.length())).sum()
        }
    }
}

/// `IntervalListScatter` over the written list (its header lines and its intervals): the header
/// every produced list carries, and the lists in order.
///
/// The lists carry the header of the preprocessed input: `sorted()` marks it coordinate-sorted, and
/// `uniqued()` (the subdivision mode's) clones the input's header as it was.
pub fn scatter(
    header: &str,
    intervals: Vec<Iv>,
    contigs: &Contigs,
    count: i32,
    mode: ScatterMode,
) -> (String, Vec<Vec<Iv>>) {
    use std::collections::VecDeque;
    // `IntervalListTools` reads the list it wrote back in: a name's final `\r` (from a last input
    // line with no newline) does not survive the line split, and a `.` left behind is no name.
    let mut intervals = intervals;
    for interval in &mut intervals {
        if let Some(name) = interval.name.as_deref().and_then(|n| n.strip_suffix('\r')) {
            interval.name = (name != ".").then(|| name.into());
        }
    }
    let length = |i: &Iv| i64::from(i.length());
    // `preprocessIntervalList`.
    let processed: Vec<Iv> = if mode == ScatterMode::IntervalSubdivision {
        ci::uniqued(contigs, intervals, true)
    } else {
        ci::sort(contigs, &mut intervals);
        intervals
    };
    // `deduceIdealSplitWeight`.
    let ideal: i64 = if mode.by_count() {
        (processed.len() as i64 / i64::from(count)).max(1)
    } else {
        let unique: i64 = ci::uniqued(contigs, processed.clone(), true)
            .iter()
            .map(length)
            .sum();
        let split = (unique / i64::from(count)).max(1);
        if mode == ScatterMode::IntervalSubdivision {
            split
        } else {
            let widest = processed.iter().map(length).max().unwrap_or(1);
            widest.max(split)
        }
    };
    // `takeSome`: the part that goes into the current list, and the part put back.
    let take_some = |interval: Iv, current: i64, remaining: f64| -> (Option<Iv>, Option<Iv>) {
        match mode {
            ScatterMode::IntervalSubdivision => {
                let amount = ideal - current;
                if amount >= length(&interval) {
                    (Some(interval), None)
                } else if amount == 0 {
                    (None, Some(interval))
                } else {
                    let cut = interval.start + amount as i32;
                    let left = Iv {
                        end: cut - 1,
                        ..interval.clone()
                    };
                    let right = Iv {
                        start: cut,
                        ..interval
                    };
                    (Some(left), Some(right))
                }
            }
            ScatterMode::BalancingWithoutIntervalSubdivision
            | ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow => {
                let projected = current + length(&interval);
                let include = projected <= ideal
                    || (mode == ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow
                        && (ideal as f64) < remaining);
                if include {
                    (Some(interval), None)
                } else {
                    (None, Some(interval))
                }
            }
            ScatterMode::IntervalCount => {
                if ideal - current > 0 {
                    (Some(interval), None)
                } else {
                    (None, Some(interval))
                }
            }
            ScatterMode::IntervalCountWithDistributedRemainder => {
                if remaining > current as f64 {
                    (Some(interval), None)
                } else {
                    (None, Some(interval))
                }
            }
        }
    };

    // The header the lists are written with.
    let mut list_header = String::new();
    for line in header.lines() {
        if line.starts_with("@HD") && mode != ScatterMode::IntervalSubdivision {
            let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            match fields.iter_mut().find(|f| f.starts_with("SO:")) {
                Some(field) => *field = "SO:coordinate".to_string(),
                None => fields.push("SO:coordinate".to_string()),
            }
            list_header.push_str(&fields.join("\t"));
        } else {
            list_header.push_str(line);
        }
        list_header.push('\n');
    }

    let mut weight_remaining = mode.list_weight(&processed);
    let mut queue: VecDeque<Iv> = processed.into_iter().collect();
    let mut returned: i64 = 0;
    let mut lists: Vec<Vec<Iv>> = Vec::new();
    while !queue.is_empty() {
        returned += 1;
        let mut running: Vec<Iv> = Vec::new();
        let mut closed = false;
        while !queue.is_empty() && returned < i64::from(count) {
            let interval = queue.pop_front().expect("the queue is not empty");
            let current = mode.list_weight(&running);
            let remaining =
                (weight_remaining - current) as f64 / (f64::from(count) - returned as f64);
            let (first, second) = take_some(interval, current, remaining);
            if let Some(second) = second {
                queue.push_front(second);
            }
            match first {
                Some(first) => running.push(first),
                None => {
                    weight_remaining -= current;
                    closed = true;
                    break;
                }
            }
        }
        if !closed {
            running.extend(queue.drain(..));
        }
        lists.push(running);
    }
    (list_header, lists)
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
