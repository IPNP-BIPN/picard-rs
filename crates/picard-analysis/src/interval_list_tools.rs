//! `IntervalListTools` (CONCAT / UNION slice).
//!
//! Ports `picard.util.IntervalListTools.doWork` at tag 3.4.0 for the two actions that need no
//! interval set-algebra: `CONCAT` (concatenate every input's intervals, then optionally sort and
//! unique) and `UNION` (`CONCAT` with `SORT` and `UNIQUE` forced on). The set actions (`INTERSECT`,
//! `SUBTRACT`, `SYMDIFF`, `OVERLAPS`) and the `INVERT` / `BREAK_BANDS_AT_MULTIPLES_OF` / `SCATTER` /
//! `PADDING` / multi-dictionary paths are separate surfaces, deferred; this slice assumes a single
//! shared sequence dictionary across the inputs (the common case) and `PADDING = 0`.
//!
//! The inputs are concatenated in file order (`openIntervalLists` reduces them with
//! `IntervalList.concatenate`), then:
//! * `SORT` (default on) sorts by the dictionary's `IntervalCoordinateComparator`;
//! * `UNIQUE` (default off) merges overlapping runs via `IntervalListTools.uniqued`, which is
//!   `getUniqueIntervals(mergeAbutting = !DONT_MERGE_ABUTTING, concatenateNames = true,
//!   enforceSameStrands = false)`: the merged interval keeps the first's start and strand, takes the
//!   max end, and joins the names of the run with `|`. With `DONT_MERGE_ABUTTING` set, only
//!   overlapping intervals merge; abutting ones (`end + 1 == next.start`) stay separate.
//! * `UNION` forces both on.
//!
//! ## Two output-header facts nailed against the oracle
//!
//! The output header's sort order is **always `SO:unsorted`** regardless of `SORT` (the concatenated
//! header carries `GroupOrder`/`SortOrder` unsorted), and the tool **always adds a `@PG`** whose `CL`
//! is the command line. The `@PG` is non-reproducible, so it is canonicalized away exactly as for
//! `DownsampleSam`: the port emits `@HD VN:1.6 SO:unsorted` + the dictionary's `@SQ` lines verbatim +
//! the interval body, with no `@PG`, and the conformance strips `@PG` from Picard's output before
//! comparing.

use htsjdk_bam::interval::{Interval, IntervalList, ParseError};

/// `IntervalListTools.Action`, restricted to the two set-algebra-free actions of this slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// `CONCAT`: concatenate all inputs, no implied sort or merge.
    Concat,
    /// `UNION`: `CONCAT` with `SORT` and `UNIQUE` forced on.
    Union,
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

/// `IntervalListTools.doWork` over the input interval lists, returning the output interval list.
///
/// All inputs are assumed to share one sequence dictionary; its `@SQ` lines (taken from the first
/// input) are emitted verbatim.
pub fn interval_list_tools(inputs: &[&str], opts: &Options) -> Result<String, ParseError> {
    let (names, sq_lines) = inputs.first().map(|s| dictionary(s)).unwrap_or_default();

    // openIntervalLists reduces the inputs with IntervalList.concatenate: all intervals in file
    // order. PADDING is 0 in this slice, so no padding is applied.
    let mut all: Vec<Interval> = Vec::new();
    for input in inputs {
        let parsed = IntervalList::parse_body(Vec::new(), input)?;
        all.extend(parsed.intervals);
    }
    let list = IntervalList {
        dictionary: names,
        intervals: all,
    };

    // UNION is CONCAT with SORT and UNIQUE forced on.
    let (sort, unique) = match opts.action {
        Action::Union => (true, true),
        Action::Concat => (opts.sort, opts.unique),
    };

    let sorted = if sort { list.sorted() } else { list };
    let final_intervals = if unique {
        uniqued(&sorted, !opts.dont_merge_abutting)
    } else {
        sorted.intervals
    };

    let mut output = String::from("@HD\tVN:1.6\tSO:unsorted\n");
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
    const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n";

    #[test]
    fn concat_sorted_orders_by_coordinate() {
        let opts = Options::default();
        let out = interval_list_tools(&[IN1, IN2], &opts).unwrap();
        assert_eq!(
            out,
            format!("{HEADER}chr1\t1\t10\t+\tA\nchr1\t5\t15\t+\tB\nchr1\t50\t60\t+\tC\nchr1\t61\t70\t+\tD\n")
        );
    }

    #[test]
    fn concat_unsorted_keeps_input_order() {
        let opts = Options {
            sort: false,
            ..Options::default()
        };
        let out = interval_list_tools(&[IN1, IN2], &opts).unwrap();
        assert_eq!(
            out,
            format!("{HEADER}chr1\t1\t10\t+\tA\nchr1\t50\t60\t+\tC\nchr1\t5\t15\t+\tB\nchr1\t61\t70\t+\tD\n")
        );
    }

    #[test]
    fn union_merges_overlaps_and_abutting_and_joins_names() {
        let opts = Options {
            action: Action::Union,
            ..Options::default()
        };
        let out = interval_list_tools(&[IN1, IN2], &opts).unwrap();
        assert_eq!(
            out,
            format!("{HEADER}chr1\t1\t15\t+\tA|B\nchr1\t50\t70\t+\tC|D\n")
        );
    }

    #[test]
    fn dont_merge_abutting_keeps_abutting_intervals_apart() {
        let opts = Options {
            unique: true,
            dont_merge_abutting: true,
            ..Options::default()
        };
        let out = interval_list_tools(&[IN1, IN2], &opts).unwrap();
        // 1-10 and 5-15 overlap -> 1-15 "A|B"; 50-60 and 61-70 only abut -> kept apart.
        assert_eq!(
            out,
            format!("{HEADER}chr1\t1\t15\t+\tA|B\nchr1\t50\t60\t+\tC\nchr1\t61\t70\t+\tD\n")
        );
    }
}
