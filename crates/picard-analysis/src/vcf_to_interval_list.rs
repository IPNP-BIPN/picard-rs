//! Port of `picard.vcf.VcfToIntervalList` (Picard 3.4.0), with the two htsjdk pieces it is made of:
//! `VCFFileReader.toIntervals` and `IntervalList.IntervalMergerIterator`.
//!
//! # The intervals are merged as they come, not sorted first
//!
//! The tool hands the reader's intervals straight to the merging iterator, in FILE order. The
//! iterator only ever compares an interval with the one it is building, so a file out of
//! coordinate order is merged pairwise along the file: two overlapping sites separated by a site
//! elsewhere are two intervals, and the output is in the file's order too.
//!
//! # A name is counted only when it is needed
//!
//! A site with no ID (`.`) is named `interval-<n>`, and `n` counts the unnamed sites that were
//! KEPT: the filter runs before the naming, so a filtered site without an ID takes a number only
//! when `INCLUDE_FILTERED` keeps it, and every later name shifts with it.
//!
//! # `CONCAT_ALL` and `USE_FIRST` differ in more than the name
//!
//! With `CONCAT_ALL` the merged interval is `IntervalList.merge` over everything folded into it:
//! the smallest start, the largest end, and the distinct names joined with `|` in first-seen
//! order. With `USE_FIRST` it is the running `MutableFeature`: the FIRST interval's start, the
//! largest end, and the first name. The starts differ only when a later interval that overlaps
//! starts before the first, which a file in coordinate order never has.

/// One interval as `Interval` holds it: always on the positive strand here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
}

/// One VCF record, reduced to what `toIntervals` reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub contig: String,
    pub start: i64,
    /// `getAttributeAsInt(END, getEnd())`.
    pub end: i64,
    /// `getID()`, which is `.` when the column was.
    pub id: String,
    pub filtered: bool,
}

/// `VCFFileReader.toIntervals(reader, includeFiltered)`.
pub fn to_intervals(sites: &[Site], include_filtered: bool) -> Vec<Interval> {
    let mut unnamed = 0;
    sites
        .iter()
        .filter(|site| include_filtered || !site.filtered)
        .map(|site| {
            let name = if site.id == "." {
                unnamed += 1;
                format!("interval-{unnamed}")
            } else {
                site.id.clone()
            };
            Interval {
                contig: site.contig.clone(),
                start: site.start,
                end: site.end,
                name,
            }
        })
        .collect()
}

/// `CoordMath.overlaps`, by way of `Locatable.withinDistanceOf`.
fn within_distance(current: &Interval, next: &Interval, distance: i64) -> bool {
    let (start, end) = (current.start, current.end);
    let (start2, end2) = (next.start - distance, next.end + distance);
    current.contig == next.contig
        && ((start2 >= start && start2 <= end)
            || (end2 >= start && end2 <= end)
            || (start >= start2 && end <= end2))
}

/// `IntervalList.merge(intervals, concatenateNames)`.
fn merge(intervals: &[Interval], concatenate_names: bool) -> Interval {
    let first = &intervals[0];
    let mut start = first.start;
    let mut end = first.end;
    let mut names: Vec<&str> = Vec::new();
    for interval in intervals {
        if !names.contains(&interval.name.as_str()) {
            names.push(&interval.name);
        }
        start = start.min(interval.start);
        end = end.max(interval.end);
    }
    let name = if concatenate_names {
        names.join("|")
    } else {
        names[0].to_string()
    };
    Interval {
        contig: first.contig.clone(),
        start,
        end,
        name,
    }
}

/// `new IntervalMergerIterator(intervals, true, false, concatenateNames)`, drained.
///
/// Abutting intervals are combined (`combineAbuttingIntervals` is true), and strands are not
/// enforced, which does not matter here: every interval is on the positive strand.
pub fn merge_intervals(intervals: &[Interval], concatenate_names: bool) -> Vec<Interval> {
    let mut out = Vec::new();
    let mut current: Option<Interval> = None;
    let mut to_be_merged: Vec<Interval> = Vec::new();
    for next in intervals {
        match &mut current {
            None => {
                current = Some(next.clone());
                to_be_merged.push(next.clone());
            }
            Some(building)
                if within_distance(building, next, 0) || within_distance(building, next, 1) =>
            {
                to_be_merged.push(next.clone());
                building.end = building.end.max(next.end);
            }
            Some(building) => {
                out.push(if concatenate_names {
                    merge(&to_be_merged, true)
                } else {
                    building.clone()
                });
                to_be_merged.clear();
                to_be_merged.push(next.clone());
                *building = next.clone();
            }
        }
    }
    if let Some(building) = current {
        out.push(if concatenate_names {
            merge(&to_be_merged, true)
        } else {
            building
        });
    }
    out
}

/// `new SAMFileHeader(dictionary)` through `SAMTextHeaderCodec.encode`, then one line per
/// interval as `IntervalListWriter.write` puts it.
///
/// The header is `@HD VN:1.6` alone, with no sort order, and one `@SQ` per contig carrying the
/// only attribute a VCF contig line gives a `SAMSequenceRecord`: its assembly, as `AS`.
pub fn render(dictionary: &[crate::vcf_io::Sequence], intervals: &[Interval]) -> String {
    let mut out = String::from("@HD\tVN:1.6\n");
    for sequence in dictionary {
        out.push_str(&format!(
            "@SQ\tSN:{}\tLN:{}",
            sequence.name, sequence.length
        ));
        if let Some(assembly) = &sequence.assembly {
            out.push_str(&format!("\tAS:{assembly}"));
        }
        out.push('\n');
    }
    for interval in intervals {
        out.push_str(&format!(
            "{}\t{}\t{}\t+\t{}\n",
            interval.contig, interval.start, interval.end, interval.name
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(contig: &str, start: i64, end: i64, id: &str, filtered: bool) -> Site {
        Site {
            contig: contig.to_string(),
            start,
            end,
            id: id.to_string(),
            filtered,
        }
    }

    #[test]
    fn a_filtered_site_takes_a_number_only_when_it_is_kept() {
        let sites = [
            site("chr1", 10, 10, ".", true),
            site("chr1", 50, 50, ".", false),
        ];
        assert_eq!(to_intervals(&sites, false)[0].name, "interval-1");
        assert_eq!(to_intervals(&sites, true)[1].name, "interval-2");
    }

    #[test]
    fn abutting_sites_merge_and_their_names_join_in_order() {
        let intervals = to_intervals(
            &[
                site("chr1", 100, 100, "rs1", false),
                site("chr1", 101, 101, ".", false),
                site("chr1", 103, 103, "rs2", false),
            ],
            false,
        );
        let merged = merge_intervals(&intervals, true);
        assert_eq!(merged.len(), 2);
        assert_eq!((merged[0].start, merged[0].end), (100, 101));
        assert_eq!(merged[0].name, "rs1|interval-1");
        assert_eq!(merge_intervals(&intervals, false)[0].name, "rs1");
    }

    #[test]
    fn the_file_order_is_kept_and_a_site_elsewhere_breaks_a_merge() {
        let intervals = to_intervals(
            &[
                site("chr1", 100, 100, "a", false),
                site("chr2", 5, 5, "b", false),
                site("chr1", 100, 100, "c", false),
            ],
            false,
        );
        let names: Vec<String> = merge_intervals(&intervals, true)
            .into_iter()
            .map(|i| i.name)
            .collect();
        assert_eq!(names, ["a", "b", "c"]);
    }
}
