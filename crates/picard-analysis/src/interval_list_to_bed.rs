//! `IntervalListToBed`.
//!
//! Ports `picard.util.IntervalListToBed.doWork` at tag 3.4.0: read an interval list and write each
//! interval as a BED line. BED is 0-based half-open, so an interval `start..=end` (1-based inclusive)
//! becomes `contig<TAB>start-1<TAB>end<TAB>name<TAB>SCORE<TAB>strand`, where `SCORE` is a constant
//! (default 500) and `strand` is `+` or `-`. When `SORT` is on (the default) the intervals are first
//! put in coordinate order by `IntervalCoordinateComparator` (sequence index from the header
//! dictionary, then start, end, positive-strand-first, then name); otherwise they are written in file
//! order.
//!
//! The header dictionary gives the contig order the sort keys off, so it is read from the `@SQ`
//! `SN:` fields of the list header. `generateOutput` writes each field with `String.valueOf`, so a
//! null name would print `null`; the interval codec never produces a null name (it takes the literal
//! 5th column), and a `.` column reads back through [`htsjdk_bam::interval`] as an absent name that
//! this port prints back as `.`.

use htsjdk_bam::interval::{IntervalList, ParseError};

/// `IntervalListToBed`'s options.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `SCORE`, the constant score written for every interval.
    pub score: i32,
    /// `SORT`: put the intervals in coordinate order before writing.
    pub sort: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            score: 500,
            sort: true,
        }
    }
}

/// The `@SQ` `SN:` contig names, in header order: the dictionary `IntervalCoordinateComparator` uses.
fn dictionary_from_header(interval_list: &str) -> Vec<String> {
    let mut dict = Vec::new();
    for line in interval_list.lines() {
        if !line.starts_with("@SQ") {
            if line.starts_with('@') {
                continue;
            }
            // The header is contiguous and precedes every interval line.
            break;
        }
        if let Some(sn) = line.split('\t').find_map(|field| field.strip_prefix("SN:")) {
            dict.push(sn.to_string());
        }
    }
    dict
}

/// `IntervalListToBed.doWork` over an interval list, returning the BED text.
pub fn interval_list_to_bed(interval_list: &str, opts: &Options) -> Result<String, ParseError> {
    let dictionary = dictionary_from_header(interval_list);
    let list = IntervalList::parse_body(dictionary, interval_list)?;
    let list = if opts.sort { list.sorted() } else { list };

    let mut out = String::new();
    for interval in &list.intervals {
        let strand = if interval.negative_strand { "-" } else { "+" };
        let name = interval.name.as_deref().unwrap_or(".");
        // BED start is 0-based: start - 1.
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            interval.contig,
            interval.start - 1,
            interval.end,
            name,
            opts.score,
            strand
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn list(body: &str) -> String {
        format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:100\n{body}")
    }

    #[test]
    fn each_interval_becomes_a_zero_based_bed_line() {
        let il = list("chr1\t10\t20\t+\ta\n");
        let out = interval_list_to_bed(&il, &Options::default()).unwrap();
        assert_eq!(out, "chr1\t9\t20\ta\t500\t+\n");
    }

    #[test]
    fn sort_orders_by_dictionary_index_not_string() {
        // chr2 appears first in the file but chr1 is dictionary index 0, so it sorts first.
        let il = list("chr2\t5\t8\t+\tb\nchr1\t1\t4\t-\ta\n");
        let out = interval_list_to_bed(&il, &Options::default()).unwrap();
        assert_eq!(out, "chr1\t0\t4\ta\t500\t-\nchr2\t4\t8\tb\t500\t+\n");
    }

    #[test]
    fn sort_off_keeps_file_order() {
        let il = list("chr2\t5\t8\t+\tb\nchr1\t1\t4\t-\ta\n");
        let opts = Options {
            sort: false,
            ..Options::default()
        };
        let out = interval_list_to_bed(&il, &opts).unwrap();
        assert_eq!(out, "chr2\t4\t8\tb\t500\t+\nchr1\t0\t4\ta\t500\t-\n");
    }

    #[test]
    fn score_is_configurable() {
        let il = list("chr1\t1\t1\t+\tp\n");
        let opts = Options {
            score: 1000,
            ..Options::default()
        };
        let out = interval_list_to_bed(&il, &opts).unwrap();
        assert_eq!(out, "chr1\t0\t1\tp\t1000\t+\n");
    }
}
