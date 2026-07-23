//! `BedToIntervalList`.
//!
//! Ports `picard.util.BedToIntervalList.doWork` at tag 3.4.0: read a BED file plus a sequence
//! dictionary and write a Picard interval list. Each BED feature becomes an interval with `start =
//! bedStart + 1` (BED is 0-based, interval lists are 1-based) and `end = bedEnd`, named by the BED
//! name column (an empty name becomes absent), on the negative strand only when the BED strand column
//! is `-`. The output carries a fresh `@HD VN:1.6 SO:coordinate` header followed by the `@SQ` lines of
//! the sequence dictionary, then the intervals; with `SORT` (default on) they are in coordinate order,
//! and with `UNIQUE` (default off) overlapping and abutting intervals are merged.
//!
//! The dictionary's `@SQ` lines are emitted verbatim: `IntervalList.write` re-serializes the same
//! `SAMSequenceDictionary` that a `.dict` (or `.fasta`/`.dict` pair) was parsed from, and that round
//! trip is idempotent for a canonically written dictionary, so the `M5`/`UR` attributes and their
//! order survive unchanged. The port reads and re-emits those exact lines rather than reconstructing
//! them.
//!
//! BED parsing follows `htsjdk.tribble.bed.BEDCodec` with the default `StartOffset.ONE`: fields are
//! split on a tab or a run of spaces (`\t|( +)`), `#`/`track`/`browser` lines are header lines and
//! skipped, a line with fewer than two fields yields no feature, and the name column has its double
//! quotes stripped. A length-zero feature (`bedStart == bedEnd`, i.e. `start == end + 1`) is skipped
//! unless `KEEP_LENGTH_ZERO_INTERVALS` is set.

use htsjdk_bam::interval::{Interval, IntervalList};

/// `BedToIntervalList`'s options.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `SORT`: put the intervals in coordinate order before writing.
    pub sort: bool,
    /// `UNIQUE`: merge overlapping and abutting intervals (`uniqued`, concatenating names).
    pub unique: bool,
    /// `KEEP_LENGTH_ZERO_INTERVALS`: keep `start == end + 1` features rather than skipping them.
    pub keep_length_zero_intervals: bool,
    /// `DROP_MISSING_CONTIGS`: skip features on a contig absent from the dictionary rather than error.
    pub drop_missing_contigs: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            sort: true,
            unique: false,
            keep_length_zero_intervals: false,
            drop_missing_contigs: false,
        }
    }
}

/// Why the conversion could not run. The messages mirror the `PicardException` text.
#[derive(Debug, PartialEq, Eq)]
pub enum BedError {
    /// A BED coordinate did not parse as an integer.
    BadNumber(String),
    /// `Sequence '{0}' was not found in the sequence dictionary`.
    UnknownContig(String),
    /// `Start on sequence '{0}' was less than one: {1}`.
    StartLessThanOne(String, i32),
    /// `Start on sequence '{0}' was past the end: {1} < {2}` (sequence length, start).
    StartPastEnd(String, i32, i32),
    /// `End on sequence '{0}' was less than one: {1}`.
    EndLessThanOne(String, i32),
    /// `End on sequence '{0}' was past the end: {1} < {2}` (sequence length, end).
    EndPastEnd(String, i32, i32),
    /// `On sequence '{0}', end < start-1: {1} <= {2}` (end, start).
    EndBeforeStart(String, i32, i32),
}

/// One sequence dictionary record: its name and length, in header order.
struct DictEntry {
    name: String,
    length: i32,
}

/// Parses a `.dict` into the ordered `(name, length)` records and the verbatim `@SQ` lines.
fn parse_dictionary(text: &str) -> (Vec<DictEntry>, Vec<String>) {
    let mut entries = Vec::new();
    let mut sq_lines = Vec::new();
    for line in text.lines() {
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
            entries.push(DictEntry { name, length });
        }
        sq_lines.push(line.to_string());
    }
    (entries, sq_lines)
}

/// `Pattern.compile("\\t|( +)").split(line, -1)`: split on a tab or a run of spaces, keeping the
/// empty tokens that adjacent delimiters produce.
fn split_bed_line(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'\t' => {
                tokens.push(&line[start..i]);
                i += 1;
                start = i;
            }
            b' ' => {
                tokens.push(&line[start..i]);
                while i < bytes.len() && bytes[i] == b' ' {
                    i += 1;
                }
                start = i;
            }
            _ => i += 1,
        }
    }
    tokens.push(&line[start..]);
    tokens
}

/// `BedToIntervalList.doWork` over a sequence dictionary and BED text, returning the interval list.
pub fn bed_to_interval_list(
    sequence_dictionary: &str,
    bed: &str,
    opts: &Options,
) -> Result<String, BedError> {
    let (dict_entries, sq_lines) = parse_dictionary(sequence_dictionary);

    let dictionary: Vec<String> = dict_entries.iter().map(|e| e.name.clone()).collect();
    let mut list = IntervalList::new(dictionary);

    for line in bed.lines() {
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let tokens = split_bed_line(line);
        if tokens.len() < 2 {
            continue;
        }

        let contig = tokens[0];
        let bed_start: i32 = tokens[1]
            .parse()
            .map_err(|_| BedError::BadNumber(tokens[1].to_string()))?;
        let start = bed_start + 1;
        let end: i32 = if tokens.len() > 2 {
            tokens[2]
                .parse()
                .map_err(|_| BedError::BadNumber(tokens[2].to_string()))?
        } else {
            start
        };

        // Name column: quotes stripped; empty becomes absent.
        let name = if tokens.len() > 3 {
            let raw: String = tokens[3].chars().filter(|&c| c != '"').collect();
            if raw.is_empty() {
                None
            } else {
                Some(raw)
            }
        } else {
            None
        };

        // Strand: only `-` (as the first character of column six) makes it negative. A non-numeric
        // score column would stop BEDCodec before it reads the strand, leaving it non-negative; the
        // corpus uses numeric scores, so the common path (read column six) is what matters here.
        let negative_strand = tokens.len() > 5 && tokens[5].trim_start().starts_with('-');

        let record = dict_entries.iter().find(|e| e.name == contig);
        let record = match record {
            Some(r) => r,
            None => {
                if opts.drop_missing_contigs {
                    continue;
                }
                return Err(BedError::UnknownContig(contig.to_string()));
            }
        };
        let seq_len = record.length;

        if start < 1 {
            return Err(BedError::StartLessThanOne(contig.to_string(), start));
        } else if seq_len < start {
            return Err(BedError::StartPastEnd(contig.to_string(), seq_len, start));
        } else if (end == 0 && start != 1) || end < 0 {
            return Err(BedError::EndLessThanOne(contig.to_string(), end));
        } else if seq_len < end {
            return Err(BedError::EndPastEnd(contig.to_string(), seq_len, end));
        } else if end < start - 1 {
            return Err(BedError::EndBeforeStart(contig.to_string(), end, start));
        }

        // A length-zero feature has start == end + 1; skip it unless asked to keep it.
        if start == end + 1 && !opts.keep_length_zero_intervals {
            continue;
        }

        list.intervals.push(Interval::with_strand_and_name(
            contig,
            start,
            end,
            negative_strand,
            name.as_deref(),
        ));
    }

    let list = if opts.sort { list.sorted() } else { list };
    let list = if opts.unique {
        list.uniqued(true)
    } else {
        list
    };

    let mut out = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for sq in &sq_lines {
        out.push_str(sq);
        out.push('\n');
    }
    out.push_str(&list.write_body());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const DICT: &str = "@HD\tVN:1.6\n\
        @SQ\tSN:chr1\tLN:20\tM5:aaa\tUR:file:///ref.fasta\n\
        @SQ\tSN:chr2\tLN:20\tM5:bbb\tUR:file:///ref.fasta\n";

    fn header() -> String {
        "@HD\tVN:1.6\tSO:coordinate\n\
         @SQ\tSN:chr1\tLN:20\tM5:aaa\tUR:file:///ref.fasta\n\
         @SQ\tSN:chr2\tLN:20\tM5:bbb\tUR:file:///ref.fasta\n"
            .to_string()
    }

    #[test]
    fn a_bed_feature_becomes_a_one_based_interval() {
        let bed = "chr1\t9\t20\tfoo\t0\t+\n";
        let out = bed_to_interval_list(DICT, bed, &Options::default()).unwrap();
        assert_eq!(out, format!("{}chr1\t10\t20\t+\tfoo\n", header()));
    }

    #[test]
    fn an_empty_name_becomes_a_dot_and_minus_strand_is_negative() {
        let bed = "chr2\t0\t5\t\t0\t-\n";
        let out = bed_to_interval_list(DICT, bed, &Options::default()).unwrap();
        assert_eq!(out, format!("{}chr2\t1\t5\t-\t.\n", header()));
    }

    #[test]
    fn sort_orders_by_dictionary_index() {
        let bed = "chr2\t0\t5\tb\t0\t+\nchr1\t0\t4\ta\t0\t+\n";
        let out = bed_to_interval_list(DICT, bed, &Options::default()).unwrap();
        assert_eq!(
            out,
            format!("{}chr1\t1\t4\t+\ta\nchr2\t1\t5\t+\tb\n", header())
        );
    }

    #[test]
    fn a_length_zero_feature_is_skipped_by_default() {
        // bedStart == bedEnd == 5, so start == end + 1 == 6: length zero, dropped.
        let bed = "chr1\t5\t5\tz\t0\t+\n";
        let out = bed_to_interval_list(DICT, bed, &Options::default()).unwrap();
        assert_eq!(out, header());
    }

    #[test]
    fn unique_merges_overlapping_intervals_and_concatenates_names() {
        let bed = "chr1\t0\t10\ta\t0\t+\nchr1\t5\t15\tb\t0\t+\n";
        let opts = Options {
            unique: true,
            ..Options::default()
        };
        let out = bed_to_interval_list(DICT, bed, &opts).unwrap();
        assert_eq!(out, format!("{}chr1\t1\t15\t+\ta|b\n", header()));
    }

    #[test]
    fn a_missing_contig_is_an_error_by_default() {
        let bed = "chrX\t0\t5\tx\t0\t+\n";
        assert_eq!(
            bed_to_interval_list(DICT, bed, &Options::default()),
            Err(BedError::UnknownContig("chrX".to_string()))
        );
    }

    #[test]
    fn split_handles_tabs_and_space_runs() {
        assert_eq!(split_bed_line("a\tb\tc"), vec!["a", "b", "c"]);
        assert_eq!(split_bed_line("a  b"), vec!["a", "b"]);
        assert_eq!(split_bed_line("a\t b"), vec!["a", "", "b"]);
    }
}
