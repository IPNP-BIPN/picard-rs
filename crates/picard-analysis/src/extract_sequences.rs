//! `ExtractSequences`.
//!
//! Ports `picard.reference.ExtractSequences.doWork` at tag 3.4.0: read an interval list, pull the
//! matching sub-sequence out of the reference FASTA for each interval, reverse-complement it when the
//! interval is on the negative strand, and write every sub-sequence as its own FASTA record named by
//! the interval, wrapped at `LINE_LENGTH` (default 80).
//!
//! `getSubsequenceAt(contig, start, end)` uses 1-based inclusive coordinates, so an interval
//! `start..=end` selects `bases[start - 1 ..= end - 1]`. The bases are taken with their case
//! preserved (as [`htsjdk_bam::fasta`] documents), and `SequenceUtil.reverseComplement`
//! complements in place, which [`htsjdk_bam::sequence::reverse_complement`] reproduces including its
//! case handling. Each record is `>name`, a newline, then the bases with a newline inserted every
//! `LINE_LENGTH` bases (so a multiple-of-`LINE_LENGTH` length does not get a trailing blank line),
//! then a closing newline.
//!
//! Scope of this slice: the extraction itself. `assertSequenceDictionariesEqual` (interval-list dict
//! vs reference dict) is a validation gate that passes for well-formed inputs and does not affect the
//! bytes written; it is not reproduced here.

use std::collections::HashMap;

use htsjdk_bam::fasta::{read_fasta, FastaError};
use htsjdk_bam::interval::{IntervalList, ParseError};
use htsjdk_bam::sequence::reverse_complement;

/// `ExtractSequences`'s options.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `LINE_LENGTH`, the wrapped line width.
    pub line_length: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options { line_length: 80 }
    }
}

/// Why the extraction could not run.
#[derive(Debug)]
pub enum ExtractError {
    /// The reference FASTA could not be read.
    Fasta(FastaError),
    /// The interval list could not be parsed.
    Interval(ParseError),
    /// An interval named a contig absent from the reference (`getSubsequenceAt` would throw).
    UnknownContig(String),
    /// An interval ran past the end of its contig (`getSubsequenceAt` would throw).
    IntervalOutOfRange {
        contig: String,
        start: i32,
        end: i32,
    },
}

impl From<FastaError> for ExtractError {
    fn from(e: FastaError) -> Self {
        ExtractError::Fasta(e)
    }
}

impl From<ParseError> for ExtractError {
    fn from(e: ParseError) -> Self {
        ExtractError::Interval(e)
    }
}

/// `ExtractSequences.doWork` over an interval list and a reference FASTA, returning the FASTA text.
pub fn extract_sequences(
    interval_list: &str,
    reference_fasta: &str,
    opts: &Options,
) -> Result<String, ExtractError> {
    let mut bases_by_contig: HashMap<String, Vec<u8>> = HashMap::new();
    for sequence in read_fasta(reference_fasta.as_bytes())? {
        bases_by_contig.insert(sequence.name, sequence.bases);
    }

    // The dictionary is only used by htsjdk to validate the list against the reference; extraction
    // reads the interval fields directly, so an empty dictionary is enough here.
    let intervals = IntervalList::parse_body(Vec::new(), interval_list)?;

    let mut out = String::new();
    for interval in &intervals.intervals {
        let bases = bases_by_contig
            .get(&interval.contig)
            .ok_or_else(|| ExtractError::UnknownContig(interval.contig.clone()))?;

        // 1-based inclusive; `getSubsequenceAt` returns `bases[start - 1 ..= end - 1]`.
        let start = interval.start;
        let end = interval.end;
        if start < 1 || (end as usize) > bases.len() {
            return Err(ExtractError::IntervalOutOfRange {
                contig: interval.contig.clone(),
                start,
                end,
            });
        }
        let mut sub = bases[(start as usize - 1)..(end as usize)].to_vec();
        if interval.negative_strand {
            reverse_complement(&mut sub);
        }

        out.push('>');
        // A null name reads back from `.` in the file; write it back as `.`.
        out.push_str(interval.name.as_deref().unwrap_or("."));
        out.push('\n');
        for (i, &base) in sub.iter().enumerate() {
            if i > 0 && i % opts.line_length == 0 {
                out.push('\n');
            }
            out.push(base as char);
        }
        out.push('\n');
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    // chr1 is 40 bases; chr2 is 12 bases.
    const REFERENCE: &str =
        ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nAAAACCCCGGGG\n";

    fn interval_list(body: &str) -> String {
        format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n@SQ\tSN:chr2\tLN:12\n{body}")
    }

    #[test]
    fn positive_strand_extracts_verbatim_and_wraps() {
        let il = interval_list("chr2\t1\t12\t+\ttwelve\n");
        let opts = Options { line_length: 4 };
        let out = extract_sequences(&il, REFERENCE, &opts).unwrap();
        assert_eq!(out, ">twelve\nAAAA\nCCCC\nGGGG\n");
    }

    #[test]
    fn negative_strand_reverse_complements() {
        // chr2[1..=12] = AAAACCCCGGGG; reverse complement = CCCCGGGGTTTT.
        let il = interval_list("chr2\t1\t12\t-\trc\n");
        let out = extract_sequences(&il, REFERENCE, &Options::default()).unwrap();
        assert_eq!(out, ">rc\nCCCCGGGGTTTT\n");
    }

    #[test]
    fn a_length_that_is_a_multiple_of_line_length_has_no_trailing_blank_line() {
        // 8 bases at LINE_LENGTH 4 wraps once, with no newline after the second full line.
        let il = interval_list("chr1\t1\t8\t+\teight\n");
        let opts = Options { line_length: 4 };
        let out = extract_sequences(&il, REFERENCE, &opts).unwrap();
        assert_eq!(out, ">eight\nACGT\nACGT\n");
    }

    #[test]
    fn several_intervals_become_several_records() {
        let il = interval_list("chr1\t1\t4\t+\ta\nchr2\t5\t8\t+\tb\n");
        let out = extract_sequences(&il, REFERENCE, &Options::default()).unwrap();
        assert_eq!(out, ">a\nACGT\n>b\nCCCC\n");
    }

    #[test]
    fn an_unknown_contig_is_an_error() {
        let il = interval_list("chrX\t1\t4\t+\tx\n");
        assert!(matches!(
            extract_sequences(&il, REFERENCE, &Options::default()),
            Err(ExtractError::UnknownContig(_))
        ));
    }
}
