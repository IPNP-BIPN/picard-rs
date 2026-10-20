//! `NormalizeFasta`.
//!
//! Ports `picard.reference.NormalizeFasta.doWork` at tag 3.4.0: read a FASTA and rewrite it so every
//! sequence's lines are the same length (`LINE_LENGTH`, default 100), except the last line of each
//! sequence. Optionally the record names are truncated at the first whitespace
//! (`TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE`, default off, so `>chr1 human` stays `chr1 human`).
//!
//! Parsing follows `FastaSequenceFile` (as [`htsjdk_bam::fasta`] documents): bases are concatenated
//! with each line's trailing whitespace trimmed and **case preserved** (soft-masked lowercase stays
//! lowercase), blank lines are skipped, and the name is the header line trimmed (then truncated at
//! whitespace only when requested). Each sequence is written as `>name`, then its bases wrapped at
//! `LINE_LENGTH`, then a newline; a zero-base sequence writes only its `>name` line.

use crate::fasta_records::records;
use htsjdk_bam::fasta::FastaError;
use std::io::{Read, Write};

/// `NormalizeFasta`'s options.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `LINE_LENGTH`, the wrapped line width.
    pub line_length: usize,
    /// `TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE`.
    pub truncate_names_at_whitespace: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            line_length: 100,
            truncate_names_at_whitespace: false,
        }
    }
}

/// Why the FASTA could not be normalized.
#[derive(Debug, PartialEq, Eq)]
pub enum NormalizeError {
    /// A base line appeared before any `>` header (`Expected > but saw ...`).
    ExpectedHeader,
    /// A header carried no name (`Missing sequence name in FASTA`).
    MissingName,
    /// The file could not be read or the output written.
    Io(String),
}

impl From<FastaError> for NormalizeError {
    fn from(error: FastaError) -> Self {
        match error {
            FastaError::ExpectedHeader(_) => NormalizeError::ExpectedHeader,
            FastaError::MissingName => NormalizeError::MissingName,
            FastaError::Io(e) => NormalizeError::Io(e.to_string()),
        }
    }
}

/// `NormalizeFasta.doWork` over FASTA text, returning the normalized FASTA text.
pub fn normalize_fasta(input: &str, opts: &Options) -> Result<String, NormalizeError> {
    let mut out = Vec::new();
    normalize_fasta_to(input.as_bytes(), opts, &mut out)?;
    Ok(String::from_utf8(out).expect("FASTA bases are ASCII"))
}

/// [`normalize_fasta`] one sequence at a time, from a reader to a writer.
pub fn normalize_fasta_to(
    input: impl Read,
    opts: &Options,
    out: &mut impl Write,
) -> Result<(), NormalizeError> {
    let io = |e: std::io::Error| NormalizeError::Io(e.to_string());
    for seq in records(input, opts.truncate_names_at_whitespace) {
        let seq = seq?;
        writeln!(out, ">{}", seq.name).map_err(io)?;
        if !seq.bases.is_empty() {
            for chunk in seq.bases.chunks(opts.line_length) {
                out.write_all(chunk).map_err(io)?;
                out.write_all(b"\n").map_err(io)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = ">chr1 a description here\nACGTacgtACGTNNNNacgt\nACGTACGT\n>chr2\nTTTTTTTTTTGGGGGGGGGGCCCCC\n>empty\n";

    #[test]
    fn default_keeps_the_full_name_preserves_case_and_wraps_at_100() {
        let out = normalize_fasta(INPUT, &Options::default()).unwrap();
        assert_eq!(
            out,
            ">chr1 a description here\nACGTacgtACGTNNNNacgtACGTACGT\n\
             >chr2\nTTTTTTTTTTGGGGGGGGGGCCCCC\n\
             >empty\n"
        );
    }

    #[test]
    fn a_short_line_length_wraps_each_sequence() {
        let opts = Options {
            line_length: 10,
            ..Options::default()
        };
        let out = normalize_fasta(INPUT, &opts).unwrap();
        assert_eq!(
            out,
            ">chr1 a description here\nACGTacgtAC\nGTNNNNacgt\nACGTACGT\n\
             >chr2\nTTTTTTTTTT\nGGGGGGGGGG\nCCCCC\n\
             >empty\n"
        );
    }

    #[test]
    fn truncating_names_stops_at_the_first_whitespace() {
        let opts = Options {
            truncate_names_at_whitespace: true,
            ..Options::default()
        };
        let out = normalize_fasta(INPUT, &opts).unwrap();
        assert!(out.starts_with(">chr1\n"), "{out}");
    }

    #[test]
    fn a_base_before_a_header_is_an_error() {
        assert_eq!(
            normalize_fasta("ACGT\n", &Options::default()),
            Err(NormalizeError::ExpectedHeader)
        );
    }
}
