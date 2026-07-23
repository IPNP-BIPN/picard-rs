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
}

/// `SAMSequenceRecord.truncateSequenceName`: everything up to the first whitespace.
fn truncate_at_whitespace(name: &str) -> &str {
    match name.find(char::is_whitespace) {
        Some(i) => &name[..i],
        None => name,
    }
}

struct Sequence {
    name: String,
    bases: Vec<u8>,
}

/// `NormalizeFasta.doWork` over FASTA text, returning the normalized FASTA text.
pub fn normalize_fasta(input: &str, opts: &Options) -> Result<String, NormalizeError> {
    let mut sequences: Vec<Sequence> = Vec::new();
    for line in input.lines() {
        // Blank lines between records are skipped (`skipNewlines`), not treated as bases.
        if line.trim().is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            let trimmed = header.trim();
            let name = if opts.truncate_names_at_whitespace {
                truncate_at_whitespace(trimmed)
            } else {
                trimmed
            };
            if name.is_empty() {
                return Err(NormalizeError::MissingName);
            }
            sequences.push(Sequence {
                name: name.to_string(),
                bases: Vec::new(),
            });
        } else {
            match sequences.last_mut() {
                Some(seq) => seq.bases.extend_from_slice(line.trim_end().as_bytes()),
                None => return Err(NormalizeError::ExpectedHeader),
            }
        }
    }

    let mut out = String::new();
    for seq in &sequences {
        out.push('>');
        out.push_str(&seq.name);
        out.push('\n');
        if !seq.bases.is_empty() {
            for (i, chunk) in seq.bases.chunks(opts.line_length).enumerate() {
                if i > 0 {
                    out.push('\n');
                }
                // Bases are ASCII, so the chunk is valid UTF-8.
                out.push_str(std::str::from_utf8(chunk).expect("FASTA bases are ASCII"));
            }
            out.push('\n');
        }
    }
    Ok(out)
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
