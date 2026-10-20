//! A FASTA one sequence at a time: `htsjdk_bam::fasta::read_fasta`'s parse, holding only the
//! sequence being read, for tools that visit each contig once.

use htsjdk_bam::fasta::{FastaError, ReferenceSequence};
use std::io::{BufRead, BufReader, Read};

pub struct Records<R> {
    reader: BufReader<R>,
    line: String,
    name: Option<String>,
    index: i32,
    truncate: bool,
    failed: bool,
}

/// The sequences of `source` in file order; `truncate` cuts names at the first whitespace, as
/// `read_fasta` always does. After an error the iterator ends.
pub fn records<R: Read>(source: R, truncate: bool) -> Records<R> {
    Records {
        reader: BufReader::new(source),
        line: String::new(),
        name: None,
        index: -1,
        truncate,
        failed: false,
    }
}

/// The error a whole-file `read_to_string` of `path` gave, which came before any parse error:
/// read again only once the streamed parse has failed.
pub fn unreadable(path: &str) -> Option<std::io::Error> {
    std::fs::read_to_string(path).err()
}

impl<R: Read> Records<R> {
    fn fail(&mut self, error: FastaError) -> Option<Result<ReferenceSequence, FastaError>> {
        self.failed = true;
        Some(Err(error))
    }

    fn sequence(&self, name: String, bases: Vec<u8>) -> ReferenceSequence {
        ReferenceSequence {
            name,
            index: self.index,
            bases,
        }
    }
}

impl<R: Read> Iterator for Records<R> {
    type Item = Result<ReferenceSequence, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        let mut bases = Vec::new();
        loop {
            // `BufRead::lines`, into one buffer: the newline, and a carriage return before it, go.
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => return self.name.take().map(|name| Ok(self.sequence(name, bases))),
                Err(e) => return self.fail(FastaError::Io(e)),
                Ok(_) => {}
            }
            let mut line = self.line.as_str();
            if let Some(rest) = line.strip_suffix('\n') {
                line = rest.strip_suffix('\r').unwrap_or(rest);
            }
            // Blank lines between records are skipped by `skipNewlines`, not treated as data.
            if line.trim().is_empty() {
                continue;
            }
            if let Some(header) = line.strip_prefix('>') {
                let trimmed = header.trim();
                let name = match (self.truncate, trimmed.find(char::is_whitespace)) {
                    (true, Some(i)) => &trimmed[..i],
                    _ => trimmed,
                };
                if name.is_empty() {
                    return self.fail(FastaError::MissingName);
                }
                let previous = self.name.replace(name.to_string());
                if let Some(previous) = previous {
                    let sequence = self.sequence(previous, bases);
                    self.index += 1;
                    return Some(Ok(sequence));
                }
                self.index += 1;
            } else if self.name.is_none() {
                return self.fail(FastaError::ExpectedHeader(line.as_bytes()[0]));
            } else {
                bases.extend_from_slice(line.trim_end().as_bytes());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::fasta::read_fasta;

    #[test]
    fn yields_what_read_fasta_reads() {
        let text = b">c1 desc\nACGT\n\nac\n>c2\n>c3\nNNNN\n";
        let streamed: Vec<_> = records(&text[..], true).map(Result::unwrap).collect();
        assert_eq!(streamed, read_fasta(&text[..]).unwrap());
        let kept: Vec<_> = records(&text[..], false).map(|r| r.unwrap().name).collect();
        assert_eq!(kept, ["c1 desc", "c2", "c3"]);
    }

    #[test]
    fn stops_at_the_error_read_fasta_reports() {
        let errors: Vec<_> = [&b"ACGT\n"[..], b">c1\nA\n>\nC\n"]
            .iter()
            .map(|text| {
                records(*text, true)
                    .find_map(Result::err)
                    .map(|e| format!("{e:?}"))
            })
            .collect();
        let expected: Vec<_> = [&b"ACGT\n"[..], b">c1\nA\n>\nC\n"]
            .iter()
            .map(|text| read_fasta(*text).err().map(|e| format!("{e:?}")))
            .collect();
        assert_eq!(errors, expected);
    }
}
