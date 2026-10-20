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
//! Nothing is held whole. The FASTA is scanned once with [`htsjdk_bam::fasta::read_fasta`]'s line
//! rules, keeping only where each contig's base lines sit in the file (one [`Run`] per stretch of
//! equally long lines, so a regular FASTA costs one or two per contig), and each interval is then
//! read back with one `pread` per run it touches. The interval list is read three times: whole, so a
//! failure to read it comes before any FASTA or parse refusal as the whole-file read had it; then
//! checked against the contigs, so a refusal leaves no output; then written record by record.
//!
//! Scope of this slice: the extraction itself. `assertSequenceDictionariesEqual` (interval-list dict
//! vs reference dict) is a validation gate that passes for well-formed inputs and does not affect the
//! bytes written; it is not reproduced here.

use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Read, Write};

use htsjdk_bam::fasta::FastaError;
use htsjdk_bam::interval::{Interval, IntervalList, ParseError};
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
    /// An input could not be read, or the output written.
    Io(io::Error),
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

impl From<io::Error> for ExtractError {
    fn from(e: io::Error) -> Self {
        ExtractError::Io(e)
    }
}

/// A source read at an offset, so each run costs one `pread` rather than a seek and a read.
pub trait ReadAt {
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()>;
}

impl ReadAt for std::fs::File {
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()> {
        std::os::unix::fs::FileExt::read_exact_at(self, buf, offset)
    }
}

impl ReadAt for io::Cursor<&[u8]> {
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()> {
        let bytes = self.get_ref();
        let source = usize::try_from(offset)
            .ok()
            .and_then(|at| bytes.get(at..at.checked_add(buf.len())?))
            .ok_or(io::ErrorKind::UnexpectedEof)?;
        buf.copy_from_slice(source);
        Ok(())
    }
}

/// `count` base lines of `line` bases each, the first at byte `file`, the next ones `stride` bytes
/// apart, holding the contig's bases from `pos` on.
struct Run {
    file: u64,
    pos: u64,
    line: u64,
    stride: u64,
    count: u64,
}

/// Where a contig's bases sit in the FASTA.
#[derive(Default)]
struct Contig {
    len: u64,
    runs: Vec<Run>,
}

impl Contig {
    fn push(&mut self, file: u64, line: u64) {
        if let Some(run) = self.runs.last_mut() {
            if run.line == line && (run.count == 1 || file == run.file + run.count * run.stride) {
                if run.count == 1 {
                    run.stride = file - run.file;
                }
                run.count += 1;
                self.len += line;
                return;
            }
        }
        self.runs.push(Run {
            file,
            pos: self.len,
            line,
            stride: 0,
            count: 1,
        });
        self.len += line;
    }

    /// The bases of `[start, end)`, 0-based, read from `fasta` into `bases`.
    fn read<F: ReadAt>(
        &self,
        fasta: &F,
        start: u64,
        end: u64,
        bases: &mut Vec<u8>,
        span: &mut Vec<u8>,
    ) -> io::Result<()> {
        bases.clear();
        let mut pos = start;
        let mut r = self
            .runs
            .partition_point(|run| run.pos + run.line * run.count <= start);
        while pos < end {
            let run = &self.runs[r];
            let last = end.min(run.pos + run.line * run.count) - 1;
            let at = |p: u64| {
                run.file + (p - run.pos) / run.line * run.stride + (p - run.pos) % run.line
            };
            let from = at(pos);
            span.resize((at(last) + 1 - from) as usize, 0);
            fasta.read_exact_at(span, from)?;
            let mut col = (pos - run.pos) % run.line;
            let mut i = 0usize;
            while pos <= last {
                let take = (run.line - col).min(last + 1 - pos);
                bases.extend_from_slice(&span[i..i + take as usize]);
                pos += take;
                i += take as usize;
                if pos <= last {
                    i += (run.stride - run.line) as usize;
                }
                col = 0;
            }
            r += 1;
        }
        Ok(())
    }
}

/// `read_fasta`'s line rules over a stream, keeping each contig's layout instead of its bases. A
/// format refusal waits for the end of the file, since the whole-file read failed on an
/// unreadable or non-UTF-8 byte anywhere before parsing began.
fn scan_fasta<F: Read>(fasta: F) -> Result<HashMap<String, Contig>, ExtractError> {
    let mut reader = BufReader::with_capacity(1 << 16, fasta);
    let mut contigs = HashMap::new();
    let mut current: Option<(String, Contig)> = None;
    let mut refusal = None;
    let mut buffer = String::new();
    let mut offset = 0u64;
    loop {
        buffer.clear();
        let n = reader.read_line(&mut buffer)?;
        if n == 0 {
            break;
        }
        let start = offset;
        offset += n as u64;
        if refusal.is_some() {
            continue;
        }
        let line = match buffer.strip_suffix('\n') {
            Some(line) => line.strip_suffix('\r').unwrap_or(line),
            None => &buffer,
        };
        if line.trim().is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            let name = header
                .trim()
                .split(char::is_whitespace)
                .next()
                .unwrap_or("");
            if name.is_empty() {
                refusal = Some(FastaError::MissingName);
                continue;
            }
            // A repeated name replaces the earlier contig, as the map insert in file order did.
            if let Some((name, contig)) = current.take() {
                contigs.insert(name, contig);
            }
            current = Some((name.to_string(), Contig::default()));
        } else {
            match current.as_mut() {
                Some((_, contig)) => contig.push(start, line.trim_end().len() as u64),
                None => refusal = Some(FastaError::ExpectedHeader(line.as_bytes()[0])),
            }
        }
    }
    if let Some(refusal) = refusal {
        return Err(refusal.into());
    }
    if let Some((name, contig)) = current {
        contigs.insert(name, contig);
    }
    Ok(contigs)
}

/// The list's intervals in file order, each line parsed as `IntervalList::parse_body` parses it.
fn intervals<I: BufRead>(source: I) -> impl Iterator<Item = Result<Interval, ExtractError>> {
    source.lines().filter_map(|line| match line {
        Err(e) => Some(Err(e.into())),
        Ok(line) => match IntervalList::parse_body(Vec::new(), &line) {
            Err(e) => Some(Err(e.into())),
            Ok(list) => list.intervals.into_iter().next().map(Ok),
        },
    })
}

/// The checks `getSubsequenceAt` and the slice made, in their order, for one interval.
fn check<'a>(
    contigs: &'a HashMap<String, Contig>,
    interval: &Interval,
    opts: &Options,
) -> Result<&'a Contig, ExtractError> {
    let contig = contigs
        .get(&interval.contig)
        .ok_or_else(|| ExtractError::UnknownContig(interval.contig.clone()))?;
    let (start, end) = (interval.start, interval.end);
    if start < 1 || (end as usize as u64) > contig.len {
        return Err(ExtractError::IntervalOutOfRange {
            contig: interval.contig.clone(),
            start,
            end,
        });
    }
    // The whole-file port panicked on these before writing anything; panic before the output too.
    assert!(
        start - 1 <= end,
        "slice index starts at {} but ends at {end}",
        start - 1
    );
    assert!(
        opts.line_length > 0 || end - start < 1,
        "attempt to calculate the remainder with a divisor of zero"
    );
    Ok(contig)
}

/// Writes `bases` as the whole-file port's `String` did: one `char` per byte, so a byte from 0x80
/// up becomes two in UTF-8.
fn write_bases<W: Write>(out: &mut W, bases: &[u8]) -> io::Result<()> {
    if bases.is_ascii() {
        return out.write_all(bases);
    }
    for &base in bases {
        let mut utf8 = [0u8; 2];
        out.write_all((base as char).encode_utf8(&mut utf8).as_bytes())?;
    }
    Ok(())
}

/// `ExtractSequences.doWork`, streaming: `open_intervals` is called once per pass over the list,
/// and the output is created only once every interval has passed its checks.
pub fn extract_sequences_to<I: BufRead, F: Read + ReadAt, W: Write>(
    mut open_intervals: impl FnMut() -> io::Result<I>,
    open_fasta: impl FnOnce() -> io::Result<F>,
    create_output: impl FnOnce() -> io::Result<W>,
    opts: &Options,
) -> Result<W, ExtractError> {
    // shortcut: the list is opened three times and the FASTA is read at offsets, so a pipe
    // (`/dev/stdin`) works for neither; spill the list to a temp file if a pipe is ever needed.
    let mut parse_error = None;
    for interval in intervals(open_intervals()?) {
        match interval {
            Err(ExtractError::Interval(e)) => {
                parse_error.get_or_insert(e);
            }
            Err(e) => return Err(e),
            Ok(_) => {}
        }
    }
    let mut fasta = open_fasta()?;
    let contigs = scan_fasta(&mut fasta)?;
    if let Some(e) = parse_error {
        return Err(e.into());
    }

    for interval in intervals(open_intervals()?) {
        check(&contigs, &interval?, opts)?;
    }

    let mut out = create_output()?;
    let (mut bases, mut span) = (Vec::new(), Vec::new());
    for interval in intervals(open_intervals()?) {
        let interval = interval?;
        let contig = check(&contigs, &interval, opts)?;
        contig.read(
            &fasta,
            interval.start as u64 - 1,
            interval.end as u64,
            &mut bases,
            &mut span,
        )?;
        if interval.negative_strand {
            reverse_complement(&mut bases);
        }

        out.write_all(b">")?;
        // A null name reads back from `.` in the file; write it back as `.`.
        out.write_all(interval.name.as_deref().unwrap_or(".").as_bytes())?;
        out.write_all(b"\n")?;
        let width = opts.line_length.max(1);
        for (i, line) in bases.chunks(width).enumerate() {
            if i > 0 {
                out.write_all(b"\n")?;
            }
            write_bases(&mut out, line)?;
        }
        out.write_all(b"\n")?;
    }
    Ok(out)
}

/// `ExtractSequences.doWork` over an interval list and a reference FASTA, returning the FASTA text.
pub fn extract_sequences(
    interval_list: &str,
    reference_fasta: &str,
    opts: &Options,
) -> Result<String, ExtractError> {
    let out = extract_sequences_to(
        || Ok(interval_list.as_bytes()),
        || Ok(io::Cursor::new(reference_fasta.as_bytes())),
        || Ok(Vec::new()),
        opts,
    )?;
    Ok(String::from_utf8(out).expect("every byte is written as a char"))
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

    #[test]
    fn irregular_lines_crlf_and_blank_lines_read_as_the_whole_file_did() {
        // Lines of 3, 3, 2, a blank, 4 and 1 bases, CRLF endings and trailing spaces: every
        // interval over chr1 must match a slice of the concatenated bases.
        let reference = ">chr1 desc\r\nACG\r\nTAC  \r\nGT\r\n\r\nacgt\r\nN\n>chr2\nTT\n";
        let bases = b"ACGTACGTacgtN";
        let il = (1..=13)
            .flat_map(|s| (s - 1..=13).map(move |e| format!("chr1\t{s}\t{e}\t+\tx\n")))
            .collect::<String>();
        let out = extract_sequences(&interval_list(&il), reference, &Options::default()).unwrap();
        let want = (1..=13usize)
            .flat_map(|s| (s - 1..=13).map(move |e| (s, e)))
            .map(|(s, e)| format!(">x\n{}\n", std::str::from_utf8(&bases[s - 1..e]).unwrap()))
            .collect::<String>();
        assert_eq!(out, want);
    }

    #[test]
    fn a_parse_error_anywhere_wins_over_an_unknown_contig_before_it() {
        let il = interval_list("chrX\t1\t4\t+\tx\nchr1\t1\n");
        assert!(matches!(
            extract_sequences(&il, REFERENCE, &Options::default()),
            Err(ExtractError::Interval(ParseError::TooFewFields(_)))
        ));
    }
}
