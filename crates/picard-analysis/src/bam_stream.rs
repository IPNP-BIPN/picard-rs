//! An input read record by record, with the answers of reading it whole.
//!
//! The tools used to read the whole file, inflate all of it and decode every record into a vector
//! before doing anything, which held the inflated stream and every decoded record at once. This
//! module hands the records over one at a time instead, and keeps what that order of work
//! answered: a damaged file is refused with the message the whole-file read gave, and before any
//! record reaches the tool.
//!
//! A BAM is inflated once. A first pass walks the compressed blocks without inflating them and
//! makes the checks htsjdk's reader makes on each, which is where a truncated file or a damaged
//! block is found: such a file is answered by reading it whole, exactly as before, before any
//! record reaches the tool. A file whose blocks are sound is then inflated on several cores with
//! noodles-bgzf and its records handed over as they are decoded. Should that pass meet what the
//! whole-file read would have answered differently (a block that does not inflate, a record that
//! does not decode, an inflated length that is not the sum of the declared ones), the whole-file
//! read gives the answer: its error, raised there, or its records, taken up after the ones
//! already handed over.

use crate::record_validation::Validator;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::{BamReader, ReadError};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_bgzf::BgzfError;
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::num::NonZeroUsize;

const HEADER: usize = 18;
const FOOTER: usize = 8;
const MAX_BLOCK: usize = 64 * 1024;
/// How much more inflated stream a refill asks for at a time.
const CHUNK: usize = 256 * 1024;
/// How many inflated bytes of whole records the first pass hands a thread at a time.
const BATCH: usize = 256 * 1024;

/// Why an input could not be read, as the whole-file read raised it.
#[derive(Debug)]
pub enum InputError {
    /// Opening or reading the file.
    Io(io::Error),
    /// A SAM input that is not UTF-8.
    Utf8(std::string::FromUtf8Error),
    /// [`crate::bgzf::decompress_owned`].
    Bgzf(BgzfError),
    /// [`BamReader::new`], or a record it iterated to.
    Bam(ReadError),
    /// [`htsjdk_bam::sam_file::read_sam`].
    Sam(ParseError),
}

/// The throwable htsjdk raises for the BGZF errors a damaged BAM meets most, out of a reader that
/// `SamReaderFactory` opened on a stream with no name.
fn thrown_message(e: &InputError) -> Option<&'static str> {
    match e {
        InputError::Bgzf(BgzfError::PrematureEnd) => {
            Some("htsjdk.samtools.FileTruncatedException: Premature end of file: data stream")
        }
        InputError::Bgzf(BgzfError::InvalidGzipHeader) => {
            Some("htsjdk.samtools.SAMFormatException: Invalid GZIP header")
        }
        _ => None,
    }
}

/// A read that failed after records were handed over: the throwable htsjdk raises, or the
/// message [`crate::metrics_cli::open_input`] gives.
fn report(e: InputError) -> ! {
    match thrown_message(&e) {
        Some(message) => crate::metrics_cli::thrown(message),
        None => crate::metrics_cli::fail(&e.message()),
    }
}

impl From<InputError> for Box<dyn std::error::Error> {
    /// What the tools' own `main` returned: the I/O and UTF-8 errors as they are, the others as
    /// their `Debug` text.
    fn from(e: InputError) -> Self {
        match e {
            InputError::Io(e) => e.into(),
            InputError::Utf8(e) => e.into(),
            e @ InputError::Bgzf(BgzfError::PrematureEnd | BgzfError::InvalidGzipHeader) => {
                crate::metrics_cli::thrown(thrown_message(&e).unwrap_or_default())
            }
            InputError::Bgzf(e) => format!("{e:?}").into(),
            InputError::Bam(e) => format!("{e:?}").into(),
            InputError::Sam(e) => format!("{e:?}").into(),
        }
    }
}

impl InputError {
    /// What [`crate::metrics_cli::read_input`] failed with: `Display` for I/O and UTF-8, `Debug`
    /// for the rest.
    pub fn message(&self) -> String {
        match self {
            InputError::Io(e) => format!("{e}"),
            InputError::Utf8(e) => format!("{e}"),
            InputError::Bgzf(e) => format!("{e:?}"),
            InputError::Bam(e) => format!("{e:?}"),
            InputError::Sam(e) => format!("{e:?}"),
        }
    }
}

/// One record's line of `write_sam`, without its newline, or `None` for a record the text form
/// cannot carry: a writer can then stream the file `write_sam` builds whole.
pub fn sam_line(header: &SamHeader, record: &BamRecord) -> Option<String> {
    let name_of = |index: i32| -> &str {
        usize::try_from(index)
            .ok()
            .and_then(|index| header.sequences.get(index))
            .map_or("*", |sequence| sequence.name.as_str())
    };
    htsjdk_bam::text::write_alignment(
        record,
        name_of(record.reference_index),
        name_of(record.mate_reference_index),
    )
}

/// An input that has been checked, ready to hand over its records.
pub struct Input {
    pub header: SamHeader,
    records: Records,
}

impl Input {
    /// Reads the input at `path`, BAM or SAM, and answers as reading it whole did; the records of
    /// a BAM that passes are then decoded as they are taken.
    pub fn open(path: &str) -> Result<Input, InputError> {
        Input::open_with(path, ValidationStringency::default())
    }

    /// [`Self::open`], a SAM input parsed at `stringency`, as the tools that pass their
    /// `VALIDATION_STRINGENCY` to `SamReaderFactory` read one.
    pub fn open_with(path: &str, stringency: ValidationStringency) -> Result<Input, InputError> {
        let mut head = Vec::new();
        File::open(path)
            .and_then(|f| f.take(4).read_to_end(&mut head))
            .map_err(InputError::Io)?;
        // A CRAM is decoded whole, against the reference, as htsjdk's reader would hand it over.
        if crate::cram_input::is_cram(&head) {
            let raw = std::fs::read(path).map_err(InputError::Io)?;
            let (header, records) = crate::cram_input::decode_or_exit(&raw);
            return Ok(Input {
                header,
                records: Records::Decoded(records.into_iter(), None),
            });
        }
        head.truncate(2);
        // A BAM begins with a BGZF block, whose first two bytes are gzip's magic; anything else is
        // read as SAM, whole.
        if head != [0x1f, 0x8b] {
            let mut raw = Vec::new();
            File::open(path)
                .and_then(|mut f| f.read_to_end(&mut raw))
                .map_err(InputError::Io)?;
            let text = String::from_utf8(raw).map_err(InputError::Utf8)?;
            let (header, records) = read_sam_with(&text, stringency).map_err(InputError::Sam)?;
            // shortcut: a SAM is not validated; SAMLineParser's checks are their own port.
            return Ok(Input {
                header,
                records: Records::Decoded(records.into_iter(), None),
            });
        }
        if structure_clean(path) == Some(true) {
            stream(path)
        } else {
            whole(path)
        }
    }

    /// Whether the records are decoded from the file as they are taken, rather than held.
    pub fn streams(&self) -> bool {
        matches!(self.records, Records::Stream(..))
    }

    /// The same input again, for a tool that reads it twice, as Picard's do. A file [`Self::open`]
    /// streamed has had its blocks checked and is streamed again without that; anything else is
    /// opened as it was the first time.
    pub fn reopen(path: &str, streamed: bool) -> Result<Input, InputError> {
        if streamed {
            stream(path)
        } else {
            Input::open(path)
        }
    }

    /// The same input with its records not validated, as `BamReader::new` hands them over.
    pub fn unvalidated(mut self) -> Input {
        match &mut self.records {
            Records::Stream(_, validator) | Records::Decoded(_, validator) => *validator = None,
        }
        self
    }

    /// The header, and the records in file order.
    pub fn into_parts(self) -> (SamHeader, Records) {
        (self.header, self.records)
    }
}

impl IntoIterator for Input {
    type Item = BamRecord;
    type IntoIter = Records;

    fn into_iter(self) -> Records {
        self.records
    }
}

/// The records of an [`Input`], a BAM's validated as they are taken.
pub enum Records {
    Stream(Box<Streaming>, Option<Validator>),
    Decoded(std::vec::IntoIter<BamRecord>, Option<Validator>),
}

/// The single inflating pass over a BAM whose blocks were checked: the stream, until it ends,
/// and what is needed to fall back on the whole-file read if it ends otherwise than that read.
pub struct Streaming {
    stream: Option<RecordStream<MultithreadedInflater<Blocks<BufReader<File>>>>>,
    path: String,
    handed: usize,
}

impl Streaming {
    /// The next record, or `Err(())` where the whole-file read would have answered otherwise.
    fn next(&mut self) -> Result<Option<BamRecord>, ()> {
        let Some(stream) = self.stream.as_mut() else {
            return Ok(None);
        };
        match stream.record() {
            Ok(Some(record)) => {
                self.handed += 1;
                Ok(Some(record))
            }
            Ok(None) if self.ended_clean() => Ok(None),
            _ => Err(()),
        }
    }

    /// Whether the stream ended as the whole-file read would have: every block clean, and as many
    /// bytes inflated as the blocks declared.
    fn ended_clean(&mut self) -> bool {
        let Some(stream) = self.stream.take() else {
            return true;
        };
        let mut inflater = stream.into_inner();
        if io::copy(&mut inflater, &mut io::sink()).is_err() {
            return false;
        }
        let inflated = inflater.count;
        let Ok(mut blocks) = inflater.inner.finish() else {
            return false;
        };
        io::copy(&mut blocks, &mut io::sink()).is_ok()
            && !blocks.damaged
            && blocks.declared == inflated
    }
}

impl Iterator for Records {
    type Item = BamRecord;

    fn next(&mut self) -> Option<BamRecord> {
        self.try_next().unwrap_or_else(|e| report(e))
    }
}

impl Records {
    /// The next record, validated, or the error the whole-file read raised where the single pass
    /// met what that read answers otherwise; [`Iterator::next`] reports that error and exits.
    pub fn try_next(&mut self) -> Result<Option<BamRecord>, InputError> {
        let record = match self {
            Records::Decoded(records, _) => records.next(),
            Records::Stream(streaming, validator) => match streaming.next() {
                Ok(record) => record,
                // The whole-file read answers: its error, or its records after the ones already
                // handed over.
                Err(()) => {
                    let handed = streaming.handed;
                    let validator = validator.take();
                    let mut rest = match whole(&streaming.path)?.records {
                        Records::Decoded(rest, _) => rest,
                        Records::Stream(..) => unreachable!("the whole-file read decodes"),
                    };
                    if handed > 0 {
                        rest.nth(handed - 1);
                    }
                    *self = Records::Decoded(rest, validator);
                    return self.try_next();
                }
            },
        };
        let validator = match self {
            Records::Stream(_, validator) | Records::Decoded(_, validator) => validator.as_mut(),
        };
        if let (Some(rec), Some(validator)) = (&record, validator) {
            if let Err(message) = validator.check(rec) {
                crate::metrics_cli::thrown(&message);
            }
        }
        Ok(record)
    }
}

/// The second pass, over a file the first has vouched for.
fn stream(path: &str) -> Result<Input, InputError> {
    let mut stream = RecordStream::new(inflater(Blocks::new(BufReader::new(
        File::open(path).map_err(InputError::Io)?,
    ))));
    // A header that does not decode is answered by the whole-file read, before any record.
    let Ok(header) = stream.header() else {
        return whole(path);
    };
    // Each pass validates, as each of Picard's readers does.
    let validator = Validator::new(&header, &stream.binary_sequences);
    Ok(Input {
        header,
        records: Records::Stream(
            Box::new(Streaming {
                stream: Some(stream),
                path: path.to_string(),
                handed: 0,
            }),
            Some(validator),
        ),
    })
}

/// A clean BAM's records as batches of their bytes, cut at the lengths they declare, for a tool
/// that decodes them on every core: the single pass [`Input::open`] streams, without decoding.
pub struct Batches {
    streaming: Streaming,
    tail: bool,
    /// Whether every block passes htsjdk's checks, found on a thread of its own while the stream
    /// is read: nothing is answered from the batches before the end, so it need not come first.
    structure: Option<std::thread::JoinHandle<Option<bool>>>,
}

/// The header, the checks [`Records`] would run on each record, and the record batches of the BAM
/// at `path`; or `None` where its blocks or its header are what [`Input::open`] answers for.
pub fn batches(path: &str) -> Option<(SamHeader, Validator, Batches)> {
    let file = File::open(path).ok()?;
    let checked = path.to_string();
    let structure = Some(std::thread::spawn(move || structure_clean(&checked)));
    let mut stream = RecordStream::new(inflater(Blocks::new(BufReader::new(file))));
    let header = stream.header().ok()?;
    let validator = Validator::new(&header, &stream.binary_sequences);
    let streaming = Streaming {
        stream: Some(stream),
        path: path.to_string(),
        handed: 0,
    };
    let batches = Batches {
        streaming,
        tail: false,
        structure,
    };
    Some((header, validator, batches))
}

impl Batches {
    /// The next batches, `count` at most, and none once the stream has ended as the whole-file
    /// read would have. `None` where it ends otherwise, which [`Input::open`] answers for.
    pub fn next_round(&mut self, count: usize) -> Option<Vec<Vec<u8>>> {
        if self.tail {
            let clean = self.structure.take().and_then(|s| s.join().ok().flatten());
            return match self.streaming.next() {
                Ok(None) if clean == Some(true) => Some(Vec::new()),
                _ => None,
            };
        }
        let stream = self.streaming.stream.as_mut()?;
        let (batches, tail) = stream.cut(count.max(1)).ok()?;
        self.tail = tail;
        if batches.is_empty() {
            return self.next_round(count);
        }
        Some(batches)
    }
}

/// Whether every block of the BAM at `path` passes the checks htsjdk's reader makes, read without
/// inflating anything; `None` when the file cannot be read.
fn structure_clean(path: &str) -> Option<bool> {
    let mut blocks = Blocks::new(BufReader::new(File::open(path).ok()?));
    io::copy(&mut blocks, &mut io::sink()).ok()?;
    Some(!blocks.damaged)
}

/// The whole-file read, for a BAM the first pass could not vouch for.
fn whole(path: &str) -> Result<Input, InputError> {
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(InputError::Io)?;
    let raw = crate::cram_input::uncram(raw);
    let plain = crate::bgzf::decompress_owned(raw).map_err(InputError::Bgzf)?;
    let reader = BamReader::new(&plain).map_err(InputError::Bam)?;
    let header = reader.header.text.clone();
    let validator = Validator::new(&header, &reader.header.binary_sequences);
    let records = reader
        .collect::<Result<Vec<_>, _>>()
        .map_err(InputError::Bam)?;
    Ok(Input {
        header,
        records: Records::Decoded(records.into_iter(), Some(validator)),
    })
}

/// The header of the BAM at `path`, read as [`BamReader::new`] reads it, when the whole file
/// inflates as cleanly as [`crate::bgzf::decompress_owned`] would have inflated it; `None` for
/// anything else, which the caller answers by reading the file whole, as before. Holds no more
/// than the header and a few blocks, for a tool that only needs the header and copies the blocks.
pub fn header_if_clean(path: &str) -> Option<SamHeader> {
    let file = File::open(path).ok()?;
    let mut stream = RecordStream::new(inflater(Blocks::new(BufReader::new(file))));
    let header = stream.header().ok()?;
    let mut inflater = stream.into_inner();
    io::copy(&mut inflater, &mut io::sink()).ok()?;
    let inflated = inflater.count;
    let mut blocks = inflater.inner.finish().ok()?;
    io::copy(&mut blocks, &mut io::sink()).ok()?;
    (!blocks.damaged && blocks.declared == inflated).then_some(header)
}

/// What noodles reads: the compressed bytes, a block at a time, refusing to go on past the first
/// block that `BgzfReader::next_block` would not have read as a clean one, and adding up the
/// inflated lengths the blocks declare.
pub struct Blocks<R> {
    inner: R,
    block: Vec<u8>,
    served: usize,
    damaged: bool,
    declared: u64,
}

impl<R: Read> Blocks<R> {
    fn new(inner: R) -> Self {
        Blocks {
            inner,
            block: Vec::with_capacity(MAX_BLOCK),
            served: 0,
            damaged: false,
            declared: 0,
        }
    }

    /// The next block into `self.block`, or `false` at the end, clean or damaged.
    fn next_block(&mut self) -> io::Result<bool> {
        self.block.clear();
        self.served = 0;
        if self.damaged {
            return Ok(false);
        }
        let got = read_up_to(&mut self.inner, &mut self.block, HEADER)?;
        if got == 0 {
            return Ok(false);
        }
        if got < HEADER {
            self.damaged = true;
            return Ok(false);
        }
        let h = &self.block[..];
        // The checks `BgzfReader::next_block` makes, and no others.
        let block_length = u16::from_le_bytes([h[16], h[17]]) as usize + 1;
        if h[0] != 0x1f
            || h[1] != 0x8b
            || h[2] != 8
            || h[3] != 4
            || u16::from_le_bytes([h[10], h[11]]) != 6
            || !(HEADER + FOOTER..=MAX_BLOCK).contains(&block_length)
        {
            self.damaged = true;
            return Ok(false);
        }
        let rest = block_length - HEADER;
        if read_up_to(&mut self.inner, &mut self.block, rest)? < rest {
            self.damaged = true;
            return Ok(false);
        }
        let f = &self.block[block_length - FOOTER..];
        self.declared += u64::from(u32::from_le_bytes([f[4], f[5], f[6], f[7]]));
        Ok(true)
    }
}

impl<R: Read> Read for Blocks<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.served == self.block.len() && !self.next_block()? {
            return Ok(0);
        }
        let rest = &self.block[self.served..];
        let n = rest.len().min(buf.len());
        buf[..n].copy_from_slice(&rest[..n]);
        self.served += n;
        Ok(n)
    }
}

/// Appends up to `n` bytes of `inner` to `out`, fewer only at the end; how many it appended.
fn read_up_to(inner: &mut impl Read, out: &mut Vec<u8>, n: usize) -> io::Result<usize> {
    inner.by_ref().take(n as u64).read_to_end(out)
}

/// The header of a BAM, for a tool that reads no record, when every block inflates cleanly and
/// the header parses: then the whole-file read would have answered with it. `None` otherwise, for
/// the caller's whole-file read to give the answer it gave.
pub fn clean_header(path: &str) -> Option<SamHeader> {
    let file = File::open(path).ok()?;
    let mut stream = RecordStream::new(inflater(Blocks::new(BufReader::new(file))));
    let header = stream.header().ok()?;
    let mut inflater = stream.into_inner();
    io::copy(&mut inflater, &mut io::sink()).ok()?;
    let inflated = inflater.count;
    let mut blocks = inflater.inner.finish().ok()?;
    io::copy(&mut blocks, &mut io::sink()).ok()?;
    (!blocks.damaged && blocks.declared == inflated).then_some(header)
}

/// A BAM decoded in one pass, for a tool that can throw away what it made of the records: the
/// header to `on_header`, then each record, validated, to `on_record`. `true` when the whole-file
/// read would have handed over the same records and raised nothing; `false` at the first sign
/// otherwise, for the caller to read the file whole for the answer it gave.
pub fn one_pass(
    path: &str,
    on_header: impl FnOnce(&SamHeader),
    mut on_record: impl FnMut(BamRecord),
) -> bool {
    let Ok(file) = File::open(path) else {
        return false;
    };
    let mut stream = RecordStream::new(inflater(Blocks::new(BufReader::new(file))));
    let Ok(header) = stream.header() else {
        return false;
    };
    let mut validator = Validator::new(&header, &stream.binary_sequences);
    on_header(&header);
    loop {
        match stream.record() {
            Ok(Some(record)) if validator.check(&record).is_ok() => on_record(record),
            Ok(None) => break,
            _ => return false,
        }
    }
    let mut inflater = stream.into_inner();
    if io::copy(&mut inflater, &mut io::sink()).is_err() {
        return false;
    }
    let inflated = inflater.count;
    let Ok(mut blocks) = inflater.inner.finish() else {
        return false;
    };
    io::copy(&mut blocks, &mut io::sink()).is_ok() && !blocks.damaged && blocks.declared == inflated
}

/// noodles' multithreaded reader, counting what it inflates.
pub struct MultithreadedInflater<R> {
    inner: noodles_bgzf::io::MultithreadedReader<R>,
    count: u64,
}

/// The cores a run may keep busy at once.
fn cores() -> usize {
    std::thread::available_parallelism().map_or(1, NonZeroUsize::get)
}

/// The inflating workers' share of the cores: half, the other half going to rayon, which deflates
/// the output and runs the tools' parallel work. Measured on 16 cores against 2, 4 and 8 inflating
/// workers, half kept the wall time of two full pools and cut their CPU time by 8 to 25%.
pub fn inflate_workers() -> usize {
    (cores() / 2).max(1)
}

/// Sizes rayon's global pool to the cores the inflating workers leave, so that the two pools
/// together never ask for more threads than there are cores. An explicit `RAYON_NUM_THREADS`
/// wins, as rayon's own default lets it.
pub fn budget_threads() {
    if std::env::var_os("RAYON_NUM_THREADS").is_some() {
        return;
    }
    let rest = cores().saturating_sub(inflate_workers()).max(1);
    // A pool someone already built keeps its size.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(rest)
        .build_global();
}

fn inflater<R: Read + Send + 'static>(inner: R) -> MultithreadedInflater<R> {
    let workers = NonZeroUsize::new(inflate_workers()).unwrap_or(NonZeroUsize::MIN);
    MultithreadedInflater {
        inner: noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, inner),
        count: 0,
    }
}

impl<R: Read + Send + 'static> Read for MultithreadedInflater<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.count += n as u64;
        Ok(n)
    }
}

/// Why decoding a stream stopped.
#[derive(Debug)]
pub enum Failure {
    /// The inflated stream could not be read.
    Io(io::Error),
    /// What [`BamReader`] raised on the same bytes.
    Bam(ReadError),
}

impl From<io::Error> for Failure {
    fn from(e: io::Error) -> Self {
        Failure::Io(e)
    }
}

/// [`BamReader`] over an inflated stream rather than a slice: the same header, the same records
/// and the same errors, holding no more of the stream than the record being decoded.
pub struct RecordStream<R> {
    inner: R,
    buf: Vec<u8>,
    start: usize,
    end_of_stream: bool,
    /// The header's binary dictionary, once [`Self::header`] has read it.
    pub binary_sequences: Vec<(String, i32)>,
}

impl<R: Read> RecordStream<R> {
    pub fn new(inner: R) -> Self {
        RecordStream {
            inner,
            buf: Vec::new(),
            start: 0,
            end_of_stream: false,
            binary_sequences: Vec::new(),
        }
    }

    /// Reads until `n` bytes past `start` are held, or the stream ends.
    fn fill(&mut self, n: usize) -> io::Result<()> {
        while self.buf.len() - self.start < n && !self.end_of_stream {
            if self.start > 0 && self.start >= self.buf.len() / 2 {
                self.buf.drain(..self.start);
                self.start = 0;
            }
            let len = self.buf.len();
            self.buf.resize(len + CHUNK, 0);
            let got = loop {
                match self.inner.read(&mut self.buf[len..]) {
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    other => break other,
                }
            };
            self.buf.truncate(len + *got.as_ref().unwrap_or(&0));
            if got? == 0 {
                self.end_of_stream = true;
            }
        }
        Ok(())
    }

    /// [`BamReader::new`]. A `Truncated` answer is only final once the whole stream is held,
    /// which is when it reports the length the whole-file read would have.
    pub fn header(&mut self) -> Result<SamHeader, Failure> {
        loop {
            match BamReader::new(&self.buf) {
                Ok(reader) => {
                    self.start = reader.records_start();
                    self.binary_sequences = reader.header.binary_sequences;
                    return Ok(reader.header.text);
                }
                Err(ReadError::Truncated { need, .. }) if !self.end_of_stream => self.fill(need)?,
                Err(e) => return Err(Failure::Bam(e)),
            }
        }
    }

    /// [`BamReader`]'s next record. `BamRecord::decode` reads no further than the length the
    /// record declares, and only reports how much it was given when that is short, by which time
    /// the rest of the stream is held, so it answers here as it did on the whole stream.
    pub fn record(&mut self) -> Result<Option<BamRecord>, Failure> {
        self.fill(4)?;
        let held = &self.buf[self.start..];
        if held.is_empty() {
            return Ok(None);
        }
        if held.len() >= 4 {
            let length = i32::from_le_bytes(held[..4].try_into().unwrap());
            if length >= 0 {
                self.fill(4 + length as usize)?;
            }
        }
        match BamRecord::decode(&self.buf[self.start..]) {
            Ok(None) => {
                self.start = self.buf.len();
                Ok(None)
            }
            Ok(Some((record, used))) => {
                self.start += used;
                Ok(Some(record))
            }
            Err(e) => {
                self.start = self.buf.len();
                Err(Failure::Bam(ReadError::Record(e)))
            }
        }
    }

    /// Up to `count` batches of whole records, each about [`BATCH`] bytes, cut at the lengths
    /// they declare; and whether the cut stopped at the stream's tail: its end, or a record whose
    /// length runs past it or cannot be a length, which is left to [`Self::record`].
    fn cut(&mut self, count: usize) -> Result<(Vec<Vec<u8>>, bool), Failure> {
        let mut batches: Vec<Vec<u8>> = Vec::with_capacity(count);
        let mut tail = false;
        while batches.len() < count && !tail {
            // What is held starts the batch, and the stream is read straight into it after that.
            let mut batch = Vec::with_capacity(BATCH + CHUNK);
            batch.extend_from_slice(&self.buf[self.start..]);
            self.buf.clear();
            self.start = 0;
            // How many bytes of whole records the batch begins with.
            let mut whole = 0;
            loop {
                while batch.len() - whole >= 4 {
                    let length = i32::from_le_bytes(batch[whole..whole + 4].try_into().unwrap());
                    if length < 0 {
                        tail = true;
                        break;
                    }
                    let total = 4 + length as usize;
                    if batch.len() - whole < total {
                        break;
                    }
                    whole += total;
                }
                if tail || whole >= BATCH {
                    break;
                }
                if self.end_of_stream {
                    tail = true;
                    break;
                }
                let len = batch.len();
                batch.resize(len + CHUNK, 0);
                let got = loop {
                    match self.inner.read(&mut batch[len..]) {
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                        other => break other,
                    }
                };
                batch.truncate(len + *got.as_ref().unwrap_or(&0));
                if got? == 0 {
                    self.end_of_stream = true;
                }
            }
            // The bytes past the whole records are held, for the next batch or for `record`.
            self.buf.extend_from_slice(&batch[whole..]);
            batch.truncate(whole);
            if !batch.is_empty() {
                batches.push(batch);
            }
        }
        Ok((batches, tail))
    }

    /// Decodes every remaining record and keeps none, answering with the first error
    /// [`Self::record`] would have met. The records are cut at the lengths they declare and
    /// decoded on every core, a batch at a time; a record whose length runs past the end of the
    /// stream, or cannot be a length, is left to [`Self::record`], after the batches before it.
    pub fn check_records(&mut self) -> Result<(), Failure> {
        let threads = rayon::current_num_threads().max(1);
        loop {
            let (batches, tail) = self.cut(threads)?;
            use rayon::prelude::*;
            let first_error = batches
                .par_iter()
                .map(|batch| {
                    let mut p = 0;
                    while p < batch.len() {
                        match BamRecord::decode(&batch[p..]) {
                            Ok(Some((_, used))) => p += used,
                            Ok(None) => break,
                            Err(e) => return Some(e),
                        }
                    }
                    None
                })
                .collect::<Vec<_>>()
                .into_iter()
                .flatten()
                .next();
            if let Some(e) = first_error {
                self.start = self.buf.len();
                return Err(Failure::Bam(ReadError::Record(e)));
            }
            if tail {
                while self.record()?.is_some() {}
                return Ok(());
            }
        }
    }

    pub fn into_inner(self) -> R {
        self.inner
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// What the tools did before: the whole file, inflated, every record decoded.
    fn whole_file(path: &str) -> Result<(SamHeader, Vec<BamRecord>), String> {
        let raw = std::fs::read(path).map_err(|e| format!("{e}"))?;
        if raw.starts_with(&[0x1f, 0x8b]) {
            let plain = crate::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
            let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
            let header = reader.header.text.clone();
            let records = reader
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("{e:?}"))?;
            Ok((header, records))
        } else {
            let text = String::from_utf8(raw).map_err(|e| format!("{e}"))?;
            htsjdk_bam::sam_file::read_sam(&text).map_err(|e| format!("{e:?}"))
        }
    }

    /// Opens and reads to the end, as a tool does, with the error a mid-stream fallback raises
    /// returned rather than reported.
    fn streamed(path: &str) -> (Result<(SamHeader, Vec<BamRecord>), String>, bool) {
        match Input::open(path) {
            Ok(input) => {
                let (header, mut records) = input.into_parts();
                let streaming = matches!(records, Records::Stream(..));
                let mut all = Vec::new();
                loop {
                    match records.try_next() {
                        Ok(Some(record)) => all.push(record),
                        Ok(None) => break,
                        Err(e) => return (Err(e.message()), streaming),
                    }
                }
                (Ok((header, all)), streaming)
            }
            Err(e) => (Err(e.message()), false),
        }
    }

    fn sam_text(n: usize) -> String {
        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n");
        for i in 0..n {
            let seq: String = (0..50).map(|j| b"ACGT"[(i + j) % 4] as char).collect();
            text.push_str(&format!(
                "read{i}\t0\tchr1\t{}\t60\t50M\t*\t0\t0\t{seq}\t{}\tNM:i:{}\n",
                i + 1,
                "I".repeat(50),
                i % 3
            ));
        }
        text
    }

    #[test]
    fn records_one_line_at_a_time_write_the_sam_the_whole_slice_does() {
        let (header, records) = htsjdk_bam::sam_file::read_sam(&sam_text(50)).unwrap();
        let mut text = header.encode();
        for record in &records {
            text.push_str(&sam_line(&header, record).unwrap());
            text.push('\n');
        }
        assert_eq!(
            Some(text),
            htsjdk_bam::sam_file::write_sam(&header, &records)
        );
    }

    /// The plain BAM stream of `n` records, and where its records begin.
    fn plain_bam(n: usize) -> (Vec<u8>, usize) {
        let (header, records) = htsjdk_bam::sam_file::read_sam(&sam_text(n)).unwrap();
        let mut writer = htsjdk_bam::writer::BamWriter::new(Vec::new(), &header).unwrap();
        for record in &records {
            writer.write(record).unwrap();
        }
        let plain = htsjdk_bgzf::decompress_all(&writer.finish().unwrap()).unwrap();
        let start = BamReader::new(&plain).unwrap().records_start();
        (plain, start)
    }

    fn bgzf(plain: &[u8]) -> Vec<u8> {
        let mut w = htsjdk_bgzf::BgzfWriter::new(Vec::new());
        w.write_all(plain).unwrap();
        w.into_inner().unwrap()
    }

    /// Offset of the `k`th record.
    fn record_at(plain: &[u8], start: usize, k: usize) -> usize {
        let mut p = start;
        for _ in 0..k {
            p += 4 + i32::from_le_bytes(plain[p..p + 4].try_into().unwrap()) as usize;
        }
        p
    }

    fn block_end(file: &[u8], at: usize) -> usize {
        at + u16::from_le_bytes([file[at + 16], file[at + 17]]) as usize + 1
    }

    #[test]
    fn answers_as_the_whole_file_read_did() {
        let (plain, start) = plain_bam(20_000);
        let clean = bgzf(&plain);
        let second = block_end(&clean, 0);
        assert!(block_end(&clean, second) < clean.len(), "several blocks");

        // Whether the single pass streams: every file whose blocks are sound and whose header
        // decodes. A damaged block, a truncated file, or a header that does not decode is
        // answered by the whole-file read at open, before any record.
        let mut cases: Vec<(&str, Vec<u8>, bool)> = vec![("clean", clean.clone(), true)];
        cases.push(("truncated", clean[..clean.len() - 30].to_vec(), false));
        let mut corrupt = clean.clone();
        corrupt[second + HEADER + 40] ^= 0xff;
        cases.push(("corrupt payload", corrupt, true));
        let mut lying = clean.clone();
        lying[second - 4] = lying[second - 4].wrapping_add(1);
        cases.push(("lying ISIZE", lying, false));
        let mut crc = clean.clone();
        crc[second - 8] ^= 0xff;
        cases.push(("bad CRC, which htsjdk does not check", crc, false));
        let mut trailing = clean.clone();
        trailing.extend_from_slice(b"xyz");
        cases.push(("trailing bytes", trailing, false));
        let mut short = plain.clone();
        short.truncate(short.len() - 10);
        cases.push(("truncated last record", bgzf(&short), true));
        let mut bad_length = plain.clone();
        let r = record_at(&plain, start, 2);
        bad_length[r..r + 4].copy_from_slice(&5i32.to_le_bytes());
        cases.push(("bad record length", bgzf(&bad_length), true));
        let mut huge_length = plain.clone();
        huge_length[r..r + 4].copy_from_slice(&i32::MAX.to_le_bytes());
        cases.push(("record longer than the file", bgzf(&huge_length), true));
        let mut bad_cigar = plain.clone();
        let r = record_at(&plain, start, 15_000);
        let cigar = r + 36 + bad_cigar[r + 12] as usize;
        bad_cigar[cigar] |= 0x0f;
        let bad_cigar = bgzf(&bad_cigar);
        cases.push(("unknown CIGAR operator", bad_cigar.clone(), true));
        cases.push((
            "unknown CIGAR operator, then a truncated file",
            bad_cigar[..bad_cigar.len() - 30].to_vec(),
            false,
        ));
        cases.push(("truncated header", bgzf(&plain[..20]), false));
        cases.push(("not a BAM", bgzf(b"SAM\x01 and more"), false));
        cases.push(("SAM", sam_text(20).into_bytes(), false));
        cases.push(("empty", Vec::new(), false));

        let dir = std::env::temp_dir().join(format!("bam-stream-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for (i, (name, bytes, streams)) in cases.iter().enumerate() {
            let path = dir.join(format!("{i}.bam"));
            std::fs::write(&path, bytes).unwrap();
            let path = path.to_str().unwrap();
            let (ours, streaming) = streamed(path);
            assert_eq!(
                format!("{ours:?}"),
                format!("{:?}", whole_file(path)),
                "{name}"
            );
            assert_eq!(streaming, *streams, "{name}");
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
