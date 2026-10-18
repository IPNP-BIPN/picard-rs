//! `SortSam` as a runnable binary: the covering array's port side.
//!
//! The array varies the sort order over `coordinate`, `queryname` and `duplicate`, and the third
//! is the interesting one. A duplicate sort compares records by their MATE's unclipped coordinate,
//! which it reads out of the `MC` tag, so a file whose paired records carry no mate cigar cannot
//! be sorted that way at all: the reference refuses it before writing a record, naming the first
//! record it asked. Every fixture in this corpus is such a file, so all three of the array's
//! duplicate rows are refusals rather than outputs.
//!
//! The output is a BAM unless it is named `.sam`, as `SAMFileWriterFactory` decides it, with a BAI
//! beside it when `CREATE_INDEX=true`; the input is a BAM or a SAM, told apart by its first bytes.

use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::{BamReader, ReadError};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text::write_alignment;
use htsjdk_bgzf::read::BgzfReader;
use picard_analysis::args::arg;
use picard_analysis::record_validation::Validator;
use picard_analysis::sort_sam::{
    mate_cigar_refusal, MateCigarCheck, RequestedOrder, SortingCollection,
    DEFAULT_MAX_RECORDS_IN_RAM,
};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;

/// The name htsjdk gives the index it writes beside an output.
fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
    }
}

/// The last reference position an alignment covers, which the refusal's message names.
fn alignment_end(record: &BamRecord) -> i32 {
    record.alignment_start + reference_span(&record.cigar) - 1
}

fn reference_span(cigar: &Cigar) -> i32 {
    cigar
        .elements
        .iter()
        .filter(|element| element.op.consumes_reference_bases())
        .map(|element| element.length as i32)
        .sum()
}

fn check_of(record: &BamRecord, header: &SamHeader) -> MateCigarCheck {
    let reference = if record.reference_index >= 0 {
        header
            .sequences
            .get(record.reference_index as usize)
            .map(|sequence| sequence.name.clone())
    } else {
        None
    };
    MateCigarCheck {
        read_name: record.read_name.clone(),
        paired: record.flags & READ_PAIRED != 0,
        first_of_pair: record.flags & FIRST_OF_PAIR != 0,
        unmapped: record.flags & READ_UNMAPPED != 0,
        mate_unmapped: record.flags & MATE_UNMAPPED != 0,
        read_length: record.read_bases.len(),
        reference,
        start: record.alignment_start,
        end: alignment_end(record),
        mate_cigar: match record.tags.get(Tag::new(b"MC")) {
            Some(TagValue::Str(text)) => Some(text.to_string()),
            _ => None,
        },
        sam_string: picard_analysis::view_sam::sam_string(header, record),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SortSam");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let requested = arg(&args, "SORT_ORDER=")
        .or_else(|| arg(&args, "SO="))
        .ok_or("SORT_ORDER= is required")?;
    let order = RequestedOrder::parse(&requested)
        .ok_or_else(|| format!("unknown SORT_ORDER: {requested}"))?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);
    // Read so that a row naming it is not refused for naming it; the corpus's fixtures are valid
    // under all three levels.
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let max_records_in_ram = match arg(&args, "MAX_RECORDS_IN_RAM=") {
        Some(value) => value
            .parse()
            .map_err(|_| format!("MAX_RECORDS_IN_RAM is not an integer: {value}"))?,
        None => DEFAULT_MAX_RECORDS_IN_RAM,
    };
    // shortcut: the first TMP_DIR only, where htsjdk picks the one with the most free space.
    let tmp_dir = arg(&args, "TMP_DIR=")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);

    let mut head = Vec::new();
    std::fs::File::open(&input)?
        .take(4)
        .read_to_end(&mut head)?;
    let (mut header, records): (
        SamHeader,
        Box<dyn Iterator<Item = Result<BamRecord, String>>>,
    ) = if picard_analysis::cram_input::is_cram(&head) {
        // shortcut: a CRAM is decoded whole, then validated as the BAM it decodes to.
        let raw = picard_analysis::cram_input::uncram(std::fs::read(&input)?);
        let plain = picard_analysis::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
        let reader =
            picard_analysis::record_validation::reader(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Vec<_>>();
        (header, Box::new(records.into_iter()))
    } else if head.starts_with(&[0x1f, 0x8b]) {
        let (header, stream) = BamStream::open(&input)?;
        (header, Box::new(stream))
    } else {
        // shortcut: a SAM input is read whole; stream it if a SAM larger than RAM turns up.
        let text = String::from_utf8(std::fs::read(&input)?)?;
        let (header, records) = read_sam(&text).map_err(|e| format!("{e:?}"))?;
        (header, Box::new(records.into_iter().map(Ok)))
    };

    let Some(sort_order) = order.sort_order() else {
        let records = records.collect::<Result<Vec<_>, _>>()?;
        // The duplicate order, which this corpus never reaches: the sort asks the first record it
        // touches for its mate's unclipped coordinate, and a record with no mate cigar has none.
        let checks: Vec<MateCigarCheck> = records
            .iter()
            .map(|record| check_of(record, &header))
            .collect();
        if let Some(message) = mate_cigar_refusal(&checks) {
            // The reference's own line, printed as the reference prints it: an uncaught exception
            // on stderr and a status of one. Wrapping it in Rust's `Error:` would be a different
            // message for the same refusal.
            eprintln!("Exception in thread \"main\" htsjdk.samtools.SAMException: {message}");
            std::process::exit(1);
        }
        eprintln!("sorting by duplicate order is a separate surface");
        std::process::exit(1);
    };

    header.set_sort_order(sort_order.name());
    let mut sorting = SortingCollection::new(sort_order.comparator(), max_records_in_ram, tmp_dir);
    // Every record is read before the output is created, so a damaged input leaves no output.
    for record in records {
        sorting.add(record?)?;
    }
    let sorted = sorting.into_sorted()?;

    if output.ends_with(".sam") {
        let name_of = |index: i32| -> &str {
            usize::try_from(index)
                .ok()
                .and_then(|i| header.sequences.get(i))
                .map_or("*", |s| s.name.as_str())
        };
        let mut out = BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(header.encode().as_bytes())?;
        for record in sorted {
            let record = record?;
            let line = write_alignment(
                &record,
                name_of(record.reference_index),
                name_of(record.mate_reference_index),
            )
            .ok_or("records failed to re-encode as SAM")?;
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
        out.flush()?;
        return Ok(());
    }

    if picard_analysis::cram_output::is_cram(&output) {
        let mut writer = picard_analysis::output::bam_writer(Vec::new(), &header)?;
        for record in sorted {
            writer.write(&record?).map_err(|e| format!("{e:?}"))?;
        }
        let bam = writer.finish()?;
        picard_analysis::cram_output::write_from_bam(&output, &bam, create_index)?;
        return Ok(());
    }

    let out = BufWriter::new(std::fs::File::create(&output)?);
    let mut writer = picard_analysis::output::bam_writer(out, &header)?;
    if create_index {
        writer = writer.with_index();
    }
    for record in sorted {
        writer.write(&record?).map_err(|e| format!("{e:?}"))?;
    }
    if create_index {
        let (out, index) = writer.finish_with_index()?;
        out.into_inner().map_err(|e| e.into_error())?;
        std::fs::write(index_path(&output), index)?;
    } else {
        writer.finish()?.into_inner().map_err(|e| e.into_error())?;
    }
    Ok(())
}

/// A BAM inflated block by block, answering as inflating it whole and then decoding did: a BGZF
/// error anywhere outranks a header or record error, because a decode that fails pulls the rest of
/// the stream before it answers.
struct BamStream {
    bgzf: BgzfReader<BufReader<std::fs::File>>,
    buf: Vec<u8>,
    pos: usize,
    eof: bool,
    /// The checks htsjdk's reader runs on each record it returns; set once the header is read.
    validator: Option<Validator>,
}

impl BamStream {
    fn open(path: &str) -> Result<(SamHeader, BamStream), Box<dyn std::error::Error>> {
        let mut stream = BamStream {
            bgzf: BgzfReader::new(BufReader::new(std::fs::File::open(path)?)),
            buf: Vec::new(),
            pos: 0,
            eof: false,
            validator: None,
        };
        loop {
            match BamReader::new(&stream.buf) {
                Ok(reader) => {
                    let header = reader.header.text.clone();
                    let validator = Validator::new(&header, &reader.header.binary_sequences);
                    stream.pos = reader.records_start();
                    stream.validator = Some(validator);
                    return Ok((header, stream));
                }
                Err(e) => {
                    if !stream.pull()? {
                        return Err(format!("{e:?}").into());
                    }
                }
            }
        }
    }

    /// Inflates one more block, or answers false at the end of the stream.
    fn pull(&mut self) -> Result<bool, String> {
        if self.eof {
            return Ok(false);
        }
        if self.pos * 2 >= self.buf.len() {
            self.buf.drain(..self.pos);
            self.pos = 0;
        }
        match self.bgzf.next_block().map_err(|e| format!("{e:?}"))? {
            Some(block) => {
                self.buf.extend_from_slice(&block.data);
                Ok(true)
            }
            None => {
                self.eof = true;
                Ok(false)
            }
        }
    }
}

impl Iterator for BamStream {
    type Item = Result<BamRecord, String>;

    fn next(&mut self) -> Option<Result<BamRecord, String>> {
        loop {
            let answer = BamRecord::decode(&self.buf[self.pos..]);
            if let Ok(Some((record, used))) = answer {
                self.pos += used;
                if let Some(Err(message)) = self.validator.as_mut().map(|v| v.check(&record)) {
                    picard_analysis::metrics_cli::thrown(&message);
                }
                return Some(Ok(record));
            }
            match self.pull() {
                Err(e) => return Some(Err(e)),
                Ok(true) => {}
                Ok(false) => {
                    return match answer {
                        Ok(_) => None,
                        Err(e) => Some(Err(format!("{:?}", ReadError::Record(e)))),
                    }
                }
            }
        }
    }
}
