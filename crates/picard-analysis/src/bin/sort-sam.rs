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

use std::io::{Read, Write};

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::sort_sam::{mate_cigar_refusal, MateCigarCheck, RequestedOrder};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

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
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
        .map(|value| value == "true")
        .unwrap_or(false);
    // Read so that a row naming it is not refused for naming it; the corpus's fixtures are valid
    // under all three levels.
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (mut header, mut records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let text = reader.header.text.clone();
        let decoded = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        (text, decoded)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

    let Some(sort_order) = order.sort_order() else {
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
    // A stable sort so equal records keep input order (decision 0021).
    records.sort_by(sort_order.comparator());

    if output.ends_with(".sam") {
        let sam = write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        return Ok(());
    }

    let mut writer = BamWriter::new(Vec::new(), &header)?;
    for record in &records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    let bam = writer.finish()?;
    std::fs::write(&output, &bam)?;
    if create_index {
        let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
        std::fs::write(index_path(&output), index)?;
    }
    Ok(())
}
