//! `CleanSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CleanSam.doWork` at tag 3.4.0 down to the two decisions the array varies
//! around the per-record cleaning, which lives in `picard_analysis::clean_sam`.
//!
//! * `VALIDATION_STRINGENCY`: `doWork` opens its reader at LENIENT when the caller asked for
//!   STRICT, because the records CleanSam exists to fix are exactly the ones STRICT refuses to
//!   read. LENIENT and SILENT are passed through as given. All three rows therefore read the same
//!   corpus, and the argument is covered rather than merely accepted.
//! * `CREATE_INDEX`: the writer indexes only a coordinate-sorted BAM. `initializeBAMWriter` calls
//!   `enableBamIndexConstruction` under `createIndex && sortOrder == coordinate`, so a
//!   queryname-sorted input with `CREATE_INDEX=true` is not an error and not an index: it is
//!   silently the same output as `CREATE_INDEX=false`, which is behaviour a port that threw would
//!   get wrong.
//!
//! The output format is `SAMFileWriterFactory.makeWriter`'s decision, not the tool's: `.sam` is
//! SAM text, `.cram` is CRAM (not ported), and *anything else* is a BAM, including the
//! `output.txt` the arrays pass. The index beside it is named as `createBamIndex` names it, which
//! drops only a `.bam` suffix, so `output.txt` is indexed as `output.txt.bai`.
//!
//! `REFERENCE_SEQUENCE` and `TMP_DIR` are accepted and unused: the reader needs a reference only
//! for CRAM, and nothing here spills to disk. Accepting them is what lets a row that names them
//! measure the tool rather than the argument parser.

use htsjdk_bam::text_parse::ValidationStringency;
use picard_analysis::args::arg;
use picard_analysis::clean_sam::{clean_records, CleanError};
use picard_analysis::java_writer;
use picard_analysis::stream_output::Output;

/// How many records are cleaned in parallel at a time, between reading and writing.
const CHUNK: usize = 16 * 1024;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CleanSam");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);
    // `SamReaderFactory` at the stringency `doWork` asks for: STRICT becomes LENIENT, so that the
    // reader accepts the unmapped read with a nonzero MAPQ that the tool is there to zero.
    let stringency = match arg(&args, "VALIDATION_STRINGENCY=").as_deref() {
        None | Some("STRICT") | Some("LENIENT") => ValidationStringency::Lenient,
        Some("SILENT") => ValidationStringency::Silent,
        Some(other) => return Err(format!("unknown VALIDATION_STRINGENCY: {other}").into()),
    };
    picard_analysis::record_validation::set_stringency(
        if matches!(stringency, ValidationStringency::Silent) {
            picard_analysis::record_validation::Stringency::Silent
        } else {
            picard_analysis::record_validation::Stringency::Lenient
        },
    );

    let (header, records) =
        picard_analysis::bam_stream::Input::open_with(&input, stringency)?.into_parts();

    let filename = if output.ends_with(".sam") {
        "null".to_string()
    } else {
        java_writer::bam_filename(&output)
    };
    // The header is written back unchanged: CleanSam hands the reader's header straight to the
    // writer, adds no `@PG` of its own, and does not re-sort, so `SO` still describes the records.
    // Nothing reached the disk before every record had been read, so a reader that throws leaves
    // no file.
    let mut out = Output::create(&output, &header, false)?;
    let mut records = records.peekable();
    // Each chunk after the first starts with the last record of the one before, already written,
    // for the order check across the two.
    let mut chunk = Vec::with_capacity(CHUNK + 1);
    while records.peek().is_some() {
        let carried = chunk.len();
        chunk.extend(records.by_ref().take(CHUNK));
        // Per record, `doWork` cleans and then hands the record to a presorted writer, which
        // checks the order: whichever of the two throws first, in file order, ends the run. The
        // records before it were written, and what of them reached the disk is what the writer's
        // buffers let through.
        let thrown_at = match clean_records(&header, &mut chunk[carried..]) {
            Ok(()) => None,
            Err(CleanError::Thrown { index, thrown }) => Some((carried + index, thrown)),
            Err(CleanError::Parse(e)) => return Err(format!("{e:?}").into()),
        };
        let out_of_order = java_writer::first_out_of_order(&header, &chunk, &filename);
        let failure = match (thrown_at, out_of_order) {
            (Some(clean), Some(order)) => Some(if clean.0 <= order.0 { clean } else { order }),
            (clean, order) => clean.or(order),
        };
        let written = failure.as_ref().map_or(chunk.len(), |(index, _)| *index);
        for record in &chunk[carried..written] {
            out.write(record)?;
        }
        if let Some((_, thrown)) = failure {
            out.abandon_unclosed()?;
            eprintln!("Exception in thread \"main\" {thrown}");
            std::process::exit(1);
        }
        chunk.drain(..chunk.len() - 1);
    }
    // Only a coordinate-sorted file is indexed, and asking for an index on any other order is not
    // an error: `initializeBAMWriter` simply does not enable the indexer.
    let coordinate = header.attributes.get("SO") == Some("coordinate");
    out.finish(create_index, create_index && coordinate)?;
    Ok(())
}
