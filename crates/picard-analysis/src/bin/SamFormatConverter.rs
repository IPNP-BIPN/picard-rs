//! `SamFormatConverter` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.SamFormatConverter.convert` at tag 3.4.0. The tool applies no transform: the
//! records pass through in file order, and the header is the reader's own, so the conversion is a
//! decode and a re-encode. What the array varies around that is the refusal.
//!
//! `convert` creates the writer **before** it checks anything, then throws
//! `Can't CREATE_INDEX unless sort order is coordinate` when `CREATE_INDEX` was asked for and the
//! header does not say `coordinate`. So a queryname-sorted input with `CREATE_INDEX=true` is a
//! refusal here, where the same pair is silently accepted by `CleanSam`: `CleanSam` leaves the
//! decision to `SAMFileWriterFactory`, which simply does not enable the indexer, and this tool
//! makes it an error of its own. Two tools, one argument, opposite answers; the array is what
//! makes the difference visible.
//!
//! Two consequences of the order of operations that the port reproduces rather than tidies up:
//! the check happens after the output file has been created, so the refused row still leaves an
//! empty file behind, and it happens before any record is written, so nothing partial is in it.
//! The port therefore creates the file and then refuses.
//!
//! Input and output formats are the file's, not the argument's: the reader sniffs BAM by its gzip
//! magic, and `makeSAMOrBAMWriter` writes SAM only for a `.sam` name and BAM for anything else,
//! including the `output.txt` the arrays pass.

use std::io::Write;

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::text_parse::ValidationStringency;
use picard_analysis::args::arg;

/// `BAMFileWriter.createBamIndex`: the `.bam` suffix is dropped, anything else is kept.
fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SamFormatConverter");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);
    // `SamReaderFactory.setDefaultValidationStringency(VALIDATION_STRINGENCY)`, which
    // `CommandLineProgram` calls before `doWork`: this tool asks for no downgrade, unlike
    // `CleanSam`, so the value is honoured as given.
    let stringency = match arg(&args, "VALIDATION_STRINGENCY=").as_deref() {
        None | Some("STRICT") => ValidationStringency::Strict,
        Some("LENIENT") => ValidationStringency::Lenient,
        Some("SILENT") => ValidationStringency::Silent,
        Some(other) => return Err(format!("unknown VALIDATION_STRINGENCY: {other}").into()),
    };

    let (header, records) =
        picard_analysis::bam_stream::Input::open_with(&input, stringency)?.into_parts();

    // `new SAMFileWriterFactory().makeWriter(...)` opens the output before the check below runs,
    // which is why the refused row leaves an empty file rather than no file.
    std::fs::File::create(&output)?;

    if create_index && header.attributes.get("SO") != Some("coordinate") {
        // The reference's own line: an uncaught PicardException on stderr and a status of one.
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: \
             Can't CREATE_INDEX unless sort order is coordinate"
        );
        std::process::exit(1);
    }

    if output.ends_with(".sam") {
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(header.encode().as_bytes())?;
        for record in records {
            // `write_sam` refused the file for one record it could not encode, which left the
            // empty file opened above: written line by line, the same file is emptied again.
            let Some(line) = picard_analysis::bam_stream::sam_line(&header, &record) else {
                drop(out);
                std::fs::File::create(&output)?;
                return Err("records failed to re-encode as SAM".into());
            };
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
        out.flush()?;
        return Ok(());
    }

    if picard_analysis::cram_output::is_cram(&output) {
        let mut writer = picard_analysis::output::bam_writer(Vec::new(), &header)?;
        for record in records {
            writer.write(&record).map_err(|e| format!("{e:?}"))?;
        }
        let bam = writer.finish()?;
        picard_analysis::cram_output::write_from_bam(&output, &bam, create_index)?;
        return Ok(());
    }
    let file = std::io::BufWriter::new(std::fs::File::create(&output)?);
    let mut writer = picard_analysis::output::bam_writer(file, &header)?;
    for record in records {
        if let Err(e) = writer.write(&record) {
            // The BAM was built whole before it was written, so a record the encoder refused left
            // the empty file opened above.
            drop(writer);
            std::fs::File::create(&output)?;
            return Err(format!("{e:?}").into());
        }
    }
    writer.finish()?.flush()?;
    if create_index {
        let bam = std::fs::read(&output)?;
        let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
        std::fs::write(index_path(&output), index)?;
    }
    Ok(())
}
