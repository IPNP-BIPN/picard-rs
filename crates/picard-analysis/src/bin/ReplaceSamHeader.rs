//! `ReplaceSamHeader` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.ReplaceSamHeader.doWork` at tag 3.4.0, which is a fork the module does not
//! take:
//!
//! ```java
//! if (BamFileIoUtils.isBamFile(INPUT)) blockCopyReheader(replacementHeader);
//! else                                 standardReheader(replacementHeader);
//! ```
//!
//! The two branches do not agree about what is legal, nor about what they check first.
//!
//! `standardReheader` refuses when the two files declare different sort orders, as a
//! `PicardException` naming both. `blockCopyReheader` is `BamFileIoUtils.reheaderBamFile`, which
//! asserts the INPUT is **writable** before it copies anything -- the block copy is written for the
//! in-place case too -- and then applies a *different* sort-order rule, in htsjdk: an `unsorted`
//! new header is accepted against any original, and the message names what the header would have
//! to be instead.
//!
//! The writability assert is why six of this array's nine rows are refusals: the corpus is mounted
//! read-only, as the oracle contract mounts it, so a BAM cannot be reheadered from it at all.
//!
//! What the block copy means for the port: the records are the input's own bytes, so writing them
//! back through the encoder has to reproduce them exactly. It does, and the covering array is where
//! that is checked rather than assumed.
//!
//! `rec.setHeader(replacementHeader)` keeps each record's integer reference index and re-resolves
//! RNAME/RNEXT against the new dictionary at that index; when the two dictionaries share their
//! `@SQ` block, as here, the reads are untouched.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::text_parse::ValidationStringency;
use picard_analysis::args::arg;
use picard_analysis::bam_stream::Input;
use picard_analysis::stream_output::Output;

/// Read a SAM or BAM, told apart by the gzip magic as `SamReaderFactory` tells them apart: whether
/// it is a BAM, and the input, its records not validated.
fn read_any(path: &str) -> Result<(Input, bool), String> {
    let mut magic = Vec::new();
    std::fs::File::open(path)
        .and_then(|f| f.take(2).read_to_end(&mut magic))
        .map_err(|e| e.to_string())?;
    // `BamFileIoUtils.isBamFile` decides the block copy, and a CRAM is not one: it is decoded and
    // takes the standard path.
    let is_bam = magic == [0x1f, 0x8b];
    // `standardReheader` opens INPUT at SILENT whatever VALIDATION_STRINGENCY says, and a BAM's
    // records are read by `BamReader::new`, unvalidated.
    let input = Input::open_with(path, ValidationStringency::Silent)
        .map_err(|e| e.message())?
        .unvalidated();
    Ok((input, is_bam))
}

/// `SAMFileHeader.getSortOrder()`: an absent `SO` is `unsorted`, and the name is what the refusal
/// prints.
fn sort_order(header: &SamHeader) -> &str {
    match header.attributes.get("SO") {
        Some(order @ ("coordinate" | "queryname" | "unsorted" | "duplicate")) => order,
        Some(_) => "unknown",
        None => "unsorted",
    }
}

/// `reheaderBamFile`'s first refusal: the INPUT must be writable.
fn check_writable(input: &str) {
    if std::fs::OpenOptions::new().write(true).open(input).is_err() {
        // `File.getAbsolutePath()`: the working directory joined to a relative path, with no
        // symlink resolution. Canonicalizing instead would rewrite /var to /private/var on
        // macOS and print a path the reference never would.
        let absolute = if std::path::Path::new(input).is_absolute() {
            input.to_string()
        } else {
            std::env::current_dir()
                .map(|cwd| cwd.join(input).display().to_string())
                .unwrap_or_else(|_| input.to_string())
        };
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMException: \
             File exists but is not writable: {absolute}"
        );
        std::process::exit(1);
    }
}

/// `reheaderBamFile`'s sort-order refusal, on the two headers.
fn check_sort_order(input_header: &SamHeader, replacement: &SamHeader) {
    // The block copy has a sort-order check of its own, in htsjdk rather than in Picard, and
    // it is not the same rule as the text path's: an `unsorted` new header is accepted against
    // any original, and the message names what the header would have to be.
    if sort_order(replacement) != "unsorted" && sort_order(replacement) != sort_order(input_header)
    {
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMException: \
             Sort order of new header does not match the original file, needs to be {}",
            sort_order(input_header)
        );
        std::process::exit(1);
    }
}

/// Whether the INPUT is a BGZF stream, by its gzip magic, as `SamReaderFactory` tells them apart.
fn is_bam(path: &str) -> std::io::Result<bool> {
    let mut magic = Vec::new();
    std::fs::File::open(path)?.take(2).read_to_end(&mut magic)?;
    Ok(magic == [0x1f, 0x8b])
}

/// The last block's worth of the file, all `checkTermination` reads.
fn tail(path: &str) -> std::io::Result<Vec<u8>> {
    use std::io::{Seek, SeekFrom};
    let mut file = std::fs::File::open(path)?;
    let len = file.metadata()?.len();
    file.seek(SeekFrom::Start(
        len.saturating_sub(htsjdk_bgzf::MAX_COMPRESSED_BLOCK_SIZE as u64),
    ))?;
    let mut data = Vec::new();
    file.read_to_end(&mut data)?;
    Ok(data)
}

/// The INPUT's header, inflating no further than it reaches, as `findVirtualOffsetOfFirstRecordInBam`
/// reads it.
fn header_only(path: &str) -> Result<SamHeader, picard_analysis::bam_stream::Failure> {
    let file = std::io::BufReader::new(std::fs::File::open(path)?);
    picard_analysis::bam_stream::RecordStream::new(htsjdk_bgzf::BgzfReader::new(file)).header()
}

/// HEADER's header, as `SamReaderFactory.getFileHeader` reads it: a BAM no further than its header,
/// so a large or damaged file past the header costs nothing and refuses nothing. A header that does
/// not read is read whole, which throws what htsjdk throws.
fn replacement_header(path: &str) -> Result<SamHeader, Box<dyn std::error::Error>> {
    if is_bam(path)? {
        if let Ok(header) = header_only(path) {
            return Ok(header);
        }
        return Ok(Input::open_with(path, ValidationStringency::Silent)?.header);
    }
    Ok(read_any(path)?.0.header)
}

/// The block copy: the new header's blocks, then the input's record blocks as they are.
fn block_copy(
    input: &str,
    output: &str,
    replacement: &SamHeader,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::io::Write;
    let (level, deflater) = picard_analysis::output::current();
    let copied = (|| -> Result<(), Box<dyn std::error::Error>> {
        let mut source = std::fs::File::open(input)?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(output)?);
        out.write_all(&htsjdk_bam::writer::write_bam_header_block_with(
            replacement,
            level,
            deflater,
        )?)?;
        htsjdk_bam::reheader::block_copy_to(&mut source, &mut out, true, false, level, deflater)
            .map_err(|e| format!("{e:?}"))?;
        out.flush()?;
        Ok(())
    })();
    if copied.is_err() {
        let _ = std::fs::remove_file(output);
    }
    copied
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("ReplaceSamHeader");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let header_path = arg(&args, "HEADER=").ok_or("HEADER= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);

    // `doWork` reads HEADER first; a BAM INPUT then goes to `reheaderBamFile`, which asserts it is
    // writable and that its last block is a valid one, compares the sort orders on the header
    // alone, and copies the record blocks without inflating them. Anything else, and a header that
    // does not read, takes the whole read below.
    // shortcut: CREATE_INDEX=true still re-encodes, to index as it writes; index the copy instead
    // once there is a streaming indexer.
    if is_bam(&input)? {
        let replacement = replacement_header(&header_path)?;
        check_writable(&input);
        if htsjdk_bgzf::check_termination(&tail(&input)?) == htsjdk_bgzf::FileTermination::Defective
        {
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: {} does not have a \
                 valid GZIP block at the end of the file.",
                picard_analysis::fingerprint::file_uri(&input)
            );
            std::process::exit(1);
        }
        if !create_index {
            if let Ok(input_header) = header_only(&input) {
                check_sort_order(&input_header, &replacement);
                return block_copy(&input, &output, &replacement);
            }
        }
    }

    let (input_file, input_is_bam) = read_any(&input)?;
    let (input_header, records) = input_file.into_parts();
    let replacement = replacement_header(&header_path)?;

    if input_is_bam {
        check_writable(&input);
        check_sort_order(&input_header, &replacement);
    } else if sort_order(&replacement) != sort_order(&input_header) {
        // The text path's own rule, in Picard: the two orders must be equal, `unsorted` included.
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: Sort orders of INPUT ({}) and \
             HEADER ({}) do not agree.",
            sort_order(&input_header),
            sort_order(&replacement)
        );
        std::process::exit(1);
    }

    let mut out = Output::create(&output, &replacement, false)?;
    for record in records {
        out.write(&record)?;
    }
    // `reheaderBamFile` is handed CREATE_INDEX, and the writer factory indexes a coordinate-sorted
    // BAM only.
    out.finish(
        create_index,
        create_index && sort_order(&replacement) == "coordinate",
    )?;
    Ok(())
}
