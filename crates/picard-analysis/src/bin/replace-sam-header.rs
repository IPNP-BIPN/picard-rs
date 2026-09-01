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

use std::io::{Read, Write};

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::writer::BamWriter;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
    }
}

/// Read a SAM or BAM, told apart by the gzip magic as `SamReaderFactory` tells them apart.
fn read_any(path: &str) -> Result<(SamHeader, Vec<BamRecord>, bool), String> {
    let mut raw = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| e.to_string())?
        .read_to_end(&mut raw)
        .map_err(|e| e.to_string())?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((header, records, true))
    } else {
        let text = String::from_utf8(raw).map_err(|e| e.to_string())?;
        // `standardReheader` opens INPUT at SILENT whatever VALIDATION_STRINGENCY says.
        let (header, records) =
            read_sam_with(&text, ValidationStringency::Silent).map_err(|e| format!("{e:?}"))?;
        Ok((header, records, false))
    }
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let header_path = arg(&args, "HEADER=").ok_or("HEADER= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(false);

    let (input_header, records, input_is_bam) = read_any(&input)?;
    let (replacement, _, _) = read_any(&header_path)?;

    if input_is_bam {
        // `reheaderBamFile` asserts the INPUT is readable and then WRITABLE, before it copies a
        // block, because the block copy is written for the in-place case as well. A corpus mounted
        // read-only therefore refuses every BAM row, and the message names the file.
        if std::fs::OpenOptions::new()
            .write(true)
            .open(&input)
            .is_err()
        {
            // `File.getAbsolutePath()`: the working directory joined to a relative path, with no
            // symlink resolution. Canonicalizing instead would rewrite /var to /private/var on
            // macOS and print a path the reference never would.
            let absolute = if std::path::Path::new(&input).is_absolute() {
                input.clone()
            } else {
                std::env::current_dir()
                    .map(|cwd| cwd.join(&input).display().to_string())
                    .unwrap_or_else(|_| input.clone())
            };
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: \
                 File exists but is not writable: {absolute}"
            );
            std::process::exit(1);
        }
        // The block copy has a sort-order check of its own, in htsjdk rather than in Picard, and
        // it is not the same rule as the text path's: an `unsorted` new header is accepted against
        // any original, and the message names what the header would have to be.
        if sort_order(&replacement) != "unsorted"
            && sort_order(&replacement) != sort_order(&input_header)
        {
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: \
                 Sort order of new header does not match the original file, needs to be {}",
                sort_order(&input_header)
            );
            std::process::exit(1);
        }
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

    if output.ends_with(".sam") {
        let sam = write_sam(&replacement, &records).ok_or("records failed to re-encode as SAM")?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        return Ok(());
    }

    let mut writer = BamWriter::new(Vec::new(), &replacement)?;
    for record in &records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    let bam = writer.finish()?;
    std::fs::write(&output, &bam)?;
    // `reheaderBamFile` is handed CREATE_INDEX, and the writer factory indexes a coordinate-sorted
    // BAM only.
    if create_index && sort_order(&replacement) == "coordinate" {
        let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
        std::fs::write(index_path(&output), index)?;
    }
    Ok(())
}
