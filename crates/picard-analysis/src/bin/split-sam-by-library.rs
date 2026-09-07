//! `SplitSamByLibrary` as a runnable binary: the covering array's port side.
//!
//! One output per `@RG`'s `LB`, named after the library, in a directory -- plus `unknown.bam` for
//! records whose read group names no library. Each output carries the INPUT's whole header, read
//! groups of other libraries included, which is Picard's own choice and not a simplification here.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use picard_analysis::split_sam_by_library::{split_sam_by_library, split_sam_by_library_to_bam};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    // The outputs take the INPUT's format, as `reader.type().fileExtension()` gives it.
    let bam_input = raw.starts_with(&[0x1f, 0x8b]);
    let text = if bam_input {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header: SamHeader = reader.header.text.clone();
        let records: Vec<BamRecord> = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        write_sam(&header, &records).ok_or("records that parsed re-encode as SAM text")?
    } else {
        let text = String::from_utf8(raw)?;
        let (header, records) = read_sam(&text).map_err(|e| format!("{e:?}"))?;
        write_sam(&header, &records).ok_or("records that parsed re-encode as SAM text")?
    };

    std::fs::create_dir_all(&output)?;
    if bam_input {
        for (name, bytes) in split_sam_by_library_to_bam(&text).map_err(|e| format!("{e:?}"))? {
            std::fs::write(
                std::path::Path::new(&output).join(format!("{name}.bam")),
                bytes,
            )?;
        }
    } else {
        for (name, text) in split_sam_by_library(&text).map_err(|e| format!("{e:?}"))? {
            std::fs::write(
                std::path::Path::new(&output).join(format!("{name}.sam")),
                text,
            )?;
        }
    }
    Ok(())
}
