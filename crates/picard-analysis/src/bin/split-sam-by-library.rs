//! `SplitSamByLibrary` as a runnable binary: the covering array's port side.
//!
//! One output per `@RG`'s `LB`, named after the library, in a directory -- plus `unknown.bam` for
//! records whose read group names no library. Each output carries the INPUT's whole header, read
//! groups of other libraries included, which is Picard's own choice and not a simplification here.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use picard_analysis::args::arg;
use picard_analysis::split_sam_by_library::{split_sam_by_library, split_sam_by_library_to_bam};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SplitSamByLibrary");
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
    let cram_input = picard_analysis::cram_input::is_cram(&raw);
    let raw = picard_analysis::cram_input::uncram(raw);
    // The outputs take the INPUT's format, as `reader.type().fileExtension()` gives it.
    let bam_input = raw.starts_with(&[0x1f, 0x8b]);
    let text = if bam_input {
        let plain = picard_analysis::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
        let reader =
            picard_analysis::record_validation::reader(&plain).map_err(|e| format!("{e:?}"))?;
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
    // A CRAM in is CRAMs out, `makeWriter` with the reference, and `CREATE_INDEX` is honoured by
    // `CRAMBAIIndexer`, which refuses an output that is not coordinate-sorted.
    if cram_input {
        let create_index = arg(&args, "CREATE_INDEX=").as_deref() == Some("true");
        for (name, header, records) in picard_analysis::split_sam_by_library::split_outputs(&text)
            .map_err(|e| format!("{e:?}"))?
        {
            let path = std::path::Path::new(&output).join(format!("{name}.cram"));
            picard_analysis::cram_output::write_records_indexed(
                &path.display().to_string(),
                &header,
                &records,
                create_index,
            )?;
        }
        return Ok(());
    }
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
