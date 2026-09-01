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

use std::io::{Read, Write};

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::clean_sam::clean_records;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `BAMFileWriter.createBamIndex`: the `.bam` suffix is dropped, anything else is kept.
fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
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
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(false);
    // `SamReaderFactory` at the stringency `doWork` asks for: STRICT becomes LENIENT, so that the
    // reader accepts the unmapped read with a nonzero MAPQ that the tool is there to zero.
    let stringency = match arg(&args, "VALIDATION_STRINGENCY=").as_deref() {
        None | Some("STRICT") | Some("LENIENT") => ValidationStringency::Lenient,
        Some("SILENT") => ValidationStringency::Silent,
        Some(other) => return Err(format!("unknown VALIDATION_STRINGENCY: {other}").into()),
    };

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    // BAM and SAM are told apart by the gzip magic, as `SamReaderFactory` tells them apart by
    // sniffing the stream rather than by the file's name.
    let (header, mut records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let text = reader.header.text.clone();
        let decoded = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        (text, decoded)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam_with(&text, stringency).map_err(|e| format!("{e:?}"))?
    };

    clean_records(&header, &mut records);

    // The header is written back unchanged: CleanSam hands the reader's header straight to the
    // writer, adds no `@PG` of its own, and does not re-sort, so `SO` still describes the records.
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
    // Only a coordinate-sorted file is indexed, and asking for an index on any other order is not
    // an error: `initializeBAMWriter` simply does not enable the indexer.
    if create_index && header.attributes.get("SO") == Some("coordinate") {
        let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
        std::fs::write(index_path(&output), index)?;
    }
    Ok(())
}
