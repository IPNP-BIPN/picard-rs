//! `GatherBamFiles` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.GatherBamFiles.doWork` at tag 3.4.0. Two paths, chosen by the inputs' NAMES:
//!
//! * every input ends in `.bam`: `gatherWithBlockCopying`, the library's block copy, which keeps the
//!   first file's header and the others' compressed blocks as they are;
//! * any other: `gatherNormally`, every record read and written again under the first file's
//!   header by `SAMFileWriterFactory.makeWriter`, which picks the format from the OUTPUT's name: SAM
//!   text for `.sam`, CRAM for `.cram`, and BAM for anything else, `output.txt` included.
//!
//! `IOUtil.unrollFiles` expands an input that is not a SAM, BAM or CRAM by name into the files it
//! lists, one per line.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::gather_bam_files::gather_bam_files;

/// `IOUtil.unrollFiles(inputs, ".bam", ".sam", ".cram")`.
fn unroll(inputs: &[String]) -> std::io::Result<Vec<String>> {
    let mut out = Vec::new();
    for input in inputs {
        if [".bam", ".sam", ".cram"]
            .iter()
            .any(|ext| input.ends_with(ext))
        {
            out.push(input.clone());
        } else {
            let text = std::fs::read_to_string(input)?;
            let listed: Vec<String> = text
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty() && !line.starts_with('#'))
                .map(str::to_string)
                .collect();
            out.extend(unroll(&listed)?);
        }
    }
    Ok(out)
}

fn read_any(path: &str) -> Result<(SamHeader, Vec<BamRecord>), Box<dyn std::error::Error>> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((header, records))
    } else {
        Ok(read_sam(&String::from_utf8(raw)?).map_err(|e| format!("{e:?}"))?)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let inputs: Vec<String> = args
        .iter()
        .filter_map(|a| {
            a.strip_prefix("INPUT=")
                .or_else(|| a.strip_prefix("I="))
                .map(str::to_string)
        })
        .collect();
    let output = args
        .iter()
        .find_map(|a| a.strip_prefix("OUTPUT=").or_else(|| a.strip_prefix("O=")))
        .ok_or("OUTPUT= is required")?
        .to_string();
    let inputs = unroll(&inputs)?;
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }

    if inputs.iter().all(|input| input.ends_with(".bam")) {
        let raw: Vec<Vec<u8>> = inputs.iter().map(std::fs::read).collect::<Result<_, _>>()?;
        let borrowed: Vec<&[u8]> = raw.iter().map(Vec::as_slice).collect();
        let gathered = gather_bam_files(&borrowed).map_err(|e| format!("{e:?}"))?;
        std::fs::write(&output, gathered)?;
        return Ok(());
    }

    let (header, _) = read_any(&inputs[0])?;
    let mut records = Vec::new();
    for input in &inputs {
        records.extend(read_any(input)?.1);
    }
    if output.ends_with(".cram") {
        return Err("a CRAM output is not ported".into());
    }
    if output.ends_with(".sam") {
        let sam = write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?;
        std::fs::write(&output, sam)?;
    } else {
        let mut writer = BamWriter::new(Vec::new(), &header)?;
        for record in &records {
            writer.write(record).map_err(|e| format!("{e:?}"))?;
        }
        std::fs::write(&output, writer.finish()?)?;
    }
    Ok(())
}
