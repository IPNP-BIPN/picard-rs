//! `DownsampleSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.DownsampleSam.doWork` at tag 3.4.0 for the default `ConstantMemory` strategy,
//! which is what `picard_analysis::downsample_sam` implements.
//!
//! The two arguments that decide the output are `PROBABILITY` and `RANDOM_SEED`, and the second is
//! not decoration: the keep decision is a `Murmur3` hash of the READ NAME seeded with it, so two
//! seeds give two different subsets of the same size distribution, and both mates of a template
//! share a decision because they share a name.
//!
//! `STRATEGY` is held at `ConstantMemory` by the array rather than varied. `HighAccuracy` and
//! `Chained` buffer and make a second pass, and this port does not have them: a row that named one
//! would measure the argument parser refusing rather than the tool downsampling. The hold is
//! declared in the fixtures file with that reason, which is the same rule every other hold there
//! follows.
//!
//! The `@PG` record `doWork` adds carries the command line, temporary paths and all, and is
//! canonicalized away by the comparison exactly as the metrics tools' header is. This binary does
//! not write one, and the claim is over the surviving records and the rest of the header.

use std::io::{Read, Write};

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::sam_file::write_sam;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::downsample_sam::{downsample_sam, DEFAULT_SEED};

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

    if let Some(strategy) = arg(&args, "STRATEGY=") {
        if strategy != "ConstantMemory" {
            return Err(format!("STRATEGY={strategy} is not ported").into());
        }
    }

    let probability: f64 = arg(&args, "PROBABILITY=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(1.0);
    let seed: i32 = arg(&args, "RANDOM_SEED=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(DEFAULT_SEED);

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let text = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?
    } else {
        String::from_utf8(raw)?
    };

    let sam = downsample_sam(&text, probability, seed).map_err(|e| format!("{e:?}"))?;

    if output.ends_with(".sam") {
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        return Ok(());
    }
    // `SAMFileWriterFactory.makeWriter` again: anything that is not `.sam` is a BAM.
    let (header, records) = htsjdk_bam::sam_file::read_sam_with(
        &sam,
        htsjdk_bam::text_parse::ValidationStringency::Lenient,
    )
    .map_err(|e| format!("{e:?}"))?;
    let mut writer = BamWriter::new(Vec::new(), &header).map_err(|e| format!("{e:?}"))?;
    for record in &records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;
    Ok(())
}
