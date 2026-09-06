//! `FixMateInformation` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.FixMateInformation.doWork` at tag 3.4.0. The mate fixing, the supplemental
//! records' mate info and the two sorts live in `picard_analysis::fix_mate_information`.
//!
//! `SORT_ORDER`'s default is the trap: it is `null`, not a sort order, and a null means "the input
//! header's order" rather than "coordinate". The tool sorts into queryname order to group
//! templates whatever the caller asked for, fixes each template, and re-sorts into the output
//! order.
//!
//! `ADD_MATE_CIGAR` decides whether `MC` is written on the primaries and copied onto the
//! supplementary records, or cleared on both. `IGNORE_MISSING_MATES` decides whether a template
//! with one end is passed through or refused.

use std::io::{Read, Write};

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::sam_file::write_sam;
use picard_analysis::fix_mate_information::{
    fix_mate_information, fix_mate_information_to_bam, SortOrder,
};

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

    let sort_order = match arg(&args, "SORT_ORDER=").as_deref() {
        None => None,
        Some("coordinate") => Some(SortOrder::Coordinate),
        Some("queryname") => Some(SortOrder::Queryname),
        Some(other) => return Err(format!("unknown SORT_ORDER: {other}").into()),
    };

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

    if output.ends_with(".sam") {
        let sam = fix_mate_information(&text, sort_order).map_err(|e| format!("{e:?}"))?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        return Ok(());
    }
    let bam = fix_mate_information_to_bam(&text, sort_order).map_err(|e| format!("{e:?}"))?;
    std::fs::write(&output, &bam)?;
    Ok(())
}
