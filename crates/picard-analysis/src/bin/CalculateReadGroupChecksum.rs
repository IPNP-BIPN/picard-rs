//! `CalculateReadGroupChecksum` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CalculateReadGroupChecksum.doWork` at tag 3.4.0. The digest itself is
//! `picard_analysis::calculate_read_group_checksum`, over the header's read groups; the tool is a
//! wrapper around it, and the two things worth reproducing here are both about the wrapper.
//!
//! The output is the 32 hex characters and nothing else: `FileWriter.write(hashText)` with no
//! newline, so a port that used `println!` would differ from the reference by one byte on every
//! row.
//!
//! `OUTPUT` is optional, and when it is absent the file goes *next to the input*, named
//! `<input file name>.read_group_md5`. The arrays always pass one, but the port implements the
//! default because that is the path a reader of the source would expect to work, and because on a
//! read-only corpus it is the difference between a file and a refusal.
//!
//! `SAMUtils.calculateReadGroupRecordChecksum` opens the file itself, so only the header is read
//! and the records are never decoded: a BAM and the SAM it came from give the same digest.

use std::io::{Read, Write};

use htsjdk_bam::read_group_checksum::calculate_read_group_record_checksum;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::ValidationStringency;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    // `new File(INPUT.getParentFile(), INPUT.getName() + ".read_group_md5")`.
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .unwrap_or_else(|| format!("{input}.read_group_md5"));

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let read_groups = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        reader.header.text.read_groups.clone()
    } else {
        let text = String::from_utf8(raw)?;
        let (header, _) =
            read_sam_with(&text, ValidationStringency::Lenient).map_err(|e| format!("{e:?}"))?;
        header.read_groups
    };

    let digest = calculate_read_group_record_checksum(&read_groups);
    let mut out = std::fs::File::create(&output)?;
    // No newline: the reference writes the digest and closes the writer.
    out.write_all(digest.as_bytes())?;
    out.flush()?;
    Ok(())
}
