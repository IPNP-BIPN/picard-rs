//! `BuildBamIndex` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.BuildBamIndex.doWork` at tag 3.4.0. The index itself is
//! `picard_analysis::build_bam_index`, which is htsjdk's `BAMIndexer.createIndex`; this is the
//! file around it and the two refusals the tool makes before it indexes anything, in their order:
//! a SAM text file is refused by TYPE, and a BAM that is not coordinate-sorted by its header's
//! `SO`, whatever order its records are actually in.
//!
//! With no `OUTPUT` the index goes to the working directory, named after the input's FILE NAME
//! rather than beside it: `x.bai` for `x.bam`, `<name>.bai` otherwise.

use picard_analysis::build_bam_index::build_bam_index;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn refuse(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" htsjdk.samtools.SAMException: {message}");
    std::process::exit(1);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .filter(|value| value != "null")
        .unwrap_or_else(|| {
            let name = std::path::Path::new(&input)
                .file_name()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_default();
            match name.strip_suffix(".bam") {
                Some(stem) => format!("{stem}.bai"),
                None => format!("{name}.bai"),
            }
        });

    let raw = std::fs::read(&input)?;
    if !raw.starts_with(&[0x1f, 0x8b]) {
        refuse("Input file must be bam file, not sam file.");
    }
    let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
    let reader = htsjdk_bam::reader::BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
    if reader.header.text.attributes.get("SO") != Some("coordinate") {
        refuse("Input bam file must be sorted by coordinate");
    }
    let bai = build_bam_index(&raw).map_err(|e| format!("{e:?}"))?;
    std::fs::write(&output, bai)?;
    Ok(())
}
