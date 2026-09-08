//! `SamToFastqWithTags` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.SamToFastqWithTags.doWork` at tag 3.4.0, which is `SamToFastq`'s with one
//! addition: beside the read FASTQ it writes one FASTQ per `SEQUENCE_TAG_GROUP`, whose reads are
//! built from TAG VALUES rather than from bases. Those files are named after the group, with the
//! commas turned into underscores, and land beside the `FASTQ` the caller named.
//!
//! The array compares a tag file, which is what this tool adds. What it mostly measures is
//! refusals, and the three are worth naming:
//!
//!  * `INTERLEAVE` with `SECOND_END_FASTQ` is refused by Barclay before the tool runs;
//!  * a read missing a tag a group names is a `PicardException`, per read;
//!  * and a file whose reads did not all pair up ends in a `MATE_NOT_FOUND` validation error --
//!    which `VALIDATION_STRINGENCY` decides the fate of. Under `STRICT` it throws; under `LENIENT`
//!    or `SILENT` the tool writes its files and returns zero, so the same corpus is an answer or a
//!    refusal depending on an argument that touches nothing else.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use picard_analysis::sam_to_fastq::{
    sam_to_fastq_paired_counting_leftovers, sam_to_fastq_unpaired, Options,
};
use picard_analysis::sam_to_fastq_with_tags::{
    sam_to_fastq_with_tags_paired_with, sam_to_fastq_with_tags_unpaired_with, TagGroup,
};

const READ_PAIRED: u16 = 0x1;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn args_of(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let fastq = arg(&args, "FASTQ=")
        .or_else(|| arg(&args, "F="))
        .ok_or("FASTQ= is required")?;
    let second = arg(&args, "SECOND_END_FASTQ=").or_else(|| arg(&args, "F2="));
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);

    if flag("INTERLEAVE=", false) && second.is_some() {
        eprintln!("Cannot set INTERLEAVE to true and pass in a SECOND_END_FASTQ");
        std::process::exit(1);
    }

    let stringency = arg(&args, "VALIDATION_STRINGENCY=").unwrap_or_else(|| "STRICT".to_string());
    if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
        return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
    }

    let sequence_groups = args_of(&args, "SEQUENCE_TAG_GROUP=");
    let quality_groups = args_of(&args, "QUALITY_TAG_GROUP=");
    let separators = args_of(&args, "TAG_GROUP_SEPERATOR=");
    let groups: Vec<TagGroup> = sequence_groups
        .iter()
        .enumerate()
        .map(|(index, spec)| {
            let mut group = TagGroup::new(spec);
            if let Some(quality) = quality_groups.get(index) {
                group = group.with_quality(quality);
            }
            if let Some(separator) = separators.get(index) {
                group = group.with_separator(separator);
            }
            group
        })
        .collect();

    let options = Options {
        re_reverse: flag("RE_REVERSE=", true),
        include_non_primary: flag("INCLUDE_NON_PRIMARY_ALIGNMENTS=", false),
        include_non_pf: flag("INCLUDE_NON_PF_READS=", false),
    };

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let records: Vec<BamRecord> = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        BamReader::new(&plain)
            .map_err(|e| format!("{e:?}"))?
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?.1
    };

    // The tag files land beside the read FASTQ, which is what `new File(FASTQ.getParent(), name)`
    // means when no per-read-group output was asked for.
    let directory = std::path::Path::new(&fastq)
        .parent()
        .map(std::path::Path::to_path_buf)
        .unwrap_or_default();

    let paired = records.iter().any(|rec| rec.flags & READ_PAIRED != 0);
    let leftovers = if paired {
        let (first, second_text, leftovers) =
            sam_to_fastq_paired_counting_leftovers(&records, &options);
        let tag_files = match sam_to_fastq_with_tags_paired_with(&records, &groups, &options) {
            Ok(files) => files,
            Err(missing) => {
                eprintln!(
                    "Exception in thread \"main\" picard.PicardException: Record: {} does have a value for tag: {}",
                    missing.read_name, missing.tag
                );
                std::process::exit(1);
            }
        };
        std::fs::write(&fastq, first)?;
        if let Some(path) = second {
            std::fs::write(path, second_text)?;
        }
        for (name, text) in tag_files {
            std::fs::write(directory.join(name), text)?;
        }
        leftovers
    } else {
        let tag_files = match sam_to_fastq_with_tags_unpaired_with(&records, &groups, &options) {
            Ok(files) => files,
            Err(missing) => {
                eprintln!(
                    "Exception in thread \"main\" picard.PicardException: Record: {} does have a value for tag: {}",
                    missing.read_name, missing.tag
                );
                std::process::exit(1);
            }
        };
        std::fs::write(&fastq, sam_to_fastq_unpaired(&records, &options))?;
        for (name, text) in tag_files {
            std::fs::write(directory.join(name), text)?;
        }
        0
    };

    // `SAMUtils.processValidationError`, which is where the stringency finally decides something:
    // STRICT throws after everything has been written, and the other two do not.
    if leftovers > 0 && stringency == "STRICT" {
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMFormatException: SAM validation error: ERROR::MATE_NOT_FOUND:Found {leftovers} unpaired mates"
        );
        std::process::exit(1);
    }
    Ok(())
}
