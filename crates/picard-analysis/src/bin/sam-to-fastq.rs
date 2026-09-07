//! `SamToFastq` as a runnable binary: the covering array's port side.
//!
//! What the array measures here is mostly which reads survive: `INCLUDE_NON_PRIMARY_ALIGNMENTS`
//! and `INCLUDE_NON_PF_READS` decide that, and `RE_REVERSE` decides whether a negative-strand read
//! is put back in sequencing orientation before it is written -- which changes every base and
//! every quality of that read, not its presence.
//!
//! `INTERLEAVE` and `SECOND_END_FASTQ` are mutually exclusive, and the refusal comes from Barclay
//! rather than the tool: `Cannot set INTERLEAVE to true and pass in a SECOND_END_FASTQ`.
//!
//! A file of paired reads writes two files; the array compares the first, which is the one
//! `--FASTQ` names.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use picard_analysis::sam_to_fastq::{sam_to_fastq_paired, sam_to_fastq_unpaired, Options};

const READ_PAIRED: u16 = 0x1;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
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

    // `customCommandLineValidation`, which Barclay prints after the usage block.
    if flag("INTERLEAVE=", false) && second.is_some() {
        eprintln!("Cannot set INTERLEAVE to true and pass in a SECOND_END_FASTQ");
        std::process::exit(1);
    }

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

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

    // A file whose reads are paired writes two files; the tool decides per record, and this corpus
    // is one or the other throughout.
    let paired = records.iter().any(|rec| rec.flags & READ_PAIRED != 0);
    if paired {
        let (first, second_text) = sam_to_fastq_paired(&records, &options);
        std::fs::write(&fastq, first)?;
        if let Some(path) = second {
            std::fs::write(path, second_text)?;
        }
    } else {
        std::fs::write(&fastq, sam_to_fastq_unpaired(&records, &options))?;
    }
    Ok(())
}
