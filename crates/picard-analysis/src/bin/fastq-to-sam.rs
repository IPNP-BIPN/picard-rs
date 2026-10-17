//! `FastqToSam` as a runnable binary: the covering array's port side.
//!
//! The whole tool is `picard_analysis::fastq_to_sam::fastq_to_bam_full`: one FASTQ or a pair, the
//! read group built from every field the command line gives, the comments, the quality range
//! `MIN_Q..MAX_Q` checked read by read, and the sort order. Before any of that, two refusals:
//!
//! * `QUALITY_FORMAT` is a CLAIM about the file, and htsjdk checks it. The corpus is
//!   Standard-encoded, so claiming Illumina or Solexa is refused by `FastqQualityFormat`'s range
//!   check before a record is written -- the message names the format that was claimed, not the
//!   one that was found;
//! * `USE_SEQUENTIAL_FASTQS` reads the FASTQ's NAME rather than its content: the file must end in
//!   `_001` and one of five extensions so that the next file in the series can be found. The
//!   corpus's `reads_1.fastq` does not, and the refusal quotes the whole list it expected.

use picard_analysis::add_or_replace_read_groups::iso8601;
use picard_analysis::fastq_to_sam::{fastq_to_bam_full, FullOptions};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `FastqQualityFormat`'s range check, as `QualityEncodingDetector` applies it: a Standard file
/// read as Illumina or Solexa has quality bytes below those encodings' floor.
fn quality_out_of_range(text: &str, format: &str) -> bool {
    // Standard is 33-based, Illumina and Solexa 64-based, so any character below `;` (59) proves
    // the file is not in a 64-based encoding.
    if format == "Standard" {
        return false;
    }
    text.lines()
        .skip(3)
        .step_by(4)
        .flat_map(str::bytes)
        .any(|byte| byte < b';')
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let fastq = arg(&args, "FASTQ=")
        .or_else(|| arg(&args, "F1="))
        .ok_or("FASTQ= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let sample = arg(&args, "SAMPLE_NAME=").ok_or("SAMPLE_NAME= is required")?;
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let int = |key: &str, default: i32| -> Result<i32, Box<dyn std::error::Error>> {
        Ok(arg(&args, key)
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(default))
    };
    let min_q = int("MIN_Q=", 0)?;
    let max_q = int("MAX_Q=", 93)?;
    // `customCommandLineValidation`, one message and the first that applies.
    if min_q < 0 {
        eprintln!("MIN_Q must be >= 0");
        std::process::exit(1);
    }
    if max_q > 93 {
        eprintln!("MAX_Q must be <= 93");
        std::process::exit(1);
    }
    let text = std::fs::read_to_string(&fastq)?;
    let second = arg(&args, "FASTQ2=").or_else(|| arg(&args, "F2="));
    let second_text = match &second {
        Some(path) => Some(std::fs::read_to_string(path)?),
        None => None,
    };
    let format = arg(&args, "QUALITY_FORMAT=").unwrap_or_else(|| "Standard".to_string());
    if !matches!(format.as_str(), "Standard" | "Illumina" | "Solexa") {
        return Err(format!("unknown QUALITY_FORMAT: {format}").into());
    }
    let empty = text.trim().is_empty();
    if empty && !flag("ALLOW_EMPTY_FASTQ=", false) {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: Input fastq is empty. Set ALLOW_EMPTY_FASTQ if you still want to write out a file with no reads."
        );
        std::process::exit(1);
    }
    let both = format!("{text}{}", second_text.as_deref().unwrap_or(""));
    if !empty && quality_out_of_range(&both, &format) {
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMException: The quality values do not \
             fall in the range appropriate for the expected quality of {format}."
        );
        std::process::exit(1);
    }

    // `getSequentialFileList`, which reads the NAME rather than the content. It is checked
    // AFTER the quality format, which is why a row that gets both wrong is refused for the
    // quality and not for the name.
    if flag("USE_SEQUENTIAL_FASTQS=", false) {
        let name = fastq
            .rsplit('/')
            .next()
            .unwrap_or(&fastq)
            .to_ascii_uppercase();
        let sequential = [
            "_001.FASTQ",
            "_001.FASTQ.GZ",
            "_001.FQ",
            "_001.FQ.GZ",
            "_001.BFQ",
        ]
        .iter()
        .any(|suffix| name.ends_with(suffix));
        if !sequential {
            eprintln!(
                "Exception in thread \"main\" picard.PicardException: Could not parse the FASTQ \
                 extension (expected '_001' + '[FASTQ, FASTQ_GZ, FQ, FQ_GZ, BFQ]'): {fastq}"
            );
            std::process::exit(1);
        }
    }

    let sort_order = arg(&args, "SORT_ORDER=")
        .or_else(|| arg(&args, "SO="))
        .unwrap_or_else(|| "queryname".to_string());
    if !matches!(sort_order.as_str(), "queryname" | "coordinate" | "unsorted") {
        return Err(format!("unknown SORT_ORDER: {sort_order}").into());
    }
    let run_date = match arg(&args, "RUN_DATE=").or_else(|| arg(&args, "DT=")) {
        Some(value) => Some(iso8601(&value).ok_or(format!("unreadable RUN_DATE: {value}"))?),
        None => None,
    };
    let options = FullOptions {
        read_group_name: arg(&args, "READ_GROUP_NAME=")
            .or_else(|| arg(&args, "RG="))
            .unwrap_or_else(|| "A".to_string()),
        sample_name: sample,
        sort_order,
        library_name: arg(&args, "LIBRARY_NAME=").or_else(|| arg(&args, "LB=")),
        platform: arg(&args, "PLATFORM=").or_else(|| arg(&args, "PL=")),
        platform_unit: arg(&args, "PLATFORM_UNIT=").or_else(|| arg(&args, "PU=")),
        sequencing_center: arg(&args, "SEQUENCING_CENTER=").or_else(|| arg(&args, "CN=")),
        predicted_insert_size: arg(&args, "PREDICTED_INSERT_SIZE=")
            .or_else(|| arg(&args, "PI="))
            .map(|v| v.parse())
            .transpose()?,
        description: arg(&args, "DESCRIPTION=").or_else(|| arg(&args, "DS=")),
        run_date,
        platform_model: arg(&args, "PLATFORM_MODEL=").or_else(|| arg(&args, "PM=")),
        program_group: arg(&args, "PROGRAM_GROUP=").or_else(|| arg(&args, "PG=")),
        comments: args
            .iter()
            .filter_map(|a| a.strip_prefix("COMMENT=").or_else(|| a.strip_prefix("CO=")))
            .map(str::to_string)
            .collect(),
        min_q,
        max_q,
        quality_format: format,
        allow_and_ignore_empty_lines: flag("ALLOW_AND_IGNORE_EMPTY_LINES=", false),
    };
    match fastq_to_bam_full(&text, second_text.as_deref(), &options) {
        Ok(bam) => picard_analysis::cram_output::write(&output, bam)?,
        Err(thrown) => {
            eprintln!("Exception in thread \"main\" {thrown}");
            std::process::exit(1);
        }
    }
    Ok(())
}
