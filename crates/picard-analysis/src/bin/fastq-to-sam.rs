//! `FastqToSam` as a runnable binary: the covering array's port side.
//!
//! Eleven of the array's twelve rows are refusals, and the two kinds say different things:
//!
//! * `QUALITY_FORMAT` is a CLAIM about the file, and htsjdk checks it. The corpus is
//!   Standard-encoded, so claiming Illumina or Solexa is refused by `FastqQualityFormat`'s range
//!   check before a record is written -- the message names the format that was claimed, not the
//!   one that was found;
//! * `USE_SEQUENTIAL_FASTQS` reads the FASTQ's NAME rather than its content: the file must end in
//!   `_001` and one of five extensions so that the next file in the series can be found. The
//!   corpus's `reads_1.fastq` does not, and the refusal quotes the whole list it expected.
//!
//! The row that runs writes an unmapped, queryname-sorted BAM with one read group and no `@PG`.

use picard_analysis::fastq_to_sam::{fastq_to_sam_unpaired_to_bam, Options};

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

    let text = std::fs::read_to_string(&fastq)?;
    let format = arg(&args, "QUALITY_FORMAT=").unwrap_or_else(|| "Standard".to_string());
    if !matches!(format.as_str(), "Standard" | "Illumina" | "Solexa") {
        return Err(format!("unknown QUALITY_FORMAT: {format}").into());
    }
    if quality_out_of_range(&text, &format) {
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

    let mut options = Options::new(&sample);
    if let Some(name) = arg(&args, "READ_GROUP_NAME=") {
        options.read_group_name = name;
    }
    if let Some(order) = arg(&args, "SORT_ORDER=").or_else(|| arg(&args, "SO=")) {
        if !matches!(order.as_str(), "queryname" | "coordinate" | "unsorted") {
            return Err(format!("unknown SORT_ORDER: {order}").into());
        }
        options.sort_order = order;
    }
    let bam = fastq_to_sam_unpaired_to_bam(&text, &options).map_err(|e| format!("{e:?}"))?;
    std::fs::write(&output, bam)?;
    Ok(())
}
