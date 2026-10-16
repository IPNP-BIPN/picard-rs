//! `VcfToIntervalList` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.VcfToIntervalList.doWork` at tag 3.4.0. The naming, the merge and the
//! rendering live in `picard_analysis::vcf_to_interval_list`.
//!
//! The output header is built from the VCF's contig lines, so a VCF without any has no dictionary
//! to give it. `new SAMFileHeader(null)` is accepted, and the refusal comes a step later, from the
//! header codec asking the missing dictionary for its sequences: a `NullPointerException`, not a
//! message of Picard's.

use picard_analysis::vcf_io::{die, header_dictionary, read_path};
use picard_analysis::vcf_to_interval_list::{merge_intervals, render, to_intervals, Site};

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
    let include_filtered = arg(&args, "INCLUDE_FILTERED=")
        .or_else(|| arg(&args, "IF="))
        .map(|value| value == "true")
        .unwrap_or(false);
    let concatenate_ids = match arg(&args, "VARIANT_ID_METHOD=").as_deref() {
        None | Some("CONCAT_ALL") => true,
        Some("USE_FIRST") => false,
        Some(other) => return Err(format!("unknown VARIANT_ID_METHOD: {other}").into()),
    };
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    let Some(dictionary) = header_dictionary(&vcf.file.header) else {
        die("java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because the return value of \
             \"htsjdk.samtools.SAMFileHeader.getSequenceDictionary()\" is null");
    };

    let sites: Vec<Site> = vcf
        .records
        .iter()
        .map(|record| {
            let variant = &record.variant;
            // `getCommonInfo().getAttributeAsInt(END, getEnd())`.
            let end = variant
                .attributes
                .iter()
                .find(|(key, _)| key == "END")
                .and_then(|(_, value)| value.format())
                .and_then(|text| text.parse().ok())
                .unwrap_or(variant.stop);
            Site {
                contig: variant.contig.clone(),
                start: variant.start,
                end,
                id: variant.id.clone(),
                filtered: variant.is_filtered(),
            }
        })
        .collect();

    let intervals = merge_intervals(&to_intervals(&sites, include_filtered), concatenate_ids);
    std::fs::write(&output, render(&dictionary, &intervals))?;
    Ok(())
}
