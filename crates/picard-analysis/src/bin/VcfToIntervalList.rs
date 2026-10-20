//! `VcfToIntervalList` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.VcfToIntervalList.doWork` at tag 3.4.0. The naming, the merge and the
//! rendering live in `picard_analysis::vcf_to_interval_list`.
//!
//! The output header is built from the VCF's contig lines, so a VCF without any has no dictionary
//! to give it. `new SAMFileHeader(null)` is accepted, and the refusal comes a step later, from the
//! header codec asking the missing dictionary for its sequences: a `NullPointerException`, not a
//! message of Picard's.

use picard_analysis::args::arg;
use picard_analysis::vcf_io::{die, header_dictionary, VcfReader};
use picard_analysis::vcf_to_interval_list::{push_interval, render_header, Merger, Site};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("VcfToIntervalList");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let include_filtered = arg(&args, "INCLUDE_FILTERED=")
        .or_else(|| arg(&args, "IF="))
        .map(|value| picard_analysis::args::is_true(&value))
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

    let mut vcf = VcfReader::open(&input).unwrap_or_else(|exception| die(&exception));
    let Some(dictionary) = header_dictionary(&vcf.file.header) else {
        vcf.fail(
            "java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because the return value of \
             \"htsjdk.samtools.SAMFileHeader.getSequenceDictionary()\" is null",
        );
    };

    // The output is held as text, not written as it comes: a read failure exits before anything
    // is written, as the whole-file read did. The sites and intervals are never held.
    let mut out = render_header(&dictionary);
    let mut merger = Merger::new(include_filtered, concatenate_ids);
    for record in vcf {
        let variant = &record.variant;
        // `getCommonInfo().getAttributeAsInt(END, getEnd())`.
        let end = variant
            .attributes
            .iter()
            .find(|(key, _)| key == "END")
            .and_then(|(_, value)| value.format())
            .and_then(|text| text.parse().ok())
            .unwrap_or(variant.stop);
        let site = Site {
            contig: variant.contig.clone(),
            start: variant.start,
            end,
            id: variant.id.clone(),
            filtered: variant.is_filtered(),
        };
        if let Some(interval) = merger.push(site) {
            push_interval(&mut out, &interval);
        }
    }
    if let Some(interval) = merger.finish() {
        push_interval(&mut out, &interval);
    }
    std::fs::write(&output, out)?;
    Ok(())
}
