//! `UpdateVcfSequenceDictionary` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.UpdateVcfSequenceDictionary.doWork` at tag 3.4.0. The dictionary extraction,
//! the header rewrite and the read and write round trip are `picard_analysis::vcf_io`.
//!
//! `VCFHeader.setSequenceDictionary` removes every contig line and adds one per sequence, and a
//! line built from a `SAMSequenceRecord` carries `ID`, `length` and -- only when the record has
//! one -- `assembly`, nothing else: a `.dict`'s `M5`, `UR` and `SP` do not survive the trip. The
//! dictionary is taken as it is, never compared with the records: a record on a contig the new
//! dictionary lacks is written anyway.
//!
//! The records go to the writer as the codec produced them, under the file's own header object,
//! so the sample columns keep their order and a file whose names were already sorted has its
//! genotype blocks copied as read.

use htsjdk_vcf::vcf_file::reject_vcf_v43_headers;
use picard_analysis::vcf_io::{
    die, extract_dictionary, read_path, set_sequence_dictionary, write_output,
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
    let sequence_dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "SD="))
        .ok_or("SEQUENCE_DICTIONARY= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(false);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let dictionary = extract_dictionary(&sequence_dictionary)?;
    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));

    // A VCF without contig lines extracts to `null`. The builder takes that happily unless it is
    // to index, and then refuses; otherwise the refusal is `setSequenceDictionary`'s.
    if create_index && dictionary.is_none() {
        die(
            "java.lang.IllegalArgumentException: A reference dictionary is required for creating \
             Tribble indices on the fly",
        );
    }
    let Some(dictionary) = dictionary else {
        die(
            "java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \"dictionary\" is null",
        );
    };

    // `fileHeader.setSequenceDictionary(samSequenceDictionary)` on the reader's own header, which
    // is why a 4.3 file keeps its version into `writeHeader` and is refused there.
    let mut header = vcf.file.header.clone();
    set_sequence_dictionary(&mut header, &dictionary);
    if let Err(message) = reject_vcf_v43_headers(vcf.file.header_version) {
        die(&format!("java.lang.IllegalArgumentException: {message}"));
    }

    let index_dictionary = if create_index {
        Some(dictionary.as_slice())
    } else {
        None
    };
    write_output(&output, &header, &vcf.records, index_dictionary)?;
    Ok(())
}
