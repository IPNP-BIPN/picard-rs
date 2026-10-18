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
use picard_analysis::args::arg;
use picard_analysis::vcf_io::{
    die, extract_dictionary, is_compressed_vcf, set_sequence_dictionary, VcfReader,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("UpdateVcfSequenceDictionary");
    picard_analysis::output::install();
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
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let dictionary = extract_dictionary(&sequence_dictionary)?;
    let mut vcf = VcfReader::open(&input).unwrap_or_else(|exception| die(&exception));

    // A VCF without contig lines extracts to `null`. The builder takes that happily unless it is
    // to index, and then refuses; otherwise the refusal is `setSequenceDictionary`'s.
    // The builder's refusal is `VCF`'s alone: a `BLOCK_COMPRESSED_VCF` gets a tabix creator with
    // or without a dictionary.
    if create_index && dictionary.is_none() && !is_compressed_vcf(&output) {
        vcf.fail(
            "java.lang.IllegalArgumentException: A reference dictionary is required for creating \
             Tribble indices on the fly",
        );
    }
    let Some(dictionary) = dictionary else {
        vcf.fail(
            "java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \"dictionary\" is null",
        );
    };

    // `fileHeader.setSequenceDictionary(samSequenceDictionary)` on the reader's own header, which
    // is why a 4.3 file keeps its version into `writeHeader` and is refused there.
    let mut header = vcf.file.header.clone();
    set_sequence_dictionary(&mut header, &dictionary);
    if let Err(message) = reject_vcf_v43_headers(vcf.file.header_version) {
        vcf.fail(&format!("java.lang.IllegalArgumentException: {message}"));
    }

    let index_dictionary = if create_index {
        Some(dictionary.as_slice())
    } else {
        None
    };
    vcf.copy_to(&output, &header, index_dictionary, Ok);
    Ok(())
}
