//! `VcfFormatConverter` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.VcfFormatConverter.doWork` at tag 3.4.0 for VCF and block-compressed VCF,
//! read and written. BCF is not ported: htsjdk-rs has no BCF codec, so a `.bcf` on either side is
//! refused here with a message of the port's own.
//!
//! The records are copied under a copy of the input's header, so the conversion is the shared
//! round trip (`picard_analysis::vcf_io`) and what varies is the container. The output's format
//! is the writer builder's guess from the name: `.vcf.gz` and `.vcf.bgz` are BGZF, `.vcf` is text,
//! and any other name is the builder's `IllegalArgumentException`.
//!
//! # `REQUIRE_INDEX` looks for the index the reader would use
//!
//! `AbstractFeatureReader.getFeatureReader` takes a block-compressed file with a `.tbi` beside it
//! as tabix; anything else goes to the Tribble reader, which with `requireIndex` looks for a
//! `.idx` and refuses without one -- so the message names `.idx` even for a `.vcf.gz`.
//!
//! # A compressed output is indexed with tabix, which checks the order
//!
//! The plain writer's Tribble index takes records in any order. The tabix creator does not: a
//! record that starts before the one written before it on the same contig is an uncaught
//! `IllegalArgumentException` naming both, their file positions BGZF virtual offsets taken before
//! each line is written. Both are `vcf_io::write_output`, which writes the `.tbi` beside a
//! `.vcf.gz` and the `.idx` beside a `.vcf`.

use std::io::Write;

use picard_analysis::args::arg;
use picard_analysis::vcf_io::{
    die, header_dictionary, is_compressed_vcf, read_path, write_output, write_records, Record,
};

/// `IOUtil.hasBlockCompressedExtension`.
fn has_block_compressed_extension(path: &str) -> bool {
    [".gz", ".gzip", ".bgz", ".bgzf"]
        .iter()
        .any(|ext| path.to_lowercase().ends_with(ext))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("VcfFormatConverter");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let require_index = arg(&args, "REQUIRE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(true);
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if input.ends_with(".bcf") || output.ends_with(".bcf") {
        return Err("BCF is not ported: htsjdk-rs has no BCF codec".into());
    }

    // `new VCFFileReader(INPUT, REQUIRE_INDEX)`: the header is read before the index is looked for.
    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    if require_index {
        let tabix = has_block_compressed_extension(&input)
            && std::path::Path::new(&format!("{input}.tbi")).exists();
        if !tabix && !std::path::Path::new(&format!("{input}.idx")).exists() {
            let absolute = std::path::absolute(&input)?;
            die(&format!(
                "htsjdk.tribble.TribbleException: An index is required, but none found with file \
                 ending .idx, for input source: file://{}",
                absolute.display()
            ));
        }
    }

    // `new VCFHeader(reader.getFileHeader())`: the same lines and samples.
    let header = vcf.file.header.clone();
    let dictionary = header_dictionary(&header);
    if create_index && dictionary.is_none() {
        die(
            "picard.PicardException: A sequence dictionary must be available in the input file \
             when creating indexed output.",
        );
    }
    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };

    if output.ends_with(".vcf") || is_compressed_vcf(&output) {
        write_output(&output, &header, &vcf.records, index_dictionary)?;
    } else {
        die(
            "java.lang.IllegalArgumentException: Output format type is not set, or could not be \
             inferred from the output path. If a path was used, does it have a valid VCF \
             extension (.vcf, .vcf.gz, .vcf.bgz, .bcf)?",
        );
    }
    Ok(())
}
