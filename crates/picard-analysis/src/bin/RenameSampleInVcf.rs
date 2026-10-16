//! `RenameSampleInVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.RenameSampleInVcf.doWork` at tag 3.4.0. The read and write round trip is
//! `picard_analysis::vcf_io`.
//!
//! The rename touches only the header. Each record goes to the writer as the codec produced it,
//! and a file of one sample has its sample names "already sorted", so every genotype block is
//! still lazy and is copied as it was read: the writer never looks up the name the new header
//! gives the column, which is the only reason a header naming a sample no genotype is keyed by
//! writes anything at all.
//!
//! Two refusals, and a third that is not one. More than one sample is refused. An
//! `OLD_SAMPLE_NAME` other than the file's sample is refused and names the one found. A file with
//! NO sample is accepted, and its records -- which have no genotypes -- come out with a `GT`
//! column of no-calls for the new sample; but give it an `OLD_SAMPLE_NAME` and the check reads
//! the first of zero samples, which is the JVM's index-out-of-bounds.

use picard_analysis::vcf_io::{die, header_dictionary, read_path, write_output};

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
    let new_sample_name = arg(&args, "NEW_SAMPLE_NAME=").ok_or("NEW_SAMPLE_NAME= is required")?;
    // Barclay reads the value `null` as no value, which is this argument's default.
    let old_sample_name = arg(&args, "OLD_SAMPLE_NAME=").filter(|value| value != "null");
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(false);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    let samples = &vcf.file.header.samples;
    if samples.len() > 1 {
        die("java.lang.IllegalArgumentException: Input VCF must be single-sample.");
    }
    if let Some(old) = &old_sample_name {
        match samples.first() {
            None => die("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"),
            Some(found) if found != old => die(&format!(
                "java.lang.IllegalArgumentException: Input VCF did not contain expected sample. \
                 Contained: {found}"
            )),
            Some(_) => {}
        }
    }

    // `new VCFHeader(header.getMetaDataInInputOrder(), CollectionUtil.makeList(NEW_SAMPLE_NAME))`.
    let mut header = vcf.file.header.clone();
    header.samples = vec![new_sample_name];

    let dictionary = header_dictionary(&header);
    if create_index && dictionary.is_none() {
        die(
            "java.lang.IllegalArgumentException: A reference dictionary is required for creating \
             Tribble indices on the fly",
        );
    }
    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    write_output(&output, &header, &vcf.records, index_dictionary)?;
    Ok(())
}
