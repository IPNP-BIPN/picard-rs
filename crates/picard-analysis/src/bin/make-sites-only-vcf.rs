//! `MakeSitesOnlyVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.MakeSitesOnlyVcf.doWork` at tag 3.4.0. The read and write round trip is
//! `picard_analysis::vcf_io`.
//!
//! The name undersells the tool: it keeps the genotypes of every sample named by `SAMPLE`, and
//! only with the default empty set is the output sites-only. The kept samples are a `TreeSet`, so
//! the output's columns are in sorted order whatever order the input had them in, and a name the
//! input does not have still gets a column: `subsetToSamples` skips it, and the writer, which walks
//! the HEADER's samples, fills the gap with a no-call at the record's ploidy.
//!
//! `subsetToSamples` reads the genotypes by name, which decodes them, so no block is ever copied
//! here: the kept columns are re-encoded (FORMAT keys sorted after `GT`, trailing missing fields
//! trimmed) even when the input's sample names were sorted.
//!
//! `CREATE_INDEX` defaults to true in this tool's constructor, and an index needs a dictionary, so
//! an input with no contig lines is refused before anything is written -- unless `CREATE_INDEX` is
//! false, in which case the same file is accepted.

use picard_analysis::vcf_io::{
    die, header_dictionary, java_compare, read_path, write_output, Record,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// A collection argument: every occurrence appends, and the value `null` empties it (Barclay).
fn collection(args: &[String], keys: &[&str]) -> Vec<String> {
    let mut values = Vec::new();
    for a in args {
        if let Some(value) = keys.iter().find_map(|key| a.strip_prefix(key)) {
            if value == "null" {
                values.clear();
            } else {
                values.push(value.to_string());
            }
        }
    }
    values
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    // `Set<String> SAMPLE = new TreeSet<String>()`: sorted and without repeats.
    let mut samples = collection(&args, &["SAMPLE=", "S="]);
    samples.sort_by(|a, b| java_compare(a, b));
    samples.dedup();

    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    let dictionary = header_dictionary(&vcf.file.header);
    if create_index && dictionary.is_none() {
        die(
            "picard.PicardException: A sequence dictionary must be available (either through the \
             input file or by setting it explicitly) when creating indexed output.",
        );
    }

    // `new VCFHeader(inputVcfHeader.getMetaDataInInputOrder(), SAMPLE)`: the input's lines, and
    // the requested samples as the columns.
    let mut header = vcf.file.header.clone();
    header.samples = samples.clone();

    let records: Vec<Record> = vcf
        .records
        .into_iter()
        .map(|mut record| {
            // `ctx.getGenotypes().subsetToSamples(samples)`: an empty set is `NO_GENOTYPES`, and
            // otherwise each named sample the record has, looked up by name.
            record.decode();
            let genotypes = std::mem::take(&mut record.variant.genotypes);
            record.variant.genotypes = samples
                .iter()
                .filter_map(|sample| genotypes.iter().find(|g| &g.sample_name == sample))
                .cloned()
                .collect();
            record
        })
        .collect();

    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    write_output(&output, &header, &records, index_dictionary)?;
    Ok(())
}
