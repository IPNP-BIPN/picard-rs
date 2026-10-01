//! `VcfToAdpc` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.VcfToAdpc.doWork` at tag 3.4.0, with `IlluminaAdpcFileWriter` and
//! `InfiniumDataFile`'s writers, whose layout is `picard_analysis::vcf_to_adpc`.
//!
//! The file is SAMPLE-major: each input is walked once per genotype sample, in the header's column
//! order, so a two-sample VCF is every locus of its first sample and then every locus of its
//! second. The record count of every walk must agree with the first one's, across files.
//!
//! Every failure after the arguments is caught by the tool itself: it logs the exception through
//! `log.error(e)` -- an `ERROR` line with an empty message and the stack trace after it -- and
//! returns 1, whatever the exception was. Only the opening checks (`unrollFiles` and the
//! readable/writable assertions) escape as an uncaught exception.

use htsjdk_vcf::variant::{Genotype, Value, VariantContext};
use picard_analysis::java_number::{parse_float, parse_int};
use picard_analysis::metrics_cli::Args;
use picard_analysis::vcf_io::{die, log_error, read_path, unroll_paths};
use picard_analysis::vcf_to_adpc::{
    write_record, IlluminaGenotype, Record, HEADER, MAX_UNSIGNED_SHORT,
};

const TOOL: &str = "VcfToAdpc";

/// `Object.toString()` of an attribute value as the codec stored it.
fn java_string(value: &Value) -> String {
    match value {
        Value::Str(s) => s.clone(),
        Value::Int(i) => i.to_string(),
        Value::Double(d) => d.to_string(),
        Value::Bool(b) => b.to_string(),
        Value::Missing => ".".to_string(),
        Value::List(items) => format!(
            "[{}]",
            items.iter().map(java_string).collect::<Vec<_>>().join(", ")
        ),
    }
}

fn info_attribute(context: &VariantContext, key: &str) -> Result<String, String> {
    context
        .attributes
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| java_string(v))
        .ok_or_else(|| {
            format!(
                "picard.PicardException: Unable to find attribute {key} in VCF.  Is this an \
                 Arrays VCF file?"
            )
        })
}

fn number_format(message: String) -> String {
    format!("java.lang.NumberFormatException: {message}")
}

/// `getUnsignedShortAttributeAsInt`.
fn unsigned_short(genotype: &Genotype, key: &str) -> Result<i32, String> {
    let value = genotype.get(key).map(java_string).ok_or_else(|| {
        format!(
            "picard.PicardException: Unable to find attribute {key} in VCF Genotype field.  Is \
             this an Arrays VCF file?"
        )
    })?;
    let value = parse_int(&value).map_err(number_format)?;
    if value < 0 {
        return Err(format!(
            "picard.PicardException: Value for key {key} ({value}) is <= 0!  Invalid value for \
             unsigned int"
        ));
    }
    if value > MAX_UNSIGNED_SHORT {
        eprintln!(
            "WARNING\t{TOOL}\tValue for key {key} ({value}) is > {MAX_UNSIGNED_SHORT} (truncating it)"
        );
        return Ok(MAX_UNSIGNED_SHORT);
    }
    Ok(value)
}

/// `getFloatAttribute(Genotype, key)`: absent or `?` is null, which the record writes as NaN.
fn optional_float(genotype: &Genotype, key: &str) -> Result<f32, String> {
    match genotype.get(key).map(java_string) {
        None => Ok(f32::NAN),
        Some(value) if value == "?" => Ok(f32::NAN),
        Some(value) => parse_float(&value).map_err(number_format),
    }
}

/// `Allele.basesMatch(String)`: a symbolic allele never matches; a no-call has no bases.
fn bases_match(genotype: &Genotype, index: usize, bases: &str) -> bool {
    let allele = &genotype.alleles[index];
    !allele.is_symbolic() && allele.display_string() == bases
}

/// `getIlluminaGenotype`.
fn illumina_genotype(
    genotype: &Genotype,
    context: &VariantContext,
) -> Result<IlluminaGenotype, String> {
    if !genotype.is_called() {
        return Ok(IlluminaGenotype::Nn);
    }
    // `StringUtils.stripEnd(value, "*")`: every trailing star, which marks the reference.
    let allele_a = info_attribute(context, "ALLELE_A")?;
    let allele_a = allele_a.trim_end_matches('*');
    let allele_b = info_attribute(context, "ALLELE_B")?;
    let allele_b = allele_b.trim_end_matches('*');
    let mismatch = || {
        format!(
            "picard.PicardException: Error matching called alleles to Illumina alleles.  \
             Context: {}:{}",
            context.contig, context.start
        )
    };
    if genotype.alleles.len() != 2 {
        return Err(format!(
            "picard.PicardException: Unexpected number of called alleles in variant context \
             {}:{}",
            context.contig, context.start
        ));
    }
    if bases_match(genotype, 0, allele_a) {
        if bases_match(genotype, 1, allele_a) {
            Ok(IlluminaGenotype::Aa)
        } else if bases_match(genotype, 1, allele_b) {
            Ok(IlluminaGenotype::Ab)
        } else {
            Err(mismatch())
        }
    } else if bases_match(genotype, 0, allele_b) {
        if bases_match(genotype, 1, allele_a) {
            Ok(IlluminaGenotype::Ab)
        } else if bases_match(genotype, 1, allele_b) {
            Ok(IlluminaGenotype::Bb)
        } else {
            Err(mismatch())
        }
    } else {
        Err(mismatch())
    }
}

/// One record, in the order `doWork` asks for its parts.
fn record(context: &VariantContext, sample: usize) -> Result<Record, String> {
    let gc_score = parse_float(&info_attribute(context, "GC_SCORE")?).map_err(number_format)?;
    let genotype = &context.genotypes[sample];
    let genotype_code = illumina_genotype(genotype, context)?;
    let a_intensity = unsigned_short(genotype, "X")?;
    let b_intensity = unsigned_short(genotype, "Y")?;
    let a_normalized = optional_float(genotype, "NORMX")?;
    let b_normalized = optional_float(genotype, "NORMY")?;
    Ok(Record {
        a_intensity: a_intensity as u16,
        b_intensity: b_intensity as u16,
        a_normalized,
        b_normalized,
        gc_score,
        genotype: genotype_code,
    })
}

/// The body of the `try`: everything it throws is the tool's own `log.error` and exit code 1.
fn convert(
    inputs: &[String],
    output: &str,
    samples_file: &str,
    markers_file: &str,
) -> Result<(), String> {
    let mut adpc: Vec<u8> = HEADER.to_vec();
    let mut sample_names: Vec<String> = Vec::new();
    let mut number_of_loci: Option<usize> = None;
    // The writer is opened before any input is read, so even a refusal leaves its header behind.
    std::fs::write(output, &adpc).map_err(|e| format!("java.io.FileNotFoundException: {e}"))?;
    for input in inputs {
        let vcf = read_path(input)?;
        for (sample_number, sample_name) in vcf.file.header.samples.iter().enumerate() {
            sample_names.push(sample_name.clone());
            let mut loci = 0usize;
            for record_in in &vcf.records {
                let record = record(&record_in.variant, sample_number)?;
                adpc.extend(write_record(&record));
                loci += 1;
            }
            if loci == 0 {
                return Err(format!(
                    "picard.PicardException: Found no records in VCF' {}'",
                    absolute(input)
                ));
            }
            match number_of_loci {
                None => number_of_loci = Some(loci),
                Some(expected) if expected != loci => {
                    return Err(
                        "picard.PicardException: VCFs have differing number of loci".to_string()
                    )
                }
                Some(_) => {}
            }
        }
    }
    let io = |e: std::io::Error| format!("java.io.IOException: {e}");
    std::fs::write(output, &adpc).map_err(io)?;
    std::fs::write(samples_file, sample_names.join("\n")).map_err(io)?;
    // `"" + numberOfLoci`: a run over files with no samples never sets it, and writes `null`.
    let markers = number_of_loci
        .map(|n| n.to_string())
        .unwrap_or_else(|| "null".to_string());
    std::fs::write(markers_file, markers).map_err(io)?;
    Ok(())
}

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

fn main() {
    let args = Args::from_env(&[
        ("O", "OUTPUT"),
        ("SF", "SAMPLES_FILE"),
        ("NMF", "NUM_MARKERS_FILE"),
    ]);
    let vcfs = args.all("VCF");
    let output = args.required("OUTPUT");
    let samples_file = args.required("SAMPLES_FILE");
    let markers_file = args.required("NUM_MARKERS_FILE");

    let inputs = unroll_paths(&vcfs).unwrap_or_else(|exception| die(&exception));
    for input in &inputs {
        if !std::path::Path::new(input).is_file() {
            die(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(input)
            ));
        }
    }

    if let Err(exception) = convert(&inputs, &output, &samples_file, &markers_file) {
        log_error(TOOL, "");
        eprintln!("{exception}");
        std::process::exit(1);
    }
}
