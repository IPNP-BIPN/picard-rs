//! `FixVcfHeader` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.FixVcfHeader.doWork` at tag 3.4.0. The read and write round trip is
//! `picard_analysis::vcf_io`.
//!
//! # Without `HEADER`, the header is rebuilt as a set, and a set keeps near-duplicates
//!
//! The records are read once to collect every FILTER, INFO key and extended FORMAT key the header
//! does not define, each given a placeholder line (`Number=.`, `Type=String`). The new header is
//! then `new HashSet<>(existing lines)`, plus the standard FORMAT lines of `Genotype.PRIMARY_KEYS`
//! whether the file uses them or not, plus the placeholders. A compound line's `equals` includes
//! its description, so a file whose `DP` says `Depth` keeps that line AND gains the standard one:
//! two `FORMAT=<ID=DP` lines, both written. The samples are the sorted names.
//!
//! # With `HEADER`, the header is replaced and nothing is checked until a record is written
//!
//! `ENFORCE_SAME_SAMPLES` compares the two files' SORTED names, so the same samples in another
//! column order pass -- and the replacement header's order is what is written above the input's
//! copied genotype blocks. Without it, the replacement's lines go over the input's sorted names.
//! The writer refuses a missing field (`ALLOW_MISSING_FIELDS_IN_HEADER` is unset) only for what it
//! encodes: a copied genotype block's FORMAT keys are never looked at.
//!
//! The writer always indexes, whatever `CREATE_INDEX` says, so the output header needs contig
//! lines.

use htsjdk_vcf::encoder::EncodeError;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::standard_header_lines::{standard_format_line, PRIMARY_KEYS};
use picard_analysis::vcf_io::{
    die, header_dictionary, read_path, sample_names_in_order, write_output, write_records,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

const MISSING: &str = "Missing description: this {} line was added by Picard's FixVCFHeader";

/// `VCFCompoundHeaderLine.equals`: the key, ID, count, type and description; the extra fields are
/// not compared.
fn same_line(a: &HeaderLine, b: &HeaderLine) -> bool {
    match (a, b) {
        (
            HeaderLine::Compound {
                key: ak,
                id: ai,
                number: an,
                line_type: at,
                description: ad,
                ..
            },
            HeaderLine::Compound {
                key: bk,
                id: bi,
                number: bn,
                line_type: bt,
                description: bd,
                ..
            },
        ) => ak == bk && ai == bi && an == bn && at == bt && ad == bd,
        _ => a == b,
    }
}

fn encode_error(error: &EncodeError) -> String {
    match error {
        EncodeError::MissingFromHeader {
            key,
            field,
            contig,
            start,
        } => format!(
            "java.lang.IllegalStateException: Key {key} found in VariantContext field {field} at \
             {contig}:{start} but this key isn't defined in the VCFHeader.  We require all VCFs \
             to have complete VCF headers by default."
        ),
        other => format!("{other:?}"),
    }
}

/// `enforceSameSamples`: the sorted names, pairwise.
fn assert_same_samples(reader: &VcfHeader, input: &VcfHeader) {
    let reader_samples = sample_names_in_order(reader);
    let input_samples = sample_names_in_order(input);
    if reader_samples.len() != input_samples.len() {
        die(
            "picard.PicardException: The input VCF had a different # of samples than the input \
             VCF header.",
        );
    }
    for (i, (a, b)) in reader_samples.iter().zip(&input_samples).enumerate() {
        if a != b {
            die(&format!(
                "picard.PicardException: Mismatch in the {i}th sample: '{a}' != '{b}'"
            ));
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let check_first_n: i64 = arg(&args, "CHECK_FIRST_N_RECORDS=")
        .or_else(|| arg(&args, "N="))
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(-1);
    let header_path = arg(&args, "HEADER=")
        .or_else(|| arg(&args, "H="))
        .filter(|value| value != "null");
    let enforce_same_samples = arg(&args, "ENFORCE_SAME_SAMPLES=")
        .map(|value| value == "true")
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if !output.ends_with(".vcf") {
        return Err("only a .vcf OUTPUT is ported".into());
    }

    // `customCommandLineValidation`: printed after the usage, and the program exits 1.
    if header_path.is_some() && check_first_n >= 0 {
        eprintln!("CHECK_FIRST_N_RECORDS should no be specified when HEADER is specified");
        std::process::exit(1);
    }

    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    let existing = &vcf.file.header;

    let out_header = match &header_path {
        Some(path) => {
            let replacement = read_path(path).unwrap_or_else(|exception| die(&exception));
            if enforce_same_samples {
                let header = replacement.file.header;
                assert_same_samples(existing, &header);
                header
            } else {
                // `new VCFHeader(inputHeader.getMetaDataInInputOrder(), getSampleNamesInOrder())`.
                VcfHeader {
                    lines: replacement.file.header.lines,
                    samples: sample_names_in_order(existing),
                }
            }
        }
        None => {
            let mut filters: Vec<HeaderLine> = Vec::new();
            let mut infos: Vec<HeaderLine> = Vec::new();
            let mut formats: Vec<HeaderLine> = Vec::new();
            let scanned = if check_first_n > 0 {
                vcf.records.len().min(check_first_n as usize)
            } else {
                vcf.records.len()
            };
            let defines = |lines: &[HeaderLine], id: &str| {
                lines.iter().any(|line| match line {
                    HeaderLine::Compound { id: i, .. } | HeaderLine::Filter { id: i, .. } => {
                        i == id
                    }
                    _ => false,
                })
            };
            for record in &vcf.records[..scanned] {
                let ctx = &record.variant;
                for filter in ctx.filters.iter().flatten() {
                    if !existing.has_filter_line(filter) && !defines(&filters, filter) {
                        filters.push(HeaderLine::filter(filter, &MISSING.replace("{}", "FILTER")));
                    }
                }
                for (id, _) in &ctx.attributes {
                    if !existing.has_info_line(id) && !defines(&infos, id) {
                        infos.push(HeaderLine::info(
                            id,
                            Cardinality::Unbounded,
                            LineType::String,
                            &MISSING.replace("{}", "INFO"),
                        ));
                    }
                }
                for genotype in &ctx.genotypes {
                    for (id, _) in &genotype.extended {
                        if !existing.has_format_line(id) && !defines(&formats, id) {
                            formats.push(HeaderLine::format(
                                id,
                                Cardinality::Unbounded,
                                LineType::String,
                                &MISSING.replace("{}", "FORMAT"),
                            ));
                        }
                    }
                }
            }

            // `new HashSet<>(getMetaDataInInputOrder())`, then the standard FORMAT lines and the
            // placeholders: a line equal to one already there is not added twice.
            let mut lines: Vec<HeaderLine> = Vec::new();
            let mut add = |line: HeaderLine| {
                if !lines.iter().any(|kept| same_line(kept, &line)) {
                    lines.push(line);
                }
            };
            existing.lines.iter().cloned().for_each(&mut add);
            PRIMARY_KEYS
                .iter()
                .filter_map(|id| standard_format_line(id))
                .for_each(&mut add);
            filters.into_iter().for_each(&mut add);
            infos.into_iter().for_each(&mut add);
            formats.into_iter().for_each(&mut add);
            VcfHeader {
                lines,
                samples: sample_names_in_order(existing),
            }
        }
    };

    // `setOption(INDEX_ON_THE_FLY)` with `setReferenceDictionary(outHeader.getSequenceDictionary())`.
    let Some(dictionary) = header_dictionary(&out_header) else {
        die(
            "java.lang.IllegalArgumentException: A reference dictionary is required for creating \
             Tribble indices on the fly",
        );
    };
    if let Err(error) = write_records(&out_header, &vcf.records) {
        die(&encode_error(&error));
    }
    write_output(&output, &out_header, &vcf.records, Some(&dictionary))?;
    Ok(())
}
