//! `SplitVcfs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.SplitVcfs.doWork` at tag 3.4.0 for plain-text outputs. The read and write
//! round trip is `picard_analysis::vcf_io`.
//!
//! Both outputs get the input's header unchanged, and each record goes to one of them by
//! `VariantContext.getType()`: `INDEL` to `INDEL_OUTPUT`, `SNP` to `SNP_OUTPUT`. Every other type
//! -- `NO_VARIATION`, `MNP`, `SYMBOLIC`, `MIXED` -- is dropped when `STRICT` is false and is an
//! uncaught `IllegalStateException` naming the type when it is true. The type is decided by the
//! alleles alone, so a multiallelic SNP is a SNP and a site whose alternates disagree (a SNP and
//! an insertion) is `MIXED`.
//!
//! `SEQUENCE_DICTIONARY` only replaces the dictionary the indexes are built against, which is
//! also why a VCF with no contig lines is refused without it only when indexing.

use htsjdk_vcf::variant::VariantContext;
use picard_analysis::vcf_io::{
    die, header_dictionary, parse_sam_dictionary, read_path, write_output,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `VariantContext.getType()`'s name.
fn type_name(vc: &VariantContext) -> &'static str {
    // `determineType`: one allele is no variation, whatever it is.
    if vc.alleles.len() <= 1 {
        return "NO_VARIATION";
    }
    let reference = &vc.alleles[0];
    let mut found: Option<&'static str> = None;
    for allele in &vc.alleles[1..] {
        // `typeOfBiallelicVariant`.
        let biallelic = if allele.is_symbolic() {
            "SYMBOLIC"
        } else if reference.len() == allele.len() {
            if allele.len() == 1 {
                "SNP"
            } else {
                "MNP"
            }
        } else {
            "INDEL"
        };
        match found {
            None => found = Some(biallelic),
            Some(previous) if previous != biallelic => return "MIXED",
            Some(_) => {}
        }
    }
    found.unwrap_or("NO_VARIATION")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let snp_output = arg(&args, "SNP_OUTPUT=").ok_or("SNP_OUTPUT= is required")?;
    let indel_output = arg(&args, "INDEL_OUTPUT=").ok_or("INDEL_OUTPUT= is required")?;
    let sequence_dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "D="))
        .filter(|value| value != "null");
    let strict = arg(&args, "STRICT=")
        .map(|value| value == "true")
        .unwrap_or(true);
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if !snp_output.ends_with(".vcf") || !indel_output.ends_with(".vcf") {
        return Err("only .vcf outputs are ported".into());
    }

    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    // Read through `SamReaderFactory`, as a SAM header: a `.dict`'s `@SQ` lines.
    let dictionary = match &sequence_dictionary {
        Some(path) => Some(parse_sam_dictionary(&std::fs::read_to_string(path)?)),
        None => header_dictionary(&vcf.file.header),
    };
    if create_index && dictionary.is_none() {
        die(
            "picard.PicardException: A sequence dictionary must be available (either through the \
             input file or by setting it explicitly) when creating indexed output.",
        );
    }

    let mut snps = Vec::new();
    let mut indels = Vec::new();
    for record in vcf.records {
        match type_name(&record.variant) {
            "INDEL" => indels.push(record),
            "SNP" => snps.push(record),
            other => {
                if strict {
                    die(&format!(
                        "java.lang.IllegalStateException: Found a record with type {other}"
                    ));
                }
            }
        }
    }

    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    write_output(&snp_output, &vcf.file.header, &snps, index_dictionary)?;
    write_output(&indel_output, &vcf.file.header, &indels, index_dictionary)?;
    Ok(())
}
