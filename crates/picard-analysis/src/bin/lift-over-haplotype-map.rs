//! `LiftOverHaplotypeMap` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.LiftOverHaplotypeMap.doWork` at tag 3.4.0. The loop between the chain
//! and the table lives in `picard_analysis::lift_over_haplotype_map`; the table itself is
//! `picard_analysis::haplotype_map`.
//!
//! A SNP that does not lift is DROPPED and the run carries on, so the answer is a shorter table
//! rather than none at all -- and a database whose every SNP failed still leaves a file holding
//! its header and its column line. The exit code for any failure is 101, not the usual 1.
//!
//! The alleles are carried over unchanged: the chain may put a SNP on the negative strand and the
//! tool still writes the bases it read, never their complements.

use picard_analysis::haplotype_map::{format_frequency, parse_haplotype_database};
use picard_analysis::lift_over_haplotype_map::lift_over_haplotype_map;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`: a FASTA is not parsed, the `.dict` beside
/// it is.
fn read_dictionary(path: &str) -> std::io::Result<String> {
    let candidate = std::path::Path::new(path);
    if let Some(extension) = candidate.extension().and_then(|e| e.to_str()) {
        if matches!(extension, "fasta" | "fa" | "fna") {
            return std::fs::read_to_string(candidate.with_extension("dict"));
        }
    }
    std::fs::read_to_string(path)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let chain = arg(&args, "CHAIN=").ok_or("CHAIN= is required")?;
    let dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "SD="))
        .ok_or("SEQUENCE_DICTIONARY= is required")?;

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let blocks = match parse_haplotype_database(&std::fs::read_to_string(&input)?) {
        Ok(blocks) => blocks,
        Err(message) => {
            eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
            std::process::exit(1);
        }
    };
    let lift = match htsjdk_bam::liftover::LiftOver::load(&std::fs::read_to_string(&chain)?) {
        Ok(lift) => lift,
        Err(error) => {
            eprintln!("Exception in thread \"main\" picard.PicardException: {error:?}");
            std::process::exit(1);
        }
    };

    // `new SAMFileHeader(dict)` written through `SAMTextHeaderCodec`: the version line the codec
    // writes is the CURRENT one, not the dictionary file's, and then the `@SQ` lines as they are.
    let dictionary_text = read_dictionary(&dictionary)?;
    let sequence_lines: Vec<&str> = dictionary_text
        .lines()
        .filter(|line| line.starts_with("@SQ"))
        .collect();
    let names: Vec<String> = sequence_lines
        .iter()
        .filter_map(|line| {
            line.split('\t')
                .find_map(|field| field.strip_prefix("SN:"))
                .map(str::to_string)
        })
        .collect();

    let result = match lift_over_haplotype_map(&blocks, &lift, &names) {
        Ok(result) => result,
        Err(missing) => {
            eprintln!("Exception in thread \"main\" picard.PicardException: {missing}");
            std::process::exit(1);
        }
    };

    let mut text = String::from("@HD\tVN:1.6\n");
    for line in &sequence_lines {
        text.push_str(line);
        text.push('\n');
    }
    text.push_str(
        "#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS\n",
    );
    for row in &result.rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.chromosome,
            row.position,
            row.name,
            row.major_allele as char,
            row.minor_allele as char,
            format_frequency(row.minor_allele_frequency),
            row.anchor.clone().unwrap_or_default(),
            row.panels.clone().unwrap_or_default(),
        ));
    }
    std::fs::write(&output, text)?;
    std::process::exit(result.return_code);
}
