//! `ConvertHaplotypeDatabaseToVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.ConvertHaplotypeDatabaseToVcf.doWork` and `HaplotypeMap.writeAsVcf`
//! at tag 3.4.0. The table and the records live in `picard_analysis::haplotype_map`; this is the
//! file around them.
//!
//! * the reference is asked about EVERY SNP, in `Snp` order (contig name, then position), before
//!   any block is built, so the SNP a refusal names is the first of them that agrees with neither
//!   allele and not the first of the first block;
//! * the header is the one `writeAsVcf` builds: the contigs come from the haplotype map's own
//!   `@SQ` lines, the `##reference` line is the FASTA's `file:` URI (the later line that tries to
//!   overwrite it with the tool's name loses, because the header keeps the first of a key), and
//!   the one sample is `HetGenotypeForPhasing`;
//! * `AF` is written the way htsjdk formats a double between 0.01 and 1: three decimals.

use std::collections::HashMap;

use picard_analysis::haplotype_map::{
    as_vcf, parse_haplotype_database, ALLELE_DISAGREEMENT_PREFIX, HET_GENOTYPE_FOR_PHASING,
    VCF_SOURCE,
};
use picard_analysis::metrics_cli::{absolute, thrown, Args};

fn fasta_bases(text: &str) -> HashMap<String, Vec<u8>> {
    let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current: Option<String> = None;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            contigs.entry(name.clone()).or_default();
            current = Some(name);
        } else if let Some(name) = &current {
            if let Some(bases) = contigs.get_mut(name) {
                bases.extend(line.trim_end().bytes());
            }
        }
    }
    contigs
}

/// The `@SQ` lines of the table's own header: the name and the length of each contig, in order.
fn dictionary(text: &str) -> Vec<(String, String)> {
    text.lines()
        .filter(|line| line.starts_with("@SQ"))
        .filter_map(|line| {
            let mut name = None;
            let mut length = None;
            for field in line.split('\t') {
                if let Some(value) = field.strip_prefix("SN:") {
                    name = Some(value.to_string());
                } else if let Some(value) = field.strip_prefix("LN:") {
                    length = Some(value.to_string());
                }
            }
            Some((name?, length?))
        })
        .collect()
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("R", "REFERENCE_SEQUENCE")]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let reference = args.required("REFERENCE_SEQUENCE");
    if let Some(stringency) = args.get("VALIDATION_STRINGENCY") {
        if !matches!(stringency, "STRICT" | "LENIENT" | "SILENT") {
            picard_analysis::metrics_cli::fail(&format!(
                "Argument 'VALIDATION_STRINGENCY' cannot be set to '{stringency}'"
            ));
        }
    }

    let table = std::fs::read_to_string(&input).unwrap_or_else(|e| thrown(&format!("{e}")));
    let blocks = match parse_haplotype_database(&table) {
        Ok(blocks) => blocks,
        Err(message) => thrown(&format!("picard.PicardException: {message}")),
    };
    let fasta = std::fs::read_to_string(&reference).unwrap_or_else(|e| thrown(&format!("{e}")));
    let contigs = fasta_bases(&fasta);
    let base_at = |contig: &str, position: i32| -> Option<u8> {
        contigs
            .get(contig)?
            .get((position - 1) as usize)
            .map(u8::to_ascii_uppercase)
    };

    // `asVcf` checks every SNP in `Snp.compareTo` order first.
    let mut snps: Vec<_> = blocks
        .iter()
        .flat_map(|block| block.sorted_snps())
        .collect();
    snps.sort_by(|a, b| {
        a.chromosome
            .cmp(&b.chromosome)
            .then(a.position.cmp(&b.position))
    });
    for snp in &snps {
        let base = base_at(&snp.chromosome, snp.position).unwrap_or_else(|| {
            thrown(&format!(
                "java.lang.IllegalArgumentException: Unknown contig {}",
                snp.chromosome
            ))
        });
        if !base.eq_ignore_ascii_case(&snp.major_allele)
            && !base.eq_ignore_ascii_case(&snp.minor_allele)
        {
            thrown(&format!(
                "java.lang.RuntimeException: {ALLELE_DISAGREEMENT_PREFIX}{}:{}",
                snp.chromosome, snp.position
            ));
        }
    }

    let dict = dictionary(&table);
    let order: Vec<String> = dict.iter().map(|(name, _)| name.clone()).collect();
    let records = as_vcf(&blocks, &order, |contig, position| {
        base_at(contig, position)
    })
    .unwrap_or_else(|message| thrown(&format!("java.lang.RuntimeException: {message}")));

    let mut text = String::from("##fileformat=VCFv4.2\n");
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    text.push_str(
        "##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase-set identifier for phased genotypes.\">\n",
    );
    text.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n");
    for (name, length) in &dict {
        text.push_str(&format!("##contig=<ID={name},length={length}>\n"));
    }
    text.push_str(&format!("##reference=file://{}\n", absolute(&reference)));
    text.push_str(&format!("##source={VCF_SOURCE}\n"));
    text.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{HET_GENOTYPE_FOR_PHASING}\n"
    ));
    for record in &records {
        let (format, sample) = match record.phase_set {
            Some(set) => ("GT:PS", format!("{}:{set}", record.genotype)),
            None => ("GT", record.genotype.clone()),
        };
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t.\t.\tAF={:.3}\t{format}\t{sample}\n",
            record.chromosome,
            record.position,
            record.id,
            record.reference,
            record.alternate,
            record.allele_frequency,
        ));
    }
    if let Err(e) = std::fs::write(&output, text) {
        thrown(&format!(
            "picard.PicardException: Problem writing haplotype map to {output}: {e}"
        ));
    }
}
