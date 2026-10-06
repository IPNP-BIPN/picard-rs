//! `ExtractFingerprint.doWork`, which `IdentifyContaminant.doWork` also runs: it builds an
//! `ExtractFingerprint` with its own defaults and the contamination flag inverted.
//!
//! Ports `picard.fingerprint.ExtractFingerprint.doWork`, `FingerprintChecker.identifyContaminant`
//! and `FingerprintUtils.writeFingerPrint` at tag 3.4.0, over `crate::fingerprint`.
//!
//! # A failure to write is logged, and the run succeeds
//!
//! `doWork` catches whatever `writeFingerPrint` throws and returns 0. The writer is built, and its
//! header written, before the records are, so a map naming a SNP twice, or a SNP neither of whose
//! alleles is the reference base, leaves a VCF with a header and no records, and exit code 0.

use crate::crosscheck_run::{merge_by, DataType};
use crate::fingerprint::{
    fingerprint_sam_file, reference_path, Evidence, HaplotypeMap, Probs, SamOptions, SharedRandom,
    Snp,
};
use crate::metrics_cli::{thrown, Args};
use crate::theoretical_sensitivity::JavaRandom;
use crate::vcf_io::{log_error, parse_sam_dictionary, set_sequence_dictionary, Record};
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, VariantContext};

/// `Double.toString`.
fn java_double_to_string(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let magnitude = value.abs();
    if magnitude == 0.0 {
        return format!("{sign}0.0");
    }
    let scientific = format!("{magnitude:e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if (1e-3..1e7).contains(&magnitude) {
        if exponent >= 0 {
            let whole_len = exponent as usize + 1;
            let mut padded = digits.clone();
            while padded.len() < whole_len {
                padded.push('0');
            }
            let (whole, fraction) = padded.split_at(whole_len);
            let fraction = if fraction.is_empty() { "0" } else { fraction };
            format!("{sign}{whole}.{fraction}")
        } else {
            let zeros = "0".repeat((-exponent - 1) as usize);
            format!("{sign}0.{zeros}{digits}")
        }
    } else {
        let (first, rest) = digits.split_at(1);
        let rest = if rest.is_empty() { "0" } else { rest };
        format!("{sign}{first}.{rest}E{exponent}")
    }
}

/// `GenotypeLikelihoods.GLsToPLs`.
fn pls(gls: [f64; 3]) -> Vec<i32> {
    let adjust = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter()
        .map(|gl| ((-10.0 * (gl - adjust)).min(f64::from(i32::MAX)) + 0.5).floor() as i32)
        .collect()
}

/// `getVariantContextFromSnp`.
fn variant(
    reference: &[(String, Vec<u8>)],
    sample: &str,
    probs: &Probs,
    snp: &Snp,
) -> Result<VariantContext, String> {
    let bases = reference
        .iter()
        .find(|(n, _)| *n == snp.chrom)
        .map(|(_, b)| b)
        .ok_or_else(|| format!("Unable to find entry for contig: {}", snp.chrom))?;
    let ref_base = bases[(snp.pos - 1) as usize].to_ascii_uppercase();
    if snp.allele1 != ref_base && snp.allele2 != ref_base {
        return Err(
            "picard.PicardException: Don't know how to deal with missing reference allele in fingerprinting map"
                .to_string(),
        );
    }
    let swap = snp.allele2 == ref_base;
    let (r, a, obs_ref, obs_alt) = if swap {
        (snp.allele2, snp.allele1, probs.obs2(), probs.obs1())
    } else {
        (snp.allele1, snp.allele2, probs.obs1(), probs.obs2())
    };
    let mut gls = probs.log_likelihoods();
    if swap {
        gls.reverse();
    }
    let alleles = vec![
        htsjdk_vcf::allele::Allele::create(&[r], true).map_err(|e| e.to_string())?,
        htsjdk_vcf::allele::Allele::create(&[a], false).map_err(|e| e.to_string())?,
    ];
    let mut genotype = Genotype::new(sample, Vec::new());
    genotype.dp = Some(probs.total_obs());
    genotype.pl = Some(pls(gls));
    genotype.ad = Some(vec![obs_ref, obs_alt]);
    let mut vc = VariantContext::new(&snp.chrom, i64::from(snp.pos), alleles);
    vc.genotypes = vec![genotype];
    Ok(vc)
}

/// `createVCSetFromFingerprint`: one record per SNP (or per block's representative), in
/// dictionary order, the first at a position kept.
fn records(
    map: &HaplotypeMap,
    fp: &crate::fingerprint::Fingerprint,
    reference: &[(String, Vec<u8>)],
    dictionary: &[String],
    sample: &str,
    representative_only: bool,
) -> Result<Vec<VariantContext>, String> {
    let snps_of = |p: &Probs| -> Vec<Snp> {
        if representative_only {
            vec![p.representative().clone()]
        } else {
            map.blocks[p.block].snps.clone()
        }
    };
    let mut names: Vec<String> = Vec::new();
    for p in fp.blocks.values() {
        for snp in snps_of(p) {
            if snp.name.is_empty() {
                continue;
            }
            if names.contains(&snp.name) {
                return Err(format!("java.lang.IllegalArgumentException: Found same SNP name twice ({}) in fingerprint. Cannot create a VCF.", snp.name));
            }
            names.push(snp.name.clone());
        }
    }
    let mut out: Vec<VariantContext> = Vec::new();
    for p in fp.blocks.values() {
        for snp in snps_of(p) {
            let vc = variant(reference, sample, p, &snp)?;
            let index = |c: &str| {
                dictionary
                    .iter()
                    .position(|d| d == c)
                    .map_or(-1, |i| i as i64)
            };
            let at = out.binary_search_by(|o| {
                index(&o.contig)
                    .cmp(&index(&vc.contig))
                    .then(o.start.cmp(&vc.start))
            });
            if let Err(at) = at {
                out.insert(at, vc);
            }
        }
    }
    Ok(out)
}

/// `ExtractFingerprint.doWork`, with the defaults and the inversion `IdentifyContaminant` applies
/// already folded into the arguments.
pub fn run(args: &Args, extract_contamination: bool, locus_max_reads: i64) {
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let reference = args.required("REFERENCE_SEQUENCE");
    let map = HaplotypeMap::load(&args.required("HAPLOTYPE_MAP")).unwrap_or_else(|e| thrown(&e));
    let sample_alias = args.get("SAMPLE_ALIAS").map(str::to_string);
    let mut contamination = args.double("CONTAMINATION", 0.0);
    if !extract_contamination {
        contamination = 1.0 - contamination;
    }
    let strict = !matches!(
        args.get("VALIDATION_STRINGENCY"),
        Some("LENIENT") | Some("SILENT")
    );
    let options = SamOptions {
        locus_max_reads: locus_max_reads.max(0) as usize,
        strict,
        default_sample: sample_alias
            .clone()
            .unwrap_or_else(|| "<UNKNOWN>".to_string()),
        ..SamOptions::default()
    };
    let mut random = SharedRandom(JavaRandom::new(42));
    let c = contamination;
    let by_group = fingerprint_sam_file(
        &input,
        &map,
        &options,
        &mut random,
        &move |m, b| {
            Probs::with(
                m,
                b,
                Evidence::Contaminator {
                    map: [[0.0; 3]; 3],
                    contamination: c,
                    obs1: 0,
                    obs2: 0,
                    other: 0,
                },
            )
        },
        &|p, snp, base, qual| p.add_contaminator_base(snp, base, qual),
    )
    .unwrap_or_else(|e| thrown(&e));
    let by_sample = merge_by(&by_group, DataType::Sample);
    if by_sample.len() != 1 {
        log_error(
            "ExtractFingerprint",
            &format!("Expected exactly 1 fingerprint, found {}", by_sample.len()),
        );
        thrown(&format!(
            "java.lang.IllegalArgumentException: Expected exactly 1 fingerprint in Input file, found {}",
            by_sample.len()
        ));
    }
    let (details, fp) = &by_sample[0];
    let sample = match &sample_alias {
        Some(s) => s.clone(),
        None => format!(
            "{}{}",
            details.sample.as_deref().unwrap_or("null"),
            if extract_contamination {
                "-contaminant"
            } else {
                ""
            }
        ),
    };

    // `getVariantContextWriter`: the header, written before any record is built.
    let fasta = std::fs::read(&reference).unwrap_or_else(|e| thrown(&format!("{e}")));
    let sequences: Vec<(String, Vec<u8>)> = htsjdk_bam::fasta::read_fasta(&fasta[..])
        .unwrap_or_else(|e| thrown(&format!("{e:?}")))
        .into_iter()
        .map(|s| (s.name, s.bases))
        .collect();
    let dict_path = match reference.rsplit_once('.') {
        Some((stem, _)) => format!("{stem}.dict"),
        None => format!("{reference}.dict"),
    };
    let dictionary = std::fs::read_to_string(&dict_path)
        .map(|t| parse_sam_dictionary(&t))
        .unwrap_or_default();
    let mut header = VcfHeader::new();
    header.lines.push(HeaderLine::Unstructured {
        key: "reference".to_string(),
        value: reference_path(&reference),
    });
    header.lines.push(HeaderLine::Unstructured {
        key: "source".to_string(),
        value: format!(
            "PLs derived from {input} using an assumed contamination of {}",
            java_double_to_string(contamination)
        ),
    });
    header.lines.push(HeaderLine::Unstructured {
        key: "fileDate".to_string(),
        value: "<now>".to_string(),
    });
    header.lines.push(HeaderLine::format(
        "PL",
        Cardinality::G,
        LineType::Integer,
        "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification",
    ));
    header.lines.push(HeaderLine::format(
        "AD",
        Cardinality::R,
        LineType::Integer,
        "Allelic depths for the ref and alt alleles in the order listed",
    ));
    header.lines.push(HeaderLine::format(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
    ));
    header.samples = vec![sample.clone()];
    set_sequence_dictionary(&mut header, &dictionary);
    let names: Vec<String> = dictionary.iter().map(|d| d.name.clone()).collect();
    let representative_only = !args.bool("EXTRACT_NON_REPRESENTATIVES_TOO", false);
    let body = match records(&map, fp, &sequences, &names, &sample, representative_only) {
        Ok(vcs) => vcs,
        Err(e) => {
            log_error("ExtractFingerprint", &e);
            Vec::new()
        }
    };
    let rows: Vec<Record> = body
        .into_iter()
        .map(|variant| Record {
            variant,
            lazy_genotypes: None,
        })
        .collect();
    let text =
        crate::vcf_io::write_records(&header, &rows).unwrap_or_else(|e| thrown(&format!("{e:?}")));
    if let Err(e) = std::fs::write(&output, text) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
