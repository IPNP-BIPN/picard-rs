//! `ExtractFingerprint` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.ExtractFingerprint.doWork` at tag 3.4.0, through
//! `FingerprintChecker.identifyContaminant` and `FingerprintUtils.writeFingerPrint`:
//!
//! * `CONTAMINATION` is flipped to `1 - CONTAMINATION` unless `EXTRACT_CONTAMINATION`, and it is
//!   the FLIPPED value the `##source` line reports;
//! * one `HaplotypeProbabilitiesFromContaminatorSequence` per block and read group, fed by the
//!   pileup at every site (`LOCUS_MAX_READS` sampled through the checker's static
//!   `Random(42)`), then merged by sample; anything but exactly one sample is refused;
//! * one record per block's representative SNP (the block's first by position), or per SNP with
//!   `EXTRACT_NON_REPRESENTATIVES_TOO`, in a `TreeSet` on the dictionary's contig order and the
//!   position; REF is whichever allele the reference carries and the PLs and ADs swap with it;
//! * a failure while building the records is LOGGED and the run still exits 0, with the header
//!   already written and no records: a map whose alleles miss the reference leaves exactly that.
//!
//! The `##fileDate` line is the wall clock, which the harness strips.

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::genotype_likelihoods::gls_to_pls;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, VariantContext};
use picard_analysis::fingerprinting::{
    absolute_path, fingerprint_sam, log, merge_fingerprints_by, reference_view, uri_of, DataType,
    HaplotypeMap, Probabilities, SamOptions, Stringency,
};
use picard_analysis::metrics_cli::{fail, read_input, thrown, Args};
use picard_analysis::theoretical_sensitivity::JavaRandom;

const TOOL: &str = "ExtractFingerprint";

/// `Double.toString`.
fn java_double(value: f64) -> String {
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0.0".into()
        } else {
            "0.0".into()
        };
    }
    let abs = value.abs();
    if (1e-3..1e7).contains(&abs) {
        let text = format!("{value}");
        if text.contains('.') {
            text
        } else {
            format!("{text}.0")
        }
    } else {
        let text = format!("{value:e}");
        let (mantissa, exponent) = text.split_once('e').expect("exponent form");
        let mantissa = if mantissa.contains('.') {
            mantissa.to_string()
        } else {
            format!("{mantissa}.0")
        };
        format!("{mantissa}E{exponent}")
    }
}

/// `new Date().toString()` in UTC: `EEE MMM dd HH:mm:ss zzz yyyy`.
fn java_date_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let weekday = ["Thu", "Fri", "Sat", "Sun", "Mon", "Tue", "Wed"][days.rem_euclid(7) as usize];
    // Civil from days (Howard Hinnant).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    let month = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ][(m - 1) as usize];
    format!(
        "{weekday} {month} {d:02} {:02}:{:02}:{:02} UTC {y}",
        rem / 3600,
        (rem / 60) % 60,
        rem % 60
    )
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("C", "CONTAMINATION"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let haplotype_map = args.required("HAPLOTYPE_MAP");
    let mut contamination = args.double("CONTAMINATION", 0.0);
    let sample_alias = args.get("SAMPLE_ALIAS").map(str::to_string);
    let locus_max_reads = args.int("LOCUS_MAX_READS", 50) as i32;
    let extract_contamination = args.bool("EXTRACT_CONTAMINATION", false);
    let _test_readability = args.bool("TEST_INPUT_READABILITY", true);
    let non_representatives = args.bool("EXTRACT_NON_REPRESENTATIVES_TOO", false);
    let reference = args.required("REFERENCE_SEQUENCE");
    let stringency = match args.get("VALIDATION_STRINGENCY").unwrap_or("STRICT") {
        "STRICT" => Stringency::Strict,
        "LENIENT" => Stringency::Lenient,
        "SILENT" => Stringency::Silent,
        other => fail(&format!(
            "Argument 'VALIDATION_STRINGENCY' cannot be set to '{other}'"
        )),
    };
    if !(0.0..=1.0).contains(&contamination) {
        fail(&format!(
            "Argument 'CONTAMINATION' has a value {contamination} outside its range"
        ));
    }

    let map_text = std::fs::read_to_string(&haplotype_map).unwrap_or_else(|e| fail(&e.to_string()));
    let map = HaplotypeMap::from_database(&map_text, &absolute_path(&haplotype_map))
        .unwrap_or_else(|(class, message)| thrown(&format!("{class}: {message}")));

    if !extract_contamination {
        contamination = 1.0 - contamination;
    }
    let options = SamOptions {
        locus_max_reads,
        stringency,
        default_sample: sample_alias
            .clone()
            .unwrap_or_else(|| "<UNKNOWN>".to_string()),
        ..SamOptions::default()
    };

    let (header, records) = read_input(&input);
    let uri = uri_of(&input);
    let mut random = JavaRandom::new(42);
    let c = contamination;
    let by_group = fingerprint_sam(
        &header,
        &records,
        &input,
        &uri,
        &map,
        &options,
        &mut random,
        &|b| Probabilities::contaminator(b, c),
    )
    .unwrap_or_else(|(class, message)| thrown(&format!("{class}: {message}")));
    let by_sample = merge_fingerprints_by(&by_group, DataType::Sample, &map)
        .unwrap_or_else(|(class, message)| thrown(&format!("{class}: {message}")));
    if by_sample.len() != 1 {
        log(
            "ERROR",
            TOOL,
            &format!("Expected exactly 1 fingerprint, found {}", by_sample.len()),
        );
        thrown(&format!(
            "java.lang.IllegalArgumentException: Expected exactly 1 fingerprint in Input file, found {}",
            by_sample.len()
        ));
    }
    let (id, fingerprint) = by_sample.iter().next().expect("one entry");
    let fp_sample = id.sample.clone().unwrap_or_else(|| "null".into());
    let sample = match &sample_alias {
        Some(alias) => alias.clone(),
        None if extract_contamination => format!("{fp_sample}-contaminant"),
        None => fp_sample,
    };

    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| fail(&format!("{e:?}")));
    let source = format!(
        "PLs derived from {} using an assumed contamination of {}",
        reference_view(&input),
        java_double(contamination)
    );
    let mut vcf_header = VcfHeader::new();
    vcf_header.lines.push(HeaderLine::Unstructured {
        key: "reference".into(),
        value: reference_view(&absolute_path(&reference)),
    });
    vcf_header.lines.push(HeaderLine::Unstructured {
        key: "source".into(),
        value: source,
    });
    vcf_header.lines.push(HeaderLine::Unstructured {
        key: "fileDate".into(),
        value: java_date_now(),
    });
    vcf_header.lines.push(HeaderLine::format(
        "PL",
        Cardinality::G,
        LineType::Integer,
        "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification",
    ));
    vcf_header.lines.push(HeaderLine::format(
        "AD",
        Cardinality::R,
        LineType::Integer,
        "Allelic depths for the ref and alt alleles in the order listed",
    ));
    vcf_header.lines.push(HeaderLine::format(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
    ));
    for (index, contig) in contigs.iter().enumerate() {
        vcf_header.lines.push(HeaderLine::contig(
            &contig.name,
            contig.bases.len() as i64,
            index as i32,
        ));
    }
    vcf_header.samples = vec![sample.clone()];
    let mut text = vcf_header.write();

    // `createVCSetFromFingerprint`: every record is built before any is written, so a failure
    // leaves the header alone.
    match records_of(&map, fingerprint, &contigs, &sample, !non_representatives) {
        Ok(vcs) => {
            let encoder = VcfEncoder::new(&vcf_header);
            for vc in &vcs {
                encoder
                    .encode_into(vc, &mut text)
                    .unwrap_or_else(|e| fail(&format!("{e:?}")));
                text.push('\n');
            }
        }
        Err(message) => log("ERROR", TOOL, &message),
    }
    std::fs::write(&output, text).unwrap_or_else(|e| fail(&e.to_string()));
}

/// `createVCSetFromFingerprint` and `getVariantContextFromSnp`.
fn records_of(
    map: &HaplotypeMap,
    fingerprint: &picard_analysis::fingerprinting::Fingerprint,
    contigs: &[htsjdk_bam::fasta::ReferenceSequence],
    sample: &str,
    representative_only: bool,
) -> Result<Vec<VariantContext>, String> {
    // The same SNP name twice is refused before anything is built.
    let mut names: Vec<&str> = Vec::new();
    for hp in fingerprint.map.values() {
        let snps: Vec<usize> = if representative_only {
            vec![hp.representative_snp(map)]
        } else {
            map.blocks[hp.block].snps.clone()
        };
        for s in snps {
            let name = map.snps[s].name.as_str();
            if name.is_empty() {
                continue;
            }
            if names.contains(&name) {
                return Err(format!(
                    "Found same SNP name twice ({name}) in fingerprint. Cannot create a VCF."
                ));
            }
            names.push(name);
        }
    }
    let mut set: Vec<(i32, i32, VariantContext)> = Vec::new();
    for hp in fingerprint.map.values() {
        let snps: Vec<usize> = if representative_only {
            vec![hp.representative_snp(map)]
        } else {
            map.blocks[hp.block].snps.clone()
        };
        for s in snps {
            let snp = &map.snps[s];
            let contig_index = contigs.iter().position(|c| c.name == snp.chrom);
            let Some(contig_index) = contig_index else {
                return Err(format!("Unknown contig {}", snp.chrom));
            };
            let base = contigs[contig_index]
                .bases
                .get((snp.pos - 1) as usize)
                .copied()
                .unwrap_or(b'N')
                .to_ascii_uppercase();
            if snp.allele1 != base && snp.allele2 != base {
                return Err(
                    "Don't know how to deal with missing reference allele in fingerprinting map"
                        .into(),
                );
            }
            let swap = snp.allele2 == base;
            let (reference, alternate, obs_ref, obs_alt) = if swap {
                (snp.allele2, snp.allele1, hp.obs_allele2(), hp.obs_allele1())
            } else {
                (snp.allele1, snp.allele2, hp.obs_allele1(), hp.obs_allele2())
            };
            let mut gls = hp.log_likelihoods(map);
            if swap {
                gls.reverse();
            }
            let alleles = vec![
                Allele::create(&[reference], true).map_err(|e| format!("{e:?}"))?,
                Allele::create(&[alternate], false).map_err(|e| format!("{e:?}"))?,
            ];
            let mut genotype = Genotype::new(sample, Vec::new());
            genotype.dp = Some(hp.total_obs());
            genotype.pl = Some(gls_to_pls(&gls));
            genotype.ad = Some(vec![obs_ref, obs_alt]);
            let mut vc = VariantContext::new(&snp.chrom, i64::from(snp.pos), alleles);
            vc.genotypes = vec![genotype];
            let key = (contig_index as i32, snp.pos);
            if !set.iter().any(|(c, p, _)| (*c, *p) == key) {
                set.push((key.0, key.1, vc));
            }
        }
    }
    set.sort_by_key(|(c, p, _)| (*c, *p));
    Ok(set.into_iter().map(|(_, _, vc)| vc).collect())
}
