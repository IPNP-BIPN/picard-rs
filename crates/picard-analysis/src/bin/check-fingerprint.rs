//! `CheckFingerprint` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CheckFingerprint.doWork` at tag 3.4.0, through
//! `FingerprintChecker.checkFingerprints` (reads) and `checkFingerprintsFromPaths` (a VCF):
//!
//! * `customCommandLineValidation`: `IGNORE_READ_GROUPS` only with reads and
//!   `OBSERVED_SAMPLE_ALIAS` only with a VCF, told apart by the file's extension;
//! * the observed sample from the BAM's read groups (one sample only) or the VCF's header, the
//!   expected one defaulting to it, and the expected sample's absence answered by
//!   `EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND` with nothing written;
//! * the expected fingerprints are the named sample's AND an empty one for every other sample in
//!   the genotypes file, and each result row keeps only the first of its `TreeSet` of matches,
//!   the highest LOD: an observed sample that matches nobody well reports an empty fingerprint's
//!   zeros under the expected sample's name, with no detail rows;
//! * one row per read group in the `HashMap` order of their identities (or one for the file with
//!   `IGNORE_READ_GROUPS`, or one per sample of an observed VCF), the detail rows in SNP order;
//! * `EXIT_CODE_WHEN_NO_VALID_CHECKS` when every LOD is zero, after both files are written.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprinting::{
    absolute_path, calculate_match_results, compare_match_results, fingerprint_sam,
    load_fingerprints, log, read_genotype_file, uri_of, Fingerprint, HaplotypeMap, MatchResults,
    SamOptions,
};
use picard_analysis::metrics_cli::{fail, read_input, refuse_validation, thrown, Args};
use picard_analysis::theoretical_sensitivity::JavaRandom;

const TOOL: &str = "CheckFingerprint";

struct Summary {
    read_group: Option<String>,
    sample: String,
    ll_expected: f64,
    ll_random: f64,
    lod: f64,
    with_genotypes: i64,
    checked: i64,
    matching: i64,
    het_as_hom: i64,
    hom_as_het: i64,
    hom_as_other_hom: i64,
}

impl MetricBean for Summary {
    fn class_name(&self) -> &str {
        "picard.analysis.FingerprintingSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "READ_GROUP",
            "SAMPLE",
            "LL_EXPECTED_SAMPLE",
            "LL_RANDOM_SAMPLE",
            "LOD_EXPECTED_SAMPLE",
            "HAPLOTYPES_WITH_GENOTYPES",
            "HAPLOTYPES_CONFIDENTLY_CHECKED",
            "HAPLOTYPES_CONFIDENTLY_MATCHING",
            "HET_AS_HOM",
            "HOM_AS_HET",
            "HOM_AS_OTHER_HOM",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            opt(&self.read_group),
            Value::Str(self.sample.clone()),
            Value::Double(self.ll_expected),
            Value::Double(self.ll_random),
            Value::Double(self.lod),
            Value::Long(self.with_genotypes),
            Value::Long(self.checked),
            Value::Long(self.matching),
            Value::Long(self.het_as_hom),
            Value::Long(self.hom_as_het),
            Value::Long(self.hom_as_other_hom),
        ]
    }
}

struct Detail {
    read_group: Option<String>,
    sample: String,
    snp: String,
    alleles: String,
    chrom: String,
    position: i64,
    expected: String,
    observed: String,
    lod: f64,
    obs_a: i64,
    obs_b: i64,
}

impl MetricBean for Detail {
    fn class_name(&self) -> &str {
        "picard.analysis.FingerprintingDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "READ_GROUP",
            "SAMPLE",
            "SNP",
            "SNP_ALLELES",
            "CHROM",
            "POSITION",
            "EXPECTED_GENOTYPE",
            "OBSERVED_GENOTYPE",
            "LOD",
            "OBS_A",
            "OBS_B",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            opt(&self.read_group),
            Value::Str(self.sample.clone()),
            Value::Str(self.snp.clone()),
            Value::Str(self.alleles.clone()),
            Value::Str(self.chrom.clone()),
            Value::Long(self.position),
            Value::Str(self.expected.clone()),
            Value::Str(self.observed.clone()),
            Value::Double(self.lod),
            Value::Long(self.obs_a),
            Value::Long(self.obs_b),
        ]
    }
}

fn opt(value: &Option<String>) -> Value {
    match value {
        Some(s) => Value::Str(s.clone()),
        None => Value::Null,
    }
}

/// `CheckFingerprint.fileContainsReads`: by the URI's path, so by extension.
fn file_contains_reads(path: &str) -> bool {
    path.ends_with(".bam") || path.ends_with(".sam") || path.ends_with(".cram")
}

fn throw(thrown_: (String, String)) -> ! {
    thrown(&format!("{}: {}", thrown_.0, thrown_.1))
}

fn read_genotypes(path: &str) -> picard_analysis::fingerprinting::GenotypeFile {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&e.to_string()));
    read_genotype_file(&text).unwrap_or_else(|t| throw(t))
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("S", "SUMMARY_OUTPUT"),
        ("D", "DETAIL_OUTPUT"),
        ("G", "GENOTYPES"),
        ("SAMPLE_ALIAS", "EXPECTED_SAMPLE_ALIAS"),
        ("H", "HAPLOTYPE_MAP"),
        ("LOD", "GENOTYPE_LOD_THRESHOLD"),
        ("IGNORE_RG", "IGNORE_READ_GROUPS"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let input = args.required("INPUT");
    let observed_alias = args.get("OBSERVED_SAMPLE_ALIAS").map(str::to_string);
    let output = args.get("OUTPUT").map(str::to_string);
    let summary_output = args.get("SUMMARY_OUTPUT").map(str::to_string);
    let detail_output = args.get("DETAIL_OUTPUT").map(str::to_string);
    let genotypes = args.required("GENOTYPES");
    let mut expected_alias = args.get("EXPECTED_SAMPLE_ALIAS").map(str::to_string);
    let haplotype_map = args.required("HAPLOTYPE_MAP");
    let lod_threshold = args.double("GENOTYPE_LOD_THRESHOLD", 5.0);
    let ignore_read_groups = args.bool("IGNORE_READ_GROUPS", false);
    let exit_not_found = args.int("EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND", 1) as i32;
    let exit_no_valid = args.int("EXIT_CODE_WHEN_NO_VALID_CHECKS", 2) as i32;

    // customCommandLineValidation.
    let reads = file_contains_reads(&input);
    if !reads && ignore_read_groups {
        refuse_validation(
            TOOL,
            &[
                "The parameter IGNORE_READ_GROUPS can only be used with BAM/SAM/CRAM inputs."
                    .to_string(),
            ],
        );
    }
    if reads && observed_alias.is_some() {
        refuse_validation(
            TOOL,
            &[
                "The parameter OBSERVED_SAMPLE_ALIAS can only be used with a VCF input."
                    .to_string(),
            ],
        );
    }
    if args.get("REFERENCE_SEQUENCE").is_none() && input.ends_with(".cram") {
        refuse_validation(
            TOOL,
            &["REFERENCE must be provided when using CRAM as input.".to_string()],
        );
    }

    let (detail_path, summary_path) = match &output {
        Some(o) => (
            format!("{o}.fingerprinting_detail_metrics"),
            format!("{o}.fingerprinting_summary_metrics"),
        ),
        None => (
            detail_output.unwrap_or_else(|| fail("DETAIL_OUTPUT is required")),
            summary_output.unwrap_or_else(|| fail("SUMMARY_OUTPUT is required")),
        ),
    };

    let map_text = std::fs::read_to_string(&haplotype_map).unwrap_or_else(|e| fail(&e.to_string()));
    let map = HaplotypeMap::from_database(&map_text, &absolute_path(&haplotype_map))
        .unwrap_or_else(|t| throw(t));

    let genotype_file = read_genotypes(&genotypes);
    let genotypes_uri = uri_of(&genotypes);

    // The observed sample, and the reads or genotypes it comes from.
    let mut sam: Option<(
        htsjdk_bam::header::SamHeader,
        Vec<htsjdk_bam::record::BamRecord>,
    )> = None;
    let mut observed_vcf = None;
    let observed_sample: Option<String> = if reads {
        let (header, records) = read_input(&input);
        let mut sample: Option<String> = None;
        for rg in &header.read_groups {
            let s = rg.attributes.get("SM").map(str::to_string);
            match &sample {
                None => sample = s,
                Some(existing) => {
                    if Some(existing) != s.as_ref() {
                        thrown(
                            "picard.PicardException: inputPath SAM/BAM file must not contain data from multiple samples.",
                        );
                    }
                }
            }
        }
        sam = Some((header, records));
        sample
    } else {
        let file = read_genotypes(&input);
        if file.samples.is_empty() {
            thrown("picard.PicardException: inputPath VCF file must contain at least one sample.");
        }
        if file.samples.len() > 1 && observed_alias.is_none() {
            thrown(
                "picard.PicardException: inputPath VCF file contains multiple samples and yet the OBSERVED_SAMPLE_ALIAS parameter is not set.",
            );
        }
        let sample = observed_alias
            .clone()
            .unwrap_or_else(|| file.samples[0].clone());
        if !file.samples.contains(&sample) {
            thrown(&format!(
                "picard.PicardException: inputPath VCF file does not contain OBSERVED_SAMPLE_ALIAS: {sample}"
            ));
        }
        observed_vcf = Some(file);
        Some(sample)
    };
    if expected_alias.is_none() {
        expected_alias = observed_sample.clone();
    }
    let contains = expected_alias
        .as_ref()
        .is_some_and(|s| genotype_file.samples.contains(s));
    if !contains {
        log(
            "WARN",
            TOOL,
            &format!(
                "Sample {} where fingerprint was expected not found in {}",
                expected_alias.as_deref().unwrap_or("null"),
                genotypes
            ),
        );
        std::process::exit(exit_not_found);
    }
    let expected_sample = expected_alias.clone().expect("checked above");

    // The expected fingerprints: the named sample's and every other sample's empty one.
    let expected: Vec<Fingerprint> = load_fingerprints(
        &genotype_file,
        &genotypes_uri,
        &map,
        Some(&expected_sample),
        0.01,
    )
    .unwrap_or_else(|t| throw(t))
    .iter()
    .map(|(_, fp)| fp.clone())
    .collect();
    if expected.is_empty() {
        thrown(&format!(
            "java.lang.IllegalStateException: Could not find any fingerprints in: [{}]",
            genotypes
        ));
    }

    // (read group, observed fingerprint) per result row.
    let mut observed: Vec<(Option<String>, Fingerprint)> = Vec::new();
    if let Some((header, records)) = &sam {
        let uri = uri_of(&input);
        let mut random = JavaRandom::new(42);
        let by_group = fingerprint_sam(
            header,
            records,
            &input,
            &uri,
            &map,
            &SamOptions::default(),
            &mut random,
            &picard_analysis::fingerprinting::Probabilities::sequence,
        )
        .unwrap_or_else(|t| throw(t));
        if ignore_read_groups {
            let mut combined = Fingerprint::new(Some(expected_sample.clone()), None, None);
            for (_, fp) in by_group.iter() {
                combined.merge(fp, &map).unwrap_or_else(|t| throw(t));
            }
            observed.push((None, combined));
        } else {
            for (id, fp) in by_group.iter() {
                observed.push((id.platform_unit.clone(), fp.clone()));
            }
        }
    } else if let Some(file) = &observed_vcf {
        let by_sample = load_fingerprints(
            file,
            &uri_of(&input),
            &map,
            observed_sample.as_deref(),
            0.01,
        )
        .unwrap_or_else(|t| throw(t));
        if by_sample.is_empty() {
            thrown(&format!(
                "java.lang.IllegalStateException: Found no fingerprints in observed genotypes file: [{input}]"
            ));
        }
        for (_, fp) in by_sample.iter() {
            observed.push((None, fp.clone()));
        }
    }

    let mut summary_file = MetricsFile::new();
    let mut detail_file = MetricsFile::new();
    for file in [&mut summary_file, &mut detail_file] {
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
    }
    let mut all_zero = true;
    for (read_group, fp) in &observed {
        let mut results: Vec<MatchResults> = Vec::new();
        for e in &expected {
            let r =
                calculate_match_results(fp, e, &map, 0.0, true, true).unwrap_or_else(|t| throw(t));
            // `TreeSet.add`: a result comparing equal to one already held is dropped.
            if !results
                .iter()
                .any(|held| compare_match_results(held, &r).is_eq())
            {
                results.push(r);
            }
        }
        results.sort_by(compare_match_results);
        let mr = &results[0];
        let mut s = Summary {
            read_group: read_group.clone(),
            sample: expected_sample.clone(),
            ll_expected: mr.sample_likelihood,
            ll_random: mr.population_likelihood,
            lod: mr.lod,
            with_genotypes: 0,
            checked: 0,
            matching: 0,
            het_as_hom: 0,
            hom_as_het: 0,
            hom_as_other_hom: 0,
        };
        for lr in &mr.locus_results {
            let expected_gt = lr.expected_genotype;
            let observed_gt = lr.most_likely_genotype;
            s.with_genotypes += 1;
            if lr.lod_genotype >= lod_threshold {
                s.checked += 1;
                if expected_gt == observed_gt {
                    s.matching += 1;
                }
                if expected_gt.is_heterozygous() && observed_gt.is_homozygous() {
                    s.het_as_hom += 1;
                }
                if expected_gt.is_homozygous() && observed_gt.is_heterozygous() {
                    s.hom_as_het += 1;
                }
                if expected_gt.is_homozygous()
                    && observed_gt.is_homozygous()
                    && expected_gt != observed_gt
                {
                    s.hom_as_other_hom += 1;
                }
            }
            let snp = &map.snps[lr.snp];
            detail_file.add_metric(&Detail {
                read_group: read_group.clone(),
                sample: expected_sample.clone(),
                snp: snp.name.clone(),
                alleles: snp.allele_string(),
                chrom: snp.chrom.clone(),
                position: i64::from(snp.pos),
                expected: expected_gt.name().to_string(),
                observed: observed_gt.name().to_string(),
                lod: lr.lod_genotype,
                obs_a: i64::from(lr.allele1_count),
                obs_b: i64::from(lr.allele2_count),
            });
        }
        all_zero &= s.lod == 0.0;
        summary_file.add_metric(&s);
    }
    std::fs::write(&summary_path, summary_file.write()).unwrap_or_else(|e| fail(&e.to_string()));
    std::fs::write(&detail_path, detail_file.write()).unwrap_or_else(|e| fail(&e.to_string()));
    if all_zero {
        log(
            "ERROR",
            TOOL,
            "No non-zero results found. This is likely an error. Probable cause: EXPECTED_SAMPLE \
             (if provided) or the sample name from INPUT (if EXPECTED_SAMPLE isn't provided)isn't \
             a sample in GENOTYPES file.",
        );
        std::process::exit(exit_no_valid);
    }
}
