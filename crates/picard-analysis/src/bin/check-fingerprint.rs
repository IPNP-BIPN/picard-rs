//! `CheckFingerprint` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CheckFingerprint.doWork` at tag 3.4.0, over
//! `picard_analysis::fingerprint`.
//!
//! # The best match is not always the expected sample's
//!
//! `loadFingerprints` names a sample to read, and then adds an EMPTY fingerprint for every other
//! sample in the file's header (`computeIfAbsent`). Each is compared, the results sit in a
//! `TreeSet` ordered by LOD, and the summary row takes the first. An empty fingerprint scores 0,
//! so when the expected sample's own LOD is negative the row reports one of the empty ones: no
//! haplotypes, a LOD of 0, and with every row at 0 the run exits `EXIT_CODE_WHEN_NO_VALID_CHECKS`.
//! The expected sample's name is written whichever fingerprint won.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprint::{
    calculate_match_results, fingerprint_sam_file, java_compare, load_vcf_fingerprints,
    Fingerprint, HaplotypeMap, MatchResults, Probs, SamOptions, SharedRandom,
};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};
use picard_analysis::theoretical_sensitivity::JavaRandom;

const TOOL: &str = "CheckFingerprint";

struct Row {
    class: &'static str,
    columns: &'static [&'static str],
    values: Vec<Value>,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        self.class
    }
    fn columns(&self) -> &[&'static str] {
        self.columns
    }
    fn values(&self) -> Vec<Value> {
        self.values.clone()
    }
}

const SUMMARY: &[&str] = &[
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
];

const DETAIL: &[&str] = &[
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
];

/// `SamReader.Type`'s extensions, tested on the path as `fileContainsReads` tests the URI.
fn contains_reads(path: &str) -> bool {
    path.ends_with(".bam") || path.ends_with(".sam") || path.ends_with(".cram")
}

fn text(value: Option<&str>) -> Value {
    value.map_or(Value::Null, |v| Value::Str(v.to_string()))
}

fn is_het(genotype: &str) -> bool {
    genotype.as_bytes()[0] != genotype.as_bytes()[1]
}

/// `FingerprintResults`: one input fingerprint against every expected one, best first.
struct Results {
    read_group: Option<String>,
    matches: Vec<MatchResults>,
}

fn compare_all(
    map: &HaplotypeMap,
    observed: &Fingerprint,
    expected: &[(String, Fingerprint)],
    read_group: Option<String>,
) -> Results {
    let mut matches: Vec<MatchResults> = expected
        .iter()
        .map(|(_, fp)| calculate_match_results(map, observed, fp, 0.0, true, true))
        .collect();
    // `MatchResults.compareTo`: LOD descending, then sample; equal ones collapse in the TreeSet.
    matches.sort_by(|a, b| {
        b.lod
            .partial_cmp(&a.lod)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                java_compare(
                    a.sample.as_deref().unwrap_or(""),
                    b.sample.as_deref().unwrap_or(""),
                )
            })
    });
    matches.dedup_by(|a, b| a.lod == b.lod && a.sample == b.sample);
    Results {
        read_group,
        matches,
    }
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
    let genotypes = args.required("GENOTYPES");
    let map_path = args.required("HAPLOTYPE_MAP");
    let observed_alias = args.get("OBSERVED_SAMPLE_ALIAS").map(str::to_string);
    let ignore_read_groups = args.bool("IGNORE_READ_GROUPS", false);
    let lod_threshold = args.double("GENOTYPE_LOD_THRESHOLD", 5.0);
    let not_found_code = args.int("EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND", 1) as i32;
    let no_valid_code = args.int("EXIT_CODE_WHEN_NO_VALID_CHECKS", 2) as i32;
    let strict = !matches!(
        args.get("VALIDATION_STRINGENCY"),
        Some("LENIENT") | Some("SILENT")
    );

    let reads = contains_reads(&input);
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

    let (detail_path, summary_path) = match args.get("OUTPUT") {
        Some(prefix) => (
            format!("{prefix}.fingerprinting_detail_metrics"),
            format!("{prefix}.fingerprinting_summary_metrics"),
        ),
        None => (
            args.required("DETAIL_OUTPUT"),
            args.required("SUMMARY_OUTPUT"),
        ),
    };

    let map = HaplotypeMap::load(&map_path).unwrap_or_else(|e| thrown(&e));

    // `extractObservedSampleName`.
    let observed_sample: Option<String> = if reads {
        let (header, _) = picard_analysis::metrics_cli::read_input(&input);
        let mut observed: Option<String> = None;
        for rg in &header.read_groups {
            let sample = rg.attributes.get("SM").map(str::to_string);
            match &observed {
                None => observed = sample,
                Some(o) if Some(o) != sample.as_ref() => thrown(
                    "picard.PicardException: inputPath SAM/BAM file must not contain data from multiple samples.",
                ),
                _ => {}
            }
        }
        observed
    } else {
        let vcf = picard_analysis::vcf_io::read_path(&input).unwrap_or_else(|e| thrown(&e));
        let samples = &vcf.file.header.samples;
        if samples.is_empty() {
            thrown("picard.PicardException: inputPath VCF file must contain at least one sample.");
        }
        if samples.len() > 1 && observed_alias.is_none() {
            thrown("picard.PicardException: inputPath VCF file contains multiple samples and yet the OBSERVED_SAMPLE_ALIAS parameter is not set.");
        }
        let observed = observed_alias.clone().unwrap_or_else(|| samples[0].clone());
        if !samples.contains(&observed) {
            thrown(&format!(
                "picard.PicardException: inputPath VCF file does not contain OBSERVED_SAMPLE_ALIAS: {observed}"
            ));
        }
        Some(observed)
    };
    let expected_sample = args
        .get("EXPECTED_SAMPLE_ALIAS")
        .map(str::to_string)
        .or_else(|| observed_sample.clone());

    let genotype_vcf =
        picard_analysis::vcf_io::read_path(&genotypes).unwrap_or_else(|e| thrown(&e));
    if !expected_sample
        .as_ref()
        .is_some_and(|s| genotype_vcf.file.header.samples.contains(s))
    {
        eprintln!(
            "WARN\t{TOOL}\tSample {} where fingerprint was expected not found in {genotypes}",
            expected_sample.as_deref().unwrap_or("null")
        );
        std::process::exit(not_found_code);
    }

    let expected = load_vcf_fingerprints(&genotypes, &map, expected_sample.as_deref(), 0.01)
        .unwrap_or_else(|e| thrown(&e));
    if expected.is_empty() {
        thrown(&format!(
            "java.lang.IllegalStateException: Could not find any fingerprints in: [{genotypes}]"
        ));
    }

    let mut random = SharedRandom(JavaRandom::new(42));
    let results: Vec<Results> = if reads {
        let options = SamOptions {
            strict,
            ..SamOptions::default()
        };
        let by_group = fingerprint_sam_file(
            &input,
            &map,
            &options,
            &mut random,
            &|m, b| Probs::sequence(m, b),
            &|p, snp, base, qual| p.add_base(snp, base, qual),
        )
        .unwrap_or_else(|e| thrown(&e));
        if ignore_read_groups {
            let mut combined = Fingerprint::new(expected_sample.clone(), Some(input.clone()), None);
            for (_, fp) in &by_group {
                combined.merge(fp);
            }
            vec![compare_all(&map, &combined, &expected, None)]
        } else {
            by_group
                .iter()
                .map(|(details, fp)| {
                    compare_all(&map, fp, &expected, details.platform_unit.clone())
                })
                .collect()
        }
    } else {
        let observed = load_vcf_fingerprints(&input, &map, observed_sample.as_deref(), 0.01)
            .unwrap_or_else(|e| thrown(&e));
        observed
            .iter()
            .map(|(_, fp)| compare_all(&map, fp, &expected, None))
            .collect()
    };

    let mut summary = MetricsFile::new();
    summary.add_header(&format!("{TOOL} <command line>"));
    summary.add_header("Started on: <timestamp>");
    let mut details = MetricsFile::new();
    details.add_header(&format!("{TOOL} <command line>"));
    details.add_header("Started on: <timestamp>");
    let mut all_zero = true;
    let sample = expected_sample.as_deref();
    for result in &results {
        let mr = &result.matches[0];
        let (mut with, mut checked, mut matching, mut het_hom, mut hom_het, mut hom_other) =
            (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
        for lr in &mr.locus_results {
            let (e, o) = (&lr.expected_genotype, &lr.most_likely_genotype);
            with += 1;
            if lr.lod_genotype >= lod_threshold {
                checked += 1;
                if e == o {
                    matching += 1;
                }
                if is_het(e) && !is_het(o) {
                    het_hom += 1;
                }
                if !is_het(e) && is_het(o) {
                    hom_het += 1;
                }
                if !is_het(e) && !is_het(o) && e != o {
                    hom_other += 1;
                }
            }
            details.add_metric(&Row {
                class: "picard.analysis.FingerprintingDetailMetrics",
                columns: DETAIL,
                values: vec![
                    text(result.read_group.as_deref()),
                    text(sample),
                    Value::Str(lr.snp.name.clone()),
                    Value::Str(lr.snp.allele_string()),
                    Value::Str(lr.snp.chrom.clone()),
                    Value::Long(i64::from(lr.snp.pos)),
                    Value::Str(e.clone()),
                    Value::Str(o.clone()),
                    Value::Double(lr.lod_genotype),
                    Value::Long(i64::from(lr.allele1_count)),
                    Value::Long(i64::from(lr.allele2_count)),
                ],
            });
        }
        summary.add_metric(&Row {
            class: "picard.analysis.FingerprintingSummaryMetrics",
            columns: SUMMARY,
            values: vec![
                text(result.read_group.as_deref()),
                text(sample),
                Value::Double(mr.sample_likelihood),
                Value::Double(mr.population_likelihood),
                Value::Double(mr.lod),
                Value::Long(with),
                Value::Long(checked),
                Value::Long(matching),
                Value::Long(het_hom),
                Value::Long(hom_het),
                Value::Long(hom_other),
            ],
        });
        all_zero &= mr.lod == 0.0;
    }
    for (path, file) in [(&summary_path, &summary), (&detail_path, &details)] {
        if let Err(e) = std::fs::write(path, file.write()) {
            thrown(&format!("htsjdk.samtools.SAMException: {e}"));
        }
    }
    if all_zero {
        eprintln!("ERROR\t{TOOL}\tNo non-zero results found. This is likely an error. Probable cause: EXPECTED_SAMPLE (if provided) or the sample name from INPUT (if EXPECTED_SAMPLE isn't provided)isn't a sample in GENOTYPES file.");
        std::process::exit(no_valid_code);
    }
}
