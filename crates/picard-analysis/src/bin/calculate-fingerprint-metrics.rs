//! `CalculateFingerprintMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CalculateFingerprintMetrics.doWork` at tag 3.4.0, over
//! `picard_analysis::fingerprint` and the loader in `picard_analysis::crosscheck_run`.
//!
//! `GENOTYPE_LOD_THRESHOLD` and `NUMBER_OF_SAMPLING` are declared `final` with constant
//! initialisers, so javac folds 3 and 100 into every use and the values the command line sets them
//! to are never read: they change nothing here either. The checker is never given the run's
//! stringency, so a read without a read group is refused whatever `VALIDATION_STRINGENCY` says.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::calculate_fingerprint_metrics::kl_divergence;
use picard_analysis::crosscheck_run::{fingerprint_files_concurrent, merge_by, DataType};
use picard_analysis::fingerprint::{
    calculate_match_results, file_uri, normalized_log_likelihoods, Evidence, Fingerprint,
    HaplotypeMap, Probs, SamOptions, SharedRandom,
};
use picard_analysis::math3::{chi_square_test, next_permutation, MersenneTwister};
use picard_analysis::metrics_cli::{thrown, Args};
use picard_analysis::theoretical_sensitivity::JavaRandom;

const TOOL: &str = "CalculateFingerprintMetrics";

const COLUMNS: &[&str] = &[
    "SAMPLE_ALIAS",
    "SOURCE",
    "INFO",
    "HAPLOTYPES",
    "HAPLOTYPES_WITH_EVIDENCE",
    "DEFINITE_GENOTYPES",
    "NUM_HOM_ALLELE1",
    "NUM_HOM_ALLELE2",
    "NUM_HOM_ANY",
    "NUM_HET",
    "EXPECTED_HOM_ALLELE1",
    "EXPECTED_HOM_ALLELE2",
    "EXPECTED_HET",
    "CHI_SQUARED_PVALUE",
    "LOG10_CHI_SQUARED_PVALUE",
    "CROSS_ENTROPY_LOD",
    "HET_CHI_SQUARED_PVALUE",
    "LOG10_HET_CHI_SQUARED_PVALUE",
    "HET_CROSS_ENTROPY_LOD",
    "HOM_CHI_SQUARED_PVALUE",
    "LOG10_HOM_CHI_SQUARED_PVALUE",
    "HOM_CROSS_ENTROPY_LOD",
    "LOD_SELF_CHECK",
    "DISCRIMINATORY_POWER",
];

struct Row(Vec<Value>);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.fingerprint.FingerprintMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        self.0.clone()
    }
}

/// `Math.round` over an array.
fn round(values: &[f64]) -> Vec<i64> {
    values.iter().map(|v| (v + 0.5).floor() as i64).collect()
}

/// `randomizeFingerprint`: each block's log-likelihoods permuted, as genotype likelihoods.
fn randomized(fp: &Fingerprint, rng: &mut MersenneTwister) -> Fingerprint {
    let mut out = Fingerprint::new(None, None, None);
    for (key, hp) in &fp.blocks {
        let ll = hp.log_likelihoods();
        let permutation = next_permutation(rng, 3);
        let mut permuted = [0.0; 3];
        for (i, source) in permutation.iter().enumerate() {
            permuted[i] = ll[*source];
        }
        out.blocks.insert(
            key.clone(),
            Probs {
                block: hp.block,
                priors: hp.priors,
                first: hp.first.clone(),
                evidence: Evidence::GenotypeLikelihoods {
                    ll: normalized_log_likelihoods(permuted),
                },
            },
        );
    }
    out
}

/// `getFingerprintMetrics`.
fn metrics(map: &HaplotypeMap, fp: &Fingerprint) -> Row {
    let mut rng = MersenneTwister::new(42);
    let probs: Vec<&Probs> = fp.blocks.values().collect();
    let sum3 = |f: &dyn Fn(&Probs) -> [f64; 3]| -> [f64; 3] {
        let mut total = [0.0; 3];
        for (i, p) in probs.iter().enumerate() {
            let v = f(p);
            if i == 0 {
                total = v;
            } else {
                for g in 0..3 {
                    total[g] += v[g];
                }
            }
        }
        total
    };
    let counts = sum3(&|p| p.posterior_probabilities());
    let expected = sum3(&|p| p.priors);
    let hom_vs_het_expect = [expected[0] + expected[2], expected[1]];
    let hom_vs_het_counts = [counts[0] + counts[2], counts[1]];
    let hom_expect = [expected[0], expected[2]];
    let hom_counts = [counts[0], counts[2]];
    let rounded_hom = round(&hom_counts);
    let rounded_hom_vs_het = round(&hom_vs_het_counts);
    let rounded = round(&counts);
    let chi = chi_square_test(&expected, &rounded);
    let het_chi = chi_square_test(&hom_vs_het_expect, &rounded_hom_vs_het);
    let hom_chi = chi_square_test(&hom_expect, &rounded_hom);
    // `MathUtils.RunningStat`: Knuth's running mean, not a sum divided.
    let mut mean = 0.0;
    for n in 1..=100 {
        let lod = calculate_match_results(map, fp, &randomized(fp, &mut rng), 0.0, true, true).lod;
        mean = if n == 1 {
            lod
        } else {
            mean + (lod - mean) / f64::from(n)
        };
    }
    let self_check = calculate_match_results(map, fp, fp, 0.0, true, true).lod;
    let text = |v: &Option<String>| v.clone().map_or(Value::Null, Value::Str);
    Row(vec![
        text(&fp.sample),
        Value::Str(fp.source.as_deref().map(file_uri).unwrap_or_default()),
        text(&fp.info),
        Value::Long(probs.len() as i64),
        Value::Long(probs.iter().filter(|p| p.has_evidence()).count() as i64),
        Value::Long(
            probs
                .iter()
                .filter(|p| p.lod_most_probable_genotype() >= 3.0)
                .count() as i64,
        ),
        Value::Long(rounded[0]),
        Value::Long(rounded[2]),
        Value::Long(rounded_hom_vs_het[1]),
        Value::Long(rounded[1]),
        Value::Double(expected[0]),
        Value::Double(expected[2]),
        Value::Double(expected[1]),
        Value::Double(chi),
        Value::Double(chi.log10()),
        Value::Double(kl_divergence(&counts, &expected)),
        Value::Double(het_chi),
        Value::Double(het_chi.log10()),
        Value::Double(kl_divergence(&hom_vs_het_counts, &hom_vs_het_expect)),
        Value::Double(hom_chi),
        Value::Double(hom_chi.log10()),
        Value::Double(kl_divergence(&hom_counts, &hom_expect)),
        Value::Double(self_check),
        Value::Double(self_check - mean),
    ])
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let inputs = args.collection("INPUT", &[]);
    let output = args.required("OUTPUT");
    let by = DataType::parse(args.get("CALCULATE_BY").unwrap_or("READGROUP"))
        .unwrap_or(DataType::ReadGroup);
    let map = HaplotypeMap::load(&args.required("HAPLOTYPE_MAP")).unwrap_or_else(|e| thrown(&e));
    let mut random = SharedRandom(JavaRandom::new(42));
    let fps =
        fingerprint_files_concurrent(&inputs, &map, &SamOptions::default(), false, &mut random)
            .unwrap_or_else(|e| thrown(&e));
    let merged = merge_by(&fps, by);
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for (_, fp) in &merged {
        file.add_metric(&metrics(&map, fp));
    }
    if let Err(e) = std::fs::write(&output, file.write()) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
