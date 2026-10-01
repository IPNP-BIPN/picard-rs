//! `CalculateFingerprintMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CalculateFingerprintMetrics.doWork` at tag 3.4.0:
//!
//! * every input fingerprinted in order (`fingerprintFiles` with one thread), merged by
//!   `CALCULATE_BY`, and one row per merged fingerprint in the `HashMap` order of the merged
//!   identities;
//! * the expected and observed genotype counts summed over the blocks in their `TreeMap` order,
//!   the three chi-squared tests (Commons Math, which refuses an expectation of zero), the three
//!   cross-entropies, and `NUM_HOM_ANY` read from the hom-versus-het pair's SECOND slot, which is
//!   the het count;
//! * `DISCRIMINATORY_POWER`: the self-check LOD less the running mean of a hundred LODs against
//!   permutations of the fingerprint, each block's log-likelihoods shuffled by a Mersenne Twister
//!   seeded with 42 afresh for every fingerprint;
//! * `GENOTYPE_LOD_THRESHOLD` and `NUMBER_OF_SAMPLING` are `final` fields initialised with
//!   constants, which javac inlines at their uses: the parser sets them and nothing reads them, so
//!   the port reads neither.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprinting::{
    absolute_path, calculate_match_results, fingerprint_files, merge_entries_by,
    p_normalize_vector, DataType, Fingerprint, HaplotypeMap, Probabilities, SamOptions,
};
use picard_analysis::math3::{chi_square_test, next_permutation, MersenneTwister};
use picard_analysis::metrics_cli::{fail, thrown, Args};

const TOOL: &str = "CalculateFingerprintMetrics";
const GENOTYPE_LOD_THRESHOLD: f64 = 3.0;
const NUMBER_OF_SAMPLING: usize = 100;
const RANDOM_SEED: i32 = 42;

struct Row {
    sample: Option<String>,
    source: String,
    info: Option<String>,
    longs: [i64; 7],
    doubles: [f64; 14],
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.fingerprint.FingerprintMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
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
        ]
    }
    fn values(&self) -> Vec<Value> {
        let mut out = vec![
            self.sample.clone().map_or(Value::Null, Value::Str),
            Value::Str(self.source.clone()),
            self.info.clone().map_or(Value::Null, Value::Str),
        ];
        out.extend(self.longs.iter().map(|v| Value::Long(*v)));
        out.extend(self.doubles.iter().map(|v| Value::Double(*v)));
        out
    }
}

fn throw(t: (String, String)) -> ! {
    thrown(&format!("{}: {}", t.0, t.1))
}

/// `ChiSquareTest.chiSquareTest`, which first checks every expectation is strictly positive.
fn chi_square(expected: &[f64], observed: &[i64]) -> f64 {
    for e in expected {
        if *e <= 0.0 {
            thrown(&format!(
                "org.apache.commons.math3.exception.NotStrictlyPositiveException: {} is smaller than, or equal to, the minimum (0)",
                java_double(*e)
            ));
        }
    }
    chi_square_test(expected, observed)
}

/// `Double.toString` for the values a refusal quotes.
fn java_double(value: f64) -> String {
    let text = format!("{value}");
    if text.contains('.') || text.contains('e') || text.contains("inf") || text.contains("NaN") {
        text
    } else {
        format!("{text}.0")
    }
}

/// `MathUtil.klDivergance`.
fn kl_divergence(measured: &[f64], distribution: &[f64]) -> f64 {
    let m = p_normalize_vector(measured);
    let d = p_normalize_vector(distribution);
    let mut sum = 0.0;
    for i in 0..m.len() {
        sum += m[i] * jmath::math::log(d[i] / m[i]);
    }
    -sum
}

fn sum3(acc: &mut [f64; 3], v: &[f64; 3]) {
    for i in 0..3 {
        acc[i] += v[i];
    }
}

fn metrics(fp: &Fingerprint, map: &HaplotypeMap) -> Row {
    let mut rng = MersenneTwister::new(RANDOM_SEED);
    // `reduce(MathUtil::sum)`, or zeros for an empty fingerprint.
    let mut counts = [0.0; 3];
    let mut expected = [0.0; 3];
    for (i, hp) in fp.map.values().enumerate() {
        let post = hp.posterior_probabilities(map);
        let prior = *hp.priors(map);
        if i == 0 {
            counts = post;
            expected = prior;
        } else {
            sum3(&mut counts, &post);
            sum3(&mut expected, &prior);
        }
    }
    let hom_vs_het_expect = [expected[0] + expected[2], expected[1]];
    let hom_vs_het_counts = [counts[0] + counts[2], counts[1]];
    let hom1_vs_hom2_expect = [expected[0], expected[2]];
    let hom1_vs_hom2_counts = [counts[0], counts[2]];
    let round = |v: &[f64]| -> Vec<i64> { v.iter().map(|x| jmath::math::round(*x)).collect() };
    let rounded_hom = round(&hom1_vs_hom2_counts);
    let rounded_het = round(&hom_vs_het_counts);
    let rounded = round(&counts);

    let haplotypes = fp.map.len() as i64;
    let with_evidence = fp.map.values().filter(|h| h.has_evidence(map)).count() as i64;
    let definite = fp
        .map
        .values()
        .filter(|h| h.lod_most_probable_genotype(map) >= GENOTYPE_LOD_THRESHOLD)
        .count() as i64;

    let chi = chi_square(&expected, &rounded);
    let cross = kl_divergence(&counts, &expected);
    let het_chi = chi_square(&hom_vs_het_expect, &rounded_het);
    let het_cross = kl_divergence(&hom_vs_het_counts, &hom_vs_het_expect);
    let hom_chi = chi_square(&hom1_vs_hom2_expect, &rounded_hom);
    let hom_cross = kl_divergence(&hom1_vs_hom2_counts, &hom1_vs_hom2_expect);

    // `RunningStat`: Knuth's running mean.
    let mut n = 0i64;
    let mut old_mean = 0.0;
    let mut new_mean = 0.0;
    for _ in 0..NUMBER_OF_SAMPLING {
        let randomized = randomize(fp, map, &mut rng);
        let lod = calculate_match_results(fp, &randomized, map, 0.0, true, true)
            .unwrap_or_else(|t| throw(t))
            .lod;
        n += 1;
        if n == 1 {
            old_mean = lod;
            new_mean = lod;
        } else {
            new_mean = old_mean + (lod - old_mean) / n as f64;
            old_mean = new_mean;
        }
    }
    let mean = if n > 0 { new_mean } else { 0.0 };
    let self_check = calculate_match_results(fp, fp, map, 0.0, true, true)
        .unwrap_or_else(|t| throw(t))
        .lod;

    Row {
        sample: fp.sample.clone(),
        source: fp.source.clone().unwrap_or_default(),
        info: fp.info.clone(),
        longs: [
            haplotypes,
            with_evidence,
            definite,
            rounded[0],
            rounded[2],
            rounded_het[1],
            rounded[1],
        ],
        doubles: [
            expected[0],
            expected[2],
            expected[1],
            chi,
            jmath::math::log10(chi),
            cross,
            het_chi,
            jmath::math::log10(het_chi),
            het_cross,
            hom_chi,
            jmath::math::log10(hom_chi),
            hom_cross,
            self_check,
            self_check - mean,
        ],
    }
}

/// `randomizeFingerprint`: each block's log-likelihoods permuted, into genotype-likelihood
/// evidence on the block's representative SNP.
fn randomize(fp: &Fingerprint, map: &HaplotypeMap, rng: &mut MersenneTwister) -> Fingerprint {
    let mut out = Fingerprint::new(None, None, None);
    for hp in fp.map.values() {
        let ll = hp.log_likelihoods(map);
        let permutation = next_permutation(rng, 3);
        let permuted = [ll[permutation[0]], ll[permutation[1]], ll[permutation[2]]];
        let snp = &map.snps[hp.representative_snp(map)];
        let mut p = Probabilities::genotype_likelihoods(hp.block);
        p.add_to_log_likelihoods(snp, (snp.allele1, snp.allele2), permuted);
        out.add(map, p);
    }
    out
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("H", "HAPLOTYPE_MAP")]);
    let inputs = args.all("INPUT");
    if inputs.is_empty() {
        fail("Argument 'INPUT' is required");
    }
    let output = args.required("OUTPUT");
    let haplotype_map = args.required("HAPLOTYPE_MAP");
    let by = DataType::parse(args.get("CALCULATE_BY").unwrap_or("READGROUP"))
        .unwrap_or_else(|| fail("Argument 'CALCULATE_BY' has an invalid value"));
    // Parsed for their refusals only: javac inlined both constants.
    let _ = args.double("GENOTYPE_LOD_THRESHOLD", GENOTYPE_LOD_THRESHOLD);
    let _ = args.int("NUMBER_OF_SAMPLING", NUMBER_OF_SAMPLING as i64);

    let map_text = std::fs::read_to_string(&haplotype_map).unwrap_or_else(|e| fail(&e.to_string()));
    let map = HaplotypeMap::from_database(&map_text, &absolute_path(&haplotype_map))
        .unwrap_or_else(|t| throw(t));
    let fingerprints = fingerprint_files(&inputs, &map, &SamOptions::default(), false)
        .unwrap_or_else(|t| throw(t));
    let merged = merge_entries_by(&fingerprints, by, &map).unwrap_or_else(|t| throw(t));

    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for (_, fp) in merged.iter() {
        file.add_metric(&metrics(fp, &map));
    }
    std::fs::write(&output, file.write()).unwrap_or_else(|e| fail(&e.to_string()));
}
