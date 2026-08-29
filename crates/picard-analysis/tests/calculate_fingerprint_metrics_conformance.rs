//! Conformance for `CalculateFingerprintMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/checkfingerprint-conformance/`: six runs over a three-site haplotype map,
//! with the reads carrying one allele or the other.
//!
//! # What this suite is for
//!
//!  * **the counts being arithmetic over one fingerprint**, with no comparison involved;
//!  * **the chi-squared p-values being the reference library's**, to the last decimal the file
//!    prints;
//!  * **the cross-entropy being measured against the panel's own frequencies**;
//!  * **`DISCRIMINATORY_POWER` coming out of a hundred permutations of a seeded generator**, so a
//!    port either reproduces the generator or reproduces nothing;
//!  * **and a row being per fingerprint**, so what `--CALCULATE_BY` merges decides how many rows
//!    there are and what is in them.

use std::io::Read;

use picard_analysis::calculate_fingerprint_metrics::{
    accumulate, fingerprint_metrics, kl_divergence, randomized_lods, render_row, Haplotype,
    Observation, GENOTYPE_LOD_THRESHOLD, NUMBER_OF_SAMPLING, RANDOM_SEED,
};
use picard_analysis::extract_fingerprint::haplotype_frequencies;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/fingerprint_metrics.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

/// The rows of one case's table, the header dropped.
fn rows(text: &str, case: &str) -> Vec<String> {
    field(text, "metrics", case)
        .expect("the golden")
        .lines()
        .skip(1)
        .map(str::to_string)
        .collect()
}

/// The fixture's two blocks: `rs1` and `rs2` together at a minor allele frequency of four tenths,
/// and `rs3` alone at three.
fn blocks() -> [Haplotype; 2] {
    [
        Haplotype::new(haplotype_frequencies(0.4)),
        Haplotype::new(haplotype_frequencies(0.3)),
    ]
}

/// The reads of one read group, as the fixture writes them.
///
/// Thirty reads start at 99 and cover the first block's two sites; thirty start at 199 and cover
/// the second block's one site. A read group that was given only one of the two sets sees only
/// one block.
fn fingerprint(allele: u8, first_block: bool, second_block: bool) -> Vec<Haplotype> {
    let mut haplotypes = blocks().to_vec();
    let mut observations = Vec::new();
    for copy in 0..30 {
        if first_block {
            // One read over both of the block's sites, which is one read either way.
            for _site in 0..2 {
                observations.push(Observation {
                    read_name: format!("r{copy}"),
                    block: 0,
                    base: allele,
                    allele1: b'A',
                    allele2: b'C',
                    quality: 40,
                });
            }
        }
        if second_block {
            observations.push(Observation {
                read_name: format!("s{copy}"),
                block: 1,
                base: allele,
                allele1: b'A',
                allele2: b'C',
                quality: 40,
            });
        }
    }
    accumulate(&mut haplotypes, &observations);
    haplotypes
}

/// One file of one read group sees every site, and the row is the whole fingerprint.
#[test]
fn a_row_is_one_fingerprint() {
    let text = corpus();
    let whole = fingerprint(b'C', true, true);
    let row = fingerprint_metrics("sample1", "file://<dir>/in.bam", "unit1", &whole);
    assert_eq!(vec![render_row(&row)], rows(&text, "one-read-group"));

    // Both blocks have evidence, both are called definitely, and both are called homozygous for
    // the minor allele.
    assert_eq!(row.haplotypes, 2);
    assert_eq!(row.haplotypes_with_evidence, 2);
    assert_eq!(row.definite_genotypes, 2);
    assert_eq!(row.num_hom_allele2, 2);
    assert_eq!(row.num_het, 0);
}

/// A read group that saw one block of two is a row about a half-covered fingerprint.
#[test]
fn a_read_group_is_a_row_of_its_own() {
    let text = corpus();
    let first = fingerprint_metrics(
        "sample1",
        "file://<dir>/in.bam",
        "unit1",
        &fingerprint(b'C', true, false),
    );
    let second = fingerprint_metrics(
        "sample1",
        "file://<dir>/in.bam",
        "unit2",
        &fingerprint(b'C', false, true),
    );
    assert_eq!(
        vec![render_row(&first), render_row(&second)],
        rows(&text, "two-read-groups")
    );

    // The two rows are not the same row: the blocks differ in their allele frequencies, so the
    // same evidence over a different block is a different fingerprint.
    assert_ne!(first.lod_self_check, second.lod_self_check);
    assert_eq!(first.haplotypes_with_evidence, 1);
    assert_eq!(second.haplotypes_with_evidence, 1);

    // Two read groups of two samples are the same two rows under different names, because a row
    // is the fingerprint and not the sample.
    let two_samples = rows(&text, "two-samples");
    assert_eq!(
        two_samples[1].replace("sample2", "sample1"),
        rows(&text, "two-read-groups")[1]
    );
}

/// Rolling up puts the read groups back together, and the whole fingerprint comes back.
#[test]
fn rolling_up_merges_the_evidence() {
    let text = corpus();
    let whole = fingerprint(b'C', true, true);
    // By sample, the source is dropped and the INFO is the sample; by file, the INFO is the file
    // and the sample together. The numbers are the one-read-group row either way, because the
    // evidence is the same evidence.
    let by_sample = fingerprint_metrics("sample1", "", "sample1", &whole);
    assert_eq!(
        vec![render_row(&by_sample)],
        rows(&text, "two-read-groups-by-sample")
    );
    let by_file = fingerprint_metrics("sample1", "", "file://<dir>/in.bam::sample1", &whole);
    assert_eq!(
        vec![render_row(&by_file)],
        rows(&text, "two-read-groups-by-file")
    );
}

/// The other allele is a different fingerprint, and every column moves.
#[test]
fn the_major_allele_is_a_different_fingerprint() {
    let text = corpus();
    let major = fingerprint_metrics(
        "sample1",
        "file://<dir>/in.bam",
        "unit1",
        &fingerprint(b'A', true, true),
    );
    assert_eq!(vec![render_row(&major)], rows(&text, "the-major-allele"));

    // Homozygous for the major allele this time, and a fingerprint the panel finds unremarkable:
    // the p-value is a quarter rather than a thousandth.
    assert_eq!(major.num_hom_allele1, 2);
    assert_eq!(major.num_hom_allele2, 0);
    assert!(major.chi_squared_pvalue > 0.25);
}

/// The sampled column is a hundred permutations of a generator seeded with a constant.
#[test]
fn the_discriminatory_power_is_sampled() {
    let text = corpus();
    let whole = fingerprint(b'C', true, true);
    let row = fingerprint_metrics("sample1", "file://<dir>/in.bam", "unit1", &whole);

    let trials = randomized_lods(&whole, NUMBER_OF_SAMPLING, RANDOM_SEED);
    assert_eq!(trials.len(), 100);
    let mean: f64 = trials.iter().sum::<f64>() / trials.len() as f64;
    assert!((row.discriminatory_power - (row.lod_self_check - mean)).abs() < 1e-12);

    // The number in the golden is the one that seed produces, and a different seed produces a
    // different one, which is what makes the constant part of the answer.
    let golden = rows(&text, "one-read-group");
    let recorded: Vec<&str> = golden[0].split('\t').collect();
    assert_eq!(recorded[23], "13.633269");
    let other = randomized_lods(&whole, NUMBER_OF_SAMPLING, 43);
    let other_mean: f64 = other.iter().sum::<f64>() / other.len() as f64;
    assert_ne!(mean, other_mean);
}

/// The cross-entropy is measured against the panel's frequencies, not against the counts.
#[test]
fn the_cross_entropy_is_against_the_panel() {
    let text = corpus();
    let row = &rows(&text, "one-read-group")[0];
    let columns: Vec<&str> = row.split('\t').collect();
    // Two homozygous minor calls out of an expectation of a quarter: the cross entropy is the log
    // of the ratio, which is the log of eight.
    assert_eq!(columns[15], "2.079441");
    assert!((kl_divergence(&[0.0, 0.0, 2.0], &[0.85, 0.9, 0.25]) - 8f64.ln()).abs() < 1e-9);
    // And a definite genotype is one three logs better than the next.
    assert_eq!(GENOTYPE_LOD_THRESHOLD, 3.0);
}

/// A read is used once, whatever it covers.
#[test]
fn a_read_is_counted_once() {
    // The fixture's first block has two sites and its reads cover both, and the block ends up
    // with exactly the evidence the one-site block has: thirty reads, thirty observations.
    let both = fingerprint(b'C', true, true);
    let counted = |haplotype: &Haplotype| haplotype.log_likelihoods[0];
    assert_eq!(counted(&both[0]), counted(&both[1]));

    // Counting the same read twice would double the evidence, which is what the rule is there to
    // stop: the two sites of one block are in linkage, so one molecule is one observation.
    let mut doubled = blocks().to_vec();
    accumulate(
        &mut doubled,
        &(0..60)
            .map(|index| Observation {
                read_name: format!("r{index}"),
                block: 0,
                base: b'C',
                allele1: b'A',
                allele2: b'C',
                quality: 40,
            })
            .collect::<Vec<_>>(),
    );
    assert_eq!(counted(&doubled[0]), 2.0 * counted(&both[0]));
}
