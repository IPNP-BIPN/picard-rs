//! Conformance for `CheckFingerprint` against Picard 3.4.0.
//!
//! Golden from `tools/checkfingerprint-conformance/CheckFingerprintDump.java`, eleven runs over
//! the `extractfingerprint` suite's own haplotype map.
//!
//! # What this suite is for
//!
//!  * **the LOD being the difference of two log-likelihoods, so agreement is positive**;
//!  * **an uncovered block keeping its row, its likelihoods being the priors**;
//!  * **the two exit codes, one of which writes its files anyway**;
//!  * **and the two file names one basename builds.**

use std::io::Read;

use picard_analysis::check_fingerprint::{
    block_contribution, file_names, match_result, posterior_likelihoods,
    shifted_log_evidence_probability, BlockContribution, EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND,
    EXIT_CODE_WHEN_NO_VALID_CHECKS,
};
use picard_analysis::extract_fingerprint::haplotype_frequencies;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/check_fingerprint.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

/// One table's rows, by column name.
fn rows(text: &str, kind: &str, case: &str) -> Vec<Vec<(String, String)>> {
    let body = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = body.split('\n').filter(|line| !line.is_empty());
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .map(str::to_string)
        .collect();
    lines
        .map(|line| {
            header
                .iter()
                .cloned()
                .zip(line.split('\t').map(str::to_string))
                .collect()
        })
        .collect()
}

fn number(text: &str, kind: &str, case: &str, column: &str) -> f64 {
    rows(text, kind, case)[0]
        .iter()
        .find(|(name, _)| name == column)
        .map(|(_, value)| value.parse::<f64>().expect("a number"))
        .unwrap_or_else(|| panic!("{kind}/{case}/{column}"))
}

/// The LOD is the difference of the two models' log-likelihoods.
#[test]
fn the_lod_is_the_difference_of_the_two_models() {
    let text = corpus();
    for case in ["agreeing", "disagreeing", "heterozygous-genotypes"] {
        let expected = number(&text, "summary", case, "LL_EXPECTED_SAMPLE");
        let random = number(&text, "summary", case, "LL_RANDOM_SAMPLE");
        let lod = number(&text, "summary", case, "LOD_EXPECTED_SAMPLE");
        // The file rounds to six places, so the comparison is against that rounding.
        assert!((expected - random - lod).abs() < 1e-6, "{case}");
        let result = match_result(&[BlockContribution {
            no_swap: expected,
            swap: random,
        }]);
        assert!((result.lod_expected_sample - lod).abs() < 1e-6, "{case}");
    }
    // Agreement is positive and disagreement is not, which is the whole of the claim.
    assert!(number(&text, "summary", "agreeing", "LOD_EXPECTED_SAMPLE") > 0.0);
    assert!(number(&text, "summary", "disagreeing", "LOD_EXPECTED_SAMPLE") < 0.0);
    // And reads that agree with a homozygous-alternate genotype are as positive as the reference
    // case, the score being about the agreement and not about the allele.
    assert!(
        number(
            &text,
            "summary",
            "homozygous-alternate-genotypes",
            "LOD_EXPECTED_SAMPLE"
        ) > 0.0
    );
}

/// The arithmetic itself: a block's two numbers, from likelihoods and frequencies.
#[test]
fn a_blocks_contribution_is_two_dot_products() {
    // The map's first block has a minor allele frequency of 0.4.
    let priors = haplotype_frequencies(0.4);
    // A fingerprint certain of the reference homozygote, on both sides.
    let certain = [1.0, 0.0, 0.0];
    let block = block_contribution(certain, certain, priors);
    // Under one sample the evidence is that genotype's prior ONCE, the other fingerprint's
    // posterior standing in for it; under two it is that prior TWICE, once per sample. The
    // difference is what the LOD is, and for a certain agreement it is the prior's own log.
    assert!((block.no_swap - priors[0].log10()).abs() < 1e-12);
    assert!((block.swap - 2.0 * priors[0].log10()).abs() < 1e-12);
    assert!(block.no_swap > block.swap);
    // Two fingerprints that disagree give a no-swap term of zero probability, which is what makes
    // the LOD go negative.
    let other = [0.0, 0.0, 1.0];
    let disagreeing = block_contribution(certain, other, priors);
    assert!(disagreeing.no_swap.is_infinite() && disagreeing.no_swap.is_sign_negative());
    // The posterior likelihoods are the numerator of Bayes' rule and are not normalised.
    let posterior = posterior_likelihoods([0.5, 0.25, 0.25], priors);
    assert!((posterior.iter().sum::<f64>() - 1.0).abs() > 1e-6);
    assert_eq!(
        shifted_log_evidence_probability([1.0, 1.0, 1.0], priors),
        priors.iter().sum::<f64>().log10()
    );
}

/// A block counts only when both fingerprints have evidence for it.
#[test]
fn an_uncovered_block_keeps_its_row_and_its_priors() {
    let text = corpus();
    // The map has three sites in two blocks, and a run over both blocks reports two detail rows.
    assert_eq!(rows(&text, "detail", "agreeing").len(), 2);
    assert_eq!(
        number(&text, "summary", "agreeing", "HAPLOTYPES_WITH_GENOTYPES"),
        2.0
    );
    // A file covering ONE block writes the same two rows: the uncovered block's observed
    // likelihoods are the priors, so it keeps its row, reads zero observations, and still counts
    // among the haplotypes with genotypes.
    let uncovered = rows(&text, "detail", "one-block-covered");
    assert_eq!(uncovered.len(), 2);
    assert_eq!(
        number(
            &text,
            "summary",
            "one-block-covered",
            "HAPLOTYPES_WITH_GENOTYPES"
        ),
        2.0
    );
    let second = &uncovered[1];
    let value = |column: &str| {
        second
            .iter()
            .find(|(name, _)| name == column)
            .map(|(_, value)| value.clone())
            .unwrap_or_else(|| panic!("{column}"))
    };
    assert_eq!(
        (value("OBS_A"), value("OBS_B")),
        ("0".to_string(), "0".to_string())
    );
    // Its own LOD is small and positive, which is the priors agreeing with themselves.
    let lod: f64 = value("LOD").parse().expect("a number");
    assert!(lod > 0.0 && lod < 0.5, "{lod}");
    // And the run's total is smaller than the one over a file that covers both.
    assert!(
        number(&text, "summary", "one-block-covered", "LOD_EXPECTED_SAMPLE")
            < number(&text, "summary", "agreeing", "LOD_EXPECTED_SAMPLE")
    );
    // And the depth is in the detail row, so a shallower pileup is the same rows with smaller
    // numbers rather than fewer of them.
    assert_eq!(rows(&text, "detail", "one-read-each").len(), 2);
    let deep: f64 = rows(&text, "detail", "agreeing")[0]
        .iter()
        .find(|(name, _)| name == "OBS_A")
        .map(|(_, value)| value.parse().expect("a number"))
        .expect("OBS_A");
    let shallow: f64 = rows(&text, "detail", "one-read-each")[0]
        .iter()
        .find(|(name, _)| name == "OBS_A")
        .map(|(_, value)| value.parse().expect("a number"))
        .expect("OBS_A");
    assert_eq!((deep, shallow), (10.0, 1.0));
}

/// The two exit codes, and the files one basename builds.
#[test]
fn the_codes_carry_answers_the_files_do_not() {
    let text = corpus();
    // A sample the genotypes do not carry: code one, and nothing written.
    let refusal = field(&text, "error", "a-sample-the-genotypes-do-not-have").expect("the code");
    assert_eq!(
        refusal,
        format!("exit {EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND}")
    );
    assert_eq!(
        field(&text, "files", "a-sample-the-genotypes-do-not-have").as_deref(),
        Some("")
    );
    // A run with nothing to check: code two, and its files written anyway.
    for case in ["no-reads", "the-other-sample-named"] {
        assert_eq!(
            field(&text, "error", case).unwrap_or_else(|| panic!("{case}")),
            format!("exit {EXIT_CODE_WHEN_NO_VALID_CHECKS}"),
            "{case}"
        );
        let written = field(&text, "files", case).unwrap_or_else(|| panic!("{case}"));
        assert_eq!(written.split(' ').count(), 2, "{case}");
    }
    // The names are one basename and two suffixes.
    let mut expected = file_names("check").to_vec();
    expected.sort();
    let mut written: Vec<String> = field(&text, "files", "agreeing")
        .expect("the files")
        .split(' ')
        .map(str::to_string)
        .collect();
    written.sort();
    assert_eq!(written, expected);
    // Naming the expected sample explicitly finds it where the default did not.
    assert_eq!(field(&text, "error", "an-expected-sample-named"), None);
    assert!(
        number(
            &text,
            "summary",
            "an-expected-sample-named",
            "LOD_EXPECTED_SAMPLE"
        ) > 0.0
    );
}
