//! Conformance for `CollectIndependentReplicateMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/replicatemetrics-conformance/CollectIndependentReplicateMetricsDump.java`,
//! twenty-three runs over one heterozygous site.
//!
//! # What this suite is for
//!
//!  * **a set that disagrees on the allele being what the estimate is built on**;
//!  * **the four levels the filters sit at, and what each of them leaves behind**;
//!  * **a set of four being neither a double nor a triple**;
//!  * **and the counter a run that examined nothing still increments.**

use std::io::Read;

use picard_analysis::collect_independent_replicate_metrics::{
    barcodes_are_usable, base_is_counted, classify_set, edit_distance, read_is_used, set_size,
    site_is_used, three_allele_sites_from_the_tail, SetClassification, SetSize,
    DEFAULT_MINIMUM_BARCODE_BQ, DEFAULT_MINIMUM_BQ, DEFAULT_MINIMUM_GQ, DEFAULT_MINIMUM_MQ,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/collect_independent_replicate_metrics.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

/// One case's metrics row, by column name.
fn metrics(text: &str, case: &str) -> Vec<(String, String)> {
    let prefix = format!("metrics\t{case}\t");
    let table = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
        .unwrap_or_else(|| panic!("metrics/{case}"));
    let mut lines = table.split('\n');
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let values: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    header
        .into_iter()
        .zip(values)
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

fn count(text: &str, case: &str, column: &str) -> i64 {
    metrics(text, case)
        .into_iter()
        .find(|(name, _)| name == column)
        .map(|(_, value)| value)
        .unwrap_or_else(|| panic!("{case}/{column}"))
        .parse()
        .unwrap_or_else(|_| panic!("{case}/{column} is a number"))
}

fn refusal(text: &str, case: &str) -> String {
    let prefix = format!("error\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .unwrap_or_else(|| panic!("error/{case}"))
}

/// A set that disagrees is what the estimate is built on, and one that agrees is filed by allele.
#[test]
fn a_set_that_disagrees_is_an_independent_replicate() {
    let text = corpus();
    assert_eq!(classify_set(1, 1, 0), SetClassification::DifferentAlleles);
    assert_eq!(
        count(&text, "doubleton-disagreeing", "nDifferentAllelesBiDups"),
        1
    );
    assert_eq!(
        count(&text, "doubleton-disagreeing", "nReferenceAllelesBiDups"),
        0
    );
    assert_eq!(classify_set(2, 0, 0), SetClassification::ReferenceAllele);
    assert_eq!(
        count(
            &text,
            "doubleton-agreeing-on-the-reference",
            "nReferenceAllelesBiDups"
        ),
        1
    );
    assert_eq!(classify_set(0, 2, 0), SetClassification::AlternateAllele);
    assert_eq!(
        count(
            &text,
            "doubleton-agreeing-on-the-alternate",
            "nAlternateAllelesBiDups"
        ),
        1
    );
    // A third allele makes the set mismatching whatever else it carries, the test coming first.
    assert_eq!(classify_set(1, 1, 1), SetClassification::MismatchingAllele);
    // And it takes the whole SITE out, which is a different counter again.
    assert_eq!(count(&text, "a-third-allele", "nSites"), 0);
    assert_eq!(count(&text, "a-third-allele", "nThreeAllelesSites"), 1);
}

/// The set sizes are counted apart, and a set of four is neither.
#[test]
fn a_set_of_four_is_neither_a_double_nor_a_triple() {
    let text = corpus();
    assert_eq!(set_size(2), SetSize::Doubleton);
    assert_eq!(set_size(3), SetSize::Tripleton);
    assert_eq!(set_size(4), SetSize::Big);
    assert_eq!(set_size(1), SetSize::Singleton);
    assert_eq!(count(&text, "doubleton-disagreeing", "nExactlyDouble"), 1);
    assert_eq!(count(&text, "doubleton-disagreeing", "nExactlyTriple"), 0);
    assert_eq!(count(&text, "tripleton", "nExactlyTriple"), 1);
    assert_eq!(count(&text, "tripleton", "nExactlyDouble"), 0);
    let four = metrics(&text, "a-set-of-four");
    let value = |column: &str| count(&text, "a-set-of-four", column);
    assert_eq!(value("nExactlyDouble"), 0);
    assert_eq!(value("nExactlyTriple"), 0);
    assert_eq!(value("nReadsInBigSets"), 4);
    assert_eq!(value("nTotalReads"), 4);
    // Its reads are counted, and none of the allele counters is: the set is one duplicate set.
    assert_eq!(value("nDuplicateSets"), 1);
    assert!(four
        .iter()
        .filter(|(name, _)| name.contains("AllelesBiDups") || name.contains("AllelesTriDups"))
        .all(|(_, value)| value == "0"));
}

/// The four levels the filters sit at, and what each leaves behind.
#[test]
fn the_filters_sit_at_four_levels() {
    let text = corpus();
    // The SITE: a homozygous genotype and one under the quality floor are both no site at all.
    assert!(!site_is_used(false, 99, DEFAULT_MINIMUM_GQ));
    assert!(!site_is_used(true, 10, DEFAULT_MINIMUM_GQ));
    assert!(site_is_used(true, 10, 5));
    assert_eq!(count(&text, "a-homozygous-site", "nSites"), 0);
    assert_eq!(count(&text, "a-low-quality-site", "nSites"), 0);
    assert_eq!(count(&text, "a-low-quality-site-allowed", "nSites"), 1);
    // The READ: a low mapping quality takes it out of the set before the set is built, so the
    // set that is left is a singleton and no doubleton is counted.
    assert!(!read_is_used(10, true, DEFAULT_MINIMUM_MQ, true));
    assert_eq!(count(&text, "a-low-mapping-quality-read", "nTotalReads"), 1);
    assert_eq!(
        count(&text, "a-low-mapping-quality-read", "nDuplicateSets"),
        0
    );
    // And an unpaired read is filtered by default, which leaves no site examined at all.
    assert!(!read_is_used(60, false, DEFAULT_MINIMUM_MQ, true));
    assert!(read_is_used(60, false, DEFAULT_MINIMUM_MQ, false));
    assert_eq!(count(&text, "unpaired-reads", "nTotalReads"), 0);
    assert_eq!(
        count(&text, "unpaired-reads-kept", "nDifferentAllelesBiDups"),
        1
    );
    // The BASE: the read stays in the set and its allele is gone, so the doubleton is still
    // counted and it is no longer a disagreeing one.
    assert!(!base_is_counted(3, DEFAULT_MINIMUM_BQ));
    assert!(base_is_counted(18, DEFAULT_MINIMUM_BQ));
    // The bound is strictly greater, so a base exactly at the floor is skipped.
    assert!(!base_is_counted(DEFAULT_MINIMUM_BQ, DEFAULT_MINIMUM_BQ));
    assert_eq!(count(&text, "a-low-base-quality", "nExactlyDouble"), 1);
    assert_eq!(
        count(&text, "a-low-base-quality", "nDifferentAllelesBiDups"),
        0
    );
    assert_eq!(
        count(
            &text,
            "a-low-base-quality-allowed",
            "nDifferentAllelesBiDups"
        ),
        1
    );
}

/// The barcode decides the UMI counters and nothing else.
#[test]
fn the_barcode_decides_only_the_umi_counters() {
    let text = corpus();
    assert!(barcodes_are_usable(
        &["IIII", "IIII"],
        DEFAULT_MINIMUM_BARCODE_BQ
    ));
    assert!(!barcodes_are_usable(
        &["!!!!", "IIII"],
        DEFAULT_MINIMUM_BARCODE_BQ
    ));
    // A read with no quality tag contributes an empty string, which has no base under the floor.
    assert!(barcodes_are_usable(&["", ""], DEFAULT_MINIMUM_BARCODE_BQ));
    // A bad barcode leaves the allele counters alone and moves the barcode counters.
    assert_eq!(
        count(&text, "a-low-quality-barcode", "nDifferentAllelesBiDups"),
        1
    );
    assert_eq!(count(&text, "a-low-quality-barcode", "nBadBarcodes"), 1);
    assert_eq!(count(&text, "doubleton-disagreeing", "nGoodBarcodes"), 1);
    assert_eq!(count(&text, "doubleton-disagreeing", "nBadBarcodes"), 0);
    // Reading the barcode from a tag nothing carries is a set of two empty barcodes, which are
    // usable and identical.
    assert_eq!(count(&text, "another-barcode-tag", "nGoodBarcodes"), 1);
    assert_eq!(
        count(&text, "another-barcode-tag", "nMatchingUMIsInDiffBiDups"),
        1
    );
    // The distance is a Hamming one over equal lengths, and the counters split on it being zero.
    assert_eq!(edit_distance("AAAA", "CCCC"), Some(4));
    assert_eq!(edit_distance("AAAA", "AAAA"), Some(0));
    assert_eq!(edit_distance("AAAA", "AAA"), None);
    assert_eq!(
        count(
            &text,
            "doubleton-disagreeing",
            "nMismatchingUMIsInDiffBiDups"
        ),
        1
    );
    assert_eq!(
        count(
            &text,
            "doubleton-disagreeing-same-umi",
            "nMatchingUMIsInDiffBiDups"
        ),
        1
    );
}

/// A run that examined nothing still reports one three-allele site.
#[test]
fn a_run_that_examined_nothing_reports_one_three_allele_site() {
    let text = corpus();
    assert_eq!(three_allele_sites_from_the_tail(false), 1);
    assert_eq!(three_allele_sites_from_the_tail(true), 0);
    for case in ["a-homozygous-site", "a-low-quality-site", "unpaired-reads"] {
        assert_eq!(count(&text, case, "nSites"), 0, "{case}");
        assert_eq!(count(&text, case, "nTotalReads"), 0, "{case}");
        assert_eq!(count(&text, case, "nThreeAllelesSites"), 1, "{case}");
    }
    // A run that did examine a site reports none of them, which is what says the counter is the
    // tail's and not the site's.
    assert_eq!(
        count(&text, "doubleton-disagreeing", "nThreeAllelesSites"),
        0
    );
}

/// The sample rules, and the short name that is not `S`.
#[test]
fn the_sample_may_be_omitted_only_when_there_is_one() {
    let text = corpus();
    let two = refusal(&text, "two-samples-without-a-name");
    assert!(
        two.contains("VCF must have exactly 1 sample. found 2"),
        "{two}"
    );
    // Named, the same VCF is accepted.
    assert_eq!(count(&text, "two-samples-with-a-name", "nSites"), 1);
    let missing = refusal(&text, "a-sample-that-is-not-there");
    assert!(
        missing.contains("Cannot find sample nobody in vcf"),
        "{missing}"
    );
    // And the short name is ALIAS: `S=` is not an option this tool has.
    let wrong = refusal(&text, "a-sample-by-the-wrong-short-name");
    assert!(wrong.contains("Unrecognized option: S"), "{wrong}");
}
