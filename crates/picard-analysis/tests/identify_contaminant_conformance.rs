//! Conformance for `IdentifyContaminant` against Picard 3.4.0.
//!
//! The tool is `ExtractFingerprint` with one argument negated, so the model ported for that one
//! serves here and this suite asks it the same questions under the other default.
//!
//! # What this suite is for
//!
//!  * **the default being the other tool's opposite**;
//!  * **the sample gaining `-contaminant` by default here**;
//!  * **`--EXTRACT_CONTAMINATED` restoring the other tool's default**;
//!  * **the flip being a no-op at 0.5 and visible at the extremes**;
//!  * **`--LOCUS_MAX_READS` defaulting to two hundred**;
//!  * **and the model being the same one.**

use std::io::Read;

use picard_analysis::extract_fingerprint::{
    contamination_to_use, extract_contamination_for_identify, haplotype_frequencies,
    phred_likelihoods, sample_to_use, wrong_fingerprint_count_message, ContaminatorProbabilities,
    DEFAULT_LOCUS_MAX_READS, IDENTIFY_CONTAMINANT_LOCUS_MAX_READS,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("identify_contaminant.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn records(text: &str, case: &str) -> Vec<(i32, String, i32, Vec<i32>)> {
    field(text, "out", case)
        .unwrap_or_else(|| panic!("{case} wrote a VCF"))
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let values: Vec<&str> = columns[9].split(':').collect();
            (
                columns[1].parse().expect("a position"),
                values[0].to_string(),
                values[1].parse().expect("a depth"),
                values[2]
                    .split(',')
                    .map(|v| v.parse().expect("a pl"))
                    .collect(),
            )
        })
        .collect()
}

fn minor_allele_frequency(position: i32) -> f64 {
    match position {
        101 | 105 => 0.4,
        201 => 0.3,
        other => panic!("no block at {other}"),
    }
}

fn ours(base: u8, quality: u8, depth: i32, contamination: f64, position: i32) -> Vec<i32> {
    let mut model = ContaminatorProbabilities::new(contamination);
    for _ in 0..depth {
        model.add(base, b'A', b'C', quality);
    }
    phred_likelihoods(
        model.log_likelihoods(haplotype_frequencies(minor_allele_frequency(position))),
    )
    .to_vec()
}

/// The default is the other tool's opposite, and the model is the same one.
#[test]
fn the_default_is_the_other_tools_opposite() {
    let text = corpus();
    assert!(extract_contamination_for_identify(false));
    assert!(!extract_contamination_for_identify(true));
    // A default run here uses the contamination as given, not one less it.
    let used = contamination_to_use(0.5, extract_contamination_for_identify(false));
    assert_eq!(used, 0.5);
    for (position, ad, depth, pl) in records(&text, "major-default") {
        assert_eq!(depth, 10);
        assert_eq!(ad, "10,0");
        assert_eq!(pl, ours(b'A', 40, 10, used, position));
    }
    for (position, ad, _, pl) in records(&text, "minor-default") {
        assert_eq!(ad, "0,10");
        assert_eq!(pl, ours(b'C', 40, 10, used, position));
    }
}

/// The sample gains `-contaminant` by default here, and loses it under the flag.
#[test]
fn the_sample_gains_the_suffix_by_default() {
    let text = corpus();
    assert_eq!(
        field(&text, "sample", "major-default").as_deref(),
        Some("sample1-contaminant")
    );
    assert_eq!(
        field(&text, "sample", "extract-contaminated").as_deref(),
        Some("sample1")
    );
    assert_eq!(
        field(&text, "sample", "sample-alias").as_deref(),
        Some("named")
    );
    assert_eq!(
        sample_to_use("sample1", None, extract_contamination_for_identify(false)),
        "sample1-contaminant"
    );
    assert_eq!(
        sample_to_use("sample1", None, extract_contamination_for_identify(true)),
        "sample1"
    );
}

/// At 0.5 the flip is a no-op, so only the name tells the two settings apart; at the extremes the
/// PLs are each other's.
#[test]
fn the_flip_shows_only_at_the_extremes() {
    let text = corpus();
    let default = records(&text, "major-default");
    let flagged = records(&text, "extract-contaminated");
    // Same numbers at 0.5.
    assert_eq!(
        default.iter().map(|r| r.3.clone()).collect::<Vec<_>>(),
        flagged.iter().map(|r| r.3.clone()).collect::<Vec<_>>()
    );
    assert_ne!(
        field(&text, "sample", "major-default"),
        field(&text, "sample", "extract-contaminated")
    );
    // And the two extremes are each other's.
    let nought = records(&text, "contamination-nought");
    let one = records(&text, "contamination-one");
    assert_eq!(nought[0].3, vec![0, 0, 0]);
    assert_eq!(one[0].3, vec![0, 30, 400]);
    assert_eq!(nought[0].3, ours(b'A', 40, 10, 0.0, 101));
    assert_eq!(one[0].3, ours(b'A', 40, 10, 1.0, 101));
}

/// The cap defaults to two hundred here and fifty there.
#[test]
fn the_cap_defaults_to_two_hundred() {
    let text = corpus();
    assert_eq!(IDENTIFY_CONTAMINANT_LOCUS_MAX_READS, 200);
    assert_eq!(DEFAULT_LOCUS_MAX_READS, 50);
    // A hundred reads pass uncapped.
    let uncapped = records(&text, "hundred-reads");
    assert!(uncapped.iter().all(|r| r.2 == 100));
    // The same hundred under an explicit fifty are bounded by the block, not the record.
    let capped = records(&text, "hundred-reads-capped");
    let first = capped.iter().find(|r| r.0 == 101).expect("rs1");
    let second = capped.iter().find(|r| r.0 == 201).expect("rs3");
    assert_eq!(first.2, 79);
    assert_eq!(second.2, 50);
    assert_eq!(first.3, ours(b'A', 40, 79, 0.5, 101));
    assert_eq!(second.3, ours(b'A', 40, 50, 0.5, 201));
}

/// The rest of the model answers as the other tool's does.
#[test]
fn the_model_is_the_same_one() {
    let text = corpus();
    for (_, ad, depth, pl) in records(&text, "neither-allele") {
        assert_eq!(ad, "0,0");
        assert_eq!(depth, 10);
        assert_eq!(pl, vec![0, 0, 0]);
    }
    for (_, ad, depth, pl) in records(&text, "no-reads") {
        assert_eq!(ad, "0,0");
        assert_eq!(depth, 0);
        assert_eq!(pl, vec![0, 0, 0]);
    }
    let error = field(&text, "error", "two-samples").expect("its refusal");
    assert_eq!(
        error,
        format!(
            "java.lang.IllegalArgumentException:{}",
            wrong_fingerprint_count_message(2)
        )
    );
}
