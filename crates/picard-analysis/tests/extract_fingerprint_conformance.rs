//! Conformance for `ExtractFingerprint` against Picard 3.4.0.
//!
//! Each case carries the file the tool read and the VCF it wrote. The pileup is not ported, so
//! the port is given the bases the fixture put at each site and asked for the same PLs.
//!
//! # What this suite is for
//!
//!  * **the PLs being the contaminator's, so an all-major pileup calls the minor genotype**;
//!  * **the contamination argument being flipped rather than the output**;
//!  * **the block's minor-allele frequency being the prior**;
//!  * **a base matching neither allele reaching DP and nothing else**;
//!  * **a base under the quality floor reaching nothing at all**;
//!  * **`--LOCUS_MAX_READS` bounding the block and not the record**;
//!  * **the sample name gaining `-contaminant` unless aliased or contaminated**;
//!  * **one record per representative SNP unless every SNP is asked for**;
//!  * **and a file naming two samples being refused by a message counting them.**

use std::io::Read;

use picard_analysis::extract_fingerprint::{
    contamination_to_use, error_probability, haplotype_frequencies, phred_likelihoods,
    sample_to_use, wrong_fingerprint_count_message, ContaminatorProbabilities,
    DEFAULT_LOCUS_MAX_READS, IDENTIFY_CONTAMINANT_LOCUS_MAX_READS,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("extract_fingerprint.txt.gz");
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

/// One case's records, as (position, AD, DP, PL).
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
                    .map(|v| v.parse().expect("a likelihood"))
                    .collect(),
            )
        })
        .collect()
}

/// The two blocks of the fixture's haplotype map, by their representative SNP's position.
fn minor_allele_frequency(position: i32) -> f64 {
    match position {
        101 | 105 => 0.4,
        201 => 0.3,
        other => panic!("no block at {other}"),
    }
}

/// The port's PLs for a pileup of one base repeated, at one site.
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

/// The PLs are the contaminator's: an all-major pileup at half contamination is not a flat call.
#[test]
fn the_pls_are_the_contaminators() {
    let text = corpus();
    let major = records(&text, "major-half-contaminated");
    let minor = records(&text, "minor-half-contaminated");
    assert_eq!(major.len(), 2);
    for (position, ad, depth, pl) in &major {
        assert_eq!(*depth, 10, "{position}");
        assert_eq!(ad, "10,0", "{position}");
        assert_eq!(
            *pl,
            ours(b'A', 40, 10, contamination_to_use(0.5, false), *position),
            "{position}"
        );
    }
    for (position, ad, depth, pl) in &minor {
        assert_eq!(*depth, 10, "{position}");
        assert_eq!(ad, "0,10", "{position}");
        assert_eq!(
            *pl,
            ours(b'C', 40, 10, contamination_to_use(0.5, false), *position),
            "{position}"
        );
    }
    // The two are mirror images: the major pileup's best contaminant genotype is the minor's
    // worst and the other way about.
    assert_eq!(major[0].3[0], 0);
    assert_eq!(minor[0].3[2], 0);
}

/// The contamination argument is flipped rather than the output, so the same number means
/// opposite things under the two settings.
#[test]
fn the_contamination_argument_is_flipped() {
    assert_eq!(contamination_to_use(0.5, true), 0.5);
    assert_eq!(contamination_to_use(0.5, false), 0.5);
    assert_eq!(contamination_to_use(0.0, true), 0.0);
    assert_eq!(contamination_to_use(0.0, false), 1.0);
    let text = corpus();
    // A contamination of nought under the default setting, and of nought under the other, are the
    // two ends of the same scale.
    let nought = records(&text, "contamination-nought");
    let flipped = records(&text, "extract-contamination-nought");
    assert_ne!(nought[0].3, flipped[0].3);
    // The default run is EXTRACT_CONTAMINATION=false, which flips the argument to one.
    assert_eq!(
        nought[0].3,
        ours(b'A', 40, 10, contamination_to_use(0.0, false), 101)
    );
    assert_eq!(
        flipped[0].3,
        ours(b'A', 40, 10, contamination_to_use(0.0, true), 101)
    );
}

/// The block's minor-allele frequency is the prior, so two blocks of the same pileup differ.
#[test]
fn the_frequency_is_the_prior() {
    let text = corpus();
    let rows = records(&text, "middling-base-quality");
    assert_eq!(rows.len(), 2);
    let first = rows.iter().find(|r| r.0 == 101).expect("rs1");
    let second = rows.iter().find(|r| r.0 == 201).expect("rs3");
    assert_ne!(first.3, second.3);
    assert_eq!(first.3, vec![0, 13, 30]);
    assert_eq!(second.3, vec![0, 12, 30]);
    assert_eq!(first.3, ours(b'A', 20, 10, 0.5, 101));
    assert_eq!(second.3, ours(b'A', 20, 10, 0.5, 201));
    // Hardy-Weinberg, which is where the difference comes from.
    let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-12);
    assert!(close(haplotype_frequencies(0.4), [0.36, 0.48, 0.16]));
    assert!(close(haplotype_frequencies(0.3), [0.49, 0.42, 0.09]));
}

/// A base matching neither allele reaches DP and nothing else, and one under the quality floor
/// reaches nothing at all.
#[test]
fn the_two_kinds_of_ignored_base() {
    let text = corpus();
    let neither = records(&text, "neither-allele");
    for (_, ad, depth, pl) in &neither {
        assert_eq!(ad, "0,0");
        assert_eq!(*depth, 10);
        assert_eq!(*pl, vec![0, 0, 0]);
    }
    let low = records(&text, "base-quality-under-the-floor");
    for (_, ad, depth, pl) in &low {
        assert_eq!(ad, "0,0");
        // Not even the depth: the base never reached the pileup.
        assert_eq!(*depth, 0);
        assert_eq!(*pl, vec![0, 0, 0]);
    }
    // The port counts the third kind apart and moves no likelihood for it.
    let mut model = ContaminatorProbabilities::new(0.5);
    model.add(b'G', b'A', b'C', 40);
    assert_eq!(model.observed_other, 1);
    assert_eq!(model.log_likelihoods, [[0.0; 3]; 3]);
}

/// The cap bounds the block and not the record.
#[test]
fn the_cap_bounds_the_block() {
    let text = corpus();
    let uncapped = records(&text, "deep-uncapped");
    assert!(uncapped.iter().all(|r| r.2 == 40));
    let capped = records(&text, "deep-capped");
    let first = capped.iter().find(|r| r.0 == 101).expect("rs1");
    let second = capped.iter().find(|r| r.0 == 201).expect("rs3");
    // Ten was asked for; the block of two SNPs reports sixteen and the block of one reports ten.
    assert_eq!(first.2, 16);
    assert_eq!(second.2, 10);
    assert_eq!(first.3, ours(b'A', 40, 16, 0.5, 101));
    assert_eq!(second.3, ours(b'A', 40, 10, 0.5, 201));
    // The two tools' defaults are different numbers.
    assert_eq!(DEFAULT_LOCUS_MAX_READS, 50);
    assert_eq!(IDENTIFY_CONTAMINANT_LOCUS_MAX_READS, 200);
}

/// The sample name gains `-contaminant` unless aliased or contaminated.
#[test]
fn the_sample_name_says_which_it_is() {
    let text = corpus();
    // The default run is EXTRACT_CONTAMINATION=false, so no suffix.
    assert_eq!(
        field(&text, "sample", "major-half-contaminated").as_deref(),
        Some("sample1")
    );
    assert_eq!(
        field(&text, "sample", "sample-alias").as_deref(),
        Some("named")
    );
    assert_eq!(
        field(&text, "sample", "extract-contamination").as_deref(),
        Some("sample1-contaminant")
    );
    // Another tool's argument named here is an exit code and not a warning.
    assert_eq!(
        field(&text, "error", "unknown-argument").as_deref(),
        Some("exit 1")
    );
    assert_eq!(sample_to_use("sample1", None, true), "sample1-contaminant");
    assert_eq!(sample_to_use("sample1", None, false), "sample1");
    assert_eq!(sample_to_use("sample1", Some("named"), true), "named");
    // The alias replaces rather than adds.
    assert_eq!(sample_to_use("sample1", Some("named"), false), "named");
}

/// One record per representative SNP, unless every SNP is asked for.
#[test]
fn one_record_per_block_unless_told_otherwise() {
    let text = corpus();
    let representatives = records(&text, "major-half-contaminated");
    assert_eq!(
        representatives.iter().map(|r| r.0).collect::<Vec<_>>(),
        vec![101, 201]
    );
    let all = records(&text, "all-snps");
    assert_eq!(
        all.iter().map(|r| r.0).collect::<Vec<_>>(),
        vec![101, 105, 201]
    );
    // The two SNPs of one block carry the same PLs, being one model read twice.
    assert_eq!(all[0].3, all[1].3);
}

/// A file naming two samples is refused by a message counting the fingerprints.
#[test]
fn two_samples_are_refused() {
    let text = corpus();
    let error = field(&text, "error", "two-samples").expect("its refusal");
    assert_eq!(
        error,
        format!(
            "java.lang.IllegalArgumentException:{}",
            wrong_fingerprint_count_message(2)
        )
    );
}

/// A site with no reads still gets a record, and its PLs are nought across the board.
#[test]
fn a_site_with_no_reads_still_gets_a_record() {
    let text = corpus();
    let rows = records(&text, "no-reads");
    assert_eq!(rows.len(), 2);
    for (_, ad, depth, pl) in &rows {
        assert_eq!(ad, "0,0");
        assert_eq!(*depth, 0);
        assert_eq!(*pl, vec![0, 0, 0]);
    }
    // Which the port reaches too: no base moves any likelihood, so the three are equal and the
    // phred shift takes them all to nought.
    assert_eq!(ours(b'A', 40, 0, 0.5, 101), vec![0, 0, 0]);
}

/// The phred conversion, which the PLs rest on.
#[test]
fn the_phred_conversion_shifts_to_the_best() {
    assert_eq!(phred_likelihoods([0.0, -1.0, -2.0]), [0, 10, 20]);
    assert_eq!(phred_likelihoods([-2.0, -1.0, 0.0]), [20, 10, 0]);
    assert_eq!(phred_likelihoods([-5.0, -5.0, -5.0]), [0, 0, 0]);
    assert!((error_probability(40) - 1e-4).abs() < 1e-12);
    assert!((error_probability(20) - 1e-2).abs() < 1e-12);
    assert!((error_probability(0) - 1.0).abs() < 1e-12);
}
