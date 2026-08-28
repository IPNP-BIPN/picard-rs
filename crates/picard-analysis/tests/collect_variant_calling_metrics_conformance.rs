//! Conformance for `CollectVariantCallingMetrics` against Picard 3.4.0.
//!
//! Each case carries the detail and summary tables the tool wrote. The golden prints the input VCF
//! and the dbSNP one, so the variants are read back from the reference's own files.
//!
//! # What this suite is for
//!
//!  * **the output argument being a prefix**;
//!  * **known and novel being counted apart, and summing to the total**;
//!  * **the two TI/TV ratios being counted separately**;
//!  * **a filtered variant reaching neither tally**;
//!  * **indels and multiallelic SNPs having columns of their own**;
//!  * **the summary not being the detail rows' sum**;
//!  * **`--TARGET_INTERVALS` dropping what falls outside**;
//!  * **and an empty VCF writing two kinds of empty ratio.**

use std::io::Read;

use picard_analysis::accumulate_variant_calling_metrics::{DETAIL_EXTENSION, SUMMARY_EXTENSION};
use picard_analysis::collect_variant_calling_metrics::{
    accumulate, collect, file_names, is_transition, Counts, Variant,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_variant_calling_metrics.txt.gz");
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

fn rows(text: &str, kind: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let table = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = table.lines().filter(|line| !line.is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

fn number(row: &std::collections::HashMap<String, String>, name: &str) -> f64 {
    let value = row.get(name).unwrap_or_else(|| panic!("{name}"));
    if value == "?" {
        f64::NAN
    } else {
        value.parse().unwrap_or_else(|_| panic!("{name}={value}"))
    }
}

/// The variants the golden's own input VCF holds, against its own dbSNP one.
fn variants(text: &str) -> (Vec<Variant>, Vec<String>) {
    let vcf = field(text, "in", "mixed").expect("the vcf");
    let db = field(text, "in", "dbsnp").expect("the dbsnp");
    let known: Vec<i64> = db
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            line.split('\t')
                .nth(1)
                .expect("a position")
                .parse()
                .expect("a number")
        })
        .collect();
    let mut samples = Vec::new();
    let mut out = Vec::new();
    for line in vcf.lines() {
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let position: i64 = c[1].parse().expect("a position");
        out.push(Variant {
            reference: c[3].to_string(),
            alternates: c[4].split(',').map(str::to_string).collect(),
            filtered: c[6] != "PASS" && c[6] != ".",
            in_db_snp: known.contains(&position),
            genotypes: c[9..]
                .iter()
                .map(|genotype| {
                    let gt = genotype.split(':').next().expect("a GT");
                    let mut alleles: Vec<Option<usize>> = gt
                        .split(['/', '|'])
                        .map(|a| a.parse::<usize>().ok())
                        .collect();
                    alleles.resize(2, None);
                    [alleles[0], alleles[1]]
                })
                .collect(),
        });
    }
    (out, samples)
}

/// The output argument is a prefix, and the two extensions are the accumulator's own.
#[test]
fn the_output_argument_is_a_prefix() {
    let (detail, summary) = file_names("out");
    assert_eq!(detail, format!("out.{DETAIL_EXTENSION}"));
    assert_eq!(summary, format!("out.{SUMMARY_EXTENSION}"));
}

/// The golden's own input, counted by the port, is what the tool wrote.
#[test]
fn the_printed_input_counts_to_the_written_row() {
    let text = corpus();
    let (variants, samples) = variants(&text);
    let (details, summary) = collect(&variants, &samples);
    let theirs = &rows(&text, "detail", "known-and-novel")[0];
    let ours = &details["s1"];
    assert_eq!(ours.total_snps.to_string(), theirs["TOTAL_SNPS"]);
    assert_eq!(ours.num_in_db_snp.to_string(), theirs["NUM_IN_DB_SNP"]);
    assert_eq!(ours.novel_snps.to_string(), theirs["NOVEL_SNPS"]);
    assert_eq!(ours.num_singletons.to_string(), theirs["NUM_SINGLETONS"]);
    assert_eq!(ours.db_snp_titv(), number(theirs, "DBSNP_TITV"));
    assert_eq!(ours.novel_titv(), number(theirs, "NOVEL_TITV"));
    assert_eq!(ours.pct_db_snp(), number(theirs, "PCT_DBSNP"));
    // And the summary row of the same file.
    let their_summary = &rows(&text, "summary", "known-and-novel")[0];
    assert_eq!(summary.total_snps.to_string(), their_summary["TOTAL_SNPS"]);
}

/// Known and novel are counted apart and sum to the total.
#[test]
fn known_and_novel_sum_to_the_total() {
    let text = corpus();
    for case in ["known-and-novel", "nothing-known", "everything-known"] {
        let row = &rows(&text, "detail", case)[0];
        let total: i64 = row["TOTAL_SNPS"].parse().expect("a count");
        let known: i64 = row["NUM_IN_DB_SNP"].parse().expect("a count");
        let novel: i64 = row["NOVEL_SNPS"].parse().expect("a count");
        assert_eq!(known + novel, total, "{case}");
    }
    assert_eq!(
        rows(&text, "detail", "nothing-known")[0]["NUM_IN_DB_SNP"],
        "0"
    );
    assert_eq!(
        rows(&text, "detail", "everything-known")[0]["NOVEL_SNPS"],
        "0"
    );
}

/// The two TI/TV ratios are counted separately, and the fixture makes them differ.
#[test]
fn the_two_ratios_are_counted_apart() {
    let text = corpus();
    let row = &rows(&text, "detail", "known-and-novel")[0];
    assert_eq!(number(row, "DBSNP_TITV"), 2.0);
    assert_eq!(number(row, "NOVEL_TITV"), 0.5);
    // Which is what the port's own classification gives.
    assert!(is_transition(b'A', b'G'));
    assert!(is_transition(b'C', b'T'));
    assert!(is_transition(b'T', b'C'));
    assert!(!is_transition(b'A', b'C'));
    assert!(!is_transition(b'G', b'T'));
    assert!(!is_transition(b'A', b'T'));
}

/// A filtered variant is counted as filtered and nowhere else.
#[test]
fn a_filtered_variant_reaches_neither_tally() {
    let text = corpus();
    let row = &rows(&text, "detail", "filtered")[0];
    assert_eq!(row["FILTERED_SNPS"], "1");
    assert_eq!(row["TOTAL_SNPS"], "1");
    assert_eq!(row["NUM_IN_DB_SNP"], "1");
    assert_eq!(row["NOVEL_SNPS"], "0");
    // The port keeps the same books: the filtered one raises the total and the filtered count.
    let mut counts = Counts::default();
    accumulate(
        &mut counts,
        &Variant {
            reference: "C".to_string(),
            alternates: vec!["T".to_string()],
            filtered: true,
            in_db_snp: true,
            genotypes: vec![[Some(0), Some(1)]],
        },
    );
    assert_eq!(counts.total_snps, 1);
    assert_eq!(counts.filtered_snps, 1);
    assert_eq!(counts.num_in_db_snp, 0);
    assert_eq!(counts.novel_snps, 0);
}

/// Indels and multiallelic SNPs have columns of their own.
#[test]
fn indels_and_multiallelics_are_counted_apart() {
    let text = corpus();
    let indels = &rows(&text, "detail", "indels")[0];
    assert_eq!(indels["TOTAL_INDELS"], "2");
    assert_eq!(indels["TOTAL_SNPS"], "1");
    let multi = &rows(&text, "detail", "multiallelic")[0];
    assert_eq!(multi["TOTAL_MULTIALLELIC_SNPS"], "1");
    assert_eq!(multi["TOTAL_SNPS"], "1");
    // A multiallelic SNP is not among the plain ones.
    let variant = Variant {
        reference: "C".to_string(),
        alternates: vec!["T".to_string(), "A".to_string()],
        filtered: false,
        in_db_snp: true,
        genotypes: vec![[Some(0), Some(1)]],
    };
    assert!(variant.is_multiallelic());
    let mut counts = Counts::default();
    accumulate(&mut counts, &variant);
    assert_eq!(counts.total_multiallelic_snps, 1);
    assert_eq!(counts.total_snps, 0);
}

/// The summary is not the detail rows' sum.
#[test]
fn the_summary_is_not_the_details_sum() {
    let text = corpus();
    let details = rows(&text, "detail", "two-samples");
    let summary = &rows(&text, "summary", "two-samples")[0];
    assert_eq!(details.len(), 2);
    let s1: i64 = details
        .iter()
        .find(|r| r["SAMPLE_ALIAS"] == "s1")
        .expect("s1")["TOTAL_SNPS"]
        .parse()
        .expect("a count");
    let s2: i64 = details
        .iter()
        .find(|r| r["SAMPLE_ALIAS"] == "s2")
        .expect("s2")["TOTAL_SNPS"]
        .parse()
        .expect("a count");
    let total: i64 = summary["TOTAL_SNPS"].parse().expect("a count");
    assert_eq!((s1, s2, total), (4, 2, 4));
    assert_ne!(s1 + s2, total);
    // The summary row is the file's, so it names no sample.
    assert!(!summary.contains_key("SAMPLE_ALIAS"));
}

/// The target intervals drop what falls outside them.
#[test]
fn the_target_intervals_drop_the_rest() {
    let text = corpus();
    let targeted = &rows(&text, "detail", "targeted")[0];
    let all = &rows(&text, "detail", "known-and-novel")[0];
    assert_eq!(all["TOTAL_SNPS"], "6");
    assert_eq!(targeted["TOTAL_SNPS"], "2");
    // The two that survive are the two known ones inside the interval.
    assert_eq!(targeted["NUM_IN_DB_SNP"], "2");
    assert_eq!(targeted["NOVEL_SNPS"], "0");
}

/// An empty VCF writes both files, with two different kinds of empty ratio.
#[test]
fn an_empty_vcf_writes_two_kinds_of_empty() {
    let text = corpus();
    let row = &rows(&text, "detail", "empty")[0];
    assert_eq!(row["TOTAL_SNPS"], "0");
    // Nought over nought is NaN, which the writer renders as `?`.
    assert_eq!(row["PCT_DBSNP"], "?");
    assert!(number(row, "PCT_DBSNP").is_nan());
    // Where the TI/TV columns come out as a plain zero.
    assert_eq!(row["DBSNP_TITV"], "0");
    assert_eq!(row["NOVEL_TITV"], "0");
    // The port's own empty counts reach the NaN the same way.
    let counts = Counts::default();
    assert!(counts.pct_db_snp().is_nan());
    assert!(rows(&text, "summary", "empty").len() == 1);
}
