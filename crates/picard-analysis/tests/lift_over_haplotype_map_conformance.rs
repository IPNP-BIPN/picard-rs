//! Conformance for `LiftOverHaplotypeMap` against Picard 3.4.0.
//!
//! Each case carries the haplotype database and the chain the tool read, the lifted table it
//! wrote and the exit code it returned. The port reads the same table with
//! [`picard_analysis::haplotype_map`], lifts it with [`htsjdk_bam::liftover`], and must write the
//! same rows.
//!
//! # What this suite is for
//!
//!  * **the exit code for a failed liftover being 101 and not 1**;
//!  * **a SNP that does not lift being dropped while the file is still written**;
//!  * **the alleles never being complemented across a negative-strand chain**;
//!  * **the frequency staying the minor allele's, reformatted on the way out**;
//!  * **the anchor column being rewritten rather than carried**;
//!  * **the panels being carried, comma-separated**;
//!  * **a dictionary missing a contig the chain lifts to being refused before any SNP**;
//!  * **and a contig the chain does not cover not being refused at all.**

use std::io::Read;

use htsjdk_bam::liftover::LiftOver;
use picard_analysis::haplotype_map::{format_frequency, parse_haplotype_database, Row};
use picard_analysis::lift_over_haplotype_map::{
    lift_over_haplotype_map, LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("lift_over_haplotype_map.txt.gz");
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

/// The rows of a written table, without its header or its column line.
fn written(text: &str, case: &str) -> Vec<String> {
    field(text, "out", case)
        .unwrap_or_else(|| panic!("{case} wrote a table"))
        .lines()
        .filter(|line| !line.starts_with('@') && !line.starts_with('#') && !line.is_empty())
        .map(str::to_string)
        .collect()
}

fn code(text: &str, case: &str) -> i32 {
    field(text, "code", case)
        .unwrap_or_else(|| panic!("{case} returned a code"))
        .trim()
        .parse()
        .expect("a code")
}

/// `writeToFile`'s line: the eight columns, the last two possibly empty.
fn render(row: &Row) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.chromosome,
        row.position,
        row.name,
        row.major_allele as char,
        row.minor_allele as char,
        format_frequency(row.minor_allele_frequency),
        row.anchor.clone().unwrap_or_default(),
        row.panels.clone().unwrap_or_default()
    )
}

fn run(text: &str, database: &str, chain: &str, contigs: &[&str]) -> (Vec<String>, i32) {
    let blocks = parse_haplotype_database(&field(text, "db", database).expect("a database"))
        .expect("the database parses");
    let lift = LiftOver::load(&field(text, "chain", chain).expect("a chain")).expect("the chain");
    let order: Vec<String> = contigs.iter().map(|name| (*name).to_string()).collect();
    let result = lift_over_haplotype_map(&blocks, &lift, &order).expect("the dictionary is whole");
    (result.rows.iter().map(render).collect(), result.return_code)
}

/// The cases that write a table, with the database, the chain and the dictionary each was given.
const WRITTEN: &[(&str, &str, &str, &[&str])] = &[
    ("all-lift", "liftable", "two-blocks", &["chrA", "chrB"]),
    ("one-fails", "partly", "two-blocks", &["chrA", "chrB"]),
    (
        "whole-block-fails",
        "whole-block-fails",
        "two-blocks",
        &["chrA", "chrB"],
    ),
    (
        "negative-strand",
        "reversed",
        "two-blocks",
        &["chrA", "chrB"],
    ),
    ("no-chain-for-contig", "reversed", "one-block", &["chrA"]),
];

/// Every case that writes a table writes the rows the port produces, and returns the same code.
#[test]
fn every_case_writes_the_same_rows() {
    let text = corpus();
    for (case, database, chain, contigs) in WRITTEN {
        let (rows, return_code) = run(&text, database, chain, contigs);
        assert_eq!(rows, written(&text, case), "{case}");
        assert_eq!(return_code, code(&text, case), "{case}");
    }
}

/// The exit code for a failed liftover is 101 and not 1, whatever the count of failures.
#[test]
fn a_failed_liftover_is_a_hundred_and_one() {
    let text = corpus();
    assert_eq!(code(&text, "all-lift"), 0);
    for case in [
        "one-fails",
        "whole-block-fails",
        "no-chain-for-contig",
        "everything-fails",
    ] {
        assert_eq!(
            code(&text, case),
            LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS,
            "{case}"
        );
    }
    assert_eq!(LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS, 101);
}

/// A SNP that does not lift is dropped and the file is still written: a database whose every SNP
/// fails leaves a table holding nothing.
#[test]
fn a_dropped_snp_still_leaves_a_file() {
    let text = corpus();
    assert_eq!(written(&text, "one-fails").len(), 1);
    assert_eq!(written(&text, "whole-block-fails").len(), 1);
    assert!(written(&text, "everything-fails").is_empty());
    assert!(written(&text, "no-chain-for-contig").is_empty());
    // The file exists all the same: it carries its header and its column line.
    let empty = field(&text, "out", "everything-fails").expect("a file");
    assert!(empty.contains("#CHROMOSOME\tPOSITION\tNAME"), "{empty}");
}

/// The alleles are carried over unchanged across the negative-strand chain, and the frequency
/// with them: `chr2:50 A/C` becomes `chrB:851 A/C`.
#[test]
fn the_alleles_are_not_complemented() {
    let text = corpus();
    let rows = written(&text, "negative-strand");
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0], "chrB\t851\trs1\tA\tC\t0.1\t\t");
    // The input said chr2:50 A/C, and the complement of A/C would be T/G.
    let database = field(&text, "db", "reversed").expect("the database");
    assert!(database.contains("chr2\t50\trs1\tA\tC\t0.10"), "{database}");
}

/// The frequency is reformatted on the way out: `0.10` comes back as `0.1`.
#[test]
fn the_frequency_is_reformatted() {
    let text = corpus();
    for row in written(&text, "all-lift") {
        let frequency = row.split('\t').nth(5).expect("a frequency");
        assert!(!frequency.ends_with('0') || frequency == "0", "{row}");
    }
    assert_eq!(format_frequency(0.10), "0.1");
    assert_eq!(format_frequency(0.20), "0.2");
}

/// The anchor column is rewritten: the first row of a block by position gets an EMPTY anchor and
/// the later ones get that first row's name, whatever the input named.
#[test]
fn the_anchor_column_is_rewritten() {
    let text = corpus();
    let rows = written(&text, "all-lift");
    assert_eq!(rows.len(), 2);
    let anchors: Vec<&str> = rows
        .iter()
        .map(|row| row.split('\t').nth(6).expect("an anchor"))
        .collect();
    assert_eq!(anchors, vec!["", "rs1"]);
    // The input named rs2, the SECOND row, as the anchor of both.
    let database = field(&text, "db", "liftable").expect("the database");
    assert!(database.contains("rs1\tA\tC\t0.10\trs2"), "{database}");
    // And the panels are carried, comma-separated.
    let panels: Vec<&str> = rows
        .iter()
        .map(|row| row.split('\t').nth(7).expect("panels"))
        .collect();
    assert_eq!(panels, vec!["panelA", "panelA,panelB"]);
}

/// A dictionary that does not name a contig the chain lifts to is refused before any SNP is
/// looked at, while a contig the chain does not cover is not refused at all.
#[test]
fn a_missing_target_contig_is_refused() {
    let text = corpus();
    let error = field(&text, "error", "dictionary-missing-contig").expect("the refusal");
    assert!(
        error.contains("Sequence chrB from chain file is not found in sequence dictionary."),
        "{error}"
    );
    let blocks = parse_haplotype_database(&field(&text, "db", "liftable").expect("a database"))
        .expect("the database parses");
    let lift = LiftOver::load(&field(&text, "chain", "two-blocks").expect("a chain")).expect("ok");
    assert!(lift_over_haplotype_map(&blocks, &lift, &["chrA".to_string()]).is_err());
    // The one-block chain over the same dictionary is fine, and its SNP simply fails.
    let one = LiftOver::load(&field(&text, "chain", "one-block").expect("a chain")).expect("ok");
    let reversed = parse_haplotype_database(&field(&text, "db", "reversed").expect("a database"))
        .expect("the database parses");
    let result =
        lift_over_haplotype_map(&reversed, &one, &["chrA".to_string()]).expect("not refused");
    assert!(result.rows.is_empty());
    assert_eq!(result.return_code, LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS);
}

/// A block spanning two contigs never reaches the liftover: the reader refuses it.
#[test]
fn a_block_spanning_two_contigs_is_refused_by_the_reader() {
    let text = corpus();
    let error = field(&text, "error", "block-across-contigs").expect("the refusal");
    assert!(
        error.contains("does not agree with chromosome of existing snp(s)"),
        "{error}"
    );
    assert!(
        parse_haplotype_database(&field(&text, "db", "two-contigs").expect("a database")).is_err()
    );
}
