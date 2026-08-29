//! Conformance for `LiftoverVcf` against Picard 3.4.0.
//!
//! Golden from `tools/liftovervcf-conformance/`: eleven runs over a chain with one forward block
//! and one reversed one.
//!
//! # What this suite is for
//!
//!  * **a variant inside a block being renumbered by the block's offset**;
//!  * **a variant outside every block being rejected rather than dropped**;
//!  * **a reference allele the target does not carry being a different rejection**, carrying the
//!    locus and the alleles that were attempted;
//!  * **`--RECOVER_SWAPPED_REF_ALT` turning one of those into a lift**, genotypes and allele
//!    frequency and all;
//!  * **a reversed block complementing the alleles and counting from the other end**;
//!  * **the original position and alleles being recorded only where they are worth recording**;
//!  * **and the output being sorted by the target's coordinates rather than the input's.**

use std::io::Read;

use picard_analysis::liftover_vcf::{
    format_vcf_double, lift_over, liftover, parse_chains, render, run, Attribute, Interval,
    Options, Outcome, Variant,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/liftover_vcf.txt.gz");
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

const LENGTH: usize = 400;

/// The target reference: `ACGT` repeating, so a base is known by arithmetic.
fn reference() -> Vec<u8> {
    (0..LENGTH).map(|index| b"ACGT"[index % 4]).collect()
}

/// One forward block and one reversed one, as the fixture writes them.
fn chains() -> Vec<picard_analysis::liftover_vcf::Chain> {
    parse_chains(&format!(
        "chain 100 chr1 {LENGTH} + 0 100 chrT {LENGTH} + 12 112 1\n100\n\n\
         chain 100 chr1 {LENGTH} + 200 300 chrT {LENGTH} - 100 200 2\n100\n\n"
    ))
}

/// One record of the source VCF: a position, a reference, an alternate and a genotype.
fn variant(position: i32, reference: &str, alternate: &str, genotype: &str) -> Variant {
    Variant {
        contig: "chr1".to_string(),
        position,
        id: ".".to_string(),
        reference: reference.to_string(),
        alternates: vec![alternate.to_string()],
        quality: "100".to_string(),
        filters: Vec::new(),
        attributes: vec![Attribute {
            key: "AF".to_string(),
            value: Some("0.25".to_string()),
        }],
        format: vec!["GT".to_string()],
        samples: vec![vec![genotype.to_string()]],
    }
}

/// Run one case and render both files.
fn case(variants: &[Variant], options: &Options) -> (String, String) {
    let (lifted, rejected) = run(variants, &chains(), &reference(), options);
    let render_all =
        |records: Vec<Variant>| records.iter().map(render).collect::<Vec<_>>().join("\n");
    (render_all(lifted), render_all(rejected))
}

/// A variant inside a block moves by the block's offset; one between the blocks does not move at
/// all, and is written out rather than dropped.
#[test]
fn a_block_renumbers_and_a_gap_rejects() {
    let text = corpus();
    let options = Options::default();

    let (lifted, rejected) = case(&[variant(21, "A", "C", "0/1")], &options);
    assert_eq!(
        lifted,
        field(&text, "lifted", "inside-a-block").expect("the golden")
    );
    assert_eq!(
        rejected,
        field(&text, "rejected", "inside-a-block").expect("the golden")
    );

    let (lifted, rejected) = case(&[variant(150, "C", "A", "0/1")], &options);
    assert_eq!(
        lifted,
        field(&text, "lifted", "between-the-blocks").expect("the golden")
    );
    assert_eq!(
        rejected,
        field(&text, "rejected", "between-the-blocks").expect("the golden")
    );

    // The offset is the chain's own twelve, and the interval keeps its length and its strand.
    let target = lift_over(
        &chains(),
        &Interval {
            contig: "chr1".to_string(),
            start: 21,
            end: 21,
            negative_strand: false,
        },
        1.0,
    )
    .expect("a target");
    assert_eq!(target.start, 33);
    assert!(!target.negative_strand);
}

/// A reversed block complements the alleles and counts the position from the far end.
#[test]
fn a_reversed_block_flips_the_strand() {
    let text = corpus();
    let (lifted, rejected) = case(&[variant(250, "C", "A", "0/1")], &Options::default());
    assert_eq!(
        lifted,
        field(&text, "lifted", "inside-the-reversed-block").expect("the golden")
    );
    assert_eq!(
        rejected,
        field(&text, "rejected", "inside-the-reversed-block").expect("the golden")
    );
    // The alleles are complemented, and the record says so.
    assert!(lifted.contains("\tG\tT\t"));
    assert!(lifted.contains("ReverseComplementedAlleles"));
}

/// A reference allele the target does not carry is rejected, and the rejection says what was
/// attempted.
#[test]
fn a_reference_mismatch_is_rejected_with_its_attempt() {
    let text = corpus();
    let options = Options::default();

    let (lifted, rejected) = case(&[variant(21, "C", "G", "0/1")], &options);
    assert_eq!(
        lifted,
        field(&text, "lifted", "a-reference-mismatch").expect("the golden")
    );
    assert_eq!(
        rejected,
        field(&text, "rejected", "a-reference-mismatch").expect("the golden")
    );
    // The record written out is the SOURCE at its own coordinates; the target's are in the INFO.
    assert!(rejected.starts_with("chr1\t21\t"));
    assert!(rejected.contains("AttemptedLocus=chrT:33-33"));
    assert!(rejected.contains("AttemptedAlleles=C*->G"));

    // The alleles the other way round are recoverable, and are still rejected until they are
    // asked to be recovered.
    let (_, rejected) = case(&[variant(21, "C", "A", "0/1")], &options);
    assert_eq!(
        rejected,
        field(&text, "rejected", "a-swapped-ref-and-alt").expect("the golden")
    );
}

/// Recovering a swap rewrites the alleles, the genotypes and the frequency.
#[test]
fn recovering_a_swap_changes_more_than_the_alleles() {
    let text = corpus();
    let recovering = Options {
        recover_swapped_ref_alt: true,
        ..Options::default()
    };
    let (lifted, rejected) = case(&[variant(21, "C", "A", "0/1")], &recovering);
    assert_eq!(
        lifted,
        field(&text, "lifted", "a-swapped-ref-and-alt-recovered").expect("the golden")
    );
    assert_eq!(
        rejected,
        field(&text, "rejected", "a-swapped-ref-and-alt-recovered").expect("the golden")
    );

    // The genotype follows the alleles rather than being re-sorted, so it is 1/0 and not 0/1.
    assert!(lifted.ends_with("\tGT\t1/0"));
    // The frequency is the complement, and it is a computed number: it reaches the file through
    // the encoder's own formatting, which is why 0.75 is written with three decimal places while
    // the 0.25 it came from was copied through as text.
    assert!(lifted.contains("AF=0.750"));
    assert_eq!(format_vcf_double(0.75), "0.750");
    assert_eq!(format_vcf_double(1.0 - 0.25), "0.750");
}

/// The original position is recorded when asked; the original alleles only when they changed.
#[test]
fn the_original_is_recorded_where_it_is_worth_recording() {
    let text = corpus();
    let inside = [variant(21, "A", "C", "0/1")];

    let (lifted, _) = case(
        &inside,
        &Options {
            write_original_position: true,
            ..Options::default()
        },
    );
    assert_eq!(
        lifted,
        field(&text, "lifted", "with-the-original-position").expect("the golden")
    );
    assert!(lifted.contains("OriginalContig=chr1;OriginalStart=21"));

    // A plain lift leaves the alleles alone, so there is nothing to record and nothing is.
    let (lifted, _) = case(
        &inside,
        &Options {
            write_original_alleles: true,
            ..Options::default()
        },
    );
    assert_eq!(
        lifted,
        field(&text, "lifted", "with-the-original-alleles").expect("the golden")
    );
    assert!(!lifted.contains("OriginalAlleles"));

    // Nor does a recovered swap record them, the record being rebuilt from the lifted alleles
    // rather than from the source's.
    let (lifted, _) = case(
        &[variant(21, "C", "A", "0/1")],
        &Options {
            recover_swapped_ref_alt: true,
            write_original_alleles: true,
            ..Options::default()
        },
    );
    assert_eq!(
        lifted,
        field(&text, "lifted", "with-the-original-alleles-after-a-swap").expect("the golden")
    );
}

/// The output is sorted by the target's coordinates, which is not the order it was read in.
#[test]
fn the_output_is_sorted_by_the_target() {
    let text = corpus();
    let (lifted, rejected) = case(
        &[
            variant(250, "C", "A", "0/1"),
            variant(21, "A", "C", "1/1"),
            variant(150, "C", "A", "0/0"),
        ],
        &Options::default(),
    );
    assert_eq!(
        lifted,
        field(&text, "lifted", "three-variants").expect("the golden")
    );
    assert_eq!(
        rejected,
        field(&text, "rejected", "three-variants").expect("the golden")
    );
    // The variant read first comes out last: 250 lifts to 251, and 21 lifts to 33.
    let positions: Vec<&str> = lifted
        .lines()
        .map(|line| line.split('\t').nth(1).expect("a position"))
        .collect();
    assert_eq!(positions, vec!["33", "251"]);
}

/// A variant that lands nowhere is filtered rather than lost.
#[test]
fn a_rejection_keeps_the_record() {
    let source = variant(150, "C", "A", "0/1");
    let Outcome::Rejected(rejected) =
        liftover(&source, &chains(), &reference(), &Options::default())
    else {
        panic!("a rejection")
    };
    assert_eq!(rejected.contig, source.contig);
    assert_eq!(rejected.position, source.position);
    assert_eq!(rejected.filters, vec!["NoTarget".to_string()]);
}
