//! Conformance for `GtcToVcf` against Picard 3.4.0.
//!
//! Golden from `tools/arrays-conformance/`: six runs over four loci on two contigs, with the calls
//! and the manifest's flags moved between them.
//!
//! # What this suite is for
//!
//!  * **a call becoming a genotype against the BUILD's alleles**, not the chip's `A` and `B`;
//!  * **a reference base that is neither chip allele making the record triallelic**;
//!  * **a no-call being written as `./.` rather than left out**;
//!  * **a flagged locus being dropped and a duplicate being kept and filtered**, which is not the
//!    same treatment;
//!  * **the record carrying the chip's own numbers**, normalized intensities and cluster included;
//!  * **and the output being sorted by the target build's coordinates.**

use std::io::Read;

use picard_analysis::gtc_to_vcf::{
    assay_alleles, format_float_for_vcf, genotype_field, header_lines, normalize,
    polar_to_euclidean, r_and_theta, records, variant_line, Call, Cluster, Flag, ManifestRecord,
    Transformation, AA_CALL, AB_CALL, BB_CALL, NO_CALL,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gtc_to_vcf.txt.gz");
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

/// The reference base at a position, the contigs being `ACGT` repeating.
fn base(position: i32) -> String {
    ["A", "C", "G", "T"][((position - 1) % 4) as usize].to_string()
}

/// The fixture's four loci, as the extended manifest holds them.
fn manifest() -> Vec<ManifestRecord> {
    [
        ("rs1", "1", 1001, "A", "G", ""),
        ("rs2", "1", 2001, "T", "C", "ACGTACGA"),
        ("rs3", "2", 3001, "A", "C", ""),
        ("rs4", "2", 4001, "A", "T", "ACGTACGA"),
    ]
    .iter()
    .map(
        |(name, chr, position, allele_a, allele_b, probe_b)| ManifestRecord {
            name: name.to_string(),
            chr: chr.to_string(),
            position: *position,
            genome_build: "37".to_string(),
            b37_chr: chr.to_string(),
            b37_pos: *position,
            ref_allele: base(*position),
            allele_a: allele_a.to_string(),
            allele_b: allele_b.to_string(),
            rs_id: name.to_string(),
            ilmn_strand: "TOP".to_string(),
            probe_a: "ACGTACGT".to_string(),
            probe_b: if probe_b.is_empty() {
                ".".to_string()
            } else {
                probe_b.to_string()
            },
            bead_set_id: 1,
            source: "source".to_string(),
            flag: Flag::Pass,
        },
    )
    .collect()
}

/// The cluster file's four entries.
fn clusters() -> Vec<Cluster> {
    (0..4)
        .map(|index| Cluster {
            total_score: 0.5 + index as f32 / 100.0,
            n: [10 + index, 20 + index, 30 + index],
            dev_r: [0.1, 0.2, 0.3],
            mean_r: [1.0, 1.1, 1.2],
            dev_theta: [0.01, 0.02, 0.03],
            mean_theta: [0.2, 0.5, 0.8],
        })
        .collect()
}

/// The call file's four records, called `AA`, `AB`, `BB` and not at all.
fn calls(genotypes: [u8; 4], scores: [f32; 4]) -> Vec<Call> {
    (0..4)
        .map(|index| Call {
            genotype: genotypes[index],
            score: scores[index],
            raw_x: 1000 * (index as i32 + 1),
            raw_y: 1000 * (index as i32 + 1) + 100,
            b_allele_freq: match genotypes[index] {
                AA_CALL => 0.0,
                AB_CALL => 0.5,
                BB_CALL => 1.0,
                _ => f32::NAN,
            },
            log_r_ratio: 0.1 * index as f32,
        })
        .collect()
}

/// The one transformation every locus goes through.
fn transformations() -> Vec<Transformation> {
    vec![
        Transformation {
            offset_x: 10.0,
            offset_y: 20.0,
            scale_x: 1.0,
            scale_y: 1.0,
            shear: 0.0,
            theta: 0.0,
        };
        4
    ]
}

fn run(manifest: &[ManifestRecord], calls: &[Call]) -> String {
    records(manifest, calls, &clusters(), &transformations()).join("\n")
}

/// A call becomes a genotype against the build's alleles, and a reference base that is neither of
/// the chip's makes the record triallelic.
#[test]
fn a_call_becomes_a_genotype_against_the_build() {
    let text = corpus();
    let called = calls([AA_CALL, AB_CALL, BB_CALL, NO_CALL], [0.7, 0.8, 0.9, 0.0]);
    assert_eq!(
        run(&manifest(), &called),
        field(&text, "vcf", "four-loci").expect("the golden")
    );

    // The second locus's reference base is neither of its alleles, so both are alternates.
    let records = manifest();
    assert_eq!(assay_alleles(&records[1]), vec!["A", "T", "C"]);
    assert_eq!(genotype_field(&records[1], &called[1]), "1/2");
    // Where the reference IS one of them, the chip's `AA` is the reference genotype.
    assert_eq!(assay_alleles(&records[0]), vec!["A", "G"]);
    assert_eq!(genotype_field(&records[0], &called[0]), "0/0");
    // And a `BB` call on such a locus is homozygous alternate.
    assert_eq!(genotype_field(&records[2], &called[2]), "1/1");
}

/// The same call at every locus is the manifest's answer rather than the chip's.
#[test]
fn the_same_call_is_not_the_same_genotype() {
    let text = corpus();
    let all_aa = calls([AA_CALL; 4], [0.9; 4]);
    assert_eq!(
        run(&manifest(), &all_aa),
        field(&text, "vcf", "all-homozygous-reference").expect("the golden")
    );
    // Three of the four come out as `0/0` and the triallelic one as `1/1`, from the same call.
    let records = manifest();
    assert_eq!(genotype_field(&records[0], &all_aa[0]), "0/0");
    assert_eq!(genotype_field(&records[1], &all_aa[1]), "1/1");
}

/// A no-call is written out rather than left out.
#[test]
fn a_no_call_is_written() {
    let text = corpus();
    let none = calls([NO_CALL; 4], [0.0; 4]);
    let written = run(&manifest(), &none);
    assert_eq!(
        written,
        field(&text, "vcf", "all-no-calls").expect("the golden")
    );
    assert_eq!(written.lines().count(), 4);
    // The genotype is `./.`, the B allele frequency is a dot, and the allele number is nought.
    assert!(written.lines().all(|line| line.contains("\t./.:.:")));
    assert!(written.lines().all(|line| line.contains(";AN=0;")));
}

/// A flagged locus is dropped; a duplicate is kept and filtered.
#[test]
fn a_flag_and_a_duplicate_are_not_the_same() {
    let text = corpus();
    let called = calls([AA_CALL, AB_CALL, BB_CALL, NO_CALL], [0.7, 0.8, 0.9, 0.0]);

    let mut flagged = manifest();
    flagged[1].flag = Flag::IlluminaFlagged;
    let written = run(&flagged, &called);
    assert_eq!(
        written,
        field(&text, "vcf", "a-flagged-locus").expect("the golden")
    );
    // Three records, and the flagged one is nowhere in the file, not even filtered.
    assert_eq!(written.lines().count(), 3);
    assert!(!written.contains("rs2"));

    let mut duped = manifest();
    duped[2].flag = Flag::Dupe;
    let written = run(&duped, &called);
    assert_eq!(
        written,
        field(&text, "vcf", "a-duplicate-locus").expect("the golden")
    );
    // Four records, and the duplicate is one of them, carrying the filter that says so.
    assert_eq!(written.lines().count(), 4);
    assert!(written
        .lines()
        .nth(2)
        .expect("the row")
        .contains("\tDUPE\t"));
    assert!(Flag::Dupe.is_dupe() && !Flag::Dupe.is_fail());
    assert!(Flag::IlluminaFlagged.is_fail());
}

/// The record carries the chip's own numbers, which are not the ones it measured.
#[test]
fn the_record_carries_the_chips_numbers() {
    let transformation = transformations()[0];
    // The offsets come off the raw intensities before anything else.
    assert_eq!(normalize(1000, 1100, &transformation), (990.0, 1080.0));
    // The total is the sum of the two, the distance being Manhattan, and the angle is the split
    // between them.
    let (r, theta) = r_and_theta(990.0, 1080.0);
    assert_eq!(r, 2070.0);
    assert_eq!(format_float_for_vcf(theta), "0.528");

    // The cluster's polar description becomes a position in the plane, deviations and all.
    let euclidean = polar_to_euclidean(1.0, 0.1, 0.2, 0.01);
    assert_eq!(format_float_for_vcf(euclidean.mean_x), "0.755");
    assert_eq!(format_float_for_vcf(euclidean.mean_y), "0.245");
    assert_eq!(format_float_for_vcf(euclidean.dev_x), "0.076");
    assert_eq!(format_float_for_vcf(euclidean.dev_y), "0.026");

    // Three decimal places at most, no grouping separator, and a dot for what is not a number.
    assert_eq!(format_float_for_vcf(1990.0), "1990");
    assert_eq!(format_float_for_vcf(f32::NAN), ".");
}

/// The output is sorted by the target build's coordinates.
#[test]
fn the_output_is_sorted_by_the_build() {
    let text = corpus();
    let called = calls([AA_CALL, AB_CALL, BB_CALL, NO_CALL], [0.7, 0.8, 0.9, 0.0]);
    // The chip lists its loci in whatever order it likes; the file comes out in the build's.
    let mut shuffled = manifest();
    shuffled.reverse();
    let mut shuffled_calls = called.clone();
    shuffled_calls.reverse();
    let mut shuffled_clusters = clusters();
    shuffled_clusters.reverse();
    let written = records(
        &shuffled,
        &shuffled_calls,
        &shuffled_clusters,
        &transformations(),
    )
    .join("\n");
    assert_eq!(
        written,
        field(&text, "vcf", "four-loci").expect("the golden")
    );
}

/// The header names the run: the sample, the pipeline and the analysis version.
#[test]
fn the_header_carries_the_runs_identity() {
    let text = corpus();
    assert_eq!(
        header_lines(Some(1), Some("1.0"), "sample1", "sample").join("\n"),
        field(&text, "header", "four-loci").expect("the golden")
    );
    // The gender the command line declares changes nothing in the records themselves.
    assert_eq!(
        field(&text, "vcf", "a-male-sample"),
        field(&text, "vcf", "four-loci")
    );
}

/// One record, built from end to end.
#[test]
fn one_record_is_the_golden_row() {
    let text = corpus();
    let called = calls([AA_CALL, AB_CALL, BB_CALL, NO_CALL], [0.7, 0.8, 0.9, 0.0]);
    let line = variant_line(
        &manifest()[0],
        &called[0],
        &clusters()[0],
        &transformations()[0],
    );
    assert_eq!(
        line,
        field(&text, "vcf", "four-loci")
            .expect("the golden")
            .lines()
            .next()
            .expect("a row")
    );
}
