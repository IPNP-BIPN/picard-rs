//! Conformance for `MergePedIntoVcf` against Picard 3.4.0.
//!
//! Each case carries the exit code, and either the merged VCF or the note that none was written.
//! The golden prints the three input files too, so the port is driven by the reference's own copies.
//!
//! # What this suite is for
//!
//!  * **every failure being swallowed, so the exit code is nought either way**;
//!  * **the thresholds map being static, so an earlier run's thresholds leak into a later one**;
//!  * **a VCF of nothing but `GT` not being processable at all**;
//!  * **the PED's alleles being `A`, `B` or `0` and looked up per record**;
//!  * **the looked-up allele being non-reference, so naming the REF fails**;
//!  * **a pair of `NA` becoming the missing value and one `NA` being refused**;
//!  * **the first six PED fields being ignored**;
//!  * **and both multi-sample refusals.**

use std::io::Read;

use picard_analysis::merge_ped_into_vcf::{
    genotype_string, parse_ped, parse_thresholds, translate_allele, zcall_alleles, MergeError,
    HALF_NA_MESSAGE, LONG_ALLELE_MESSAGE, MISSING_VALUE, MULTI_SAMPLE_PED_MESSAGE,
    MULTI_SAMPLE_VCF_MESSAGE, PED_OFFSET,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("merge_ped_into_vcf.txt.gz");
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

/// One case's records, as (id, info, format, sample).
fn records(text: &str, case: &str) -> Vec<(String, String, String, String)> {
    match field(text, "out", case) {
        None => Vec::new(),
        Some(vcf) => vcf
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(|line| {
                let c: Vec<&str> = line.split('\t').collect();
                (
                    c[2].to_string(),
                    c[7].to_string(),
                    c[8].to_string(),
                    c[9].to_string(),
                )
            })
            .collect(),
    }
}

fn code(text: &str, case: &str) -> i32 {
    field(text, "code", case)
        .unwrap_or_else(|| panic!("{case} has a code"))
        .trim()
        .parse()
        .expect("a code")
}

/// The value of one sub-field of a record's genotype column.
fn genotype_field(format: &str, sample: &str, name: &str) -> String {
    let index = format
        .split(':')
        .position(|f| f == name)
        .unwrap_or_else(|| panic!("{name} is in {format}"));
    sample.split(':').nth(index).expect("a value").to_string()
}

/// Failures BEFORE the writer opens, which leave no file at all.
const NO_FILE: &[&str] = &["half-na", "long-allele", "two-vcf-samples"];

/// Failures INSIDE the record loop, which leave a file holding its header and no records: the
/// writer had already opened it and written the header when the exception was thrown.
const HEADER_ONLY: &[&str] = &["a-allele-is-the-reference", "no-extended-attributes"];

/// Every failure is swallowed: the exit code is nought whichever kind it is.
#[test]
fn every_failure_is_swallowed() {
    let text = corpus();
    for case in NO_FILE {
        assert_eq!(code(&text, case), 0, "{case}");
        assert!(field(&text, "none", case).is_some(), "{case}");
        assert!(field(&text, "out", case).is_none(), "{case}");
    }
    for case in HEADER_ONLY {
        assert_eq!(code(&text, case), 0, "{case}");
        // A file, and a header in it, but not one record.
        assert!(field(&text, "out", case).is_some(), "{case}");
        assert!(records(&text, case).is_empty(), "{case}");
    }
    for case in ["agreeing", "disagreeing", "static-map-leak"] {
        assert_eq!(code(&text, case), 0, "{case}");
        assert_eq!(records(&text, case).len(), 3, "{case}");
    }
}

/// The original genotype is kept as GTA and the zCall one added as GTZ.
#[test]
fn both_genotypes_are_kept() {
    let text = corpus();
    let rows = records(&text, "agreeing");
    for (id, _, format, _) in &rows {
        assert!(format.contains("GTA"), "{id}");
        assert!(format.contains("GTZ"), "{id}");
    }
    // rs1's VCF genotype is 0/0 and its PED calls B B, so the two disagree.
    let rs1 = rows.iter().find(|r| r.0 == "rs1").expect("rs1");
    assert_eq!(genotype_field(&rs1.2, &rs1.3, "GTA"), "0/0");
    assert_eq!(genotype_field(&rs1.2, &rs1.3, "GTZ"), "1/1");
    // And the record's own GT is now the zCall one.
    assert_eq!(genotype_field(&rs1.2, &rs1.3, "GT"), "1/1");
}

/// A `0` in the PED is a no-call, which zeroes the record's AC and AN.
#[test]
fn a_zero_allele_is_a_no_call() {
    let text = corpus();
    let rows = records(&text, "disagreeing");
    let rs2 = rows.iter().find(|r| r.0 == "rs2").expect("rs2");
    assert_eq!(genotype_field(&rs2.2, &rs2.3, "GTZ"), "./.");
    assert_eq!(genotype_field(&rs2.2, &rs2.3, "GTA"), "0/1");
    assert!(rs2.1.contains("AC=0"), "{}", rs2.1);
    assert!(rs2.1.contains("AN=0"), "{}", rs2.1);
    // Which the port reaches: `0` translates to nothing at all.
    assert_eq!(translate_allele('0', "G", "T"), Ok(None));
    assert_eq!(translate_allele('A', "G", "T"), Ok(Some("G".to_string())));
    assert_eq!(translate_allele('B', "G", "T"), Ok(Some("T".to_string())));
    assert_eq!(genotype_string(&[None, None], "T"), "./.".to_string());
    assert_eq!(
        genotype_string(&[Some("T".to_string()), Some("T".to_string())], "T"),
        "1/1".to_string()
    );
}

/// The looked-up allele is non-reference, so a PED naming the record's REF fails.
#[test]
fn naming_the_reference_fails() {
    let text = corpus();
    assert!(records(&text, "a-allele-is-the-reference").is_empty());
    // rs1 is REF=A ALT=C with ALLELE_A=A, so a PED of `A A` names the reference.
    assert_eq!(
        zcall_alleles("AA", "A", "C", "A"),
        Err(MergeError::AlleleNotInContext("A".to_string()))
    );
    // Where `B B` names the alternate and is fine.
    assert_eq!(
        zcall_alleles("BB", "A", "C", "A"),
        Ok([Some("C".to_string()), Some("C".to_string())])
    );
}

/// A pair of NA becomes the missing value; one NA of a pair is refused.
#[test]
fn the_na_pair_is_the_missing_value() {
    let text = corpus();
    let rows = records(&text, "na-thresholds");
    let rs2 = rows.iter().find(|r| r.0 == "rs2").expect("rs2");
    assert!(rs2.1.contains("zthresh_X=."), "{}", rs2.1);
    assert_eq!(
        parse_thresholds("rs2\tNA\tNA\n"),
        Ok([(
            "rs2".to_string(),
            [MISSING_VALUE.to_string(), MISSING_VALUE.to_string()]
        )]
        .into_iter()
        .collect())
    );
    assert_eq!(parse_thresholds("rs1\tNA\t0.6\n"), Err(MergeError::HalfNa));
    assert_eq!(MergeError::HalfNa.message(), HALF_NA_MESSAGE);
    // Which the tool swallows.
    assert!(records(&text, "half-na").is_empty());
}

/// The thresholds map is static, so an earlier run's thresholds leak into a later one.
#[test]
fn the_thresholds_map_leaks_between_runs() {
    let text = corpus();
    // Its own file names rs3 alone, and all three records carry thresholds.
    let rows = records(&text, "static-map-leak");
    assert_eq!(rows.len(), 3);
    for (id, info, _, _) in &rows {
        assert!(info.contains("zthresh_X="), "{id}: {info}");
    }
    // The first run's file also named rs1 alone, and only rs1 carried them.
    let first = records(&text, "agreeing");
    let carried: Vec<&String> = first
        .iter()
        .filter(|r| r.1.contains("zthresh_X="))
        .map(|r| &r.0)
        .collect();
    assert_eq!(carried, vec!["rs1"]);
}

/// A VCF of nothing but GT cannot be processed at all.
#[test]
fn a_vcf_of_only_gt_is_not_processable() {
    let text = corpus();
    assert!(records(&text, "no-extended-attributes").is_empty());
    assert_eq!(code(&text, "no-extended-attributes"), 0);
    assert_eq!(
        MergeError::NoExtendedAttributes.message(),
        "java.lang.UnsupportedOperationException"
    );
}

/// The PED and MAP are read in step, the first six PED fields ignored.
#[test]
fn the_ped_and_map_are_read_in_step() {
    let text = corpus();
    let ped = field(&text, "in", "ped").expect("the ped");
    let map = field(&text, "in", "map").expect("the map");
    let alleles = parse_ped(&ped, &map).expect("parsed");
    assert_eq!(alleles.len(), 3);
    assert_eq!(alleles["rs1"], "BB");
    assert_eq!(PED_OFFSET, 6);
    // A PED of two samples, and an allele of two characters.
    assert_eq!(
        parse_ped(
            "FAM\tA\t0\t0\t0\t-9\tB\tB\nFAM\tB\t0\t0\t0\t-9\tB\tB\n",
            "chr1\trs1\t0\t1\n"
        ),
        Err(MergeError::MultiSamplePed)
    );
    assert_eq!(
        parse_ped("FAM\tIND\t0\t0\t0\t-9\tBB\tB\n", "chr1\trs1\t0\t1\n"),
        Err(MergeError::LongAllele)
    );
    assert_eq!(
        MergeError::MultiSamplePed.message(),
        MULTI_SAMPLE_PED_MESSAGE
    );
    assert_eq!(MergeError::LongAllele.message(), LONG_ALLELE_MESSAGE);
    assert_eq!(
        MergeError::MultiSampleVcf.message(),
        MULTI_SAMPLE_VCF_MESSAGE
    );
    // Both of which the tool swallows.
    assert!(records(&text, "long-allele").is_empty());
    assert!(records(&text, "two-vcf-samples").is_empty());
}
