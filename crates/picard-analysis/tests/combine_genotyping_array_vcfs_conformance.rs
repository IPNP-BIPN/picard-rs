//! Conformance for `CombineGenotypingArrayVcfs` against Picard 3.4.0.
//!
//! Golden from `tools/combinearrays-conformance/CombineGenotypingArrayVcfsDump.java`, nineteen
//! merges of single-sample array VCFs.
//!
//! # What this suite is for
//!
//!  * **the lockstep, which refuses what a merge by position would have reordered**;
//!  * **the twelve refusals, each with its own message and its own exception class**;
//!  * **the seven attributes that may differ and the one that may not**;
//!  * **the depth, which is the one attribute the merge adds up and the one that makes the tool
//!    throw**;
//!  * **and the merged header, which drops what belonged to one sample.**

use std::io::Read;

use picard_analysis::combine_genotyping_array_vcfs::{
    attribute_must_agree, check_attributes, check_step, header_line_is_kept, sample_list, Refusal,
    Site, DEPTH_KEY, EXEMPT_ATTRIBUTES, SAMPLE_SPECIFIC_HEADER_LINES,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/combine_genotyping_array_vcfs.txt.gz");
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

fn refusal(text: &str, case: &str) -> String {
    field(text, "error", case).unwrap_or_else(|| panic!("error/{case}"))
}

/// A refusal as the golden wrote it: the class, the message, and the frame it came from.
fn carries(written: &str, refusal: &Refusal) -> bool {
    written.starts_with(refusal.class()) && written.contains(&refusal.message())
}

fn site(position: i64, id: &str, reference: &str, alternates: &[&str]) -> Site {
    Site {
        contig: "chr1".to_string(),
        start: position,
        id: id.to_string(),
        reference: reference.to_string(),
        alternates: alternates.iter().map(|value| value.to_string()).collect(),
        attributes: vec![("BEADSET".to_string(), "7".to_string())],
    }
}

/// The inputs are walked in lockstep, so the same loci in another order are a refusal.
#[test]
fn the_inputs_are_walked_in_lockstep() {
    let text = corpus();
    // A merge by position would have reordered these; this one refuses them at the first step.
    let written = refusal(&text, "loci-in-another-order");
    assert!(carries(&written, &Refusal::Locus), "{written}");
    assert_eq!(
        check_step(
            &site(100, "rs100", "A", &["C"]),
            &site(200, "rs200", "A", &["C"])
        ),
        Err(Refusal::Locus)
    );
    // A different NUMBER of variants is found out only when an iterator runs dry, so the message
    // names no locus.
    let count = refusal(&text, "a-different-number-of-variants");
    assert!(carries(&count, &Refusal::VariantCount), "{count}");
    assert!(!count.contains("chr1"));
    // And two files that do line up merge into one file with both samples, in the input's order.
    let merged = field(&text, "merged", "two-samples").expect("the merge");
    assert_eq!(merged.lines().count(), 2);
    let header = field(&text, "header", "two-samples").expect("the header");
    assert!(header.ends_with("sampleA\tsampleB"), "{header}");
    let three = field(&text, "header", "three-samples").expect("the header");
    assert!(three.ends_with("sampleC\tsampleA\tsampleB"), "{three}");
    assert_eq!(
        sample_list(&[
            ("in0.vcf".to_string(), vec!["sampleC".to_string()]),
            ("in1.vcf".to_string(), vec!["sampleA".to_string()]),
        ]),
        Ok(vec!["sampleC".to_string(), "sampleA".to_string()])
    );
}

/// Each way two variants can fail to line up has its own message.
#[test]
fn every_mismatch_has_its_own_message() {
    let text = corpus();
    let first = site(100, "rs100", "A", &["C"]);
    let cases = [
        (
            "a-different-id",
            site(100, "other", "A", &["C"]),
            Refusal::Id,
        ),
        (
            "a-different-reference-allele",
            site(100, "rs100", "T", &["C"]),
            Refusal::ReferenceAllele,
        ),
        (
            "a-different-alternate-count",
            site(100, "rs100", "A", &["C", "G"]),
            Refusal::AlternateAlleleCount,
        ),
        (
            "a-different-alternate-allele",
            site(100, "rs100", "A", &["G"]),
            Refusal::AlternateAllele {
                contig: "chr1".to_string(),
                start: 100,
            },
        ),
    ];
    for (case, other, expected) in cases {
        let written = refusal(&text, case);
        assert!(carries(&written, &expected), "{case}: {written}");
        assert_eq!(check_step(&first, &other), Err(expected), "{case}");
    }
    // The alternate allele is the one refusal that says WHERE, and it says it with a dot.
    assert!(refusal(&text, "a-different-alternate-allele").contains("for chr1.100"));
    // A repeated sample is refused with the file that repeated it, before any variant is read.
    let repeated = refusal(&text, "a-repeated-sample");
    assert!(
        repeated.contains("contains a sample entry (sampleA)"),
        "{repeated}"
    );
    assert!(
        repeated.starts_with("java.lang.IllegalArgumentException"),
        "{repeated}"
    );
}

/// Seven attributes may differ, the rest may not, and one may not even be present.
#[test]
fn seven_attributes_are_exempt() {
    let text = corpus();
    for key in EXEMPT_ATTRIBUTES {
        assert!(!attribute_must_agree(key), "{key}");
    }
    assert!(attribute_must_agree("BEADSET"));
    // An exempt attribute that disagrees is merged rather than refused, and the merge recalculates
    // it, so the output carries a value neither input had.
    assert_eq!(
        field(&text, "error", "an-exempt-attribute-that-disagrees"),
        None
    );
    let merged = field(&text, "merged", "an-exempt-attribute-that-disagrees").expect("the merge");
    assert!(merged.contains("AC=1;AF=0.250;AN=4"), "{merged}");
    // A key that is not exempt is refused only in ONE direction. The loop runs over the other
    // files' attributes and looks each up in the first file's, so a key only the first file has
    // is kept without comment and reaches the output, and one a later file has is refused.
    assert_eq!(
        field(&text, "error", "an-attribute-in-the-first-file-only"),
        None
    );
    let kept = field(&text, "merged", "an-attribute-in-the-first-file-only").expect("the merge");
    assert!(kept.contains("EXTRA=1"), "{kept}");
    let missing = refusal(&text, "an-attribute-in-a-later-file-only");
    assert!(
        carries(&missing, &Refusal::AttributeMissing("EXTRA".to_string())),
        "{missing}"
    );
    let disagrees = refusal(&text, "an-attribute-that-disagrees");
    assert!(
        carries(
            &disagrees,
            &Refusal::AttributeDisagrees("BEADSET".to_string())
        ),
        "{disagrees}"
    );
    assert_eq!(
        check_attributes(
            &site(100, "rs100", "A", &["C"]),
            &Site {
                attributes: vec![("BEADSET".to_string(), "8".to_string())],
                ..site(100, "rs100", "A", &["C"])
            }
        ),
        Err(Refusal::AttributeDisagrees("BEADSET".to_string()))
    );
}

/// The depth is the one attribute the merge adds up, and the one that makes the tool throw.
#[test]
fn a_depth_is_never_merged() {
    let text = corpus();
    // Two equal depths reach the write-back, which puts the sum into an unmodifiable map.
    let both = refusal(&text, "a-depth-in-both-files");
    assert!(
        both.starts_with("java.lang.UnsupportedOperationException"),
        "{both}"
    );
    assert!(both.contains("Collections$UnmodifiableMap.put"), "{both}");
    assert_eq!(field(&text, "merged", "a-depth-in-both-files"), None);
    // Two that differ are refused earlier, by the agreement check, since DP is not exempt.
    let differ = refusal(&text, "depths-that-differ");
    assert!(
        carries(&differ, &Refusal::AttributeDisagrees(DEPTH_KEY.to_string())),
        "{differ}"
    );
    assert!(attribute_must_agree(DEPTH_KEY));
    // And a depth in the first file only is refused for the attribute the LATER file has and the
    // first does not, which is the ordinary check running in its own direction rather than the
    // depth's.
    let one = refusal(&text, "a-depth-in-one-file-only");
    assert!(
        carries(&one, &Refusal::AttributeMissing("BEADSET".to_string())),
        "{one}"
    );
    let with_depth = Site {
        attributes: vec![
            ("BEADSET".to_string(), "7".to_string()),
            (DEPTH_KEY.to_string(), "10".to_string()),
        ],
        ..site(100, "rs100", "A", &["C"])
    };
    assert_eq!(
        check_attributes(&with_depth, &with_depth),
        Err(Refusal::DepthIsUnwritable)
    );
    assert_eq!(Refusal::DepthIsUnwritable.message(), "");
}

/// The merged header drops what belonged to one sample, and keeps what belonged to the chip.
#[test]
fn the_merged_header_drops_the_sample_specific_lines() {
    let text = corpus();
    let header = field(&text, "header", "two-samples").expect("the header");
    for line in SAMPLE_SPECIFIC_HEADER_LINES {
        assert!(!header.contains(&format!("##{line}=")), "{line}");
        assert!(!header_line_is_kept(line), "{line}");
    }
    // The chip's own line is not one sample's, so it stays.
    assert!(header.contains("##arrayType=TestArray-24v1-0_A1"));
    assert!(header_line_is_kept("arrayType"));
    // The filters of every input are unioned onto the merged variant.
    let filtered = field(&text, "merged", "filters-that-differ").expect("the merge");
    assert!(filtered.contains("\tLOW\t"), "{filtered}");
    // The index is written whether or not `CREATE_INDEX` asked for it: the writer's builder
    // indexes any file it has a sequence dictionary for, and the argument only adds an option the
    // builder already had.
    assert_eq!(
        field(&text, "files", "an-index").as_deref(),
        Some("merged.vcf merged.vcf.idx")
    );
    assert_eq!(
        field(&text, "files", "two-samples").as_deref(),
        Some("merged.vcf merged.vcf.idx")
    );
}
