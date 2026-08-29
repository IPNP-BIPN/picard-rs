//! Conformance for `VcfToAdpc` against Picard 3.4.0.
//!
//! Golden from `tools/vcftoadpc-conformance/VcfToAdpcDump.java`, fourteen runs whose binary output
//! is in the golden as hex.
//!
//! # What this suite is for
//!
//!  * **the eighteen bytes of a record, in the writer's own order**;
//!  * **the genotype being the ARRAY's and not the VCF's**;
//!  * **an intensity over an unsigned short being truncated and a negative one refused**;
//!  * **a missing normalized intensity being a NaN rather than a shorter record**;
//!  * **and the two text files beside the binary one.**

use std::io::Read;

use picard_analysis::vcf_to_adpc::{
    illumina_genotype, markers_file, raw_intensity, samples_file, write_file, IlluminaGenotype,
    Record, HEADER, MAX_UNSIGNED_SHORT, REFUSAL_EXIT_CODE,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_to_adpc.txt.gz");
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

/// One case's binary output, as bytes.
fn output(text: &str, case: &str) -> Vec<u8> {
    let hex = field(text, "adpc", case).unwrap_or_else(|| panic!("adpc/{case}"));
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).expect("a byte"))
        .collect()
}

/// The fixture's own record, with the genotype the case produced.
fn fixture(genotype: IlluminaGenotype) -> Record {
    Record {
        a_intensity: 1000,
        b_intensity: 2000,
        a_normalized: 0.5,
        b_normalized: 1.5,
        gc_score: 0.75,
        genotype,
    }
}

/// The header and the twenty bytes of a record are the golden's.
#[test]
fn a_record_is_eighteen_bytes() {
    let text = corpus();
    let written = output(&text, "homozygous-a");
    assert_eq!(&written[..16], HEADER);
    assert_eq!(written.len(), 16 + 18);
    assert_eq!(written, write_file(&[fixture(IlluminaGenotype::Aa)]));
    // Two loci are two records in a row, and so are two samples: the walk is sample-major and the
    // file carries no separator between them.
    let two_loci = output(&text, "two-loci");
    assert_eq!(
        two_loci,
        write_file(&[fixture(IlluminaGenotype::Aa), fixture(IlluminaGenotype::Ab)])
    );
    let two_samples = output(&text, "two-samples");
    assert_eq!(
        two_samples,
        write_file(&[fixture(IlluminaGenotype::Aa), fixture(IlluminaGenotype::Aa)])
    );
    // And two VCFs are appended in the order they were given.
    assert_eq!(
        output(&text, "two-vcfs"),
        write_file(&[fixture(IlluminaGenotype::Aa), fixture(IlluminaGenotype::Bb)])
    );
}

/// The genotype is the array's, not the VCF's.
#[test]
fn the_genotype_is_the_arrays() {
    let text = corpus();
    assert_eq!(
        illumina_genotype(Some(("A", "A")), "A", "C"),
        Some(IlluminaGenotype::Aa)
    );
    assert_eq!(
        illumina_genotype(Some(("A", "C")), "A", "C"),
        Some(IlluminaGenotype::Ab)
    );
    assert_eq!(
        illumina_genotype(Some(("C", "C")), "A", "C"),
        Some(IlluminaGenotype::Bb)
    );
    assert_eq!(
        illumina_genotype(None, "A", "C"),
        Some(IlluminaGenotype::Nn)
    );
    for (case, genotype) in [
        ("homozygous-a", IlluminaGenotype::Aa),
        ("heterozygous", IlluminaGenotype::Ab),
        ("homozygous-b", IlluminaGenotype::Bb),
        ("a-no-call", IlluminaGenotype::Nn),
    ] {
        assert_eq!(
            output(&text, case),
            write_file(&[fixture(genotype)]),
            "{case}"
        );
    }
    // The same `0/0` with the array's alleles the other way round is a BB and not an AA.
    assert_eq!(
        illumina_genotype(Some(("A", "A")), "C", "A"),
        Some(IlluminaGenotype::Bb)
    );
    assert_eq!(
        output(&text, "the-alleles-reversed"),
        write_file(&[fixture(IlluminaGenotype::Bb)])
    );
    // And an allele that is the reference carries a trailing star, which is stripped first.
    assert_eq!(
        illumina_genotype(Some(("A", "A")), "A*", "C"),
        Some(IlluminaGenotype::Aa)
    );
    assert_eq!(
        output(&text, "a-reference-allele-with-a-star"),
        write_file(&[fixture(IlluminaGenotype::Aa)])
    );
}

/// An intensity over an unsigned short is truncated, and a missing one is a NaN.
#[test]
fn an_intensity_is_truncated_and_never_widened() {
    let text = corpus();
    assert_eq!(raw_intensity(65535), Some(65535));
    assert_eq!(raw_intensity(70000), Some(65535));
    assert_eq!(raw_intensity(-1), None);
    assert_eq!(MAX_UNSIGNED_SHORT, 65535);
    let at_the_limit = Record {
        a_intensity: 65535,
        b_intensity: 0,
        genotype: IlluminaGenotype::Ab,
        ..fixture(IlluminaGenotype::Ab)
    };
    assert_eq!(
        output(&text, "an-intensity-at-the-limit"),
        write_file(&[at_the_limit])
    );
    // Over the limit gives the SAME bytes, which is what says it was truncated rather than
    // refused or wrapped.
    assert_eq!(
        output(&text, "an-intensity-over-the-limit"),
        output(&text, "an-intensity-at-the-limit")
    );
    // A missing normalized intensity is a NaN and the record keeps its width.
    let without = Record {
        a_normalized: f32::NAN,
        b_normalized: f32::NAN,
        ..fixture(IlluminaGenotype::Ab)
    };
    let written = output(&text, "without-the-normalized-intensities");
    assert_eq!(written.len(), 16 + 18);
    assert_eq!(written, write_file(&[without]));
}

/// The two text files beside the binary one, and what a refusal leaves.
#[test]
fn the_text_files_carry_the_samples_and_the_count() {
    let text = corpus();
    assert_eq!(
        field(&text, "text", "two-samples/samples.txt").as_deref(),
        Some("sample1\nsample2")
    );
    assert_eq!(samples_file(&["sample1", "sample2"]), "sample1\nsample2");
    // No trailing newline, which is what makes the file a list rather than a set of lines.
    assert!(!samples_file(&["sample1"]).ends_with('\n'));
    assert_eq!(
        field(&text, "text", "two-loci/markers.txt").as_deref(),
        Some("2")
    );
    assert_eq!(markers_file(2), "2");
    // A refusal leaves an exit code and nothing else: the tool logs its own exception, and the
    // log reaches no stream the dump could capture.
    for case in ["two-vcfs-of-different-lengths", "no-records"] {
        let refusal = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        assert!(
            refusal.starts_with(&format!("exit {REFUSAL_EXIT_CODE}")),
            "{refusal}"
        );
        assert_eq!(refusal.trim(), format!("exit {REFUSAL_EXIT_CODE}"));
        assert_eq!(field(&text, "adpc", case), None, "{case}");
    }
}
