//! Conformance for the header `MergeBamAlignment` builds, against Picard 3.4.0.
//!
//! Golden from `tools/mergebamalignment-conformance/`: twelve runs whose records are measured by
//! the suites next door, compared here for the header they arrive under.
//!
//! # What this suite is for
//!
//!  * **the sequences coming from the dictionary**, `M5` and canonical `UR` and all;
//!  * **the read groups coming from the unmapped bam**, whatever the aligner wrote;
//!  * **a program record being adopted only when the aligned header holds exactly one**, so two
//!    chained programs leave the output with none;
//!  * **a program record from the command line replacing the aligner's** rather than joining it;
//!  * **the comments of both inputs being dropped**;
//!  * **and the two ways a dictionary can disagree being two different refusals.**

use std::io::Read;

use picard_analysis::merge_bam_alignment_header::{
    adopted_program, check_dictionary, command_line_program, merged_header, Options, ProgramRecord,
    Refusal, SortOrder,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/merge_bam_alignment_header.txt.gz");
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

/// The dictionary the fixture's reference has, as `CreateSequenceDictionary` wrote it.
fn dictionary() -> Vec<String> {
    vec![
        "@SQ\tSN:chr1\tLN:40\tM5:9889878875bfc855a532253c415dceb6\tUR:file://<dir>/ref.fasta"
            .to_string(),
    ]
}

/// The unmapped bam's read groups.
fn read_groups() -> Vec<String> {
    vec!["@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA".to_string()]
}

/// The aligner's own program record.
fn bwa() -> ProgramRecord {
    ProgramRecord {
        id: "bwa".to_string(),
        fields: vec![
            ("PN".to_string(), "bwa".to_string()),
            ("VN".to_string(), "1.0".to_string()),
            ("CL".to_string(), "bwa mem".to_string()),
        ],
    }
}

fn header(programs: &[ProgramRecord], options: &Options) -> String {
    merged_header(&dictionary(), &read_groups(), programs, options).join("\n")
}

/// The three files the header is built from, each contributing its own lines.
#[test]
fn the_header_comes_from_three_files() {
    let text = corpus();
    assert_eq!(
        header(&[bwa()], &Options::default()),
        field(&text, "header", "the-whole-header").expect("the golden")
    );

    // The UR is the dictionary's canonical path: naming the same reference the long way round,
    // through a directory and back out of it, writes the short spelling.
    assert_eq!(
        field(&text, "header", "a-reference-named-the-long-way-round"),
        field(&text, "header", "the-whole-header")
    );

    // The read groups are the unmapped bam's, so an aligner that rewrote the sample and the
    // library is ignored.
    assert_eq!(
        field(&text, "header", "a-read-group-the-aligner-rewrote"),
        field(&text, "header", "the-whole-header")
    );
    // And an unmapped bam with none leaves the output with none rather than borrowing.
    assert_eq!(
        merged_header(&dictionary(), &[], &[bwa()], &Options::default()).join("\n"),
        field(&text, "header", "no-read-group-at-all").expect("the golden")
    );
}

/// A program record is adopted only when there is exactly one to adopt.
#[test]
fn one_program_is_adopted_and_two_are_not() {
    let text = corpus();
    let samtools = ProgramRecord {
        id: "samtools".to_string(),
        fields: vec![
            ("PN".to_string(), "samtools".to_string()),
            ("VN".to_string(), "1.9".to_string()),
            ("PP".to_string(), "bwa".to_string()),
            ("CL".to_string(), "samtools view".to_string()),
        ],
    };
    // Two programs, chained by PP, and the output carries neither.
    assert_eq!(
        header(&[bwa(), samtools.clone()], &Options::default()),
        field(&text, "header", "a-chain-of-programs").expect("the golden")
    );
    assert_eq!(adopted_program(&[bwa(), samtools], None), None);
    assert_eq!(adopted_program(&[], None), None);
    assert_eq!(adopted_program(&[bwa()], None), Some(bwa()));

    // Two aligned files declaring the same id are merged into one record before the count is
    // taken, so one program is what the output carries.
    assert_eq!(
        field(
            &text,
            "header",
            "two-aligned-files-with-the-same-program-id"
        ),
        field(&text, "header", "the-whole-header")
    );
}

/// A program record from the command line replaces the aligner's.
#[test]
fn a_program_from_the_command_line_replaces_the_aligners() {
    let text = corpus();
    let mine = command_line_program(
        Some("mine"),
        Some("miner"),
        Some("3.0"),
        Some("mine --do-it"),
    )
    .expect("a whole record")
    .expect("a record");
    let options = Options {
        program_record: Some(mine.clone()),
        ..Options::default()
    };
    assert_eq!(
        header(&[bwa()], &options),
        field(&text, "header", "a-program-from-the-command-line").expect("the golden")
    );
    // The aligner's record is not joined by the caller's; it is gone.
    assert_eq!(adopted_program(&[bwa()], Some(&mine)), Some(mine));

    // Three of the four arguments come together or none of them do.
    assert_eq!(
        command_line_program(Some("mine"), None, None, None),
        Err(Refusal::IncompleteProgramGroup)
    );
    assert_eq!(
        Refusal::IncompleteProgramGroup.message(),
        field(&text, "refusal", "half-a-program-from-the-command-line").expect("the golden")
    );
    assert_eq!(command_line_program(None, None, None, None), Ok(None));
}

/// A comment in either input is not a comment in the output.
#[test]
fn the_comments_are_dropped() {
    let text = corpus();
    assert_eq!(
        field(&text, "header", "comments-in-the-inputs"),
        field(&text, "header", "the-whole-header")
    );
    assert!(!header(&[bwa()], &Options::default()).contains("@CO"));
}

/// The sort order asked for is the one the header declares.
#[test]
fn the_sort_order_is_on_the_header() {
    let text = corpus();
    for (order, case) in [
        (SortOrder::Queryname, "a-queryname-output"),
        (SortOrder::Unsorted, "an-unsorted-output"),
    ] {
        let options = Options {
            sort_order: order,
            ..Options::default()
        };
        assert_eq!(
            header(&[bwa()], &options),
            field(&text, "header", case).expect("the golden"),
            "{case}"
        );
    }
}

/// The two ways a dictionary can disagree are two different refusals.
#[test]
fn a_dictionary_that_disagrees_is_refused() {
    let text = corpus();
    let dictionary = [("chr1".to_string(), 40)];

    assert_eq!(check_dictionary(&dictionary, &dictionary), Ok(()));

    let longer = [("chr1".to_string(), 41)];
    let refusal = check_dictionary(&dictionary, &longer).expect_err("a refusal");
    assert_eq!(
        refusal,
        Refusal::SequenceLengthsDiffer {
            name: "chr1".to_string(),
            first: 41,
            second: 40,
        }
    );
    assert_eq!(
        format!("java.lang.IllegalArgumentException:{}", refusal.message()),
        field(&text, "error", "an-aligned-header-that-disagrees").expect("the golden")
    );

    let elsewhere = [("chr2".to_string(), 40)];
    let refusal = check_dictionary(&dictionary, &elsewhere).expect_err("a refusal");
    assert_eq!(
        refusal,
        Refusal::DifferentSequences {
            aligned: vec!["chr2".to_string()],
            dictionary: vec!["chr1".to_string()],
        }
    );
    assert_eq!(
        format!("java.lang.IllegalArgumentException:{}", refusal.message()),
        field(&text, "error", "an-aligned-header-with-another-contig").expect("the golden")
    );
}
