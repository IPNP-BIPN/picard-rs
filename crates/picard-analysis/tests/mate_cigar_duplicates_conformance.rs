//! Conformance for the two mate-cigar duplicate markers against Picard 3.4.0.
//!
//! Golden from `tools/matecigardup-conformance/MateCigarDuplicatesDump.java`: eleven inputs
//! through `MarkDuplicatesWithMateCigar`, `SimpleMarkDuplicatesWithMateCigar` and
//! `MarkDuplicates`, so what is compared is the three tools' answers to the same file.
//!
//! # What this suite is for
//!
//!  * **the seven cases where all three agree, which is most of the claim**;
//!  * **a pair with no `MC` being skipped, and the skip removing it from its set so a set of two
//!    marks neither**;
//!  * **the two refusals, each with the exception the reference throws and the wording it throws
//!    it with**;
//!  * **and both tools refusing a queryname-sorted file that `MarkDuplicates` accepts.**

use std::io::Read;

use htsjdk_bam::text_parse::parse_cigar;
use picard_analysis::mark_duplicates::{mark, Options, Record};
use picard_analysis::mate_cigar_duplicates::{
    mark_with_mate_cigar, simple_mark_with_mate_cigar, MateCigarOptions, Refusal, SortOrder,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mate_cigar_duplicates.txt.gz");
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
                .replace("\\\\", "\\")
        })
}

fn record(line: &str) -> Record {
    let columns: Vec<&str> = line.split('\t').collect();
    let flags: u16 = columns[1].parse().expect("the flags");
    let tag = |name: &str| -> Option<String> {
        columns
            .iter()
            .skip(11)
            .find(|column| column.starts_with(&format!("{name}:")))
            .map(|column| column.rsplit(':').next().expect("a tag value").to_string())
    };
    Record {
        name: columns[0].to_string(),
        flags,
        reference_index: 0,
        alignment_start: columns[3].parse().expect("the position"),
        cigar: parse_cigar(columns[5]).expect("the cigar"),
        qualities: columns[10].bytes().map(|byte| byte - 33).collect(),
        mate_reference_index: if columns[6] == "*" { -1 } else { 0 },
        library: "lib1".to_string(),
        read_group: 0,
        barcode: None,
        existing_dt: None,
        mate_cigar: tag("MC").map(|text| parse_cigar(&text).expect("the mate cigar")),
    }
}

fn records(text: &str, case: &str) -> Vec<Record> {
    field(text, "sam", case)
        .unwrap_or_else(|| panic!("{case}"))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(record)
        .collect()
}

/// The output's records as the golden wrote them: the name and the flags.
fn marked(text: &str, label: &str) -> Option<Vec<(String, u16)>> {
    field(text, "marked", label).map(|body| {
        body.lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                (
                    columns[0].to_string(),
                    columns[1].parse().expect("the flags"),
                )
            })
            .collect()
    })
}

fn produced(
    records: &[Record],
    marking: &picard_analysis::mark_duplicates::Marking,
) -> Vec<(String, u16)> {
    records
        .iter()
        .enumerate()
        .filter(|(index, _)| marking.written[*index])
        .map(|(index, record)| {
            let mut flags = record.flags & !0x400;
            if marking.duplicate[index] {
                flags |= 0x400;
            }
            (record.name.clone(), flags)
        })
        .collect()
}

/// The cases that need no argument of their own, which is most of them.
const PLAIN: [&str; 7] = [
    "two-pairs",
    "a-soft-clipped-mate",
    "a-soft-clipped-first-end",
    "a-distant-mate",
    "three-singles",
    "optical-duplicates",
    "remove-duplicates",
];

fn options(case: &str) -> Options {
    let mut options = Options::default();
    if case == "remove-duplicates" {
        options.remove_duplicates = true;
    }
    options
}

/// Both tools mark what the reference marked, on every case that reaches them.
#[test]
fn both_tools_mark_what_the_reference_marked() {
    let text = corpus();
    for case in PLAIN {
        let input = records(&text, case);
        let base = options(case);
        let with = mark_with_mate_cigar(
            &input,
            SortOrder::Coordinate,
            &MateCigarOptions {
                base: base.clone(),
                ..MateCigarOptions::default()
            },
        )
        .expect("a coordinate-sorted file");
        let simple = simple_mark_with_mate_cigar(&input, SortOrder::Coordinate, &base)
            .expect("a coordinate-sorted file");
        assert_eq!(
            produced(&input, &with),
            marked(&text, &format!("{case}.withmatecigar")).expect("the golden"),
            "{case}"
        );
        assert_eq!(
            produced(&input, &simple),
            marked(&text, &format!("{case}.simple")).expect("the golden"),
            "{case}"
        );
        // And `MarkDuplicates` on the same file, which the golden recorded beside them: reading
        // the mate's cigar changes nothing when both ends are in the file.
        let plain = mark(&input, &base);
        assert_eq!(
            produced(&input, &plain),
            marked(&text, &format!("{case}.markduplicates")).expect("the golden"),
            "{case}"
        );
    }
}

/// A pair with no mate cigar leaves the run, and takes its set's other pair with it.
#[test]
fn a_pair_with_no_mate_cigar_is_skipped() {
    let text = corpus();
    let input = records(&text, "no-mate-cigar");
    let marking = mark_with_mate_cigar(&input, SortOrder::Coordinate, &MateCigarOptions::default())
        .expect("the skip");
    assert_eq!(
        produced(&input, &marking),
        marked(&text, "no-mate-cigar.withmatecigar").expect("the golden")
    );
    // Nothing is marked, where `MarkDuplicates` on the same file marks one pair.
    assert!(marking.duplicate.iter().all(|flag| !flag));
    let plain = mark(&input, &Options::default());
    assert_eq!(plain.duplicate.iter().filter(|flag| **flag).count(), 2);
    // The metrics still examine both pairs, because the writing pass walks the whole file.
    assert_eq!(marking.metrics[0].read_pairs_examined, 2);
    assert_eq!(marking.metrics[0].read_pair_duplicates, 0);
    assert_eq!(marking.metrics[0].estimated_library_size, None);
}

/// The refusals, each with the exception the reference throws and its wording.
#[test]
fn the_refusals_are_the_reference_ones() {
    let text = corpus();
    let input = records(&text, "no-mate-cigar");

    // The skip turned off is a PicardException naming the read.
    let refused = mark_with_mate_cigar(
        &input,
        SortOrder::Coordinate,
        &MateCigarOptions {
            skip_pairs_with_no_mate_cigar: false,
            ..MateCigarOptions::default()
        },
    )
    .expect_err("the refusal");
    let recorded = field(&text, "error", "no-mate-cigar-not-skipped.withmatecigar")
        .expect("the golden's refusal");
    assert_eq!(
        recorded,
        format!("{}:{}", refused.exception(), refused.message())
    );

    // The simple one refuses the same file whatever the skip says, in htsjdk's words. The golden
    // records more of the read than the port carries, so what is compared is the prefix the
    // message is built from.
    let refused = simple_mark_with_mate_cigar(&input, SortOrder::Coordinate, &Options::default())
        .expect_err("the refusal");
    let recorded = field(&text, "error", "no-mate-cigar.simple").expect("the golden's refusal");
    assert!(
        recorded.starts_with(&format!("{}:{}", refused.exception(), refused.message())),
        "{recorded}"
    );
    assert!(matches!(refused, Refusal::MateCigarNotFound { .. }));

    // And a queryname-sorted file, which both refuse and `MarkDuplicates` accepts.
    for label in [
        "a-queryname-sorted-file.withmatecigar",
        "a-queryname-sorted-file.simple",
    ] {
        let recorded = field(&text, "error", label).expect("the golden's refusal");
        assert_eq!(
            recorded,
            format!(
                "{}:{}",
                Refusal::NotCoordinateSorted.exception(),
                Refusal::NotCoordinateSorted.message()
            ),
            "{label}"
        );
    }
    let input = records(&text, "a-queryname-sorted-file");
    assert_eq!(
        mark_with_mate_cigar(&input, SortOrder::Queryname, &MateCigarOptions::default()),
        Err(Refusal::NotCoordinateSorted)
    );
    assert_eq!(
        simple_mark_with_mate_cigar(&input, SortOrder::Queryname, &Options::default()),
        Err(Refusal::NotCoordinateSorted)
    );
    assert!(marked(&text, "a-queryname-sorted-file.markduplicates").is_some());
}
