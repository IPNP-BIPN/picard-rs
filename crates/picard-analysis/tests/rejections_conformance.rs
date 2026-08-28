//! Conformance for the refusals `SinglePassSamProgram` reaches before a record is counted, against
//! Picard 3.4.0.
//!
//! Golden from `tools/rejection-conformance/RejectionDump.java`.
//!
//! The covering arrays are generated over combinations the tool accepts, so a row spent being
//! rejected is a row not spent on the tool. These four rows are what the arrays cannot hold: the
//! obsolete `FLOW_MODE`, the sort-order refusal on both tools that share the driver, and what
//! happens when `ASSUME_SORTED=true` takes the refusal away.

use std::io::Read;

use picard_analysis::single_pass_rejections::{
    check_flow_mode, check_sort_order, walk_reference, Rejection, SortOrder,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/rejections.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

/// The `THROWN=` of one recorded case.
fn thrown(text: &str, label: &str) -> String {
    let prefix = format!("reject\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .strip_prefix("THROWN=")
        .expect("a thrown exception")
        .to_string()
}

/// The fixture directory the harness masks, which is a per-run temporary path.
const FIXTURES: &str = "<FIXTURES>";

#[test]
fn every_recorded_refusal_matches_the_golden() {
    let text = corpus();
    let mut compared = 0;

    // FLOW_MODE moved to another tool, and is refused by name.
    let error = check_flow_mode(true).expect_err("an obsolete argument");
    assert_eq!(
        error.thrown(),
        thrown(&text, "qualityyield_flow_mode_obsolete")
    );
    compared += 1;

    // The same sort-order refusal from both tools that share the driver, naming the file, the
    // order that was found, and the argument that would bypass it.
    let path = format!("{FIXTURES}/queryname.bam");
    for label in [
        "qualityyield_queryname_not_assumed_sorted",
        "alignmentsummary_queryname_not_assumed_sorted",
    ] {
        let error = check_sort_order(&path, SortOrder::Queryname, false).expect_err(label);
        assert_eq!(error.thrown(), thrown(&text, label), "{label}");
        compared += 1;
    }

    // And taking the escape moves the failure into htsjdk rather than removing it.
    assert!(check_sort_order(&path, SortOrder::Queryname, true).is_ok());
    let error = walk_reference(Some(1), 0).expect_err("a walker asked to rewind");
    assert_eq!(
        error.thrown(),
        thrown(&text, "alignmentsummary_queryname_assumed_sorted")
    );
    compared += 1;

    assert_eq!(compared, 4, "the golden's recorded refusals");
    assert_eq!(
        text.lines()
            .filter(|line| line.starts_with("reject\t"))
            .count(),
        4,
        "and there are no others"
    );
}

/// The two libraries answer differently, which is the point of the fourth row.
#[test]
fn the_two_refusals_come_from_two_libraries() {
    let text = corpus();
    assert!(thrown(&text, "qualityyield_queryname_not_assumed_sorted")
        .starts_with("picard.PicardException:"));
    assert!(thrown(&text, "alignmentsummary_queryname_assumed_sorted")
        .starts_with("htsjdk.samtools.SAMException:"));
    assert_eq!(
        Rejection::EarlierReferenceSequence {
            requested: 0,
            current: 1
        }
        .java_class(),
        "htsjdk.samtools.SAMException"
    );
}

/// A coordinate-sorted input passes the check whether or not the escape is given, and the escape
/// is what the message names.
#[test]
fn the_sort_check_names_its_own_escape() {
    for assume in [false, true] {
        assert!(check_sort_order("/tmp/x.bam", SortOrder::Coordinate, assume).is_ok());
    }
    let message = check_sort_order("/tmp/x.bam", SortOrder::Unsorted, false)
        .expect_err("unsorted")
        .message();
    assert!(message.contains("ASSUME_SORTED=true"));
    // The order that was FOUND is named, not the one that was wanted.
    assert!(message.contains("the sort order is unsorted"));
    assert!(message.starts_with("File /tmp/x.bam should be coordinate sorted"));
}

/// The walker refuses only a backwards request, and remembers where it is.
#[test]
fn the_walker_refuses_only_a_rewind() {
    assert_eq!(
        walk_reference(None, 5),
        Ok(5),
        "the first request sets the position"
    );
    assert_eq!(walk_reference(Some(1), 1), Ok(1), "the same contig again");
    assert_eq!(walk_reference(Some(1), 2), Ok(2), "forwards");
    assert!(walk_reference(Some(2), 1).is_err(), "backwards");
}
