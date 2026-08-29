//! Conformance for `CheckIlluminaDirectory` against Picard 3.4.0.
//!
//! Golden from `tools/illumina-conformance/CheckIlluminaDirectoryDump.java`, twelve runs over a
//! four-cycle directory written byte by byte.
//!
//! # What this suite is for
//!
//!  * **the read structure deciding how many cycle files are wanted**;
//!  * **the data types deciding which kinds are, the positions not being in the default set**;
//!  * **the status being the COUNT of what is missing**;
//!  * **and a tile the lane does not declare being a refusal rather than a count.**

use std::io::Read;

use picard_analysis::check_illumina_directory::{
    default_data_types, failures, needed, DataType, Needed,
};
use picard_analysis::illumina_files::parse_read_structure;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/check_illumina_directory.txt.gz");
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

/// The fixture's directory: four cycles, a filter and the lane's `s.locs`, minus what a case took.
fn present(removed: &[Needed]) -> impl Fn(&Needed) -> bool + '_ {
    move |file: &Needed| {
        let there = match file {
            Needed::BaseCall { cycle } => *cycle <= 4,
            Needed::Filter | Needed::Positions => true,
        };
        there && !removed.contains(file)
    }
}

/// A complete directory is zero, and a missing file is one.
#[test]
fn the_status_is_the_count_of_what_is_missing() {
    let text = corpus();
    let four = parse_read_structure("4T").expect("a structure");
    let wanted = needed(&four, &default_data_types());
    assert_eq!(failures(&present(&[]), &wanted), 0);
    assert_eq!(
        field(&text, "code", "a-complete-directory").as_deref(),
        Some("0")
    );

    let without_a_cycle = present(&[Needed::BaseCall { cycle: 4 }]);
    assert_eq!(failures(&without_a_cycle, &wanted), 1);
    assert_eq!(
        field(&text, "code", "a-missing-cycle").as_deref(),
        Some("1")
    );

    // The same directory asked for fewer cycles wants that file no longer.
    let three = parse_read_structure("3T").expect("a structure");
    let shorter = needed(&three, &default_data_types());
    assert_eq!(failures(&without_a_cycle, &shorter), 0);
    assert_eq!(
        field(&text, "code", "a-missing-cycle-not-asked-for").as_deref(),
        Some("0")
    );

    // And more cycles than the directory has are missing files of their own.
    let six = parse_read_structure("6T").expect("a structure");
    assert_eq!(
        failures(&present(&[]), &needed(&six, &default_data_types())),
        2
    );
    assert_eq!(
        field(&text, "code", "more-cycles-than-there-are").as_deref(),
        Some("1")
    );

    // The filter is always asked for.
    assert_eq!(failures(&present(&[Needed::Filter]), &wanted), 1);
    assert_eq!(
        field(&text, "code", "a-missing-filter").as_deref(),
        Some("1")
    );
}

/// The positions are not in the default set, so their absence is a failure only once asked for.
#[test]
fn the_positions_have_to_be_asked_for() {
    let text = corpus();
    let four = parse_read_structure("4T").expect("a structure");
    assert!(!default_data_types().contains(&DataType::Position));

    let without_locs = present(&[Needed::Positions]);
    assert_eq!(
        failures(&without_locs, &needed(&four, &default_data_types())),
        0
    );
    assert_eq!(
        field(&text, "code", "a-missing-s-locs").as_deref(),
        Some("1")
    );

    let asked = needed(&four, &[DataType::Position]);
    assert_eq!(failures(&without_locs, &asked), 1);
    assert_eq!(
        field(&text, "code", "a-missing-s-locs-with-positions-asked-for").as_deref(),
        Some("1")
    );
    // With the positions there and asked for, nothing is missing.
    assert_eq!(failures(&present(&[]), &asked), 0);
    assert_eq!(
        field(
            &text,
            "code",
            "a-complete-directory-with-positions-asked-for"
        )
        .as_deref(),
        Some("0")
    );
    // The per-tile `.locs` is not what is asked for: removing it changes nothing.
    assert_eq!(field(&text, "code", "a-missing-locs").as_deref(), Some("0"));
    assert_eq!(
        field(&text, "code", "a-missing-locs-with-basecalls-only").as_deref(),
        Some("0")
    );
}

/// The basecalls and their qualities are one file, so asking for both asks for one.
#[test]
fn the_basecalls_and_the_qualities_are_one_file() {
    let four = parse_read_structure("4T").expect("a structure");
    let both = needed(&four, &[DataType::BaseCalls, DataType::QualityScores]);
    let one = needed(&four, &[DataType::BaseCalls]);
    assert_eq!(both, one);
    assert_eq!(both.len(), 4);
    // And the filter is a file of its own beside them.
    let with_filter = needed(&four, &default_data_types());
    assert_eq!(with_filter.len(), 5);
    assert!(with_filter.contains(&Needed::Filter));
}

/// A tile the lane does not declare is a refusal rather than a count.
#[test]
fn a_tile_that_is_not_there_is_a_refusal() {
    let text = corpus();
    assert_eq!(
        field(&text, "code", "a-tile-that-is-there").as_deref(),
        Some("0")
    );
    let refusal = field(&text, "error", "a-tile-that-is-not").expect("the golden's refusal");
    assert!(
        refusal.contains("0 input tiles were specified"),
        "{refusal}"
    );
    // The port counts files and does not enumerate tiles: which tiles a lane HAS is the tile
    // metrics file's answer, which `illumina_lane_metrics` reads.
}
