//! Conformance for `CheckDuplicateMarking` against Picard 3.4.0.
//!
//! Each case carries the file the tool read, as SAM without its header, the exit code it returned
//! and the lines it wrote to `--OUTPUT`. The port walks the same records and must reach the same
//! verdict.
//!
//! # What this suite is for
//!
//!  * **the exit code being the count's sign and never the count**;
//!  * **the comparison being against the first record of each name, so a name whose flags go
//!    true, false, false is bad twice**;
//!  * **a name that appears once never being bad**;
//!  * **the output holding one line per bad record**;
//!  * **a coordinate-sorted file being sorted by query name before any of it**;
//!  * **each mode skipping what its name says, before anything is compared**;
//!  * **and a skipped first record not becoming the name's reference either.**

use std::io::Read;

use picard_analysis::check_duplicate_marking::{check, is_skipped, Mode, Record};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("check_duplicate_marking.txt.gz");
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

fn field(text: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("{kind} {name} is in the corpus"));
    unescape(&line[prefix.len()..])
}

/// The records of one case, read out of its SAM body.
fn records(text: &str, case: &str) -> Vec<Record> {
    field(text, "sam", case)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("a flag word");
            Record {
                name: columns[0].to_string(),
                duplicate: flags & 0x400 != 0,
                secondary_or_supplementary: flags & 0x100 != 0 || flags & 0x800 != 0,
                unmapped: flags & 0x4 != 0,
                proper_pair: flags & 0x2 != 0,
            }
        })
        .collect()
}

/// The bad-name lines of one case.
fn bad_names(text: &str, case: &str) -> Vec<String> {
    field(text, "bad", case)
        .lines()
        .filter(|line| !line.is_empty())
        .map(str::to_string)
        .collect()
}

fn exit_code(text: &str, case: &str) -> i32 {
    field(text, "verdict", case).trim().parse().expect("a code")
}

/// The mode each case ran under, which the dump names on the command line rather than in the
/// output. `no-output-file` ran with no `--OUTPUT` at all, so its bad lines are empty whatever it
/// found, which is why its exit code is checked apart from them.
const MODES: &[(&str, Mode)] = &[
    ("all-agree", Mode::All),
    ("one-disagreement", Mode::All),
    ("two-disagreements-one-name", Mode::All),
    ("single-record", Mode::All),
    ("coordinate-sorted", Mode::All),
    ("secondary-disagrees-all", Mode::All),
    ("secondary-disagrees-primary-only", Mode::PrimaryOnly),
    ("supplementary-disagrees-all", Mode::All),
    ("supplementary-disagrees-primary-only", Mode::PrimaryOnly),
    ("unmapped-disagrees-primary-only", Mode::PrimaryOnly),
    ("unmapped-disagrees-mapped-only", Mode::PrimaryMappedOnly),
    ("improper-disagrees-mapped-only", Mode::PrimaryMappedOnly),
    (
        "improper-disagrees-proper-only",
        Mode::PrimaryProperPairOnly,
    ),
    ("skipped-first-record", Mode::PrimaryOnly),
];

/// `getSortedRecordsFromReader`, which is a stable sort by query name.
fn sorted_by_name(records: &[Record], coordinate_sorted: bool) -> Vec<Record> {
    let mut sorted = records.to_vec();
    if coordinate_sorted {
        sorted.sort_by(|a, b| a.name.cmp(&b.name));
    }
    sorted
}

/// Every case's exit code and bad-name lines are what the port reaches.
#[test]
fn every_case_reaches_the_same_verdict() {
    let text = corpus();
    for (case, mode) in MODES {
        let records = sorted_by_name(&records(&text, case), *case == "coordinate-sorted");
        let verdict = check(&records, *mode);
        assert_eq!(
            verdict.bad_names,
            bad_names(&text, case),
            "{case} bad names"
        );
        assert_eq!(
            verdict.exit_code(),
            exit_code(&text, case),
            "{case} exit code"
        );
    }
}

/// The comparison is against the first record of the name: three records whose flags go true,
/// false, false report TWO bad ones and write the name twice.
#[test]
fn the_comparison_is_against_the_first_record() {
    let text = corpus();
    let verdict = check(&records(&text, "two-disagreements-one-name"), Mode::All);
    assert_eq!(verdict.bad_names, vec!["a".to_string(), "a".to_string()]);
    assert_eq!(verdict.exit_code(), 1);
    // The code is the SIGN and not the count, which two bad records is what shows.
    assert_eq!(exit_code(&text, "two-disagreements-one-name"), 1);
    assert_eq!(bad_names(&text, "two-disagreements-one-name").len(), 2);
}

/// A skipped record is never compared and never remembered. The query-name sort puts this case's
/// secondary record LAST, so under `PRIMARY_ONLY` the two primaries agree and nothing is bad,
/// while under `ALL` the secondary disagrees with them and one record is.
#[test]
fn a_skipped_record_is_never_compared() {
    let text = corpus();
    let records = records(&text, "skipped-first-record");
    let secondary = records
        .iter()
        .position(|record| record.secondary_or_supplementary)
        .expect("a secondary record");
    assert_eq!(secondary, records.len() - 1, "the sort put it last");
    assert!(records[secondary].duplicate);
    assert!(records[..secondary].iter().all(|record| !record.duplicate));
    assert!(is_skipped(&records[secondary], Mode::PrimaryOnly));
    assert_eq!(
        check(&records, Mode::PrimaryOnly).bad_names,
        Vec::<String>::new()
    );
    assert_eq!(check(&records, Mode::All).bad_names, vec!["a".to_string()]);
}

/// Each mode skips what its name says, and only that.
#[test]
fn each_mode_skips_what_its_name_says() {
    let text = corpus();
    let secondary = &records(&text, "secondary-disagrees-all")[1];
    let supplementary = &records(&text, "supplementary-disagrees-all")[1];
    let unmapped = &records(&text, "unmapped-disagrees-primary-only")[1];
    let improper = &records(&text, "improper-disagrees-mapped-only")[1];
    for record in [secondary, supplementary] {
        assert!(!is_skipped(record, Mode::All));
        assert!(is_skipped(record, Mode::PrimaryOnly));
    }
    assert!(!is_skipped(unmapped, Mode::PrimaryOnly));
    assert!(is_skipped(unmapped, Mode::PrimaryMappedOnly));
    assert!(!is_skipped(improper, Mode::PrimaryMappedOnly));
    assert!(is_skipped(improper, Mode::PrimaryProperPairOnly));
}

/// A name that appears once is never bad, and a file where every name agrees exits zero.
#[test]
fn a_name_that_appears_once_is_never_bad() {
    let text = corpus();
    for case in ["single-record", "all-agree"] {
        let verdict = check(&records(&text, case), Mode::All);
        assert_eq!(verdict.exit_code(), 0, "{case}");
        assert!(verdict.bad_names.is_empty(), "{case}");
    }
}

/// The sort is what lets a coordinate-sorted file be checked at all: read in file order the two
/// records of each name are never adjacent, and nothing would be compared.
#[test]
fn a_coordinate_sorted_file_is_sorted_first() {
    let text = corpus();
    let records = records(&text, "coordinate-sorted");
    assert_eq!(check(&records, Mode::All).bad_names, Vec::<String>::new());
    let sorted = sorted_by_name(&records, true);
    assert_eq!(
        check(&sorted, Mode::All).bad_names,
        bad_names(&text, "coordinate-sorted")
    );
}

/// Without an output file the verdict is the same, only nothing is written.
#[test]
fn the_output_file_is_optional() {
    let text = corpus();
    let verdict = check(&records(&text, "no-output-file"), Mode::All);
    assert_eq!(verdict.exit_code(), exit_code(&text, "no-output-file"));
    assert_eq!(verdict.exit_code(), 1);
    assert_eq!(verdict.bad_names, vec!["a".to_string()]);
    assert!(bad_names(&text, "no-output-file").is_empty());
}
