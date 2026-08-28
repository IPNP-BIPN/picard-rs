//! Conformance for `CollectDuplicateMetrics` against Picard 3.4.0.
//!
//! Each case carries the file the tool read, as SAM without its header, and the metrics table it
//! wrote. The read group of every record names its library, so the port reads the flags and the
//! `RG` tag out of the SAM and must produce the same rows.
//!
//! # What this suite is for
//!
//!  * **the four counters being a chain, so a read reaches exactly one**;
//!  * **a half-mapped pair contributing one unpaired read and one unmapped one**;
//!  * **the paired counts being halved by an integer division**;
//!  * **an unmapped or secondary duplicate counting nowhere**;
//!  * **the optical count being always zero, and the library size computed from it**;
//!  * **the percentage weighting the pairs by two on both sides**;
//!  * **one row per library named in the header, used or not**;
//!  * **a read group with no library falling under `Unknown Library`**;
//!  * **the histogram being written only for a one-library file with an estimated size**;
//!  * **and an estimated size of zero making every bin NaN, written `?`.**

use std::io::Read;

use picard_analysis::collect_duplicate_metrics::{
    collect, writes_a_histogram, DuplicationMetrics, Record, UNKNOWN_LIBRARY,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_duplicate_metrics.txt.gz");
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

/// The libraries each case's header named, which the dump's own read groups fix.
fn libraries(case: &str) -> Vec<String> {
    match case {
        "two-libraries" => vec!["lib1".to_string(), "lib2".to_string()],
        "no-library" => vec![UNKNOWN_LIBRARY.to_string()],
        _ => vec!["lib1".to_string()],
    }
}

/// The library a read group belongs to, which the dump fixes the same way.
fn library_of(case: &str, group: &str) -> String {
    match (case, group) {
        ("no-library", _) => UNKNOWN_LIBRARY.to_string(),
        (_, "rg2") => "lib2".to_string(),
        _ => "lib1".to_string(),
    }
}

fn records(text: &str, case: &str) -> Vec<Record> {
    field(text, "sam", case)
        .unwrap_or_else(|| panic!("{case} has an input"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("a flag word");
            let group = columns
                .iter()
                .find_map(|column| column.strip_prefix("RG:Z:"))
                .expect("a read group");
            Record {
                library: library_of(case, group),
                duplicate: flags & 0x400 != 0,
                secondary_or_supplementary: flags & 0x100 != 0 || flags & 0x800 != 0,
                unmapped: flags & 0x4 != 0,
                paired: flags & 0x1 != 0,
                mate_unmapped: flags & 0x8 != 0,
            }
        })
        .collect()
}

/// The metrics table of one case, as its rows of values.
fn table(text: &str, case: &str) -> Vec<Vec<String>> {
    let payload = field(text, "metrics", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = payload.lines().filter(|line| !line.is_empty());
    let header = lines.next().expect("a header line");
    assert_eq!(header.split('\t').next(), Some("LIBRARY"), "{case}");
    lines
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// The reference writes a whole number without a point, a fraction with one, and a null as an
/// empty field.
fn number(value: f64) -> String {
    if value == value.trunc() {
        format!("{}", value as i64)
    } else {
        format!("{value}")
    }
}

fn as_written(row: &DuplicationMetrics) -> Vec<String> {
    vec![
        row.library.clone(),
        row.unpaired_reads_examined.to_string(),
        row.read_pairs_examined.to_string(),
        row.secondary_or_supplementary_reads.to_string(),
        row.unmapped_reads.to_string(),
        row.unpaired_read_duplicates.to_string(),
        row.read_pair_duplicates.to_string(),
        row.read_pair_optical_duplicates.to_string(),
        number(row.percent_duplication()),
        row.estimated_library_size()
            .map(|size| size.to_string())
            .unwrap_or_default(),
    ]
}

const CASES: &[&str] = &[
    "two-pairs-one-duplicate",
    "odd-paired-count",
    "unpaired",
    "mate-unmapped",
    "unmapped-duplicate",
    "secondary-and-supplementary",
    "two-libraries",
    "no-library",
    "all-duplicates",
    "only-unmapped",
    "empty",
];

fn ours(text: &str, case: &str) -> Vec<Vec<String>> {
    collect(&libraries(case), &records(text, case))
        .iter()
        .map(as_written)
        .collect()
}

/// Every case's rows are what the port produces.
#[test]
fn every_case_writes_the_same_rows() {
    let text = corpus();
    for case in CASES {
        assert_eq!(ours(&text, case), table(&text, case), "{case}");
    }
}

/// The four counters are a chain: an unmapped duplicate is counted as unmapped and nowhere else.
#[test]
fn the_counters_are_a_chain() {
    let text = corpus();
    let rows = collect(
        &libraries("unmapped-duplicate"),
        &records(&text, "unmapped-duplicate"),
    );
    assert_eq!(rows[0].unmapped_reads, 1);
    assert_eq!(rows[0].unpaired_read_duplicates, 0);
    assert_eq!(rows[0].unpaired_reads_examined, 1);
    // And a secondary or supplementary read reaches only its own counter.
    let rows = collect(
        &libraries("secondary-and-supplementary"),
        &records(&text, "secondary-and-supplementary"),
    );
    assert_eq!(rows[0].secondary_or_supplementary_reads, 2);
    assert_eq!(rows[0].read_pair_duplicates, 0);
}

/// A half-mapped pair contributes one unpaired read and one unmapped one.
#[test]
fn a_half_mapped_pair_splits_between_two_counters() {
    let text = corpus();
    let rows = collect(
        &libraries("mate-unmapped"),
        &records(&text, "mate-unmapped"),
    );
    assert_eq!(rows[0].unpaired_reads_examined, 1);
    assert_eq!(rows[0].unmapped_reads, 1);
    assert_eq!(rows[0].read_pairs_examined, 0);
}

/// The halving is an integer division: one paired read reports no pairs, and five report the same
/// two as four do.
#[test]
fn the_halving_is_an_integer_division() {
    let text = corpus();
    let one = collect(
        &libraries("secondary-and-supplementary"),
        &records(&text, "secondary-and-supplementary"),
    );
    assert_eq!(one[0].read_pairs_examined, 0);
    let four = table(&text, "two-pairs-one-duplicate");
    let five = table(&text, "odd-paired-count");
    assert_eq!(four, five);
    assert_eq!(four[0][2], "2");
}

/// The optical count is always zero, and the library size is computed from that zero.
#[test]
fn the_optical_count_is_always_zero() {
    let text = corpus();
    for case in CASES {
        for row in table(&text, case) {
            assert_eq!(row[7], "0", "{case}");
        }
    }
    let rows = collect(
        &libraries("two-pairs-one-duplicate"),
        &records(&text, "two-pairs-one-duplicate"),
    );
    assert_eq!(rows[0].read_pair_optical_duplicates, 0);
    assert_eq!(rows[0].estimated_library_size(), Some(1));
}

/// The percentage weights the pairs by two on both sides, and is zero when nothing was examined.
#[test]
fn the_percentage_weights_the_pairs_by_two() {
    let text = corpus();
    let pairs = collect(
        &libraries("two-pairs-one-duplicate"),
        &records(&text, "two-pairs-one-duplicate"),
    );
    // One duplicate pair of two pairs: (0 + 1*2) / (0 + 2*2).
    assert_eq!(pairs[0].percent_duplication(), 0.5);
    let unpaired = collect(&libraries("unpaired"), &records(&text, "unpaired"));
    // One duplicate read of two: (1 + 0) / (2 + 0), which is the same number by another route.
    assert_eq!(unpaired[0].percent_duplication(), 0.5);
    let none = collect(
        &libraries("only-unmapped"),
        &records(&text, "only-unmapped"),
    );
    assert_eq!(none[0].percent_duplication(), 0.0);
}

/// There is one row per library named in the header, whether a read ever used it or not.
#[test]
fn a_library_gets_a_row_whether_it_is_used_or_not() {
    let text = corpus();
    let empty = table(&text, "empty");
    assert_eq!(empty.len(), 1);
    assert_eq!(empty[0][0], "lib1");
    assert!(empty[0][1..8].iter().all(|value| value == "0"));
    assert!(records(&text, "empty").is_empty());
    // And two libraries give two rows, in the header's order.
    let two = table(&text, "two-libraries");
    assert_eq!(
        two.iter().map(|row| row[0].clone()).collect::<Vec<_>>(),
        vec!["lib1".to_string(), "lib2".to_string()]
    );
}

/// A read group with no library at all falls under the literal.
#[test]
fn a_read_group_with_no_library_has_a_name_of_its_own() {
    let text = corpus();
    let rows = table(&text, "no-library");
    assert_eq!(rows[0][0], UNKNOWN_LIBRARY);
    assert_eq!(UNKNOWN_LIBRARY, "Unknown Library");
}

/// The histogram is written only for a one-library file whose estimated size is not null, and an
/// estimated size of zero makes every bin NaN, which the writer renders as `?`.
#[test]
fn the_histogram_needs_one_library_and_a_size() {
    let text = corpus();
    let histogram = |case: &str| field(&text, "histogram", case).unwrap_or_default();
    // Two libraries: none, though each row has its own numbers.
    assert!(histogram("two-libraries").is_empty());
    // One library and no estimated size: none either.
    assert!(histogram("unpaired").is_empty());
    let unpaired = collect(&libraries("unpaired"), &records(&text, "unpaired"));
    assert_eq!(unpaired[0].estimated_library_size(), None);
    assert!(!writes_a_histogram(&unpaired));
    // One library with a size: written.
    let pairs = collect(
        &libraries("two-pairs-one-duplicate"),
        &records(&text, "two-pairs-one-duplicate"),
    );
    assert!(writes_a_histogram(&pairs));
    assert!(histogram("two-pairs-one-duplicate").contains("BIN\tCoverageMult"));
    // A size of zero: written, and every bin is NaN.
    let all = collect(
        &libraries("all-duplicates"),
        &records(&text, "all-duplicates"),
    );
    assert_eq!(all[0].estimated_library_size(), Some(0));
    assert!(writes_a_histogram(&all));
    let bins = histogram("all-duplicates");
    assert!(bins.contains("1.0\t?"), "{bins}");
}
