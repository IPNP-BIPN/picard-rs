//! Conformance for `MarkDuplicates` against Picard 3.4.0.
//!
//! Golden from `tools/markduplicates-conformance/MarkDuplicatesDump.java`, eighteen runs recorded
//! as the input and the marked output in SAM text plus the metrics table and its histogram. The
//! output BAM's bytes are DEFLATE and out of reach; what the tool decides is the flag, the `DT`
//! tag and the counters, and those are all here.
//!
//! # What this suite is for
//!
//!  * **which record of a set keeps `0x400` clear, under both scoring strategies**;
//!  * **the position being the unclipped 5' one, so a soft-clipped read is a duplicate of an
//!    unclipped one**;
//!  * **a pair being a unit and a fragment sharing a key with a pair always losing**;
//!  * **optical duplicates coming off the read NAME, and the regex being turnable off**;
//!  * **`REMOVE_DUPLICATES` and `REMOVE_SEQUENCING_DUPLICATES` dropping different records**;
//!  * **`TAGGING_POLICY` and `CLEAR_DT`**;
//!  * **`BARCODE_TAG` splitting one position into two sets**;
//!  * **and the metrics, including the estimated library size and the two histograms.**

use std::io::Read;

use htsjdk_bam::text_parse::parse_cigar;
use picard_analysis::mark_duplicates::{
    mark, roi_histogram, Options, Record, ScoringStrategy, TaggingPolicy,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mark_duplicates.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> String {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
                .replace("\\\\", "\\")
        })
        .unwrap_or_else(|| panic!("{kind}/{case}"))
}

/// One SAM line of the dump, as the port's reduced record.
fn record(line: &str) -> Record {
    let columns: Vec<&str> = line.split('\t').collect();
    let flags: u16 = columns[1].parse().expect("the flags");
    let unmapped = flags & 0x4 != 0;
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
        reference_index: if unmapped { -1 } else { 0 },
        alignment_start: columns[3].parse().expect("the position"),
        cigar: if columns[5] == "*" {
            parse_cigar("").expect("an empty cigar")
        } else {
            parse_cigar(columns[5]).expect("the cigar")
        },
        qualities: columns[10].bytes().map(|byte| byte - 33).collect(),
        mate_reference_index: if columns[6] == "*" { -1 } else { 0 },
        // The fixture's one read group, whose library the header names.
        library: "lib1".to_string(),
        read_group: 0,
        barcode: tag("RX"),
        existing_dt: tag("DT"),
        // `MarkDuplicates` does not read the mate's cigar; the two tools that do have their own
        // suite, and this fixture carries the tag only where they need it.
        mate_cigar: tag("MC").map(|text| parse_cigar(&text).expect("the mate cigar")),
    }
}

fn records(text: &str, case: &str) -> Vec<Record> {
    field(text, "sam", case)
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(record)
        .collect()
}

/// The output's records, as the golden wrote them: the name, the flags and the `DT` tag.
fn marked(text: &str, case: &str) -> Vec<(String, u16, Option<String>)> {
    field(text, "marked", case)
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let dt = columns
                .iter()
                .skip(11)
                .find(|column| column.starts_with("DT:"))
                .map(|column| column.rsplit(':').next().expect("a tag").to_string());
            (
                columns[0].to_string(),
                columns[1].parse().expect("the flags"),
                dt,
            )
        })
        .collect()
}

/// One row of the metrics table, by column name.
fn metrics(text: &str, case: &str) -> std::collections::HashMap<String, String> {
    let table = field(text, "metrics", case);
    let mut lines = table.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let row: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    header
        .iter()
        .zip(row)
        .map(|(name, value)| ((*name).to_string(), value.to_string()))
        .collect()
}

/// The histogram section, as bins by column.
fn histogram(text: &str, case: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let section = field(text, "histogram", case);
    if section.trim().is_empty() {
        // A run with no pairs has no histogram at all: no library size to estimate a return on,
        // and no duplicate set to count the size of.
        return (Vec::new(), Vec::new());
    }
    let mut lines = section.lines();
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .map(str::to_string)
        .collect();
    let rows = lines
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    (header, rows)
}

/// The arguments each case ran with, which the dump's own labels name.
fn options(case: &str) -> Options {
    let mut options = Options::default();
    match case {
        "optical-duplicates-without-the-regex" => options.parse_read_names = false,
        "optical-duplicates-with-a-smaller-distance" => {
            options.optical_duplicate_pixel_distance = 1
        }
        "pairs-of-different-lengths-by-quality" => {
            options.scoring = ScoringStrategy::SumOfBaseQualities
        }
        "remove-duplicates" => options.remove_duplicates = true,
        "remove-sequencing-duplicates" => options.remove_sequencing_duplicates = true,
        "tagging-policy-all" => options.tagging_policy = TaggingPolicy::All,
        "tagging-policy-optical" => options.tagging_policy = TaggingPolicy::OpticalOnly,
        "an-existing-dt-tag-kept" => options.clear_dt = false,
        "two-barcodes" => options.barcode_tag = Some("RX".to_string()),
        _ => {}
    }
    options
}

const CASES: [&str; 18] = [
    "two-pairs",
    "optical-duplicates",
    "optical-duplicates-without-the-regex",
    "optical-duplicates-with-a-smaller-distance",
    "pairs-of-different-lengths",
    "pairs-of-different-lengths-by-quality",
    "three-singles",
    "a-soft-clipped-read",
    "an-unmapped-read",
    "a-secondary-read",
    "remove-duplicates",
    "remove-sequencing-duplicates",
    "tagging-policy-all",
    "tagging-policy-optical",
    "an-existing-dt-tag",
    "an-existing-dt-tag-kept",
    "two-barcodes",
    "two-barcodes-ignored",
];

/// Every case's output, record for record: the flag, the tag, and which records survive.
#[test]
fn every_case_marks_what_the_reference_marked() {
    let text = corpus();
    for case in CASES {
        let input = records(&text, case);
        let expected = marked(&text, case);
        let marking = mark(&input, &options(case));

        let produced: Vec<(String, u16, Option<String>)> = input
            .iter()
            .enumerate()
            .filter(|(index, _)| marking.written[*index])
            .map(|(index, record)| {
                let mut flags = record.flags & !0x400;
                if marking.duplicate[index] {
                    flags |= 0x400;
                }
                (
                    record.name.clone(),
                    flags,
                    marking.duplicate_type[index].clone(),
                )
            })
            .collect();
        assert_eq!(produced, expected, "{case}");
    }
}

/// The counters, and the estimate derived from them.
#[test]
fn the_metrics_are_the_reference_ones() {
    let text = corpus();
    for case in CASES {
        let input = records(&text, case);
        let marking = mark(&input, &options(case));
        let expected = metrics(&text, case);
        assert_eq!(marking.metrics.len(), 1, "{case}");
        let row = &marking.metrics[0];
        assert_eq!(row.library, expected["LIBRARY"], "{case}");
        for (column, value) in [
            ("UNPAIRED_READS_EXAMINED", row.unpaired_reads_examined),
            ("READ_PAIRS_EXAMINED", row.read_pairs_examined),
            (
                "SECONDARY_OR_SUPPLEMENTARY_RDS",
                row.secondary_or_supplementary,
            ),
            ("UNMAPPED_READS", row.unmapped_reads),
            ("UNPAIRED_READ_DUPLICATES", row.unpaired_read_duplicates),
            ("READ_PAIR_DUPLICATES", row.read_pair_duplicates),
            (
                "READ_PAIR_OPTICAL_DUPLICATES",
                row.read_pair_optical_duplicates,
            ),
        ] {
            assert_eq!(
                expected[column].parse::<i64>().expect("a count"),
                value,
                "{case}/{column}"
            );
        }
        // The percentage is written with the metrics file's own rounding, so the comparison is
        // against the rendered value rather than against a float this test would have to guess at.
        let recorded: f64 = expected["PERCENT_DUPLICATION"].parse().expect("a fraction");
        assert!(
            (recorded - row.percent_duplication).abs() < 1e-6,
            "{case}: {recorded} vs {}",
            row.percent_duplication
        );
        match row.estimated_library_size {
            None => assert!(
                expected["ESTIMATED_LIBRARY_SIZE"].is_empty(),
                "{case}: {}",
                expected["ESTIMATED_LIBRARY_SIZE"]
            ),
            Some(size) => assert_eq!(
                expected["ESTIMATED_LIBRARY_SIZE"]
                    .parse::<i64>()
                    .expect("a size"),
                size,
                "{case}"
            ),
        }
    }
}

/// The histogram beside the table: the ROI where there is one, and the set sizes always.
#[test]
fn the_histograms_are_the_reference_ones() {
    let text = corpus();
    for case in CASES {
        let input = records(&text, case);
        let marking = mark(&input, &options(case));
        let (header, rows) = histogram(&text, case);
        // The first column is the ROI's bin where the estimate exists and the set size where it
        // does not: `calculateRoiHistogram` returns nothing when the library size cannot be
        // estimated, and the metrics file then prints the set-size histograms alone.
        let roi = roi_histogram(&marking.metrics[0]);
        if header.is_empty() {
            assert!(roi.is_none(), "{case}");
            assert!(marking.all_sets.is_empty(), "{case}");
            assert!(marking.optical_sets.is_empty(), "{case}");
            assert!(marking.non_optical_sets.is_empty(), "{case}");
            continue;
        }
        assert_eq!(header[0] == "BIN", roi.is_some(), "{case}: {header:?}");
        // A column the golden does not print is a histogram the run never filled.
        for (name, counts) in [
            ("all_sets", &marking.all_sets),
            ("optical_sets", &marking.optical_sets),
            ("non_optical_sets", &marking.non_optical_sets),
        ] {
            if !header.iter().any(|column| column == name) {
                assert!(counts.is_empty(), "{case}/{name}: {counts:?}");
            }
        }

        for row in &rows {
            let bin: f64 = row[0].parse().expect("a bin");
            for (column, name) in header.iter().enumerate().skip(1) {
                let expected = &row[column];
                match name.as_str() {
                    "CoverageMult" => {
                        let bins = roi.as_ref().expect("the roi");
                        let value = bins
                            .iter()
                            .find(|(id, _)| *id == bin)
                            .map(|(_, value)| *value)
                            .expect("the bin");
                        // The metrics file writes a double through its own formatter, which keeps
                        // six significant digits; the comparison is at that width.
                        let recorded: f64 = expected.parse().expect("a value");
                        assert!(
                            (recorded - value).abs() < 5e-7 * value.max(1.0),
                            "{case}/{bin}: {recorded} vs {value}"
                        );
                    }
                    "all_sets" | "optical_sets" | "non_optical_sets" => {
                        let counts = match name.as_str() {
                            "all_sets" => &marking.all_sets,
                            "optical_sets" => &marking.optical_sets,
                            _ => &marking.non_optical_sets,
                        };
                        let value = counts
                            .iter()
                            .find(|(id, _)| *id == bin)
                            .map(|(_, count)| *count)
                            .unwrap_or(0.0);
                        let recorded: f64 = expected.parse().expect("a count");
                        assert_eq!(recorded, value, "{case}/{name}/{bin}");
                    }
                    other => panic!("{case}: an unexpected column {other}"),
                }
            }
        }
    }
}
