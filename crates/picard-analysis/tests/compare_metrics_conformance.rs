//! Conformance for `CompareMetrics` against Picard 3.4.0.
//!
//! Each case carries the two metrics files' verdict, the report and the table the tool wrote. The
//! port rebuilds the two files from the corpus's own inputs and must reach the same verdict.
//!
//! # What this suite is for
//!
//!  * **the exit code being the verdict**;
//!  * **--METRICS_TO_IGNORE and --METRICS_NOT_REQUIRED both dropping a column outright**;
//!  * **the relative change being relative to the FIRST file's value**;
//!  * **--IGNORE_HISTOGRAM_DIFFERENCES forgiving the histogram and not the table**;
//!  * **--KEY matching rows by a column rather than by position**;
//!  * **a row one file lacks being a difference either way**;
//!  * **and the report's two words.**

use std::collections::BTreeSet;
use std::io::Read;

use picard_analysis::compare_metrics::{
    compare, exit_code, parse_allowable_relative_change, report_header, status, values_agree,
    Arguments, Difference, MetricsFile,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("compare_metrics.txt.gz");
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

fn field(corpus: &str, kind: &str, name: &str) -> Option<String> {
    corpus
        .lines()
        .find(|line| line.starts_with(&format!("{kind}\t{name}\t")))
        .map(|line| unescape(&line[format!("{kind}\t{name}\t").len()..]))
}

/// One of the corpus's input files, parsed.
fn input(corpus: &str, name: &str) -> MetricsFile {
    let text =
        field(corpus, "input", name).unwrap_or_else(|| panic!("the corpus carries input/{name}"));
    let mut metric_class = String::new();
    let mut columns = Vec::new();
    let mut rows = Vec::new();
    let mut histogram = Vec::new();
    let mut in_histogram = false;
    let mut seen_header = false;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("## METRICS CLASS\t") {
            metric_class = rest.to_string();
            continue;
        }
        if line.starts_with("## HISTOGRAM") {
            in_histogram = true;
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if in_histogram {
            histogram.push(line.to_string());
        } else if !seen_header {
            columns = line.split('\t').map(str::to_string).collect();
            seen_header = true;
        } else {
            rows.push(line.split('\t').map(str::to_string).collect());
        }
    }
    MetricsFile {
        metric_class,
        columns,
        rows,
        histogram,
    }
}

fn verdict(corpus: &str, name: &str) -> i32 {
    corpus
        .lines()
        .find(|line| line.starts_with(&format!("verdict\t{name}\t")))
        .unwrap_or_else(|| panic!("the corpus carries verdict/{name}"))
        [format!("verdict\t{name}\t").len()..]
        .parse()
        .expect("a code")
}

fn ignore(columns: &[&str]) -> BTreeSet<String> {
    columns.iter().map(|c| c.to_string()).collect()
}

/// Every case whose two files the corpus carries.
#[test]
fn every_case_reaches_the_same_verdict() {
    let corpus = corpus();
    let base = input(&corpus, "base");
    let changed = input(&corpus, "changed");
    let fewer = input(&corpus, "fewer-columns");
    let other_histogram = input(&corpus, "other-histogram");
    let reordered = input(&corpus, "reordered");
    let one_row = input(&corpus, "one-row");

    let cases: Vec<(&str, &MetricsFile, &MetricsFile, Arguments)> = vec![
        ("identical", &base, &base, Arguments::default()),
        ("one-column-differs", &base, &changed, Arguments::default()),
        (
            "ignored-column",
            &base,
            &changed,
            Arguments {
                metrics_to_ignore: ignore(&["PERCENT_DUPLICATION"]),
                ..Arguments::default()
            },
        ),
        (
            "not-required-but-present",
            &base,
            &changed,
            Arguments {
                metrics_not_required: ignore(&["PERCENT_DUPLICATION"]),
                ..Arguments::default()
            },
        ),
        (
            "relative-change-generous",
            &base,
            &changed,
            Arguments {
                allowable_relative_change: vec![("PERCENT_DUPLICATION".to_string(), 0.2)],
                ..Arguments::default()
            },
        ),
        (
            "relative-change-tight",
            &base,
            &changed,
            Arguments {
                allowable_relative_change: vec![("PERCENT_DUPLICATION".to_string(), 0.05)],
                ..Arguments::default()
            },
        ),
        (
            "relative-change-forward",
            &base,
            &changed,
            Arguments {
                allowable_relative_change: vec![("PERCENT_DUPLICATION".to_string(), 0.095)],
                ..Arguments::default()
            },
        ),
        (
            "relative-change-reversed",
            &changed,
            &base,
            Arguments {
                allowable_relative_change: vec![("PERCENT_DUPLICATION".to_string(), 0.095)],
                ..Arguments::default()
            },
        ),
        ("missing-column", &base, &fewer, Arguments::default()),
        (
            "missing-column-not-required",
            &base,
            &fewer,
            Arguments {
                metrics_not_required: ignore(&["PERCENT_DUPLICATION"]),
                ..Arguments::default()
            },
        ),
        (
            "histogram-differs",
            &base,
            &other_histogram,
            Arguments::default(),
        ),
        (
            "histogram-ignored",
            &base,
            &other_histogram,
            Arguments {
                ignore_histogram_differences: true,
                ..Arguments::default()
            },
        ),
        ("rows-reordered", &base, &reordered, Arguments::default()),
        (
            "rows-reordered-keyed",
            &base,
            &reordered,
            Arguments {
                keys: vec!["LIBRARY".to_string()],
                ..Arguments::default()
            },
        ),
        ("missing-row", &base, &one_row, Arguments::default()),
        (
            "missing-row-keyed",
            &base,
            &one_row,
            Arguments {
                keys: vec!["LIBRARY".to_string()],
                ..Arguments::default()
            },
        ),
    ];

    let mut compared = 0;
    for (name, left, right, arguments) in cases {
        let differences = compare(left, right, &arguments);
        assert_eq!(
            exit_code(&differences),
            verdict(&corpus, name),
            "{name}: {differences:?}"
        );
        compared += 1;
    }
    assert_eq!(compared, 16, "the cases the port reproduces");
}

/// Both lists drop a column outright, whether or not both files have it.
#[test]
fn both_ignore_lists_drop_a_column() {
    let corpus = corpus();
    // The two runs that name the same differing column under the two arguments agree.
    assert_eq!(verdict(&corpus, "ignored-column"), 0);
    assert_eq!(verdict(&corpus, "not-required-but-present"), 0);
    assert_eq!(verdict(&corpus, "one-column-differs"), 1);
    let ignored = Arguments {
        metrics_to_ignore: ignore(&["PERCENT_DUPLICATION"]),
        ..Arguments::default()
    };
    let not_required = Arguments {
        metrics_not_required: ignore(&["PERCENT_DUPLICATION"]),
        ..Arguments::default()
    };
    assert!(!ignored.compares("PERCENT_DUPLICATION"));
    assert!(!not_required.compares("PERCENT_DUPLICATION"));
    assert!(ignored.compares("READ_PAIRS_EXAMINED"));
    // And a column one file lacks is forgiven by not-required too.
    assert_eq!(verdict(&corpus, "missing-column"), 1);
    assert_eq!(verdict(&corpus, "missing-column-not-required"), 0);
}

/// The same absolute difference is forgiven one way round and not the other.
#[test]
fn the_relative_change_is_relative_to_the_first_file() {
    let corpus = corpus();
    // 0.1 against 0.11 is a change of 0.1; the other way round it is about 0.0909.
    assert!(!values_agree("0.1", "0.11", Some(0.095)));
    assert!(values_agree("0.11", "0.1", Some(0.095)));
    assert_eq!(verdict(&corpus, "relative-change-forward"), 1);
    assert_eq!(verdict(&corpus, "relative-change-reversed"), 0);
    // A generous tolerance forgives either ordering and a tight one neither.
    assert!(values_agree("0.1", "0.11", Some(0.2)));
    assert!(!values_agree("0.1", "0.11", Some(0.05)));
    // Equal values agree with no tolerance at all.
    assert!(values_agree("0.1", "0.1", None));
    assert!(!values_agree("0.1", "0.11", None));
    // A value that is not a number is compared as text.
    assert!(values_agree("libA", "libA", Some(1.0)));
    assert!(!values_agree("libA", "libB", Some(1.0)));
    // The argument is a colon-separated pair.
    assert_eq!(
        parse_allowable_relative_change("PERCENT_DUPLICATION:0.2"),
        Some(("PERCENT_DUPLICATION".to_string(), 0.2))
    );
    assert_eq!(parse_allowable_relative_change("no-colon"), None);
}

/// The histogram and the table are forgiven separately.
#[test]
fn the_histogram_is_forgiven_on_its_own() {
    let corpus = corpus();
    assert_eq!(verdict(&corpus, "histogram-differs"), 1);
    assert_eq!(verdict(&corpus, "histogram-ignored"), 0);
    // Ignoring the histogram does not forgive a table difference.
    let base = input(&corpus, "base");
    let changed = input(&corpus, "changed");
    let arguments = Arguments {
        ignore_histogram_differences: true,
        ..Arguments::default()
    };
    let differences = compare(&base, &changed, &arguments);
    assert!(!differences.is_empty());
    assert!(differences
        .iter()
        .all(|difference| !matches!(difference, Difference::Histogram)));
}

/// By a column rather than by position, and a missing row is a difference either way.
#[test]
fn the_key_matches_rows_rather_than_positions() {
    let corpus = corpus();
    assert_eq!(verdict(&corpus, "rows-reordered"), 1);
    assert_eq!(verdict(&corpus, "rows-reordered-keyed"), 0);
    assert_eq!(verdict(&corpus, "missing-row"), 1);
    assert_eq!(verdict(&corpus, "missing-row-keyed"), 1);
    // The keyed comparison of a file missing a row names the row it could not find.
    let base = input(&corpus, "base");
    let one_row = input(&corpus, "one-row");
    let keyed = Arguments {
        keys: vec!["LIBRARY".to_string()],
        ..Arguments::default()
    };
    let differences = compare(&base, &one_row, &keyed);
    assert!(differences.iter().any(|difference| matches!(
        difference,
        Difference::MissingRow { key } if key == &vec!["libB".to_string()]
    )));
}

/// The two words the report uses, and the header it opens with.
#[test]
fn the_report_says_equal_or_not_equal() {
    assert_eq!(status(&[]), "equal");
    assert_eq!(status(&[Difference::Histogram]), "NOT equal");
    assert_eq!(exit_code(&[]), 0);
    assert_eq!(exit_code(&[Difference::Histogram]), 1);
    let header = report_header("picard.sam.DuplicationMetrics", "a", "b", &[]);
    assert!(header
        .starts_with("Comparison of picard.sam.DuplicationMetrics metrics between files a and b"));
    assert!(header.ends_with("Metrics are equal"));
    // Which is what the corpus's own reports say.
    let corpus = corpus();
    let identical = field(&corpus, "report", "identical").expect("its report");
    assert!(identical.contains("Metrics are equal"), "{identical}");
    let differing = field(&corpus, "report", "one-column-differs").expect("its report");
    assert!(differing.contains("Metrics are NOT equal"), "{differing}");
}

/// Neither is a difference: the reader throws before the comparison begins.
#[test]
fn a_different_class_and_a_broken_file_are_exceptions() {
    let corpus = corpus();
    for name in ["different-classes", "not-a-metrics-file"] {
        let line = corpus
            .lines()
            .find(|line| line.starts_with(&format!("error\t{name}\t")))
            .unwrap_or_else(|| panic!("the corpus carries error/{name}"));
        let message = &line[format!("error\t{name}\t").len()..];
        assert!(
            message.starts_with("htsjdk.samtools.SAMException:"),
            "{message}"
        );
        // And no verdict was reached at all.
        assert!(!corpus
            .lines()
            .any(|line| line.starts_with(&format!("verdict\t{name}\t"))));
    }
    // The two classes differ in the corpus's own inputs, which is what the reader trips on.
    let base = input(&corpus, "base");
    let other = input(&corpus, "other-class");
    assert_ne!(base.metric_class, other.metric_class);
    assert_eq!(base.columns, other.columns);
}
