//! Conformance for `CreateVerifyIDIntensityContaminationMetricsFile` against Picard 3.4.0.
//!
//! Each case carries the input the tool read and either the metrics table it wrote or the
//! exception it threw. The tool does nothing but parse, so the port is given the same input and
//! must accept or refuse it the same way.
//!
//! # What this suite is for
//!
//!  * **the output argument being a basename, not a file**;
//!  * **the first two lines being fixed, each matched by its own pattern**;
//!  * **the dashes not being counted**;
//!  * **the columns splitting on runs of whitespace**;
//!  * **the fraction opening on a dot but never carrying a sign**;
//!  * **the likelihoods carrying either sign**;
//!  * **the id being unsigned**;
//!  * **an unrecognised line quoting itself and naming the input**;
//!  * **a file that ends early being a NullPointerException, not a PicardException**;
//!  * **and a file of a header and dashes writing no table at all.**

use std::io::Read;

use picard_analysis::create_verify_id_intensity_metrics::{
    is_dashes, is_header, output_name, parse, parse_row, unrecognised_line_message, Metrics,
    ParseError, FILE_EXTENSION,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("create_verify_id_intensity_metrics.txt.gz");
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

fn input(text: &str, case: &str) -> String {
    field(text, "in", case).unwrap_or_else(|| panic!("{case} has an input"))
}

/// The rows a written metrics table holds, as their four values.
fn written(text: &str, case: &str) -> Vec<Vec<String>> {
    let table = field(text, "metrics", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = table.lines().filter(|line| !line.is_empty());
    match lines.next() {
        None => Vec::new(),
        Some(header) => {
            assert_eq!(header, "ID\tPCT_MIX\tLLK\tLLK0", "{case}");
            lines
                .map(|line| line.split('\t').map(str::to_string).collect())
                .collect()
        }
    }
}

/// The reference reformats on the way out, dropping a trailing `.0`.
fn number(value: f64) -> String {
    if value == value.trunc() {
        format!("{}", value as i64)
    } else {
        format!("{value}")
    }
}

fn as_written(row: &Metrics) -> Vec<String> {
    vec![
        row.id.to_string(),
        number(row.percent_mix),
        number(row.log_likelihood),
        number(row.log_likelihood_zero),
    ]
}

const ACCEPTED: &[&str] = &[
    "one-row",
    "three-rows",
    "spaces-not-tabs",
    "one-dash",
    "leading-dot",
    "positive-likelihood",
    "no-rows",
];

const REFUSED: &[&str] = &[
    "negative-fraction",
    "negative-id",
    "short-row",
    "wrong-header",
];

/// Every accepted case parses to the rows the tool wrote.
#[test]
fn every_accepted_case_parses_to_the_same_rows() {
    let text = corpus();
    for case in ACCEPTED {
        let rows = parse(&input(&text, case)).unwrap_or_else(|e| panic!("{case}: {e:?}"));
        assert_eq!(
            rows.iter().map(as_written).collect::<Vec<_>>(),
            written(&text, case),
            "{case}"
        );
    }
}

/// The output argument is a basename: the extension is appended to it.
#[test]
fn the_output_argument_is_a_basename() {
    let text = corpus();
    assert_eq!(output_name("out"), "out.verifyidintensity_metrics");
    assert_eq!(FILE_EXTENSION, "verifyidintensity_metrics");
    assert_eq!(
        field(&text, "name", "one-row").as_deref(),
        Some("out.verifyidintensity_metrics")
    );
}

/// The dashes are not counted: one is as good as forty, and the columns split on runs of
/// whitespace so single spaces parse as well as tabs.
#[test]
fn the_patterns_are_loose_where_it_does_not_matter() {
    let text = corpus();
    assert!(is_dashes("-"));
    assert!(is_dashes("----------------------------"));
    assert!(!is_dashes(""));
    assert!(!is_dashes("--x--"));
    assert_eq!(written(&text, "one-dash"), written(&text, "one-row"));
    assert!(is_header("ID\t%Mix\tLLK\tLLK0"));
    assert!(is_header("ID   %Mix   LLK   LLK0"));
    assert_eq!(written(&text, "spaces-not-tabs"), written(&text, "one-row"));
}

/// The fraction may open on a dot and may not carry a sign; the likelihoods may carry either.
#[test]
fn the_fraction_is_unsigned_and_the_likelihoods_are_not() {
    let text = corpus();
    assert_eq!(
        parse_row("0\t.5\t-1.0\t-2.0").map(|row| row.percent_mix),
        Some(0.5)
    );
    assert_eq!(written(&text, "leading-dot")[0][1], "0.5");
    assert!(parse_row("0\t0.05\t1234.5\t2345.6").is_some());
    assert_eq!(written(&text, "positive-likelihood")[0][2], "1234.5");
    assert!(parse_row("0\t-0.05\t-1.0\t-2.0").is_none());
    // The id is unsigned too.
    assert!(parse_row("-1\t0.05\t-1.0\t-2.0").is_none());
    assert!(parse_row("0\t0.05\t-1.0").is_none());
}

/// An unrecognised line quotes itself and names the input.
#[test]
fn an_unrecognised_line_quotes_itself() {
    let text = corpus();
    for case in REFUSED {
        let error = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        assert!(
            error.starts_with("picard.PicardException:"),
            "{case}: {error}"
        );
        let refused = match parse(&input(&text, case)) {
            Err(ParseError::Unrecognised(line)) => line,
            other => panic!("{case}: {other:?}"),
        };
        // The path the message ends on is the dump's temp directory, so the comparison is on the
        // part before it: the prefix and the line, quoted as it was read.
        let expected = unrecognised_line_message(&refused, "<path>");
        let prefix = expected.split(" in <path>").next().expect("the prefix");
        assert!(error.contains(prefix), "{case}: {error}");
    }
}

/// A file that ends early is a NullPointerException and not a PicardException: the reader answers
/// null and the matcher is handed it without a guard.
#[test]
fn a_file_that_ends_early_is_a_null_pointer() {
    let text = corpus();
    for case in ["header-only", "empty"] {
        let error = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        assert!(
            error.starts_with("java.lang.NullPointerException:"),
            "{case}: {error}"
        );
        assert!(error.contains("CharSequence.length()"), "{case}: {error}");
        assert_eq!(
            parse(&input(&text, case)),
            Err(ParseError::EndedEarly),
            "{case}"
        );
    }
}

/// A file of a header and dashes writes a metrics file with no table at all, not even the column
/// line: it parses to no rows, and the writer emits its comments and stops.
#[test]
fn a_file_of_no_rows_writes_no_table() {
    let text = corpus();
    assert_eq!(parse(&input(&text, "no-rows")), Ok(Vec::new()));
    assert!(field(&text, "metrics", "no-rows")
        .unwrap_or_default()
        .trim()
        .is_empty());
    // And it is a written file all the same.
    assert_eq!(
        field(&text, "name", "no-rows").as_deref(),
        Some("out.verifyidintensity_metrics")
    );
}
