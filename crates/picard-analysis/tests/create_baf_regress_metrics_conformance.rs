//! Conformance for `CreateBafRegressMetricsFile` against Picard 3.4.0.
//!
//! Each case carries the input the tool read and either the metrics table it wrote or the
//! exception it threw. The tool parses and derives one column, and the port is given the same
//! input and must accept or refuse it the same way.
//!
//! # What this suite is for
//!
//!  * **the output argument being a basename**;
//!  * **the header being compared as a whole string, so spaces where the tabs go are refused**;
//!  * **the rows splitting on runs of whitespace all the same**;
//!  * **LOG10_PVAL being derived and not read**;
//!  * **a p-value of zero giving an infinity written `-?`**;
//!  * **negative values and exponent notation being accepted**;
//!  * **three malformed inputs raising three different classes**;
//!  * **and an empty file being a NullPointerException.**

use std::io::Read;

use picard_analysis::create_baf_regress_metrics::{
    invalid_entry_count_message, output_name, parse, unrecognised_header_message, Metrics,
    ParseError, FILE_EXTENSION, HEADER,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("create_baf_regress_metrics.txt.gz");
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

fn written(text: &str, case: &str) -> Vec<Vec<String>> {
    let table = field(text, "metrics", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = table.lines().filter(|line| !line.is_empty());
    match lines.next() {
        None => Vec::new(),
        Some(header) => {
            assert_eq!(
                header, "SAMPLE\tESTIMATE\tSTDERR\tTVAL\tPVAL\tLOG10_PVAL\tCALL_RATE\tNHOM",
                "{case}"
            );
            lines
                .map(|line| line.split('\t').map(str::to_string).collect())
                .collect()
        }
    }
}

/// `FormatUtil`: at most six fraction digits, trailing zeros dropped, and a non-finite value is a
/// question mark carrying its sign.
fn number(value: f64) -> String {
    if value.is_infinite() {
        return if value < 0.0 {
            "-?".to_string()
        } else {
            "?".to_string()
        };
    }
    let rendered = format!("{value:.6}");
    let trimmed = rendered.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

fn as_written(row: &Metrics) -> Vec<String> {
    vec![
        row.sample.clone(),
        number(row.estimate),
        number(row.standard_error),
        number(row.t_value),
        number(row.p_value),
        number(row.log10_p_value),
        number(row.call_rate),
        row.number_homozygous.to_string(),
    ]
}

const ACCEPTED: &[&str] = &[
    "one-row",
    "three-rows",
    "spaces-in-the-row",
    "zero-pvalue",
    "negative-estimate",
    "exponent-notation",
    "header-only",
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

/// The output argument is a basename.
#[test]
fn the_output_argument_is_a_basename() {
    let text = corpus();
    assert_eq!(output_name("out"), "out.bafregress_metrics");
    assert_eq!(FILE_EXTENSION, "bafregress_metrics");
    assert_eq!(
        field(&text, "name", "one-row").as_deref(),
        Some("out.bafregress_metrics")
    );
}

/// The header is compared as a whole string, so spaces where the tabs go are refused, while the
/// rows split on runs of whitespace all the same.
#[test]
fn the_header_is_a_string_and_the_rows_are_not() {
    let text = corpus();
    assert_eq!(
        HEADER,
        "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom"
    );
    assert_eq!(
        written(&text, "spaces-in-the-row"),
        written(&text, "one-row")
    );
    for case in ["spaces-in-the-header", "wrong-header"] {
        let error = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        let refused = match parse(&input(&text, case)) {
            Err(ParseError::Header(line)) => line,
            other => panic!("{case}: {other:?}"),
        };
        assert!(
            error.starts_with("picard.PicardException:"),
            "{case}: {error}"
        );
        let expected = unrecognised_header_message(&refused, "<path>");
        let prefix = expected.split(" in <path>").next().expect("the prefix");
        assert!(error.contains(prefix), "{case}: {error}");
    }
}

/// LOG10_PVAL is derived and not read, and a p-value of zero gives an infinity written `-?`.
#[test]
fn the_logarithm_is_derived() {
    let text = corpus();
    let rows = parse(&input(&text, "one-row")).expect("parsed");
    assert_eq!(rows[0].p_value, 0.001);
    assert_eq!(rows[0].log10_p_value, -3.0);
    assert_eq!(written(&text, "one-row")[0][5], "-3");
    // The input carries seven columns and the output eight: the logarithm is the extra one.
    assert_eq!(written(&text, "one-row")[0].len(), 8);
    assert_eq!(HEADER.split('\t').count(), 7);
    // A p-value of zero.
    let zero = parse(&input(&text, "zero-pvalue")).expect("parsed");
    assert!(zero[0].log10_p_value.is_infinite());
    assert!(zero[0].log10_p_value < 0.0);
    assert_eq!(written(&text, "zero-pvalue")[0][5], "-?");
}

/// Negative values and exponent notation are accepted, the columns being parsed as doubles.
#[test]
fn the_columns_are_doubles_and_not_patterns() {
    let text = corpus();
    let negative = parse(&input(&text, "negative-estimate")).expect("parsed");
    assert_eq!(negative[0].estimate, -0.05);
    assert_eq!(written(&text, "negative-estimate")[0][1], "-0.05");
    let exponent = parse(&input(&text, "exponent-notation")).expect("parsed");
    assert_eq!(exponent[0].p_value, 1e-5);
    // And rewritten on the way out.
    assert_eq!(written(&text, "exponent-notation")[0][4], "0.00001");
    assert_eq!(written(&text, "exponent-notation")[0][5], "-5");
}

/// Three malformed inputs raise three different classes, which is what says the tool checks the
/// three things in three different places.
#[test]
fn the_three_failures_are_three_classes() {
    let text = corpus();
    // The row count is an IOException.
    for case in ["short-row", "long-row"] {
        let error = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        assert!(error.starts_with("java.io.IOException:"), "{case}: {error}");
        let (count, line) = match parse(&input(&text, case)) {
            Err(ParseError::EntryCount { count, line }) => (count, line),
            other => panic!("{case}: {other:?}"),
        };
        assert!(
            error.contains(&invalid_entry_count_message(count, &line)),
            "{case}: {error}"
        );
    }
    // A fractional Nhom is a NumberFormatException, which escapes the wrapper.
    let error = field(&text, "error", "fractional-nhom").expect("its refusal");
    assert!(
        error.starts_with("java.lang.NumberFormatException:"),
        "{error}"
    );
    assert_eq!(
        parse(&input(&text, "fractional-nhom")),
        Err(ParseError::Number("1000.5".to_string()))
    );
    assert!(error.contains("1000.5"), "{error}");
    // And the header is a PicardException, which the test above pins.
    assert!(field(&text, "error", "wrong-header")
        .expect("its refusal")
        .starts_with("picard.PicardException:"));
}

/// An empty file is a NullPointerException on `String.equals`, the header comparison being called
/// on the null the reader answered.
#[test]
fn an_empty_file_is_a_null_pointer() {
    let text = corpus();
    let error = field(&text, "error", "empty").expect("its refusal");
    assert!(
        error.starts_with("java.lang.NullPointerException:"),
        "{error}"
    );
    assert!(error.contains("String.equals(Object)"), "{error}");
    assert_eq!(parse(&input(&text, "empty")), Err(ParseError::EndedEarly));
    // A file of only a header is NOT: it parses to no rows and writes a table of none.
    assert_eq!(parse(&input(&text, "header-only")), Ok(Vec::new()));
    assert!(written(&text, "header-only").is_empty());
}
