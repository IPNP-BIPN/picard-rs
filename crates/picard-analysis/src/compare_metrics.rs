//! `CompareMetrics`: two metrics files in, one verdict out.
//!
//! What is ported is what counts as a difference and what each of the four tolerance arguments
//! forgives. Reading a metrics file is not ported: the caller hands over the two tables.
//!
//! Ported from `picard.analysis.CompareMetrics` in Picard 3.4.0.

use std::collections::BTreeSet;

/// One metrics file, reduced to what the comparison reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetricsFile {
    pub metric_class: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
    /// The histogram, as its lines.
    pub histogram: Vec<String>,
}

impl MetricsFile {
    pub fn value(&self, row: usize, column: &str) -> Option<&str> {
        let index = self.columns.iter().position(|name| name == column)?;
        self.rows.get(row)?.get(index).map(String::as_str)
    }

    /// The rows keyed by the named columns, in the order they appear.
    pub fn keyed(&self, keys: &[String]) -> Vec<(Vec<String>, usize)> {
        self.rows
            .iter()
            .enumerate()
            .map(|(index, _)| {
                let key: Vec<String> = keys
                    .iter()
                    .map(|column| self.value(index, column).unwrap_or("").to_string())
                    .collect();
                (key, index)
            })
            .collect()
    }
}

/// The arguments that decide what is forgiven.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Arguments {
    /// `--METRICS_TO_IGNORE`.
    pub metrics_to_ignore: BTreeSet<String>,
    /// `--METRICS_NOT_REQUIRED`, which forgives a column present in BOTH files just as readily as
    /// one only a single file has: it is a second ignore list and not a laxer requirement.
    pub metrics_not_required: BTreeSet<String>,
    /// `--METRIC_ALLOWABLE_RELATIVE_CHANGE`, as name and tolerance.
    pub allowable_relative_change: Vec<(String, f64)>,
    pub ignore_histogram_differences: bool,
    /// `--KEY`, the columns rows are matched on.
    pub keys: Vec<String>,
}

impl Arguments {
    /// Whether a column takes part in the comparison at all.
    pub fn compares(&self, column: &str) -> bool {
        !self.metrics_to_ignore.contains(column) && !self.metrics_not_required.contains(column)
    }

    pub fn tolerance(&self, column: &str) -> Option<f64> {
        self.allowable_relative_change
            .iter()
            .find(|(name, _)| name == column)
            .map(|(_, tolerance)| *tolerance)
    }
}

/// `--METRIC_ALLOWABLE_RELATIVE_CHANGE`, which is a colon-separated pair.
pub fn parse_allowable_relative_change(spec: &str) -> Option<(String, f64)> {
    let (name, tolerance) = spec.split_once(':')?;
    Some((name.to_string(), tolerance.parse().ok()?))
}

/// Whether two values of one column agree.
///
/// The relative change is measured against the FIRST file's value, so the same absolute
/// difference is forgiven in one ordering and not in the other: 0.1 against 0.11 is a change of
/// 0.1 and 0.11 against 0.1 is about 0.0909.
pub fn values_agree(left: &str, right: &str, tolerance: Option<f64>) -> bool {
    if left == right {
        return true;
    }
    let Some(tolerance) = tolerance else {
        return false;
    };
    let (Ok(left), Ok(right)) = (left.parse::<f64>(), right.parse::<f64>()) else {
        return false;
    };
    if left == 0.0 {
        return right == 0.0;
    }
    ((left - right) / left).abs() <= tolerance
}

/// One difference the comparison found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Difference {
    /// A column one file has and the other does not.
    MissingColumn { column: String },
    /// A row one file has and the other does not.
    MissingRow { key: Vec<String> },
    /// Two values that disagree past the tolerance.
    Value {
        row: usize,
        column: String,
        left: String,
        right: String,
    },
    /// The histograms disagree.
    Histogram,
}

/// The whole comparison.
///
/// Rows are matched by position when no key is given and by the key's own columns when one is.
/// A row one file has and the other does not is a difference either way.
pub fn compare(left: &MetricsFile, right: &MetricsFile, arguments: &Arguments) -> Vec<Difference> {
    let mut differences = Vec::new();

    // A column one file lacks, in either direction.
    for column in left.columns.iter().chain(right.columns.iter()) {
        if !arguments.compares(column) {
            continue;
        }
        let in_both = left.columns.contains(column) && right.columns.contains(column);
        if !in_both
            && !differences.iter().any(|difference| {
                matches!(difference, Difference::MissingColumn { column: name } if name == column)
            })
        {
            differences.push(Difference::MissingColumn {
                column: column.clone(),
            });
        }
    }

    let shared: Vec<&String> = left
        .columns
        .iter()
        .filter(|column| right.columns.contains(column) && arguments.compares(column))
        .collect();

    if arguments.keys.is_empty() {
        for row in 0..left.rows.len().max(right.rows.len()) {
            if row >= left.rows.len() || row >= right.rows.len() {
                differences.push(Difference::MissingRow {
                    key: vec![row.to_string()],
                });
                continue;
            }
            for column in &shared {
                let (a, b) = (
                    left.value(row, column).unwrap_or(""),
                    right.value(row, column).unwrap_or(""),
                );
                if !values_agree(a, b, arguments.tolerance(column)) {
                    differences.push(Difference::Value {
                        row,
                        column: (*column).clone(),
                        left: a.to_string(),
                        right: b.to_string(),
                    });
                }
            }
        }
    } else {
        let right_keyed = right.keyed(&arguments.keys);
        for (key, row) in left.keyed(&arguments.keys) {
            let Some((_, other)) = right_keyed.iter().find(|(other, _)| *other == key) else {
                differences.push(Difference::MissingRow { key });
                continue;
            };
            for column in &shared {
                let (a, b) = (
                    left.value(row, column).unwrap_or(""),
                    right.value(*other, column).unwrap_or(""),
                );
                if !values_agree(a, b, arguments.tolerance(column)) {
                    differences.push(Difference::Value {
                        row,
                        column: (*column).clone(),
                        left: a.to_string(),
                        right: b.to_string(),
                    });
                }
            }
        }
        // A row the SECOND file has and the first does not.
        let left_keys: BTreeSet<Vec<String>> = left
            .keyed(&arguments.keys)
            .into_iter()
            .map(|(k, _)| k)
            .collect();
        for (key, _) in right_keyed {
            if !left_keys.contains(&key) {
                differences.push(Difference::MissingRow { key });
            }
        }
    }

    if !arguments.ignore_histogram_differences && left.histogram != right.histogram {
        differences.push(Difference::Histogram);
    }

    differences
}

/// The exit code the verdict becomes: zero when the two files agree and one when they do not.
pub fn exit_code(differences: &[Difference]) -> i32 {
    if differences.is_empty() {
        0
    } else {
        1
    }
}

/// The two words the report uses.
pub fn status(differences: &[Difference]) -> &'static str {
    if differences.is_empty() {
        "equal"
    } else {
        "NOT equal"
    }
}

/// The report's own header, which names the two files and the verdict.
pub fn report_header(
    metric_class: &str,
    left: &str,
    right: &str,
    differences: &[Difference],
) -> String {
    format!(
        "Comparison of {metric_class} metrics between files {left} and {right}\n\nMetrics are {}",
        status(differences)
    )
}
