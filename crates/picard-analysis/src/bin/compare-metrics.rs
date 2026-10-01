//! `CompareMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CompareMetrics` at tag 3.4.0, with the parts of htsjdk 4.2.0 it reads
//! through: `MetricsFile.read`, `FormatUtil`'s parsing and formatting, `MetricBase.equals` and
//! `Histogram.equals`.
//!
//! # The reader is typed
//!
//! `MetricsFile.read` loads the class the file names and parses every cell into its field's Java
//! type, and the comparison is over those values, not over the text: a `double` column compares
//! as a number and prints through `Double.toString` (`200` in the file is `200.0` in a
//! difference), an `int` prints as itself, an enum compares by name, and an empty cell is `null`.
//! So the reader needs each class's fields, in `getFields()` order -- the class's own, then its
//! superclasses' -- which is the order differences are reported in. [`CLASSES`] holds them.
//!
//! # What decides the verdict
//!
//! * `MetricBase.equals` compares every field by its FORMATTED value (six decimals), so two files
//!   whose doubles differ past the sixth decimal are equal before any tolerance is consulted.
//! * Then, in order: the row counts, the classes, the column sets less `METRICS_NOT_REQUIRED`, and
//!   the names in `METRICS_TO_IGNORE` and the tolerances (which must be columns of the FIRST file).
//! * Values are compared field by field: numbers by their difference, whose relative change is
//!   taken against the FIRST file's value (`Double.MAX_VALUE` when that is zero) and forgiven by a
//!   tolerance only when strictly below it; anything else by `equals`.
//! * With `KEY`, rows are matched by those fields' values, and a key is printed with
//!   `StringUtil.join`, which calls `toString()` on each element: a `null` among them is a
//!   `NullPointerException`, whose message the tool rethrows as its own.
//!
//! The verdict is the exit code: 0 for equal, 1 for not, with the report written either way.

use std::collections::{BTreeMap, HashMap};

use htsjdk_metrics::file::{MetricBean, MetricsFile as MetricsWriter, Value as CellValue};
use htsjdk_metrics::format::{format_double, format_long};
use picard_analysis::java_hash_map::string_hash_code;
use picard_analysis::metrics_cli::{absolute, refuse_validation, thrown};

const TOOL: &str = "CompareMetrics";

/// The Java type of a metrics field.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Double,
    Int,
    Long,
    Str,
    /// An enum, by the names it accepts.
    Enum(&'static [&'static str]),
}

/// The metric classes this reader knows, each with its public fields in `getFields()` order.
const CLASSES: &[(&str, &[(&str, Kind)])] = &[(
    "picard.analysis.InsertSizeMetrics",
    &[
        ("MEDIAN_INSERT_SIZE", Kind::Double),
        ("MODE_INSERT_SIZE", Kind::Double),
        ("MEDIAN_ABSOLUTE_DEVIATION", Kind::Double),
        ("MIN_INSERT_SIZE", Kind::Int),
        ("MAX_INSERT_SIZE", Kind::Int),
        ("MEAN_INSERT_SIZE", Kind::Double),
        ("STANDARD_DEVIATION", Kind::Double),
        ("READ_PAIRS", Kind::Long),
        ("PAIR_ORIENTATION", Kind::Enum(&["FR", "RF", "TANDEM"])),
        ("WIDTH_OF_10_PERCENT", Kind::Int),
        ("WIDTH_OF_20_PERCENT", Kind::Int),
        ("WIDTH_OF_30_PERCENT", Kind::Int),
        ("WIDTH_OF_40_PERCENT", Kind::Int),
        ("WIDTH_OF_50_PERCENT", Kind::Int),
        ("WIDTH_OF_60_PERCENT", Kind::Int),
        ("WIDTH_OF_70_PERCENT", Kind::Int),
        ("WIDTH_OF_80_PERCENT", Kind::Int),
        ("WIDTH_OF_90_PERCENT", Kind::Int),
        ("WIDTH_OF_95_PERCENT", Kind::Int),
        ("WIDTH_OF_99_PERCENT", Kind::Int),
        // MultilevelMetrics.
        ("SAMPLE", Kind::Str),
        ("LIBRARY", Kind::Str),
        ("READ_GROUP", Kind::Str),
    ],
)];

/// A field's value, boxed the way reflection hands it over.
#[derive(Clone, Debug, PartialEq)]
enum Value {
    Null,
    Double(f64),
    Int(i32),
    Long(i64),
    Str(String),
    Enum(String),
}

impl Value {
    /// `String.valueOf(Object)`, which is what string concatenation prints.
    fn to_java_string(&self) -> String {
        match self {
            Value::Null => "null".to_string(),
            Value::Double(v) => java_double_to_string(*v),
            Value::Int(v) => v.to_string(),
            Value::Long(v) => v.to_string(),
            Value::Str(s) | Value::Enum(s) => s.clone(),
        }
    }

    /// `FormatUtil.format(Object)`.
    fn format(&self) -> String {
        match self {
            Value::Null => String::new(),
            Value::Double(v) => format_double(*v),
            Value::Int(v) => format_long(i64::from(*v)),
            Value::Long(v) => format_long(*v),
            Value::Str(s) | Value::Enum(s) => s.clone(),
        }
    }

    fn as_number(&self) -> Option<f64> {
        match self {
            Value::Double(v) => Some(*v),
            Value::Int(v) => Some(f64::from(*v)),
            Value::Long(v) => Some(*v as f64),
            _ => None,
        }
    }

    /// `Object.equals` between two values of one field.
    fn java_equals(&self, other: &Value) -> bool {
        match (self, other) {
            // `Double.equals` compares the bits, so NaN equals NaN and 0.0 does not equal -0.0.
            (Value::Double(a), Value::Double(b)) => a.to_bits() == b.to_bits(),
            _ => self == other,
        }
    }
}

/// `Double.toString`: plain between 10^-3 and 10^7, computerized scientific notation outside.
fn java_double_to_string(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let magnitude = value.abs();
    if magnitude == 0.0 {
        return format!("{sign}0.0");
    }
    let scientific = format!("{magnitude:e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if (1e-3..1e7).contains(&magnitude) {
        if exponent >= 0 {
            let whole_len = exponent as usize + 1;
            let mut padded = digits.clone();
            while padded.len() < whole_len {
                padded.push('0');
            }
            let (whole, fraction) = padded.split_at(whole_len);
            let fraction = if fraction.is_empty() { "0" } else { fraction };
            format!("{sign}{whole}.{fraction}")
        } else {
            let zeros = "0".repeat((-exponent - 1) as usize);
            format!("{sign}0.{zeros}{digits}")
        }
    } else {
        let (first, rest) = digits.split_at(1);
        let rest = if rest.is_empty() { "0" } else { rest };
        format!("{sign}{first}.{rest}E{exponent}")
    }
}

/// `Double.parseDouble`.
fn java_parse_double(text: &str) -> Option<f64> {
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    let body = trimmed.trim_start_matches(['+', '-']);
    match body {
        "NaN" => return Some(f64::NAN),
        "Infinity" => {
            return Some(if trimmed.starts_with('-') {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            })
        }
        _ => {}
    }
    let numeric = trimmed
        .strip_suffix(['d', 'D', 'f', 'F'])
        .unwrap_or(trimmed);
    if numeric.is_empty()
        || !numeric
            .chars()
            .all(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-'))
    {
        return None;
    }
    numeric.parse::<f64>().ok()
}

/// `String.split(regex)` for a literal one-character separator: trailing empty strings dropped.
fn java_split<'a>(text: &'a str, separator: &str) -> Vec<&'a str> {
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.len() > 1 && parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    if parts.len() == 1 && parts[0].is_empty() && !text.is_empty() {
        parts.clear();
    }
    parts
}

/// A failure inside `doWork`'s try block, which it rethrows as a `PicardException` carrying the
/// cause's message.
struct Failure(String);

/// One histogram: its labels and its bins, keyed by an integer.
#[derive(Debug)]
struct Histogram {
    bin_label: String,
    value_label: String,
    bins: BTreeMap<i64, f64>,
}

impl PartialEq for Histogram {
    fn eq(&self, other: &Self) -> bool {
        // `Bin.equals`: `Double.compare(value, other) == 0`.
        self.bin_label == other.bin_label
            && self.value_label == other.value_label
            && self.bins.len() == other.bins.len()
            && self
                .bins
                .iter()
                .zip(&other.bins)
                .all(|((ka, va), (kb, vb))| ka == kb && va.total_cmp(vb).is_eq())
    }
}

/// `MetricsFile` as `read` leaves it.
struct MetricsFile {
    class: Option<&'static (&'static str, &'static [(&'static str, Kind)])>,
    column_labels: Vec<String>,
    rows: Vec<Vec<Value>>,
    histograms: Vec<Histogram>,
}

impl MetricsFile {
    fn fields(&self) -> &'static [(&'static str, Kind)] {
        self.class.map_or(&[], |c| c.1)
    }

    fn class_name(&self) -> &'static str {
        self.class.map_or("", |c| c.0)
    }
}

/// `BufferedReader.readLine` over a whole file.
fn read_lines(text: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let bytes = text.as_bytes();
    let (mut start, mut i) = (0, 0);
    while i < bytes.len() {
        if bytes[i] == b'\n' || bytes[i] == b'\r' {
            out.push(&text[start..i]);
            if bytes[i] == b'\r' && i + 1 < bytes.len() && bytes[i + 1] == b'\n' {
                i += 1;
            }
            start = i + 1;
        }
        i += 1;
    }
    if start < bytes.len() {
        out.push(&text[start..]);
    }
    out
}

/// `FormatUtil.parseObject` into a field's type.
fn parse_value(text: &str, kind: Kind) -> Result<Value, Failure> {
    let number = || Failure(format!("For input string: \"{text}\""));
    Ok(match kind {
        Kind::Double => {
            if text == "?" || text == "-?" {
                Value::Double(f64::NAN)
            } else {
                Value::Double(java_parse_double(text).ok_or_else(number)?)
            }
        }
        Kind::Int => Value::Int(text.parse().map_err(|_| number())?),
        Kind::Long => Value::Long(text.parse().map_err(|_| number())?),
        Kind::Str => Value::Str(text.to_string()),
        Kind::Enum(names) => {
            if !names.contains(&text) {
                return Err(Failure(format!("No enum constant for value {text}")));
            }
            Value::Enum(text.to_string())
        }
    })
}

/// `MetricsFile.read`.
fn read_metrics(path: &str) -> Result<MetricsFile, Failure> {
    let text = std::fs::read(path)
        .map(|bytes| String::from_utf8_lossy(&bytes).into_owned())
        .map_err(|e| Failure(format!("{path} ({e})")))?;
    let lines = read_lines(&text);
    let mut at = 0;
    let next = |at: &mut usize| -> Option<&str> {
        let line = lines.get(*at).copied();
        *at += 1;
        line
    };
    let mut file = MetricsFile {
        class: None,
        column_labels: Vec::new(),
        rows: Vec::new(),
        histograms: Vec::new(),
    };

    // The headers.
    let mut line: Option<&str>;
    let mut pending_header = false;
    loop {
        line = next(&mut at);
        let Some(current) = line else { break };
        let trimmed = current.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("## METRICS CLASS") || trimmed.starts_with("## HISTOGRAM") {
            break;
        }
        if trimmed.starts_with("## ") {
            if pending_header {
                return Err(Failure(
                    "Consecutive header class lines encountered.".to_string(),
                ));
            }
            pending_header = true;
        } else if trimmed.starts_with("# ") {
            if !pending_header {
                return Err(Failure(format!(
                    "Header class must precede header value:{trimmed}"
                )));
            }
            pending_header = false;
        } else {
            return Err(Failure(format!(
                "Illegal state. Found following string in metrics file header: {trimmed}"
            )));
        }
    }
    while let Some(current) = line {
        if current.trim().starts_with("## ") {
            break;
        }
        line = next(&mut at);
    }

    // The metrics.
    if let Some(current) = line {
        let current = current.trim();
        if current.starts_with("## METRICS CLASS") {
            let class_name = current.split('\t').nth(1).unwrap_or_default();
            let class = CLASSES
                .iter()
                .find(|(name, _)| *name == class_name)
                .ok_or_else(|| Failure(format!("Could not locate class with name {class_name}")))?;
            file.class = Some(class);
            let labels = java_split(next(&mut at).unwrap_or_default(), "\t");
            let mut kinds = Vec::new();
            for label in &labels {
                let kind = class
                    .1
                    .iter()
                    .find(|(name, _)| name == label)
                    .map(|(_, kind)| *kind)
                    .ok_or_else(|| {
                        Failure(format!(
                            "Could not get field with name {label} from class {class_name}"
                        ))
                    })?;
                kinds.push(kind);
                if !file.column_labels.iter().any(|l| l == label) {
                    file.column_labels.push(label.to_string());
                }
            }
            loop {
                line = next(&mut at);
                let Some(row) = line else { break };
                if row.trim().is_empty() {
                    break;
                }
                let values: Vec<&str> = row.split('\t').collect();
                // Fields the file does not name keep their initial value: 0, or null.
                let mut bean: Vec<Value> = class
                    .1
                    .iter()
                    .map(|(_, kind)| match kind {
                        Kind::Double => Value::Double(0.0),
                        Kind::Int => Value::Int(0),
                        Kind::Long => Value::Long(0),
                        _ => Value::Null,
                    })
                    .collect();
                for (i, (label, kind)) in labels.iter().zip(&kinds).enumerate() {
                    let cell = values.get(i).ok_or_else(|| {
                        Failure(format!(
                            "Index {i} out of bounds for length {}",
                            values.len()
                        ))
                    })?;
                    let value = if cell.is_empty() {
                        if matches!(kind, Kind::Double | Kind::Int | Kind::Long) {
                            return Err(Failure(format!(
                                "Error setting field {label} on class of type {class_name}"
                            )));
                        }
                        Value::Null
                    } else {
                        parse_value(cell, *kind)?
                    };
                    let index = class
                        .1
                        .iter()
                        .position(|(name, _)| name == label)
                        .unwrap_or(0);
                    bean[index] = value;
                }
                file.rows.push(bean);
            }
        }
    }
    while let Some(current) = line {
        if current.trim().starts_with("## ") {
            break;
        }
        line = next(&mut at);
    }

    // The histograms.
    if let Some(current) = line {
        let current = current.trim();
        if current.starts_with("## HISTOGRAM") {
            let key_class = current.split('\t').nth(1).unwrap_or_default().trim();
            if !matches!(key_class, "java.lang.Integer" | "java.lang.Long") {
                return Err(Failure(format!(
                    "Could not load class with name {key_class}"
                )));
            }
            let labels = java_split(next(&mut at).unwrap_or_default(), "\t");
            for value_label in labels.iter().skip(1) {
                file.histograms.push(Histogram {
                    bin_label: labels[0].to_string(),
                    value_label: value_label.to_string(),
                    bins: BTreeMap::new(),
                });
            }
            while let Some(row) = next(&mut at) {
                if row.is_empty() {
                    break;
                }
                let fields = java_split(row.trim(), "\t");
                let key: i64 = fields[0]
                    .parse()
                    .map_err(|_| Failure(format!("For input string: \"{}\"", fields[0])))?;
                for (i, field) in fields.iter().enumerate().skip(1) {
                    let value = if *field == "?" || *field == "-?" {
                        f64::NAN
                    } else {
                        java_parse_double(field)
                            .ok_or_else(|| Failure(format!("For input string: \"{field}\"")))?
                    };
                    let histogram = file.histograms.get_mut(i - 1).ok_or_else(|| {
                        Failure(format!(
                            "Index {} out of bounds for length {}",
                            i - 1,
                            labels.len() - 1
                        ))
                    })?;
                    *histogram.bins.entry(key).or_insert(0.0) += value;
                }
            }
        }
    }
    Ok(file)
}

/// `MetricBase.equals` over two rows of one class: field by field, by formatted value.
fn rows_equal(a: &[Value], b: &[Value]) -> bool {
    a.iter().zip(b).all(|(lhs, rhs)| match (lhs, rhs) {
        (Value::Null, Value::Null) => true,
        (Value::Null, _) | (_, Value::Null) => false,
        _ => lhs.format() == rhs.format(),
    })
}

fn metrics_equal(a: &MetricsFile, b: &MetricsFile) -> bool {
    a.rows.len() == b.rows.len()
        && (a.rows.is_empty() || a.class_name() == b.class_name())
        && a.rows.iter().zip(&b.rows).all(|(x, y)| rows_equal(x, y))
}

/// `HashMap.hash` folded into a table of `capacity` buckets.
fn bucket(name: &str, capacity: usize) -> usize {
    let h = string_hash_code(name) as u32;
    ((h ^ (h >> 16)) as usize) & (capacity - 1)
}

/// `HashMap.tableSizeFor`.
fn table_size_for(capacity: usize) -> usize {
    capacity.next_power_of_two().max(1)
}

/// One row of `OUTPUT_TABLE`.
struct Difference {
    key: String,
    metric: String,
    value1: Value,
    value2: Value,
}

impl MetricBean for Difference {
    fn class_name(&self) -> &str {
        "picard.analysis.CompareMetrics$MetricComparisonDifferences"
    }

    fn columns(&self) -> &[&'static str] {
        &["KEY", "METRIC", "VALUE1", "VALUE2"]
    }

    fn values(&self) -> Vec<CellValue> {
        vec![
            CellValue::Str(self.key.clone()),
            CellValue::Str(self.metric.clone()),
            CellValue::Str(self.value1.format()),
            CellValue::Str(self.value2.format()),
        ]
    }
}

struct Comparison<'a> {
    ignore_histogram_differences: bool,
    metrics_to_ignore: &'a [String],
    metrics_not_required: &'a [String],
    allowable: &'a HashMap<String, f64>,
    keys: &'a [String],
    differences: Vec<String>,
    value_differences: Vec<Difference>,
    metric_class_name: String,
}

impl Comparison<'_> {
    /// `compareMetricValues`: whether two values agree, and what to say when they do not.
    fn compare_values(&self, value1: &Value, value2: &Value, metric: &str) -> (bool, String) {
        if matches!(value1, Value::Null) || matches!(value2, Value::Null) {
            if value1 != value2 {
                return (false, "One of the values is null".to_string());
            }
            return (true, String::new());
        }
        if let (Some(v1), Some(v2)) = (value1.as_number(), value2.as_number()) {
            let mut absolute_change = 0.0;
            if !v1.is_nan() || !v2.is_nan() {
                absolute_change = v2 - v1;
            }
            if absolute_change != 0.0 {
                let relative_change = if v1 == 0.0 {
                    f64::MAX
                } else {
                    absolute_change / v1
                };
                let changed = format!(
                    "Changed by {} (relative change of {})",
                    java_double_to_string(absolute_change),
                    java_double_to_string(relative_change)
                );
                return match self.allowable.get(metric) {
                    Some(&allowable) => {
                        let outside = relative_change.abs() >= allowable;
                        let side = if outside { "outside" } else { "within" };
                        (
                            !outside,
                            format!(
                                "{changed} which is {side} of the allowable relative change \
                                 tolerance of {}",
                                java_double_to_string(allowable)
                            ),
                        )
                    }
                    None => (false, changed),
                };
            }
            return (true, String::new());
        }
        (value1.java_equals(value2), String::new())
    }

    /// `compareMetricsForEntry`.
    fn compare_rows(
        &mut self,
        fields: &[(&str, Kind)],
        row1: &[Value],
        row2: &[Value],
        ignore: &[String],
        key: &str,
    ) -> bool {
        let mut differ = false;
        for (i, (name, _)) in fields.iter().enumerate() {
            if ignore.iter().any(|m| m == name) {
                continue;
            }
            let (equal, description) = self.compare_values(&row1[i], &row2[i], name);
            if !equal {
                differ = true;
                self.differences.push(format!(
                    "Key: {key} Metric: {name} values differ. Value1: {} Value2: {} {description}",
                    row1[i].to_java_string(),
                    row2[i].to_java_string()
                ));
                self.value_differences.push(Difference {
                    key: key.to_string(),
                    metric: name.to_string(),
                    value1: row1[i].clone(),
                    value2: row2[i].clone(),
                });
            }
        }
        differ
    }

    /// `StringUtil.join(",", key)`, which calls `toString()` on each element.
    fn join_key(key: &[Value]) -> Result<String, Failure> {
        let mut parts = Vec::new();
        for value in key {
            if matches!(value, Value::Null) {
                return Err(Failure(
                    "Cannot invoke \"Object.toString()\" because \"obj\" is null".to_string(),
                ));
            }
            parts.push(value.to_java_string());
        }
        Ok(parts.join(","))
    }

    /// `buildMetricsMap`: a `LinkedHashMap` from each row's key values to the row.
    fn build_map(&self, file: &MetricsFile) -> Result<Vec<(Vec<Value>, usize)>, Failure> {
        let mut indices = Vec::new();
        for key in self.keys {
            let index = file
                .fields()
                .iter()
                .position(|(name, _)| name == key)
                .ok_or_else(|| Failure(key.clone()))?;
            indices.push(index);
        }
        let mut map: Vec<(Vec<Value>, usize)> = Vec::new();
        for (row_index, row) in file.rows.iter().enumerate() {
            let key: Vec<Value> = indices.iter().map(|&i| row[i].clone()).collect();
            match map.iter_mut().find(|(k, _)| keys_equal(k, &key)) {
                Some(entry) => entry.1 = row_index,
                None => map.push((key, row_index)),
            }
        }
        Ok(map)
    }

    /// `compareMetricsFiles`.
    fn compare_files(
        &mut self,
        path1: &str,
        path2: &str,
        mf1: &MetricsFile,
        mf2: &MetricsFile,
    ) -> Result<i32, Failure> {
        self.metric_class_name = if !mf1.rows.is_empty() {
            mf1.class_name().to_string()
        } else if !mf2.rows.is_empty() {
            mf2.class_name().to_string()
        } else {
            "Unknown".to_string()
        };
        let histograms_equal = mf1.histograms == mf2.histograms;
        if metrics_equal(mf1, mf2) {
            if histograms_equal {
                return Ok(0);
            }
            if self.ignore_histogram_differences {
                self.differences.push(
                    "Metrics Histograms differ, but the 'IGNORE_HISTOGRAM_DIFFERENCES' flag is \
                     set."
                        .to_string(),
                );
                return Ok(0);
            }
            self.differences
                .push("Metrics Histograms differ".to_string());
            return Ok(1);
        }
        let (abs1, abs2) = (absolute(path1), absolute(path2));
        if mf1.rows.len() != mf2.rows.len() {
            self.differences.push(format!(
                "Number of metric rows differ between {abs1} and {abs2}"
            ));
            return Ok(1);
        }
        if mf1.class_name() != mf2.class_name() {
            return Err(Failure(format!(
                "Metrics are of differing class between {abs1} and {abs2}"
            )));
        }

        let without_not_required = |labels: &[String]| -> Vec<String> {
            labels
                .iter()
                .filter(|l| !self.metrics_not_required.contains(l))
                .cloned()
                .collect()
        };
        let columns1 = without_not_required(&mf1.column_labels);
        let columns2 = without_not_required(&mf2.column_labels);
        let only1: Vec<&String> = columns1.iter().filter(|c| !columns2.contains(c)).collect();
        let only2: Vec<&String> = columns2.iter().filter(|c| !columns1.contains(c)).collect();
        if !only1.is_empty() || !only2.is_empty() {
            // A `HashSet` copied from `columns1`, so sized for it, iterated bucket by bucket.
            let capacity = table_size_for(((columns1.len() as f32 / 0.75) as usize + 1).max(16));
            let mut missing: Vec<(usize, usize, &String)> = only1
                .iter()
                .chain(&only2)
                .enumerate()
                .map(|(order, name)| (bucket(name, capacity), order, *name))
                .collect();
            missing.sort();
            let names: Vec<&str> = missing.iter().map(|(_, _, n)| n.as_str()).collect();
            self.differences.push(format!(
                "Metric columns differ between {abs1} and {abs2} ({})",
                names.join(",")
            ));
            return Ok(1);
        }

        let validate_names = |names: Vec<&String>| -> Result<(), Failure> {
            let missing: Vec<&str> = names
                .into_iter()
                .filter(|n| !mf1.column_labels.contains(n))
                .map(String::as_str)
                .collect();
            if missing.is_empty() {
                return Ok(());
            }
            let mut unique: Vec<&str> = Vec::new();
            for name in missing {
                if !unique.contains(&name) {
                    unique.push(name);
                }
            }
            unique.sort_by_key(|n| (bucket(n, 16), 0));
            Err(Failure(format!(
                "Metric(s) of the name: {} were not found in {abs1}",
                unique.join(", ")
            )))
        };
        validate_names(self.metrics_to_ignore.iter().collect())?;
        validate_names(self.allowable.keys().collect())?;

        let ignore: Vec<String> = self
            .metrics_to_ignore
            .iter()
            .chain(self.metrics_not_required)
            .cloned()
            .collect();
        let fields = mf1.fields();
        let mut result = 0;
        if self.keys.is_empty() {
            for (row, (row1, row2)) in mf1.rows.iter().zip(&mf2.rows).enumerate() {
                if self.compare_rows(fields, row1, row2, &ignore, &row.to_string()) {
                    result = 1;
                }
            }
        } else {
            let map1 = self.build_map(mf1)?;
            let mut map2 = self.build_map(mf2)?;
            for (key, index1) in &map1 {
                let found = map2.iter().position(|(k, _)| keys_equal(k, key));
                match found {
                    Some(at) => {
                        let (_, index2) = map2.remove(at);
                        let joined = Self::join_key(key)?;
                        if self.compare_rows(
                            fields,
                            &mf1.rows[*index1],
                            &mf2.rows[index2],
                            &ignore,
                            &joined,
                        ) {
                            result = 1;
                        }
                    }
                    None => {
                        let joined = Self::join_key(key)?;
                        self.differences
                            .push(format!("KEY {joined} found in {path1} but not in {path2}"));
                        result = 1;
                    }
                }
            }
            for (key, _) in &map2 {
                let joined = Self::join_key(key)?;
                self.differences
                    .push(format!("KEY {joined} found in {path2} but not in {path1}"));
                result = 1;
            }
        }

        if !self.ignore_histogram_differences {
            if !histograms_equal {
                self.differences
                    .push("Metric Histograms differ".to_string());
            }
            if result == 0 && !histograms_equal {
                result = 1;
            }
        }
        Ok(result)
    }
}

/// `List.equals` over two keys.
fn keys_equal(a: &[Value], b: &[Value]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.java_equals(y))
}

/// The command line: each `NAME=value`, short names resolved, a collection appending and `null`
/// clearing it.
struct CommandLine {
    pairs: Vec<(String, String)>,
}

impl CommandLine {
    fn from_env() -> Self {
        const ALIASES: [(&str, &str); 7] = [
            ("I", "INPUT"),
            ("O", "OUTPUT"),
            ("MI", "METRICS_TO_IGNORE"),
            ("MNR", "METRICS_NOT_REQUIRED"),
            ("MARC", "METRIC_ALLOWABLE_RELATIVE_CHANGE"),
            ("IHD", "IGNORE_HISTOGRAM_DIFFERENCES"),
            ("R", "REFERENCE_SEQUENCE"),
        ];
        let mut pairs = Vec::new();
        for raw in std::env::args().skip(1) {
            let raw = raw.trim_start_matches('-');
            if let Some((name, value)) = raw.split_once('=') {
                let long = ALIASES
                    .iter()
                    .find(|(short, _)| *short == name)
                    .map_or(name, |(_, long)| long);
                pairs.push((long.to_string(), value.to_string()));
            }
        }
        CommandLine { pairs }
    }

    fn get(&self, name: &str) -> Option<&str> {
        self.pairs
            .iter()
            .rev()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_str())
            .filter(|v| *v != "null")
    }

    fn all(&self, name: &str) -> Vec<String> {
        let mut out = Vec::new();
        for (n, v) in &self.pairs {
            if n == name {
                if v == "null" {
                    out.clear();
                } else {
                    out.push(v.clone());
                }
            }
        }
        out
    }
}

fn main() {
    let args = CommandLine::from_env();
    let inputs = args.all("INPUT");
    if inputs.len() != 2 {
        refuse_validation(
            TOOL,
            &[format!(
                "Argument 'INPUT' was specified {} times, but must be specified exactly 2 times",
                inputs.len()
            )],
        );
    }
    let output = args.get("OUTPUT").map(str::to_string);
    let output_table = args.get("OUTPUT_TABLE").map(str::to_string);
    let metrics_to_ignore = args.all("METRICS_TO_IGNORE");
    let metrics_not_required = args.all("METRICS_NOT_REQUIRED");
    let keys = args.all("KEY");
    let ignore_histogram_differences = args
        .get("IGNORE_HISTOGRAM_DIFFERENCES")
        .is_some_and(|v| v.eq_ignore_ascii_case("true"));

    // `customCommandLineValidation`.
    let mut errors = Vec::new();
    let mut allowable: HashMap<String, f64> = HashMap::new();
    for spec in args.all("METRIC_ALLOWABLE_RELATIVE_CHANGE") {
        let pair = java_split(&spec, ":");
        if pair.len() == 2 {
            match java_parse_double(pair[1]) {
                Some(value) if value > 0.0 => {
                    allowable.insert(pair[0].to_string(), value);
                }
                Some(_) => errors.push(
                    "Value for numeric component of Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE' \
                     must be > 0.0"
                        .to_string(),
                ),
                None => errors.push(
                    "Invalid value for numeric component of Argument \
                     'METRIC_ALLOWABLE_RELATIVE_CHANGE'"
                        .to_string(),
                ),
            }
        } else {
            errors
                .push("Invalid value for Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE'".to_string());
        }
    }
    if !errors.is_empty() {
        refuse_validation(TOOL, &errors);
    }

    // `doWork`.
    for input in &inputs {
        if !std::path::Path::new(input).is_file() {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(input)
            ));
        }
    }
    let mut comparison = Comparison {
        ignore_histogram_differences,
        metrics_to_ignore: &metrics_to_ignore,
        metrics_not_required: &metrics_not_required,
        allowable: &allowable,
        keys: &keys,
        differences: Vec::new(),
        value_differences: Vec::new(),
        metric_class_name: "Unknown".to_string(),
    };
    let result = read_metrics(&inputs[0]).and_then(|mf1| {
        let mf2 = read_metrics(&inputs[1])?;
        comparison.compare_files(&inputs[0], &inputs[1], &mf1, &mf2)
    });
    let result = match result {
        Ok(result) => result,
        Err(Failure(message)) => thrown(&format!("picard.PicardException: {message}")),
    };

    let status = if result == 0 { "equal" } else { "NOT equal" };
    if let Some(output) = &output {
        let text = format!(
            "Comparison of {} metrics between files {} and {}\n\nMetrics are {status}\n\n{}",
            comparison.metric_class_name,
            absolute(&inputs[0]),
            absolute(&inputs[1]),
            comparison.differences.join("\n")
        );
        if let Err(e) = std::fs::write(output, text) {
            thrown(&format!("picard.PicardException: {e}"));
        }
    }
    if let Some(table) = &output_table {
        let mut file = MetricsWriter::new();
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
        for difference in &comparison.value_differences {
            file.add_metric(difference);
        }
        if let Err(e) = std::fs::write(table, file.write()) {
            thrown(&format!("picard.PicardException: {e}"));
        }
    }
    std::process::exit(result);
}
