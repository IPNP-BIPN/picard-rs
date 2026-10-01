//! `CreateBafRegressMetricsFile` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.CreateBafRegressMetricsFile.doWork` at tag 3.4.0. The header is compared
//! as one whole string (single tabs, exactly); each row is `line.split("\\s+")`, which keeps an
//! empty first field when the line opens on whitespace and drops trailing ones, so a row ending
//! in a tab has seven fields and a row starting with a space has eight. `LOG10_PVAL` is derived,
//! `Math.log10(PVAL)`, so a p-value of zero is negative infinity, which the metrics file writes as
//! `-?`.
//!
//! The refusals are four exceptions, all uncaught: an unrecognised header (`PicardException`,
//! the line in single quotes), a wrong field count (an `IOException`, which the tool's own catch
//! wraps as "Error parsing bafRegress Output"), a field `parseDouble` or `parseInt` refuses (a
//! `NumberFormatException`, which that catch does not take), and an empty file (`equals` called on
//! the reader's null).

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::java_number::{parse_double, parse_int};
use picard_analysis::metrics_cli::Args;
use picard_analysis::vcf_io::die;

const TOOL: &str = "CreateBafRegressMetricsFile";
const FILE_EXTENSION: &str = "bafregress_metrics";
const HEADER: &str = "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom";

struct Metrics {
    sample: String,
    estimate: f64,
    stderr: f64,
    tval: f64,
    pval: f64,
    log10_pval: f64,
    call_rate: f64,
    nhom: i32,
}

impl MetricBean for Metrics {
    fn class_name(&self) -> &str {
        "picard.arrays.BafRegressMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "SAMPLE",
            "ESTIMATE",
            "STDERR",
            "TVAL",
            "PVAL",
            "LOG10_PVAL",
            "CALL_RATE",
            "NHOM",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.sample.clone()),
            Value::Double(self.estimate),
            Value::Double(self.stderr),
            Value::Double(self.tval),
            Value::Double(self.pval),
            Value::Double(self.log10_pval),
            Value::Double(self.call_rate),
            Value::Long(i64::from(self.nhom)),
        ]
    }
}

/// `BufferedReader.readLine`, line by line: a line ends at `\n`, `\r` or `\r\n`.
fn java_lines(text: &str) -> Vec<&str> {
    let bytes = text.as_bytes();
    let mut lines = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'\n' => {
                lines.push(&text[start..i]);
                start = i + 1;
            }
            b'\r' => {
                lines.push(&text[start..i]);
                if bytes.get(i + 1) == Some(&b'\n') {
                    i += 1;
                }
                start = i + 1;
            }
            _ => {}
        }
        i += 1;
    }
    if start < bytes.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// `String.split("\\s+")`: a leading empty field is kept when the line opens on whitespace, and
/// trailing empty fields are removed. An empty line is one empty field.
fn java_split_whitespace(line: &str) -> Vec<&str> {
    if line.is_empty() {
        return vec![""];
    }
    let is_space = |c: char| matches!(c, ' ' | '\t' | '\n' | '\u{0B}' | '\u{0C}' | '\r');
    let mut fields: Vec<&str> = Vec::new();
    let mut rest = line;
    loop {
        match rest.find(is_space) {
            Some(at) => {
                fields.push(&rest[..at]);
                let after = &rest[at..];
                let skip = after.len() - after.trim_start_matches(is_space).len();
                rest = &after[skip..];
            }
            None => {
                fields.push(rest);
                break;
            }
        }
    }
    while fields.len() > 1 && fields.last() == Some(&"") {
        fields.pop();
    }
    // A line that is all whitespace splits into nothing at all.
    if fields == [""] {
        return Vec::new();
    }
    fields
}

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT")]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");

    let text = match std::fs::read(&input) {
        Ok(bytes) => String::from_utf8_lossy(&bytes).into_owned(),
        Err(_) => die(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
            absolute(&input)
        )),
    };
    let metrics_path = format!("{output}.{FILE_EXTENSION}");

    let lines = java_lines(&text);
    let mut lines = lines.into_iter();
    match lines.next() {
        None => die(
            "java.lang.NullPointerException: Cannot invoke \"String.equals(Object)\" because \
             \"line\" is null",
        ),
        Some(line) if line != HEADER => die(&format!(
            "picard.PicardException: Unrecognized header line: '{line}' in {}",
            absolute(&input)
        )),
        Some(_) => {}
    }

    let number =
        |message: String| -> ! { die(&format!("java.lang.NumberFormatException: {message}")) };
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for line in lines {
        let entries = java_split_whitespace(line);
        if entries.len() != 7 {
            die(&format!(
                "picard.PicardException: Error parsing bafRegress Output\nCaused by: \
                 java.io.IOException: Invalid number of entries ({}) in line: {line}",
                entries.len()
            ));
        }
        let double = |text: &str| parse_double(text).unwrap_or_else(|m| number(m));
        // In the reference's order, which decides the message when two fields are bad.
        let estimate = double(entries[1]);
        let stderr = double(entries[2]);
        let tval = double(entries[3]);
        let pval = double(entries[4]);
        let call_rate = double(entries[5]);
        let nhom = parse_int(entries[6]).unwrap_or_else(|m| number(m));
        let metrics = Metrics {
            sample: entries[0].to_string(),
            estimate,
            stderr,
            tval,
            pval,
            log10_pval: pval.log10(),
            call_rate,
            nhom,
        };
        file.add_metric(&metrics);
    }
    if let Err(e) = std::fs::write(&metrics_path, file.write()) {
        die(&format!(
            "htsjdk.samtools.SAMException: Could not write metrics file: {e}"
        ));
    }
}
