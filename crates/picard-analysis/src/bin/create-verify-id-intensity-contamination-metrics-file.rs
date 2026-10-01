//! `CreateVerifyIDIntensityContaminationMetricsFile` as a runnable binary: the covering array's
//! port side.
//!
//! Ports `picard.arrays.CreateVerifyIDIntensityContaminationMetricsFile.doWork` at tag 3.4.0.
//! The tool reads VerifyIDIntensity's stdout with a `BufferedReader` and matches each line against
//! one of three patterns: the header, then a run of dashes, then one row per line to the end.
//!
//! The refusals are three different exceptions, and all three escape uncaught:
//!
//! * a line no pattern accepts is a `PicardException` quoting the line and the absolute input;
//! * a file that ends before its dashes hands `Pattern.matcher` the reader's null, which is a
//!   `NullPointerException` from inside `Matcher`'s constructor;
//! * a row the pattern accepts and the number parsers do not -- `1-2`, which the likelihood
//!   pattern allows because its minus sits inside the digit class, or an ID past
//!   `Integer.MAX_VALUE` -- is a `NumberFormatException`, which the `catch (IOException)` around
//!   the loop does not catch.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::java_number::{parse_double, parse_int};
use picard_analysis::metrics_cli::Args;
use picard_analysis::vcf_io::die;

const TOOL: &str = "CreateVerifyIDIntensityContaminationMetricsFile";
const FILE_EXTENSION: &str = "verifyidintensity_metrics";

struct Metrics {
    id: i32,
    pct_mix: f64,
    llk: f64,
    llk0: f64,
}

impl MetricBean for Metrics {
    fn class_name(&self) -> &str {
        "picard.arrays.VerifyIDIntensityContaminationMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &["ID", "PCT_MIX", "LLK", "LLK0"]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Long(i64::from(self.id)),
            Value::Double(self.pct_mix),
            Value::Double(self.llk),
            Value::Double(self.llk0),
        ]
    }
}

/// `BufferedReader.readLine`, line by line: a line ends at `\n`, `\r` or `\r\n`, and a last line
/// without a terminator is still a line.
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

/// `\s` in a `java.util.regex` pattern without UNICODE_CHARACTER_CLASS.
fn is_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\u{0B}' | '\u{0C}' | '\r')
}

/// The line as `token (\s+ token)* \s*`, or `None` when it opens on whitespace (every pattern is
/// anchored at a token) or is empty.
fn tokens(line: &str) -> Option<Vec<&str>> {
    if line.is_empty() || line.starts_with(is_space) {
        return None;
    }
    Some(line.split(is_space).filter(|t| !t.is_empty()).collect())
}

/// `^ID\s+%Mix\s+LLK\s+LLK0\s*$`.
fn is_header(line: &str) -> bool {
    tokens(line).is_some_and(|t| t == ["ID", "%Mix", "LLK", "LLK0"])
}

/// `^[-]+$`.
fn is_dashes(line: &str) -> bool {
    !line.is_empty() && line.bytes().all(|b| b == b'-')
}

/// `[0-9]*\.?[0-9]+` over a whole token whose characters are drawn from `prefix_class` before
/// the optional dot.
fn matches_number(token: &str, prefix_class: impl Fn(u8) -> bool) -> bool {
    let bytes = token.as_bytes();
    match token.find('.') {
        Some(dot) => {
            let (before, after) = (&bytes[..dot], &bytes[dot + 1..]);
            before.iter().all(|b| prefix_class(*b))
                && !after.is_empty()
                && after.iter().all(u8::is_ascii_digit)
        }
        None => {
            !bytes.is_empty()
                && bytes.iter().all(|b| prefix_class(*b))
                && bytes.last().is_some_and(u8::is_ascii_digit)
        }
    }
}

/// `^(\d+)\s+([0-9]*\.?[0-9]+)\s+([-0-9]*\.?[0-9]+)\s+([-0-9]*\.?[0-9]+)\s*$`, as its four
/// groups.
fn data_groups(line: &str) -> Option<[&str; 4]> {
    let t = tokens(line)?;
    if t.len() != 4 {
        return None;
    }
    let digits = |token: &str| !token.is_empty() && token.bytes().all(|b| b.is_ascii_digit());
    let signed = |b: u8| b.is_ascii_digit() || b == b'-';
    if digits(t[0])
        && matches_number(t[1], |b| b.is_ascii_digit())
        && matches_number(t[2], signed)
        && matches_number(t[3], signed)
    {
        Some([t[0], t[1], t[2], t[3]])
    } else {
        None
    }
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

    let unrecognised = |line: &str| -> ! {
        die(&format!(
            "picard.PicardException: Unrecognized line: {line} in {}",
            absolute(&input)
        ))
    };
    // `Pattern.matcher(null)`.
    let null_line = || -> ! {
        die(
            "java.lang.NullPointerException: Cannot invoke \"java.lang.CharSequence.length()\" \
             because \"this.text\" is null",
        )
    };
    let number =
        |message: String| -> ! { die(&format!("java.lang.NumberFormatException: {message}")) };

    let lines = java_lines(&text);
    let mut lines = lines.into_iter();
    match lines.next() {
        None => null_line(),
        Some(line) if !is_header(line) => unrecognised(line),
        Some(_) => {}
    }
    match lines.next() {
        None => null_line(),
        Some(line) if !is_dashes(line) => unrecognised(line),
        Some(_) => {}
    }

    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for line in lines {
        let Some([id, pct_mix, llk, llk0]) = data_groups(line) else {
            unrecognised(line)
        };
        let metrics = Metrics {
            id: parse_int(id).unwrap_or_else(|m| number(m)),
            pct_mix: parse_double(pct_mix).unwrap_or_else(|m| number(m)),
            llk: parse_double(llk).unwrap_or_else(|m| number(m)),
            llk0: parse_double(llk0).unwrap_or_else(|m| number(m)),
        };
        file.add_metric(&metrics);
    }
    if let Err(e) = std::fs::write(&metrics_path, file.write()) {
        die(&format!(
            "htsjdk.samtools.SAMException: Could not write metrics file: {e}"
        ));
    }
}
