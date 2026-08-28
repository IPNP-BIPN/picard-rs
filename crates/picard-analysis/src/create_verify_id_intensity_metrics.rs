//! `CreateVerifyIDIntensityContaminationMetricsFile`: VerifyIDIntensity's own stdout, parsed into
//! a Picard metrics file.
//!
//! The tool does nothing but parse, so all of it is here bar the metrics writer.
//!
//! Ported from `picard.arrays.CreateVerifyIDIntensityContaminationMetricsFile` in Picard 3.4.0.

/// `CreateVerifyIDIntensityContaminationMetricsFile.FILE_EXTENSION`.
///
/// The `--OUTPUT` argument is a BASENAME: this is appended to it, after a dot.
pub const FILE_EXTENSION: &str = "verifyidintensity_metrics";

/// `lineMatch`, on a line none of the three patterns accepts.
pub fn unrecognised_line_message(line: &str, input: &str) -> String {
    format!("Unrecognized line: {line} in {input}")
}

/// One row of the table.
#[derive(Debug, Clone, PartialEq)]
pub struct Metrics {
    pub id: i32,
    pub percent_mix: f64,
    pub log_likelihood: f64,
    pub log_likelihood_zero: f64,
}

/// Why a file was not parsed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    /// `lineMatch` refused the line, which is a `PicardException`.
    Unrecognised(String),
    /// The reader answered null and the matcher was handed it, which is a
    /// `NullPointerException` and NOT a `PicardException`: there is no guard between them.
    EndedEarly,
}

/// `^ID\s+%Mix\s+LLK\s+LLK0\s*$`.
pub fn is_header(line: &str) -> bool {
    let fields: Vec<&str> = line.split_whitespace().collect();
    fields == ["ID", "%Mix", "LLK", "LLK0"] && !line.starts_with(char::is_whitespace)
}

/// `^[-]+$`: one or more dashes, and the count is not looked at.
pub fn is_dashes(line: &str) -> bool {
    !line.is_empty() && line.bytes().all(|byte| byte == b'-')
}

/// `^(\d+)\s+([0-9]*\.?[0-9]+)\s+([-0-9]*\.?[0-9]+)\s+([-0-9]*\.?[0-9]+)\s*$`.
///
/// The id is unsigned, so a negative one is refused. The fraction may open on a dot and may not
/// carry a sign. The two likelihoods may carry one, and the pattern that allows it puts the minus
/// inside the digit class, so `1-2` would be accepted as well.
pub fn parse_row(line: &str) -> Option<Metrics> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() != 4 || line.starts_with(char::is_whitespace) {
        return None;
    }
    let unsigned = |text: &str| !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit());
    let fraction = |text: &str| {
        let (whole, rest) = text.split_once('.').unwrap_or((text, ""));
        !text.is_empty()
            && whole.bytes().all(|byte| byte.is_ascii_digit())
            && rest.bytes().all(|byte| byte.is_ascii_digit())
            && (text.contains('.') == !rest.is_empty() || !rest.is_empty())
            && text.bytes().any(|byte| byte.is_ascii_digit())
    };
    let signed = |text: &str| {
        let body: String = text.chars().filter(|c| *c != '-').collect();
        fraction(&body)
            && text
                .bytes()
                .all(|b| b.is_ascii_digit() || b == b'-' || b == b'.')
    };
    if !unsigned(fields[0]) || !fraction(fields[1]) || !signed(fields[2]) || !signed(fields[3]) {
        return None;
    }
    Some(Metrics {
        id: fields[0].parse().ok()?,
        percent_mix: fields[1].parse().ok()?,
        log_likelihood: fields[2].parse().ok()?,
        log_likelihood_zero: fields[3].parse().ok()?,
    })
}

/// `doWork`: a header, a run of dashes, then a row per line to the end.
///
/// A file that runs out before the first two lines is `EndedEarly`, which is the reference's
/// unguarded null. A file of a header and dashes and nothing else parses to no rows at all, and
/// the writer then emits its comments and stops: the metrics file has no column line either.
pub fn parse(text: &str) -> Result<Vec<Metrics>, ParseError> {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Err(ParseError::EndedEarly);
    };
    if !is_header(header) {
        return Err(ParseError::Unrecognised(header.to_string()));
    }
    let Some(dashes) = lines.next() else {
        return Err(ParseError::EndedEarly);
    };
    if !is_dashes(dashes) {
        return Err(ParseError::Unrecognised(dashes.to_string()));
    }
    let mut rows = Vec::new();
    for line in lines {
        match parse_row(line) {
            Some(row) => rows.push(row),
            None => return Err(ParseError::Unrecognised(line.to_string())),
        }
    }
    Ok(rows)
}

/// The file the tool writes, whose name is the basename and the extension.
pub fn output_name(basename: &str) -> String {
    format!("{basename}.{FILE_EXTENSION}")
}
