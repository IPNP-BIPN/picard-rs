//! `CreateBafRegressMetricsFile`: bafRegress' own stdout, parsed into a Picard metrics file.
//!
//! The tool does nothing but parse and derive one column, so all of it is here bar the writer.
//!
//! Ported from `picard.arrays.CreateBafRegressMetricsFile` in Picard 3.4.0.

/// `CreateBafRegressMetricsFile.FILE_EXTENSION`. The `--OUTPUT` argument is a BASENAME.
pub const FILE_EXTENSION: &str = "bafregress_metrics";

/// The header, compared as a WHOLE STRING rather than by a pattern: single tabs, exactly.
pub const HEADER: &str = "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom";

/// `doWork`, on a first line that is not the header. It quotes the line in SINGLE QUOTES, which
/// the row refusal does not.
pub fn unrecognised_header_message(line: &str, input: &str) -> String {
    format!("Unrecognized header line: '{line}' in {input}")
}

/// `doWork`, on a row with the wrong number of columns. This one counts them and quotes nothing.
pub fn invalid_entry_count_message(count: usize, line: &str) -> String {
    format!("Invalid number of entries ({count}) in line: {line}")
}

/// One row of the table, with the column the tool derives.
#[derive(Debug, Clone, PartialEq)]
pub struct Metrics {
    pub sample: String,
    pub estimate: f64,
    pub standard_error: f64,
    pub t_value: f64,
    pub p_value: f64,
    /// Derived, not read: `Math.log10(PVAL)`.
    pub log10_p_value: f64,
    pub call_rate: f64,
    pub number_homozygous: i32,
}

/// Why a file was not parsed. The three are three different Java classes, which is what says the
/// tool checks the three things in three different places.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    /// A `PicardException`.
    Header(String),
    /// An `IOException`, wrapped in a `PicardException` by the caller.
    EntryCount { count: usize, line: String },
    /// A `NumberFormatException`, which escapes that wrapper and reaches the caller by itself.
    Number(String),
    /// A `NullPointerException`: the header comparison is called on the null the reader answered.
    EndedEarly,
}

/// `doWork`: the header, then a row per line to the end.
///
/// The rows split on runs of whitespace where the header does not, so a row of spaces parses
/// under a header of tabs. `LOG10_PVAL` is derived here and not read, so a p-value of zero gives
/// negative infinity, which the writer renders as `-?`.
pub fn parse(text: &str) -> Result<Vec<Metrics>, ParseError> {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Err(ParseError::EndedEarly);
    };
    if header != HEADER {
        return Err(ParseError::Header(header.to_string()));
    }
    let mut rows = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 7 {
            return Err(ParseError::EntryCount {
                count: fields.len(),
                line: line.to_string(),
            });
        }
        let number = |text: &str| {
            text.parse::<f64>()
                .map_err(|_| ParseError::Number(text.to_string()))
        };
        let p_value = number(fields[4])?;
        rows.push(Metrics {
            sample: fields[0].to_string(),
            estimate: number(fields[1])?,
            standard_error: number(fields[2])?,
            t_value: number(fields[3])?,
            p_value,
            log10_p_value: p_value.log10(),
            call_rate: number(fields[5])?,
            number_homozygous: fields[6]
                .parse()
                .map_err(|_| ParseError::Number(fields[6].to_string()))?,
        });
    }
    Ok(rows)
}

/// The file the tool writes, whose name is the basename and the extension.
pub fn output_name(basename: &str) -> String {
    format!("{basename}.{FILE_EXTENSION}")
}
