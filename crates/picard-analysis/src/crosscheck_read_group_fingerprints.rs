//! `CrosscheckReadGroupFingerprints`: the deprecated wrapper around `CrosscheckFingerprints`.
//!
//! The tool is its parent with the arguments moved. It rolls up through two booleans of its own
//! rather than through `--CROSSCHECK_BY`, and it refuses `--CROSSCHECK_BY` outright: a combination
//! its parent's own validation would have accepted is a refusal here, because the wrapper sets
//! that argument itself and will not be told twice.
//!
//! Rolling up also moves the output. A roll-up writes a MATRIX, so the file named by `--OUTPUT`
//! becomes the matrix and the table that would have gone there is sent to `/dev/null`.
//!
//! Ported from `picard.fingerprint.CrosscheckReadGroupFingerprints` in Picard 3.4.0.

use crate::crosscheck_fingerprints::DataType;

/// What the wrapper was asked for.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Options {
    pub crosscheck_samples: bool,
    pub crosscheck_libraries: bool,
    pub expect_all_read_groups_to_match: bool,
    /// Set only by a caller that named it, which is what the refusal is about.
    pub crosscheck_by: Option<DataType>,
    pub matrix_output: Option<String>,
    pub second_input: Vec<String>,
}

/// Where a run's two files end up.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Destination {
    pub crosscheck_by: DataType,
    /// The file the table is written to.
    pub output: String,
    /// The file the matrix is written to, if there is one.
    pub matrix_output: Option<String>,
}

/// What a roll-up does to the arguments, before the parent tool sees them.
pub fn destination(output: &str, options: &Options) -> Destination {
    if options.crosscheck_libraries {
        return Destination {
            crosscheck_by: DataType::Library,
            output: "/dev/null".to_string(),
            matrix_output: Some(output.to_string()),
        };
    }
    if options.crosscheck_samples {
        return Destination {
            crosscheck_by: DataType::Sample,
            output: "/dev/null".to_string(),
            matrix_output: Some(output.to_string()),
        };
    }
    Destination {
        crosscheck_by: DataType::ReadGroup,
        output: output.to_string(),
        matrix_output: None,
    }
}

/// Whether the parent's `--EXPECT_ALL_GROUPS_TO_MATCH` ends up set.
///
/// The wrapper's own flag is copied onto the parent's when it is true, and the two are declared
/// mutually exclusive, so a command line that sets both is refused by the parser before either is
/// read.
pub fn expect_all_groups_to_match(expect_all_groups: bool, options: &Options) -> bool {
    if options.expect_all_read_groups_to_match {
        return true;
    }
    expect_all_groups
}

/// The parser's refusal for a pair of arguments declared mutually exclusive.
pub fn mutex_refusal(named: &str, others: &[&str]) -> String {
    format!(
        "ERROR: Option '{named}' cannot be used in conjunction with option(s) {}",
        others.join(" ")
    )
}

/// The refusals the wrapper writes itself, in the order it checks them.
///
/// Each names the argument, quotes what was found, and points at the parent tool. The first of
/// them lacks the newline the other two have after the tool's name, which is the sentence the
/// golden records.
pub fn refusals(options: &Options) -> Vec<String> {
    let mut errors = Vec::new();
    if let Some(data_type) = options.crosscheck_by {
        errors.push(format!(
            "When calling CrosscheckReadGroupFingerprints, please refrain from supplying a \
             CROSSCHECK_BY argument. (Found value {}\nUse CrosscheckFingerprints if you would \
             like to do that.",
            data_type.name()
        ));
    }
    if let Some(matrix) = &options.matrix_output {
        errors.push(format!(
            "When calling CrosscheckReadGroupFingerprints, please refrain from supplying a \
             MATRIX_OUTPUT argument.\n(Found value {matrix}\nUse CrosscheckFingerprints if you \
             would like to do that."
        ));
    }
    if !options.second_input.is_empty() {
        errors.push(format!(
            "When calling CrosscheckReadGroupFingerprints, please refrain from supplying a \
             SECOND_INPUT argument.\n(Found value {:?}\nUse CrosscheckFingerprints if you would \
             like to do that.",
            options.second_input
        ));
    }
    errors
}
