//! `IntervalListTools` as a runnable binary: the covering array's port side.
//!
//! Thirty-three rows, seventeen of which the reference accepts, and every one of those seventeen
//! produces a DIFFERENT output: `ACTION`, `SORT`, `UNIQUE`, `DONT_MERGE_ABUTTING` and `INVERT`
//! interact, and the array sees all of it.
//!
//! The other sixteen are one rule in three voices. `Action.takesSecondInput` splits the six
//! actions in half, and this harness passes the same arguments to every row, so the three actions
//! that require a `SECOND_INPUT` are refused for not having one. That refusal is the tool's, and
//! its message names the action, so the three read differently and count as three.
//!
//! `INCLUDE_FILTERED`, `OUTPUT_VALUE` and `SUBDIVISION_MODE` are accepted and change nothing here,
//! which is a statement about this corpus rather than about the tool: `INCLUDE_FILTERED` filters a
//! VCF input, `OUTPUT_VALUE` decides what is COUNTED on standard output rather than what is
//! written to `OUTPUT`, and `SUBDIVISION_MODE` only does anything under `SCATTER_COUNT`.

use picard_analysis::interval_list_tools::{interval_list_tools, Action, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let second_input = arg(&args, "SECOND_INPUT=").or_else(|| arg(&args, "SI="));
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);

    let action_name = arg(&args, "ACTION=").unwrap_or_else(|| "CONCAT".to_string());
    let action = match action_name.as_str() {
        "CONCAT" => Action::Concat,
        "UNION" => Action::Union,
        "INTERSECT" => Action::Intersect,
        "OVERLAPS" => Action::Overlaps,
        "SUBTRACT" => Action::Subtract,
        "SYMDIFF" => Action::Symdiff,
        other => return Err(format!("unknown ACTION: {other}").into()),
    };
    let takes_second_input = matches!(
        action,
        Action::Overlaps | Action::Subtract | Action::Symdiff
    );

    let output_value = arg(&args, "OUTPUT_VALUE=").unwrap_or_else(|| "NONE".to_string());
    let count_output = arg(&args, "COUNT_OUTPUT=");

    // `customCommandLineValidation`, in its own order: Barclay collects every message and prints
    // them after the usage block, so a row that breaks two rules prints both.
    let mut errors: Vec<String> = Vec::new();
    if second_input.is_none() && takes_second_input {
        errors.push(format!(
            "SECOND_INPUT was not provided but action {action_name} requires a second input."
        ));
    }
    if second_input.is_some() && !takes_second_input {
        errors.push(format!(
            "SECOND_INPUT was provided but action {action_name} doesn't take a second input."
        ));
    }
    if count_output.is_some() && output_value == "NONE" {
        errors.push("COUNT_OUTPUT was provided but OUTPUT_VALUE is set to NONE.".to_string());
    }
    if !errors.is_empty() {
        for message in &errors {
            eprintln!("{message}");
        }
        std::process::exit(1);
    }

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let options = Options {
        action,
        sort: flag("SORT=", true),
        unique: flag("UNIQUE=", false),
        dont_merge_abutting: flag("DONT_MERGE_ABUTTING=", false),
        invert: flag("INVERT=", false),
        padding: arg(&args, "PADDING=")
            .map(|v| v.parse::<i32>())
            .transpose()?
            .unwrap_or(0),
        break_bands_at_multiples_of: arg(&args, "BREAK_BANDS_AT_MULTIPLES_OF=")
            .map(|v| v.parse::<i32>())
            .transpose()?
            .unwrap_or(0),
    };

    let first = std::fs::read_to_string(&input)?;
    let second = second_input
        .as_ref()
        .map(std::fs::read_to_string)
        .transpose()?;
    let seconds: Vec<&str> = second.as_deref().into_iter().collect();
    let text = interval_list_tools(&[&first], &seconds, &options).map_err(|e| format!("{e:?}"))?;
    std::fs::write(&output, text)?;
    Ok(())
}
