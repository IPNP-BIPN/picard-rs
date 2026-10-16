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
//! `INCLUDE_FILTERED` is accepted and changes nothing here, which is a statement about this corpus
//! rather than about the tool: it filters a VCF input. `SCATTER_COUNT` (or `SCATTER_CONTENT`) above
//! one turns `OUTPUT` into a directory of scattered lists, cut as `SUBDIVISION_MODE` says.

use picard_analysis::interval_list_tools::{
    interval_list_tools, scatter, Action, Options, ScatterMode,
};

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
    let mut scatter_count: i32 = arg(&args, "SCATTER_COUNT=")
        .and_then(|value| value.parse().ok())
        .unwrap_or(1);
    let scatter_content: Option<i64> = arg(&args, "SCATTER_CONTENT=")
        .filter(|value| value != "null")
        .and_then(|value| value.parse().ok());
    let mode_name =
        arg(&args, "SUBDIVISION_MODE=").unwrap_or_else(|| "INTERVAL_SUBDIVISION".to_string());
    let mode = ScatterMode::parse(&mode_name)
        .ok_or_else(|| format!("unknown SUBDIVISION_MODE: {mode_name}"))?;
    if scatter_count < 1 {
        errors.push("SCATTER_COUNT must be greater than 0.".to_string());
    }
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
    // `header.addComment` for every COMMENT: the comments close the header, after any `@PG`.
    let comments: Vec<String> = args
        .iter()
        .filter_map(|a| a.strip_prefix("COMMENT=").or_else(|| a.strip_prefix("CO=")))
        .map(|c| format!("@CO\t{c}\n"))
        .collect();
    let header_end = text
        .lines()
        .take_while(|line| line.starts_with('@'))
        .map(|line| line.len() + 1)
        .sum::<usize>();
    let text = format!(
        "{}{}{}",
        &text[..header_end],
        comments.concat(),
        &text[header_end..]
    );
    // `SCATTER_CONTENT` replaces the count with the list's weight over it, rounded.
    if let Some(content) = scatter_content {
        let intervals = picard_analysis::interval_list_tools::parse_intervals(&text)
            .map_err(|e| format!("{e:?}"))?;
        let weight = mode.list_weight(&intervals) as f64;
        scatter_count = (weight / content as f64 + 0.5).floor() as i32;
    }
    if scatter_count > 1 {
        // `writeScatterIntervals`: OUTPUT is a directory of `temp_NNNN_of_N/scattered.interval_list`.
        let lists = scatter(&text, scatter_count, mode).map_err(|e| format!("{e:?}"))?;
        let (mut bases, mut intervals) = (0i64, 0i64);
        for (index, list) in lists.iter().enumerate() {
            let directory = std::path::Path::new(&output)
                .join(format!("temp_{:04}_of_{scatter_count}", index + 1));
            std::fs::create_dir_all(&directory)?;
            std::fs::write(directory.join("scattered.interval_list"), &list.text)?;
            bases += list.bases;
            intervals += list.intervals;
        }
        if let Some(path) = &count_output {
            let value = if output_value == "BASES" {
                bases
            } else {
                intervals
            };
            std::fs::write(path, format!("{value}\n"))?;
        }
        return Ok(());
    }
    std::fs::write(&output, &text)?;
    // `OUTPUT_VALUE.output`: the bases or the intervals of what was written, one line.
    if let Some(path) = &count_output {
        let intervals: Vec<(i64, i64)> = text
            .lines()
            .filter(|line| !line.starts_with('@'))
            .filter_map(|line| {
                let mut fields = line.split('\t').skip(1);
                Some((fields.next()?.parse().ok()?, fields.next()?.parse().ok()?))
            })
            .collect();
        let value = match output_value.as_str() {
            "BASES" => intervals
                .iter()
                .map(|(start, end)| end - start + 1)
                .sum::<i64>(),
            _ => intervals.len() as i64,
        };
        std::fs::write(path, format!("{value}\n"))?;
    }
    Ok(())
}
