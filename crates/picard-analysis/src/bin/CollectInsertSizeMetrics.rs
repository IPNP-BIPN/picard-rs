//! `CollectInsertSizeMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectInsertSizeMetrics` at tag 3.4.0. The collector is
//! `picard_analysis::insert_size`; this is `SinglePassSamProgram.makeItSo` around it:
//!
//! * the sort-order refusal, which `ASSUME_SORTED` (true by default) turns into a warning, and the
//!   reference walker that refuses to rewind whenever `REFERENCE_SEQUENCE` is given;
//! * the walk stops at the first unmapped record, because this program does not use them, or after
//!   `STOP_AFTER` records;
//! * the metrics file is written only when some orientation cleared `MINIMUM_PCT`; otherwise the
//!   tool warns and writes nothing, and the chart is not drawn.
//!
//! The chart itself is R's, and not ported: the port writes the metrics and no PDF. What it does
//! reproduce is R's failure when every histogram was trimmed away, which ends the run at one after
//! the metrics are written; the reference's log lines carry a timestamp the port's do not.
//! `METRIC_ACCUMULATION_LEVEL` adds levels to its `ALL_READS` default: each level's units are the
//! header's read groups, in header order, and a read whose group lacks the field goes to an
//! `unknown` unit made when the first such read arrives.

use picard_analysis::args::arg;
use picard_analysis::insert_size::{Distributors, Options};
use picard_analysis::single_pass_rejections::{check_sort_order, walk_reference, SortOrder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CollectInsertSizeMetrics");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    // Optional since Picard 3.5.0, which refuses it only where there is no R.
    picard_analysis::metrics_cli::refuse_chart_without_r(
        "CollectInsertSizeMetrics",
        arg(&args, "HISTOGRAM_FILE=").is_some() || arg(&args, "H=").is_some(),
    );
    let flag = |key: &str, default: bool| {
        arg(&args, key)
            .map(|v| v.eq_ignore_ascii_case("true"))
            .unwrap_or(default)
    };
    let optional_int = |key: &str| {
        arg(&args, key)
            .filter(|v| v != "null")
            .and_then(|v| v.parse::<i32>().ok())
    };
    let minimum_pct: f32 = arg(&args, "MINIMUM_PCT=")
        .or_else(|| arg(&args, "M="))
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.05);
    if !(0.0..=0.5).contains(&minimum_pct) {
        eprintln!(
            "MINIMUM_PCT was set to {minimum_pct:?}. It must be between 0 and 0.5 so all data categories don't get discarded."
        );
        std::process::exit(1);
    }
    // A collection argument is appended to its default, `ALL_READS`, and `null` empties it.
    let mut levels: Vec<String> = vec!["ALL_READS".to_string()];
    for level in args.iter().filter_map(|a| {
        a.strip_prefix("METRIC_ACCUMULATION_LEVEL=")
            .or_else(|| a.strip_prefix("LEVEL="))
    }) {
        if level == "null" {
            levels.clear();
        } else if !matches!(level, "ALL_READS" | "SAMPLE" | "LIBRARY" | "READ_GROUP") {
            return Err(format!("unknown METRIC_ACCUMULATION_LEVEL: {level}").into());
        } else {
            levels.push(level.to_string());
        }
    }
    let options = Options {
        minimum_pct,
        deviations: arg(&args, "DEVIATIONS=")
            .and_then(|v| v.parse().ok())
            .unwrap_or(10.0),
        histogram_width: optional_int("HISTOGRAM_WIDTH=").or_else(|| optional_int("W=")),
        min_histogram_width: optional_int("MIN_HISTOGRAM_WIDTH=").or_else(|| optional_int("MW=")),
        include_duplicates: flag("INCLUDE_DUPLICATES=", false),
    };
    let stop_after: i64 = arg(&args, "STOP_AFTER=")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);

    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();

    let found = match header.attributes.get("SO") {
        Some("coordinate") => SortOrder::Coordinate,
        Some("queryname") => SortOrder::Queryname,
        Some("duplicate") => SortOrder::Duplicate,
        Some("unsorted") => SortOrder::Unsorted,
        _ => SortOrder::Unknown,
    };
    if let Err(rejection) = check_sort_order(&input, found, flag("ASSUME_SORTED=", true)) {
        eprintln!("Exception in thread \"main\" {}", rejection.thrown());
        std::process::exit(1);
    }
    let with_reference = arg(&args, "REFERENCE_SEQUENCE=").is_some() || arg(&args, "R=").is_some();
    let mut distributors = Distributors::new(&levels, &header, options);
    let mut current: Option<i32> = None;
    let mut count = 0i64;
    for record in records {
        let record = &record;
        if with_reference && record.reference_index != -1 {
            match walk_reference(current, record.reference_index) {
                Ok(next) => current = Some(next),
                Err(rejection) => {
                    eprintln!("Exception in thread \"main\" {}", rejection.thrown());
                    std::process::exit(1);
                }
            }
        }
        if let Err(thrown) = distributors.accept(&header, record) {
            eprintln!("Exception in thread \"main\" {thrown}");
            std::process::exit(1);
        }
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
        // `usesNoRefReads()` is false: the unmapped reads at the end are not walked.
        if record.reference_index == -1 {
            break;
        }
    }

    let headers = picard_analysis::command_line::default_headers("CollectInsertSizeMetrics");
    if let Some(failure) = distributors.write(
        &output,
        &headers,
        optional_int("HISTOGRAM_WIDTH=").or_else(|| optional_int("W=")),
        optional_int("MIN_HISTOGRAM_WIDTH=").or_else(|| optional_int("MW=")),
    )? {
        eprint!("{failure}");
        std::process::exit(1);
    }
    Ok(())
}
