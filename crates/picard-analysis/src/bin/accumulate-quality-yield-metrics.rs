//! `AccumulateQualityYieldMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.AccumulateQualityYieldMetrics.doWork` at tag 3.4.0. The parse, the merge and
//! the derived read length live in `picard_analysis::accumulate_quality_yield_metrics`.
//!
//! The tool reads `CollectQualityYieldMetrics` files rather than reads, adds their ten counters,
//! and recomputes the one column that is not a counter. Its output carries NO header comments: it
//! writes a bare `new MetricsFile<>()` instead of the base class's `getMetricsFile()`, so there is
//! no command line and no start time to canonicalize away.

use picard_analysis::accumulate_quality_yield_metrics::accumulate_quality_yield_metrics;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let mut inputs: Vec<String> = args
        .iter()
        .filter_map(|a| a.strip_prefix("INPUT=").map(str::to_string))
        .collect();
    inputs.extend(
        args.iter()
            .filter_map(|a| a.strip_prefix("I=").map(str::to_string)),
    );
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let texts: Vec<String> = inputs
        .iter()
        .map(std::fs::read_to_string)
        .collect::<Result<_, _>>()?;
    let borrowed: Vec<&str> = texts.iter().map(String::as_str).collect();
    match accumulate_quality_yield_metrics(&borrowed) {
        Ok(text) => std::fs::write(&output, text)?,
        Err(error) => {
            eprintln!("Exception in thread \"main\" picard.PicardException: {error:?}");
            std::process::exit(1);
        }
    }
    Ok(())
}
