//! `LiftOverIntervalList` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.LiftOverIntervalList.doWork` at tag 3.4.0. The lift, the sort into the
//! TARGET dictionary's order and the rendering live in
//! `picard_analysis::lift_over_interval_list`; the chain itself is `htsjdk_bam::liftover`.
//!
//! The tool's answer is in two places. The intervals that lifted go to `OUTPUT`, in the target
//! dictionary's order rather than the input's, and the exit code is one when ANY interval failed
//! to lift -- which is not a failure but a count: the run wrote its file either way.
//!
//! `MIN_LIFTOVER_PCT` is what decides a partial lift. An interval whose chain covers only part of
//! it lifts when that part is at least this fraction of it, and is rejected otherwise, so the same
//! interval list against the same chain is two different files at 0.95 and at 0.3.

use picard_analysis::lift_over_interval_list::{
    lift_over_interval_list, LiftOverIntervalListError,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`, for the two shapes the corpus has.
///
/// A `.dict` is read as it lies; a FASTA is not parsed at all -- the extractor looks for the
/// `.dict` beside it and refuses when there is none, which is why every reference in this corpus
/// has one.
fn read_dictionary(path: &str) -> std::io::Result<String> {
    let candidate = std::path::Path::new(path);
    if let Some(extension) = candidate.extension().and_then(|e| e.to_str()) {
        if matches!(extension, "fasta" | "fa" | "fna") {
            let beside = candidate.with_extension("dict");
            return std::fs::read_to_string(beside);
        }
    }
    std::fs::read_to_string(path)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let chain = arg(&args, "CHAIN=").ok_or("CHAIN= is required")?;
    let dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "SD="))
        .ok_or("SEQUENCE_DICTIONARY= is required")?;
    let min_liftover_pct = arg(&args, "MIN_LIFTOVER_PCT=")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.95);

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let result = lift_over_interval_list(
        &std::fs::read_to_string(&input)?,
        &read_dictionary(&dictionary)?,
        &std::fs::read_to_string(&chain)?,
        min_liftover_pct,
    );
    let result = match result {
        Ok(result) => result,
        Err(error) => {
            let message = match error {
                LiftOverIntervalListError::Chain(chain) => format!("{chain:?}"),
                LiftOverIntervalListError::Input(parse) => format!("{parse:?}"),
                LiftOverIntervalListError::MissingToSequence(missing) => missing.to_string(),
            };
            eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
            std::process::exit(1);
        }
    };
    std::fs::write(&output, &result.output)?;
    // `return anyRejected ? 1 : 0`, after the file has been written.
    std::process::exit(result.return_code);
}
