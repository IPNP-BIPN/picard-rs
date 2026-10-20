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

use std::fs::File;
use std::io::{BufWriter, Write};

use picard_analysis::args::arg;
use picard_analysis::lift_over_interval_list::{lift, LiftOverIntervalListError};

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
    picard_analysis::argument_check::check("LiftOverIntervalList");
    picard_analysis::output::install();
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

    // The input is read line by line; the dictionary and the chain, read whole after it, are
    // reported only once the input has read cleanly, the order whole-file reads gave.
    let input = std::io::BufReader::new(std::fs::File::open(&input)?);
    let texts = read_dictionary(&dictionary)
        .and_then(|dictionary| Ok((dictionary, std::fs::read_to_string(&chain)?)));
    let (dictionary, chain) = match texts {
        Ok(texts) => texts,
        Err(error) => {
            picard_analysis::compact_interval::drain(input)?;
            return Err(error.into());
        }
    };
    let result = lift(input, &dictionary, &chain, min_liftover_pct)?;
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
    let write = |path: &str, render: &dyn Fn(&mut BufWriter<File>) -> std::io::Result<()>| {
        let mut out = BufWriter::new(File::create(path)?);
        render(&mut out)?;
        out.flush()
    };
    write(&output, &|out| result.write_output(out))?;
    if let Some(reject) = arg(&args, "REJECT=") {
        write(&reject, &|out| result.write_rejects(out))?;
    }
    // `return anyRejected ? 1 : 0`, after the file has been written.
    std::process::exit(result.return_code());
}
