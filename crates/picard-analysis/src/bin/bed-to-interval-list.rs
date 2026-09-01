//! `BedToIntervalList` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.BedToIntervalList.doWork` at tag 3.4.0 down to its arguments; the conversion
//! itself is `picard_analysis::bed_to_interval_list`.
//!
//! This is the first tool measured against a fixture that is not a BAM. Its `--INPUT` is a BED and
//! its `--SEQUENCE_DICTIONARY` a `.dict`, both declared under `per_tool` in
//! `tools/coverage/fixtures.json`, because one shared `--INPUT` list cannot serve tools that do
//! not read the same kind of file.
//!
//! What the array varies, and what each one does here:
//!
//! * `SORT` (default true) puts the intervals in coordinate order before writing; the corpus is
//!   already in that order, so the two values agree on it. That is a property of the fixture and
//!   not of the tool, and the recorded `distinct_outputs` is where it shows.
//! * `UNIQUE` (default false) merges overlapping and abutting intervals, concatenating their
//!   names. The corpus's three targets are disjoint and not abutting, so again both values agree.
//! * `KEEP_LENGTH_ZERO_INTERVALS` keeps features whose BED start equals its end; the corpus has
//!   none, and the reference says so on stderr rather than in the output.
//! * `CREATE_INDEX` is inherited from `CommandLineProgram` and reaches nothing: the output is an
//!   interval list written through `IntervalList.write`, not a `SAMFileWriter`, so no index is
//!   ever created and no refusal is raised. The port accepts it and writes the same file, which is
//!   what the row is checking.
//!
//! `VALIDATION_STRINGENCY` and `REFERENCE_SEQUENCE` reach no reader here either: the BED is parsed
//! by `BEDCodec` and the dictionary by `SAMSequenceDictionaryExtractor`, neither of which consults
//! them. Accepting them is what lets a row that names them measure the tool.

use std::io::Write;

use picard_analysis::bed_to_interval_list::{bed_to_interval_list, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// A Barclay `Boolean` argument, which is `true`/`false` and nothing else.
fn flag(args: &[String], key: &str, default: bool) -> Result<bool, String> {
    match arg(args, key).as_deref() {
        None => Ok(default),
        Some("true") => Ok(true),
        Some("false") => Ok(false),
        Some(other) => Err(format!("unknown value for {key}{other}")),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "SD="))
        .ok_or("SEQUENCE_DICTIONARY= is required")?;

    let opts = Options {
        sort: flag(&args, "SORT=", true)?,
        unique: flag(&args, "UNIQUE=", false)?,
        keep_length_zero_intervals: flag(&args, "KEEP_LENGTH_ZERO_INTERVALS=", false)?,
        drop_missing_contigs: flag(&args, "DROP_MISSING_CONTIGS=", false)?,
    };

    let dictionary_text = std::fs::read_to_string(&dictionary)?;
    let bed = std::fs::read_to_string(&input)?;
    let list = bed_to_interval_list(&dictionary_text, &bed, &opts).map_err(|e| format!("{e:?}"))?;

    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(list.as_bytes())?;
    out.flush()?;
    Ok(())
}
