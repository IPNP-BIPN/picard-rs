//! `ScatterIntervalsByNs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.ScatterIntervalsByNs.doWork` at tag 3.4.0; the segregation itself is
//! `picard_analysis::scatter_intervals_by_ns`.
//!
//! The reference is `--REFERENCE` here and not `--REFERENCE_SEQUENCE`, which is why this tool had
//! no array until the corpus declared one for that name: the argument is required, so every row was
//! the same usage error, which measures the parser.
//!
//! The output interval list carries the reference's dictionary, and `doWork` gets it the way
//! everything else does, from the `.dict` beside the FASTA. `ReferenceSequenceFileFactory` is asked
//! for the sequence file and its `getSequenceDictionary()` is what the header is written from, so a
//! reference without a dictionary is not a tool that guesses one.
//!
//! `OUTPUT_TYPE` is the argument the array is for: `N` keeps only the `Nmer` runs, `ACGT` only the
//! `ACGTmer` ones, and `BOTH` (the default) keeps every run. `MAX_TO_MERGE` is held at its default
//! by the array -- an `int` with no declared bounds -- and the port reads it anyway.

use std::io::Write;

use picard_analysis::scatter_intervals_by_ns::{scatter_intervals_by_ns, Options, OutputType};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// The `.dict` beside a reference, which is where its dictionary lives.
fn dictionary_path(reference: &str) -> String {
    match reference.rsplit_once('.') {
        Some((stem, _)) => format!("{stem}.dict"),
        None => format!("{reference}.dict"),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let reference = arg(&args, "REFERENCE=")
        .or_else(|| arg(&args, "R="))
        .ok_or("REFERENCE= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let opts = Options {
        output_type: match arg(&args, "OUTPUT_TYPE=").as_deref() {
            None | Some("BOTH") => OutputType::Both,
            Some("N") => OutputType::N,
            Some("ACGT") => OutputType::Acgt,
            Some(other) => return Err(format!("unknown OUTPUT_TYPE: {other}").into()),
        },
        max_to_merge: match arg(&args, "MAX_TO_MERGE=") {
            Some(value) => value
                .parse()
                .map_err(|_| format!("MAX_TO_MERGE: {value}"))?,
            None => 1,
        },
    };

    let fasta = std::fs::read_to_string(&reference)?;
    let dictionary = std::fs::read_to_string(dictionary_path(&reference))?;
    let list = scatter_intervals_by_ns(&fasta, &dictionary, &opts).map_err(|e| format!("{e:?}"))?;

    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(list.as_bytes())?;
    out.flush()?;
    Ok(())
}
