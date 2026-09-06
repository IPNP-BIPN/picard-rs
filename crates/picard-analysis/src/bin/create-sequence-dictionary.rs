//! `CreateSequenceDictionary` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CreateSequenceDictionary.doWork` at tag 3.4.0. The dictionary itself, and the
//! order its attributes are written in, live in `picard_analysis::create_sequence_dictionary`.
//!
//! `URI` is the argument that shows why `UR` is a parameter of the port rather than something it
//! derives: given, it is written verbatim; absent, the reference's own `file:` URI is used, which
//! is path-dependent and therefore canonicalized away by the comparison rather than matched.
//!
//! `TRUNCATE_NAMES_AT_WHITESPACE` is accepted only at its default. The reader beneath this port
//! truncates and drops the description, so the false path cannot be answered honestly; the array
//! holds the argument rather than covering it with a value the port would answer wrongly.

use std::io::Write;

use picard_analysis::create_sequence_dictionary::{create_sequence_dictionary_with, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let reference = arg(&args, "REFERENCE=")
        .or_else(|| arg(&args, "R="))
        .ok_or("REFERENCE= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    if let Some(value) = arg(&args, "TRUNCATE_NAMES_AT_WHITESPACE=") {
        if value != "true" {
            return Err("TRUNCATE_NAMES_AT_WHITESPACE=false is not ported".into());
        }
    }

    let options = Options {
        genome_assembly: arg(&args, "GENOME_ASSEMBLY="),
        species: arg(&args, "SPECIES="),
        num_sequences: arg(&args, "NUM_SEQUENCES=").and_then(|v| v.parse().ok()),
    };

    // `createUri`: the argument when given, else the reference's own absolute `file:` URI.
    let uri = match arg(&args, "URI=") {
        Some(uri) => uri,
        None => format!("file:{}", std::fs::canonicalize(&reference)?.display()),
    };

    let fasta = std::fs::read(&reference)?;
    let dict =
        create_sequence_dictionary_with(&fasta, &uri, &options).map_err(|e| format!("{e:?}"))?;
    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(dict.as_bytes())?;
    out.flush()?;
    Ok(())
}
