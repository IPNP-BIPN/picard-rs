//! `NormalizeFasta` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.reference.NormalizeFasta.doWork` at tag 3.4.0 down to its arguments; the
//! rewrapping is `picard_analysis::normalize_fasta`.
//!
//! `TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE` is not a flag the tool applies itself: it is passed to
//! `ReferenceSequenceFileFactory.getReferenceSequenceFile(INPUT, truncateNamesAtWhitespace)`, so
//! it decides what the *reader* calls each sequence, and the tool writes back whatever name it was
//! given. On `ref.fasta`, whose headers are bare contig names, there is nothing to truncate and
//! the argument is unobservable; `described.fasta` gives each header a description, one after a
//! space and one after a tab, which is what makes the row a test.
//!
//! `LINE_LENGTH` is held at its default by the array (an `int` with no declared bounds), and the
//! port reads it anyway, so that a row naming it measures the tool rather than the parser. The
//! wrapping itself is a `i % LINE_LENGTH == 0` newline before each base, so a sequence whose
//! length is an exact multiple does not end with a blank line.
//!
//! `CREATE_INDEX` and `VALIDATION_STRINGENCY` reach nothing: the output is text through
//! `IOUtil.openFileForBufferedWriting`, and the input is read by a `ReferenceSequenceFile`, which
//! consults neither.

use std::io::Write;

use picard_analysis::args::arg;
use picard_analysis::fasta_records::{records, unreadable};
use picard_analysis::normalize_fasta::{normalize_fasta_to, NormalizeError, Options};

/// A boolean argument, which the legacy parser builds as `new Boolean(value)`.
fn flag(args: &[String], key: &str, default: bool) -> bool {
    arg(args, key).map_or(default, |value| picard_analysis::args::is_true(&value))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("NormalizeFasta");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    // `doWork` refuses this before it opens anything, and says so in those words.
    if std::fs::canonicalize(&input).ok() == std::fs::canonicalize(&output).ok() {
        eprintln!(
            "Exception in thread \"main\" java.lang.IllegalArgumentException: \
             Input and output cannot be the same file."
        );
        std::process::exit(1);
    }

    let opts = Options {
        line_length: match arg(&args, "LINE_LENGTH=") {
            Some(value) => value.parse().map_err(|_| format!("LINE_LENGTH: {value}"))?,
            None => 100,
        },
        truncate_names_at_whitespace: flag(&args, "TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=", false),
    };

    // A refusal writes no file, as the whole-file parse did, so every record is checked first.
    let check = records(
        std::fs::File::open(&input)?,
        opts.truncate_names_at_whitespace,
    );
    if let Some(error) = check.into_iter().find_map(Result::err) {
        return Err(match unreadable(&input) {
            Some(io) => io.into(),
            None => format!("{:?}", NormalizeError::from(error)).into(),
        });
    }
    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    normalize_fasta_to(std::fs::File::open(&input)?, &opts, &mut out)
        .map_err(|e| format!("{e:?}"))?;
    out.flush()?;
    Ok(())
}
