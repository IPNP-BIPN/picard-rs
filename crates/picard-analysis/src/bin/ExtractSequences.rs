//! `ExtractSequences` as a runnable binary: the covering array's port side.
//!
//! The whole tool is one loop, and `LINE_LENGTH` is the only argument that changes a byte of it,
//! so the array varies that and the interval list. The reference is scanned for its line layout
//! and read back per interval, so neither input nor the output is held whole.

use picard_analysis::args::arg;
use std::io::{BufReader, BufWriter, Write};

use picard_analysis::extract_sequences::{extract_sequences_to, ExtractError, Options};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("ExtractSequences");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let interval_list = arg(&args, "INTERVAL_LIST=").ok_or("INTERVAL_LIST= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let reference = arg(&args, "REFERENCE_SEQUENCE=")
        .or_else(|| arg(&args, "R="))
        .ok_or("REFERENCE_SEQUENCE= is required")?;
    let options = Options {
        line_length: arg(&args, "LINE_LENGTH=")
            .map(|v| v.parse::<usize>())
            .transpose()?
            .unwrap_or(80),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let open = |path: &str| std::fs::File::open(path).map(|f| BufReader::with_capacity(1 << 16, f));
    match extract_sequences_to(
        || open(&interval_list),
        || std::fs::File::open(&reference),
        || std::fs::File::create(&output).map(|f| BufWriter::with_capacity(1 << 20, f)),
        &options,
    ) {
        Ok(mut out) => out.flush()?,
        Err(ExtractError::Io(e)) => return Err(e.into()),
        // htsjdk throws from `getSubsequenceAt`, and the class is the one its own callers see.
        Err(ExtractError::UnknownContig(contig)) => {
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: Unknown contig: {contig}"
            );
            std::process::exit(1);
        }
        Err(ExtractError::IntervalOutOfRange { contig, start, end }) => {
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: Query asks for data \
                 past end of contig: {contig}:{start}-{end}"
            );
            std::process::exit(1);
        }
        Err(other) => return Err(format!("{other:?}").into()),
    }
    Ok(())
}
