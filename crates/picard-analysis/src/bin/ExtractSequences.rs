//! `ExtractSequences` as a runnable binary: the covering array's port side.
//!
//! The whole tool is one loop, and `LINE_LENGTH` is the only argument that changes a byte of it,
//! so the array varies that and the interval list. The reference is read as text, which is what
//! `ReferenceSequenceFile` does for an unindexed FASTA.

use picard_analysis::extract_sequences::{extract_sequences, ExtractError, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
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

    let intervals = std::fs::read_to_string(&interval_list)?;
    let fasta = std::fs::read_to_string(&reference)?;
    match extract_sequences(&intervals, &fasta, &options) {
        Ok(text) => std::fs::write(&output, text)?,
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
