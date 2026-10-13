//! `CreateSequenceDictionary` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CreateSequenceDictionary.doWork` at tag 3.4.0. The dictionary itself, and the
//! order its attributes are written in, live in `picard_analysis::create_sequence_dictionary`.
//!
//! `URI` is the argument that shows why `UR` is a parameter of the port rather than something it
//! derives: given, it is written verbatim; absent, the reference's own `file:` URI is used, which
//! is path-dependent and therefore canonicalized away by the comparison rather than matched.
//!
//! `TRUNCATE_NAMES_AT_WHITESPACE=false` keeps the whole header line as the sequence name. A name
//! with whitespace in it is never a legal `@SQ` name, so on that path the first contig whose header
//! carries a description is refused by `SAMSequenceRecord.validateSequenceName` as the dictionary
//! is built, and a header without one names its contig exactly as the truncating path does.

use std::io::Write;

use picard_analysis::create_sequence_dictionary::{create_sequence_dictionary_with, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `SAMSequenceRecord.LEGAL_RNAME_PATTERN`, as its `toString` prints it.
const LEGAL_RNAME: &str = "[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*";

fn legal_name(name: &str) -> bool {
    const FIRST: &[u8] = b"!#$%&+./:;?@^_|~-";
    const REST: &[u8] = b"!#$%&*+./:;=?@^_|~-";
    let bytes = name.as_bytes();
    !bytes.is_empty()
        && (bytes[0].is_ascii_alphanumeric() || FIRST.contains(&bytes[0]))
        && bytes[1..]
            .iter()
            .all(|b| b.is_ascii_alphanumeric() || REST.contains(b))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let reference = arg(&args, "REFERENCE=")
        .or_else(|| arg(&args, "R="))
        .ok_or("REFERENCE= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    let truncate = arg(&args, "TRUNCATE_NAMES_AT_WHITESPACE=").is_none_or(|v| v == "true");

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
    if !truncate {
        // The contigs are read in order and at most NUM_SEQUENCES of them, so the refusal is the
        // first whole header among those that the legal-name pattern does not match.
        let limit = options.num_sequences.unwrap_or(usize::MAX);
        let names = String::from_utf8_lossy(&fasta)
            .lines()
            .filter_map(|line| line.strip_prefix('>').map(|name| name.trim().to_string()))
            .take(limit)
            .collect::<Vec<_>>();
        if let Some(name) = names.iter().find(|name| !legal_name(name)) {
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: Sequence name '{name}' doesn't match regex: '{LEGAL_RNAME}' "
            );
            std::process::exit(1);
        }
    }
    let dict =
        create_sequence_dictionary_with(&fasta, &uri, &options).map_err(|e| format!("{e:?}"))?;
    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(dict.as_bytes())?;
    out.flush()?;
    Ok(())
}
