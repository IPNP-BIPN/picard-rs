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

use picard_analysis::args::arg;
use picard_analysis::create_sequence_dictionary::{create_sequence_dictionary_with, Options};

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

/// `loadContigAliasesMap` and the `AN` attribute it feeds: each contig's alternative names, a
/// `HashSet` joined with commas in its iteration order, after the record's other attributes.
fn with_alt_names(dict: &str, alt_names: &str, path: &str) -> String {
    let refuse = |why: String| -> ! {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: Can't read alias file {}",
            picard_analysis::metrics_cli::absolute(path)
        );
        let _ = why;
        std::process::exit(1);
    };
    let mut aliases: Vec<(String, Vec<String>)> = Vec::new();
    for line in alt_names.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let Some((contig, alt)) = line.split_once('\t') else {
            refuse(format!("tabulation missing in {line}"));
        };
        if contig.trim().is_empty() || alt.trim().is_empty() {
            refuse(format!("empty field in {line}"));
        }
        if alt == contig {
            continue;
        }
        if !legal_name(alt) {
            refuse(format!(
                "Illegal alternative reference sequence name in {line}"
            ));
        }
        if aliases.iter().any(|(c, _)| c == alt) {
            refuse(format!(
                "alternate name {alt} previously defined as a contig in {line}"
            ));
        }
        if aliases
            .iter()
            .any(|(c, names)| c != contig && names.iter().any(|n| n == contig))
        {
            refuse(format!(
                "contig {contig} previously defined as an alternate name in {line}"
            ));
        }
        match aliases.iter_mut().find(|(c, _)| c == contig) {
            Some((_, names)) => {
                if !names.iter().any(|n| n == alt) {
                    names.push(alt.to_string());
                }
            }
            None => aliases.push((contig.to_string(), vec![alt.to_string()])),
        }
    }
    let java_hash = |s: &str| {
        s.encode_utf16()
            .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(c as i32))
    };
    let mut out = String::new();
    for line in dict.lines() {
        out.push_str(line);
        if let Some(name) = line
            .strip_prefix("@SQ\tSN:")
            .and_then(|rest| rest.split('\t').next())
        {
            if let Some((_, names)) = aliases.iter().find(|(c, _)| c == name) {
                let keyed: Vec<(i32, String)> =
                    names.iter().map(|n| (java_hash(n), n.clone())).collect();
                let ordered = picard_analysis::barcode_extractor::hash_set_order(keyed, 16);
                out.push_str("\tAN:");
                out.push_str(&ordered.join(","));
            }
        }
        out.push('\n');
    }
    out
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CreateSequenceDictionary");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let reference = arg(&args, "REFERENCE=")
        .or_else(|| arg(&args, "R="))
        .ok_or("REFERENCE= is required")?;
    // `FileExtensions.FASTA`: what the reference's name has to end in.
    const FASTA: [&str; 10] = [
        ".fasta",
        ".fasta.gz",
        ".fas",
        ".fas.gz",
        ".fa",
        ".fa.gz",
        ".fna",
        ".fna.gz",
        ".txt",
        ".txt.gz",
    ];
    let path = std::path::Path::new(&reference);
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let extension = FASTA.iter().find(|e| name.ends_with(*e));
    // `customCommandLineValidation`: no OUTPUT is the reference's own name with its FASTA
    // extension replaced by `.dict` (`getDefaultDictionaryForReferenceSequence`).
    let output = match arg(&args, "OUTPUT=").or_else(|| arg(&args, "O=")) {
        Some(output) => output,
        None => {
            let Some(extension) = extension else {
                eprintln!(
                    "Exception in thread \"main\" java.lang.IllegalArgumentException: File is not a supported reference file type: {}",
                    picard_analysis::fingerprint::reference_path(&reference)
                );
                std::process::exit(1);
            };
            let stem = &name[..name.len() - extension.len()];
            path.with_file_name(format!("{stem}.dict"))
                .to_string_lossy()
                .into_owned()
        }
    };
    // `doWork` refuses to overwrite a dictionary before it reads anything.
    if std::path::Path::new(&output).exists() {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: {} already exists.  Delete this file and try again, or specify a different output file.",
            picard_analysis::fingerprint::file_uri(&output)
        );
        std::process::exit(1);
    }
    // Then the reference is opened by its extension, and one that is not a FASTA's throws inside
    // the writing block, whose handler deletes the dictionary it had begun and says so.
    if extension.is_none() {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: Unknown problem. Partial dictionary file was deleted."
        );
        std::process::exit(1);
    }

    let truncate = arg(&args, "TRUNCATE_NAMES_AT_WHITESPACE=")
        .is_none_or(|v| picard_analysis::args::is_true(&v));

    let options = Options {
        genome_assembly: arg(&args, "GENOME_ASSEMBLY="),
        species: arg(&args, "SPECIES="),
        num_sequences: arg(&args, "NUM_SEQUENCES=").and_then(|v| v.parse().ok()),
    };

    // `createUri`: the argument when given, else the reference's own absolute `file:` URI.
    let uri = match arg(&args, "URI=") {
        Some(uri) => uri,
        // `getHtsPath().getURIString()`: the path made absolute against the working directory, its
        // links left as they are, as a `file://` URI.
        // shortcut: the URI is not percent-encoded; encode it once a path with a space matters.
        None => {
            std::fs::metadata(&reference)?;
            format!("file://{}", std::path::absolute(&reference)?.display())
        }
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
    let dict = match arg(&args, "ALT_NAMES=").or_else(|| arg(&args, "AN=")) {
        Some(path) => with_alt_names(&dict, &std::fs::read_to_string(&path)?, &path),
        None => dict,
    };
    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(dict.as_bytes())?;
    out.flush()?;
    Ok(())
}
