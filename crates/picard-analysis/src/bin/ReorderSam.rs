//! `ReorderSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.ReorderSam.doWork` at tag 3.4.0 down to what the array varies. The reordering
//! and its two refusals live in `picard_analysis::reorder_sam`.
//!
//! Two arguments decide whether a mismatch between the input's dictionary and the new one is an
//! abort or a rewrite, and they are independent:
//!
//! * `ALLOW_INCOMPLETE_DICT_CONCORDANCE` lets a read contig that the new dictionary does not name
//!   become unmapped, instead of aborting the run;
//! * `ALLOW_CONTIG_LENGTH_DISCORDANCE` lets a name that matches with a different LENGTH warn
//!   instead of aborting.
//!
//! `SEQUENCE_DICTIONARY` is read the way `SAMSequenceDictionaryExtractor` reads it: a FASTA does
//! not carry a dictionary, so the `.dict` beside it is what is opened. Passing `ref.fasta` and
//! reading `ref.dict` is not a shortcut here, it is the reference's own rule.

use std::io::{Read, Write};

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::sam_file::write_sam;
use picard_analysis::reorder_sam::{reorder_sam, reorder_sam_to_bam, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn flag(args: &[String], key: &str) -> bool {
    arg(args, key).map(|v| v == "true").unwrap_or(false)
}

/// The text of a SAM or BAM input, whichever it turns out to be.
fn read_as_sam(path: &str) -> Result<String, Box<dyn std::error::Error>> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?)
    } else {
        Ok(String::from_utf8(raw)?)
    }
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`: a FASTA's dictionary is the `.dict` file
/// beside it, and any other file is read for its own header.
fn read_dictionary(path: &str) -> Result<String, Box<dyn std::error::Error>> {
    for suffix in [".fasta", ".fa"] {
        if let Some(stem) = path.strip_suffix(suffix) {
            let beside = format!("{stem}.dict");
            if std::path::Path::new(&beside).exists() {
                return Ok(std::fs::read_to_string(beside)?);
            }
        }
    }
    read_as_sam(path)
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

    let options = Options {
        allow_incomplete_dict_concordance: flag(&args, "ALLOW_INCOMPLETE_DICT_CONCORDANCE="),
        allow_contig_length_discordance: flag(&args, "ALLOW_CONTIG_LENGTH_DISCORDANCE="),
    };

    let text = read_as_sam(&input)?;
    let dict = read_dictionary(&dictionary)?;

    // `SAMFileWriterFactory.makeWriter` decides by the output's name: `.sam` is text, anything
    // else is a BAM, including the `output.txt` the arrays pass.
    if output.ends_with(".sam") {
        let sam = reorder_sam(&text, &dict, &options).map_err(|e| format!("{e:?}"))?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        return Ok(());
    }
    let bam = reorder_sam_to_bam(&text, &dict, &options).map_err(|e| format!("{e:?}"))?;
    std::fs::write(&output, &bam)?;
    Ok(())
}
