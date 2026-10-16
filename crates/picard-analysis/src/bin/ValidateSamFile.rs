//! `ValidateSamFile` as a runnable binary: the covering array's port side.
//!
//! The tool's exit code is its answer, not a failure: 0 clean, 1 warnings only, 2 errors and
//! warnings, 3 errors only, and the report is written either way. The array records the code
//! together with the report for that reason (`--exit-code-is-a-result` in the runner).
//!
//! `MAX_OUTPUT` is what makes the verbose report finite. htsjdk throws
//! `MaxOutputExceededException` from inside `addError` as soon as the hundredth error is printed,
//! which ends the whole validation -- not just the printing -- and the catch prints one last line.
//! The corpus raises 421 errors and warnings, so every verbose row ends there.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use picard_analysis::validate_sam_file::{validate_records, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// Every occurrence of a repeatable argument, which is how `IGNORE` is given.
fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=").or_else(|| arg(&args, "O="));
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);
    let validate_index = flag("VALIDATE_INDEX=", true);
    let index_stringency =
        arg(&args, "INDEX_VALIDATION_STRINGENCY=").unwrap_or_else(|| "EXHAUSTIVE".to_string());

    // `customCommandLineValidation`, printed by Barclay after the usage block.
    if (!validate_index && index_stringency != "NONE")
        || (validate_index && index_stringency == "NONE")
    {
        eprintln!(
            "VALIDATE_INDEX and INDEX_VALIDATION_STRINGENCY must be consistent: VALIDATE_INDEX is \
             {validate_index} and INDEX_VALIDATION_STRINGENCY is {index_stringency}"
        );
        std::process::exit(1);
    }

    let mode = arg(&args, "MODE=")
        .or_else(|| arg(&args, "M="))
        .unwrap_or_else(|| "VERBOSE".to_string());
    let options = Options {
        verbose: match mode.as_str() {
            "VERBOSE" => true,
            "SUMMARY" => false,
            other => return Err(format!("unknown MODE: {other}").into()),
        },
        ignore: args_all(&args, "IGNORE="),
        ignore_warnings: flag("IGNORE_WARNINGS=", false),
        skip_mate_validation: flag("SKIP_MATE_VALIDATION=", false),
        max_output: arg(&args, "MAX_OUTPUT=")
            .map(|v| v.parse::<usize>())
            .transpose()?
            .unwrap_or(100),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

    // The reference is read only for the NM VALUE check; without an `NM` tag on a record there is
    // nothing to compare, which is the case for every read in this corpus.
    let reference = arg(&args, "REFERENCE_SEQUENCE=")
        .or_else(|| arg(&args, "R="))
        .map(std::fs::read)
        .transpose()?;

    let validation = validate_records(&header, &records, reference.as_deref(), &options)
        .map_err(|e| format!("{e:?}"))?;

    match &output {
        Some(path) => std::fs::write(path, &validation.report)?,
        None => print!("{}", validation.report),
    }
    std::process::exit(validation.exit_code());
}
