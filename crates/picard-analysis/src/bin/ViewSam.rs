//! `ViewSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.ViewSam.doWork` at tag 3.4.0 down to what the array varies. The record
//! filtering and the header rules live in `picard_analysis::view_sam`; this is the argument
//! surface around them.
//!
//! Three decisions the arrays reach:
//!
//! * `ALIGNMENT_STATUS` and `PF_STATUS` filter, and they filter INDEPENDENTLY: a record must pass
//!   both, so `Aligned` with `NonPf` prints the aligned records that failed vendor quality rather
//!   than nothing.
//! * `HEADER_ONLY` and `RECORDS_ONLY` are two booleans over one output, and the reference REFUSES
//!   the pair: `customCommandLineValidation` returns "Cannot specify both HEADER_ONLY=true and
//!   RECORDS_ONLY=true." and the run exits 1. The obvious reading -- both suppress a half, so both
//!   together print nothing -- is what this binary did until the array measured it, and it turned
//!   six of fourteen rows into a silent empty file where the reference refused.
//! * The output is standard output, always. `ViewSam` has no `OUTPUT` argument, which is why the
//!   runner compares stdout for this tool rather than a file.
//!
//! `REFERENCE_SEQUENCE`, `TMP_DIR` and `MAX_RECORDS_IN_RAM` are accepted and unused: nothing here
//! spills to disk and the reader needs a reference only for CRAM. Accepting them is what lets a
//! row that names them measure the tool rather than the argument parser.

use std::io::{Read, Write};

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::{parse_line_with, ValidationStringency};
use picard_analysis::args::arg;
use picard_analysis::view_sam::{keep, render, AlignmentStatus, Options, PfStatus};
use rayon::prelude::*;

/// How many records are rendered in parallel at a time.
const CHUNK: usize = 16 * 1024;

fn flag(args: &[String], key: &str) -> bool {
    arg(args, key)
        .map(|v| picard_analysis::args::is_true(&v))
        .unwrap_or(false)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("ViewSam");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;

    // `customCommandLineValidation`, which runs before `doWork` and exits 1 with this message.
    if flag(&args, "HEADER_ONLY=") && flag(&args, "RECORDS_ONLY=") {
        eprintln!("Cannot specify both HEADER_ONLY=true and RECORDS_ONLY=true.");
        std::process::exit(1);
    }

    let options = Options {
        alignment_status: match arg(&args, "ALIGNMENT_STATUS=").as_deref() {
            None | Some("All") => AlignmentStatus::All,
            Some("Aligned") => AlignmentStatus::Aligned,
            Some("Unaligned") => AlignmentStatus::Unaligned,
            Some(other) => return Err(format!("unknown ALIGNMENT_STATUS: {other}").into()),
        },
        pf_status: match arg(&args, "PF_STATUS=").as_deref() {
            None | Some("All") => PfStatus::All,
            Some("PF") => PfStatus::Pf,
            Some("NonPF") => PfStatus::NonPf,
            Some(other) => return Err(format!("unknown PF_STATUS: {other}").into()),
        },
        header_only: flag(&args, "HEADER_ONLY="),
        records_only: flag(&args, "RECORDS_ONLY="),
    };

    // `SamReaderFactory` sniffs the stream rather than trusting the name, so a `.bam` and a `.sam`
    // reach the same code by different doors. A BAM's records were written as SAM text and that
    // text parsed back at LENIENT, as `view_sam` reads a SAM: each record makes the same round
    // trip on its own.
    let mut magic = Vec::new();
    std::fs::File::open(&input)?
        .take(4)
        .read_to_end(&mut magic)?;
    let bam = magic.starts_with(&[0x1f, 0x8b]) || picard_analysis::cram_input::is_cram(&magic);
    let (read, records) =
        picard_analysis::bam_stream::Input::open_with(&input, ValidationStringency::Lenient)?
            .into_parts();
    let header = if bam {
        read_sam_with(&read.encode(), ValidationStringency::Lenient)
            .map_err(|e| format!("{e:?}"))?
            .0
    } else {
        read.clone()
    };
    let index_of = |name: &str| {
        header
            .sequences
            .iter()
            .position(|s| s.name == name)
            .map(|i| i as i32)
    };
    let line = |record: &BamRecord| -> Result<Option<String>, String> {
        let parsed;
        let record = if bam {
            let text = picard_analysis::bam_stream::sam_line(&read, record)
                .ok_or("records failed to re-encode as SAM")?;
            parsed = parse_line_with(&text, index_of, ValidationStringency::Lenient)
                .map_err(|e| format!("{e:?}"))?;
            &parsed
        } else {
            record
        };
        Ok(keep(record, &options).then(|| render(&header, record)))
    };

    // The whole text was built before any of it was printed, so a run that fails prints nothing:
    // a failure before the buffer first fills leaves it unprinted, as the whole text was.
    let mut stdout = std::io::BufWriter::with_capacity(1 << 20, std::io::stdout().lock());
    // RECORDS_ONLY suppresses the header; otherwise print it verbatim (getTextHeader).
    if !options.records_only {
        stdout.write_all(header.encode().as_bytes())?;
    }
    // HEADER_ONLY stops after the header, once every record has been read as before.
    let mut records = records.peekable();
    while records.peek().is_some() {
        let chunk: Vec<BamRecord> = records.by_ref().take(CHUNK).collect();
        for rendered in chunk.par_iter().map(line).collect::<Vec<_>>() {
            match rendered {
                Ok(Some(text)) if !options.header_only => stdout.write_all(text.as_bytes())?,
                Ok(_) => {}
                Err(e) => {
                    std::mem::forget(stdout);
                    return Err(e.into());
                }
            }
        }
    }
    stdout.flush()?;
    Ok(())
}
