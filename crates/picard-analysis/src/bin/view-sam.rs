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

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::sam_file::write_sam;
use picard_analysis::view_sam::{view_sam, AlignmentStatus, Options, PfStatus};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn flag(args: &[String], key: &str) -> bool {
    arg(args, key).map(|v| v == "true").unwrap_or(false)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    // reach the same code by different doors.
    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let text = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?
    } else {
        String::from_utf8(raw)?
    };

    let out = view_sam(&text, &options).map_err(|e| format!("{e:?}"))?;
    let mut stdout = std::io::BufWriter::new(std::io::stdout().lock());
    stdout.write_all(out.as_bytes())?;
    stdout.flush()?;
    Ok(())
}
