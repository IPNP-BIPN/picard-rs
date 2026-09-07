//! `CheckDuplicateMarking` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.markduplicates.CheckDuplicateMarking.doWork` at tag 3.4.0. The walk and what
//! each `MODE` skips live in `picard_analysis::check_duplicate_marking`.
//!
//! The tool checks that records sharing a query name agree about their duplicate flag, and NOT
//! that the marking is right. Its answer is in three places at once: the names of the offending
//! records go to `OUTPUT`, one per line, and the exit code is the SIGN of how many there were --
//! one, never the count -- while the count itself is only logged.
//!
//! A file that is not query-name sorted is sorted here rather than refused, which makes the
//! comparator part of the answer: the order it sorts into decides which record of a name the
//! others are compared against.

use std::io::{Read, Write};

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::query_name;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use picard_analysis::check_duplicate_marking::{check, Mode, Record};

const READ_PAIRED_FLAG: u16 = 0x1;
const PROPER_PAIR_FLAG: u16 = 0x2;
const READ_UNMAPPED_FLAG: u16 = 0x4;
const DUPLICATE_FLAG: u16 = 0x400;
const SECONDARY_FLAG: u16 = 0x100;
const SUPPLEMENTARY_FLAG: u16 = 0x800;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    // "Output file into which bad querynames will be placed (if not null)": with no OUTPUT the
    // names go to a null stream and only the exit code survives.
    let output = arg(&args, "OUTPUT=").or_else(|| arg(&args, "O="));
    let mode = match arg(&args, "MODE=").as_deref() {
        None | Some("ALL") => Mode::All,
        Some("PRIMARY_ONLY") => Mode::PrimaryOnly,
        Some("PRIMARY_MAPPED_ONLY") => Mode::PrimaryMappedOnly,
        Some("PRIMARY_PROPER_PAIR_ONLY") => Mode::PrimaryProperPairOnly,
        Some(other) => {
            eprintln!("Argument MODE has bad value: '{other}' is not a valid value");
            std::process::exit(1);
        }
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, mut records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
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

    // `getSortedRecordsFromReader`: a query-name sorted file is walked as it lies, and anything
    // else goes through a `SortingCollection` with `SAMRecordQueryNameComparator` -- the full
    // comparator with its tie-breaks, not the file-order one, so records of one name come out in
    // an order the flags decide. Both sorts are stable, so records the comparator calls equal keep
    // the order the file had.
    if header.attributes.get("SO") != Some("queryname") {
        records.sort_by(query_name::compare);
    }

    let walked: Vec<Record> = records
        .iter()
        .map(|record| Record {
            name: record.read_name.clone(),
            duplicate: record.flags & DUPLICATE_FLAG != 0,
            secondary_or_supplementary: record.flags & (SECONDARY_FLAG | SUPPLEMENTARY_FLAG) != 0,
            unmapped: record.flags & READ_UNMAPPED_FLAG != 0,
            // `getProperPairFlag()` throws on a read that is not paired, so an unpaired record
            // reaching `PRIMARY_PROPER_PAIR_ONLY` is an exception rather than a skip. Every record
            // of this corpus is paired; the flag is read as false for an unpaired one here, which
            // is the value the mode would skip on anyway.
            proper_pair: record.flags & READ_PAIRED_FLAG != 0
                && record.flags & PROPER_PAIR_FLAG != 0,
        })
        .collect();

    let verdict = check(&walked, mode);

    if let Some(output) = output {
        let mut file = std::io::BufWriter::new(std::fs::File::create(&output)?);
        for name in &verdict.bad_names {
            // `PrintWriter.println`, whose line separator on the reference's Linux is a newline.
            writeln!(file, "{name}")?;
        }
        file.flush()?;
    }

    // `return numBadRecords > 0 ? 1 : 0`, which is the count's sign and never the count.
    std::process::exit(verdict.exit_code());
}
