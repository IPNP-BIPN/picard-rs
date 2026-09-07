//! `RevertSam` as a runnable binary: the covering array's port side.
//!
//! The library reverted on the default path only. This adds the binary and the arguments the array
//! varies, which for this tool are mostly refusals: four of its seventeen rows produce a file, and
//! the other thirteen are the tool saying no in four different ways. Reproducing a refusal is
//! reproducing the tool, so the messages matter as much as the bytes.
//!
//! Two of them come from Barclay's validation, which collects every failure and prints them after
//! the usage block rather than throwing; the other two are `PicardException`s thrown from
//! `doWork`, one before any record is read and one after the whole file has been written.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::revert_sam::{revert_with, validate, Options, RevertError, SortOrder};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn write_bam(path: &str, header: &SamHeader, records: &[BamRecord]) -> Result<(), String> {
    let mut writer = BamWriter::new(Vec::new(), header).map_err(|e| format!("{e:?}"))?;
    for record in records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    let bytes = writer.finish().map_err(|e| format!("{e:?}"))?;
    std::fs::write(path, bytes).map_err(|e| format!("{e:?}"))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);
    let sort_order = match arg(&args, "SORT_ORDER=").as_deref() {
        None | Some("queryname") => SortOrder::Queryname,
        Some("coordinate") => SortOrder::Coordinate,
        Some(other) => return Err(format!("unknown SORT_ORDER: {other}").into()),
    };
    let options = Options {
        sort_order,
        restore_original_qualities: flag("RESTORE_ORIGINAL_QUALITIES=", true),
        remove_duplicate_information: flag("REMOVE_DUPLICATE_INFORMATION=", true),
        remove_alignment_information: flag("REMOVE_ALIGNMENT_INFORMATION=", true),
        restore_hardclips: flag("RESTORE_HARDCLIPS=", false),
        sanitize: flag("SANITIZE=", false),
        keep_first_duplicate: flag("KEEP_FIRST_DUPLICATE=", false),
        max_discard_fraction: arg(&args, "MAX_DISCARD_FRACTION=")
            .map(|v| v.parse::<f64>())
            .transpose()?
            .unwrap_or(0.01),
    };

    // Barclay validates the command line before the tool is constructed, prints the usage, and
    // then every message it collected. The harness records that trailing block, so the messages
    // are printed here in the same order and nothing else is written.
    let errors = validate(
        &options,
        flag("OUTPUT_BY_READGROUP=", false),
        std::path::Path::new(&output).is_dir(),
        &output,
    );
    if !errors.is_empty() {
        for message in &errors {
            eprintln!("{message}");
        }
        std::process::exit(1);
    }

    // This one is a PicardException from doWork, thrown after the reader is opened and before any
    // record is written, so a refused row leaves no output behind.
    if options.restore_hardclips && !options.remove_alignment_information {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: {}",
            RevertError::HardclipsWithoutRemovingAlignment.message()
        );
        std::process::exit(1);
    }

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

    // The discard-fraction refusal happens after `out.close()`, so the output exists and then the
    // tool throws. Writing before the exit is not tidiness: a row that stopped short of the write
    // would differ from the reference on the file as well as on the message.
    let reverted = revert_with(&header, records, &options);
    write_bam(&output, &reverted.header, &reverted.records)?;
    if let Some(error) = reverted.error {
        eprintln!(
            "Exception in thread \"main\" {}: {}",
            error.java_class(),
            error.message()
        );
        std::process::exit(1);
    }
    Ok(())
}
