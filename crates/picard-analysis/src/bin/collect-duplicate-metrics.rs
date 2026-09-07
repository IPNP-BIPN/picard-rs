//! `CollectDuplicateMetrics` as a runnable binary: the covering array's port side.
//!
//! It counts what `MarkDuplicates` would have counted without deciding anything: the duplicate
//! flags are read off the file rather than computed, so a file nobody has marked reports no
//! duplicates at all. That makes the tool a pure tally, and its only real argument is the driver's
//! -- the `SinglePassSamProgram` sort check, and the reference walker that refuses to rewind.
//!
//! One row per library named by the HEADER, whether a read ever used it or not, plus a row for any
//! library only the reads name.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::collect_duplicate_metrics::{collect, Record, UNKNOWN_LIBRARY};
use picard_analysis::single_pass_rejections::{
    check_sort_order, walk_reference, Rejection, SortOrder,
};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let metrics_file = arg(&args, "METRICS_FILE=")
        .or_else(|| arg(&args, "M="))
        .ok_or("METRICS_FILE= is required")?;
    let assume_sorted = arg(&args, "ASSUME_SORTED=")
        .map(|v| v == "true")
        .unwrap_or(false);

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, records) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

    let refuse = |rejection: Rejection| -> ! {
        eprintln!(
            "Exception in thread \"main\" {}: {}",
            rejection.java_class(),
            rejection.message()
        );
        std::process::exit(1);
    };

    let found = match header.attributes.get("SO") {
        Some("coordinate") => SortOrder::Coordinate,
        Some("queryname") => SortOrder::Queryname,
        Some("duplicate") => SortOrder::Duplicate,
        Some("unsorted") => SortOrder::Unsorted,
        _ => SortOrder::Unknown,
    };
    if let Err(rejection) = check_sort_order(&input, found, assume_sorted) {
        refuse(rejection);
    }

    // The driver builds a reference walker whenever REFERENCE_SEQUENCE is given, and asks it for
    // every mapped record's contig -- so assuming a queryname file is sorted moves the refusal
    // here rather than removing it.
    let with_reference = arg(&args, "REFERENCE_SEQUENCE=").is_some() || arg(&args, "R=").is_some();
    let mut current: Option<i32> = None;
    if with_reference {
        for record in &records {
            if record.reference_index == -1 {
                continue;
            }
            match walk_reference(current, record.reference_index) {
                Ok(index) => current = Some(index),
                Err(rejection) => refuse(rejection),
            }
        }
    }

    // `LibraryIdGenerator`: the libraries the header names, in header order.
    let header_libraries: Vec<String> = header
        .read_groups
        .iter()
        .map(|group| {
            group
                .attributes
                .get("LB")
                .unwrap_or(UNKNOWN_LIBRARY)
                .to_string()
        })
        .collect();

    let library_of = |record: &BamRecord| -> String {
        match record.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(id)) => header
                .read_groups
                .iter()
                .find(|group| group.id == *id)
                .and_then(|group| group.attributes.get("LB"))
                .unwrap_or(UNKNOWN_LIBRARY)
                .to_string(),
            _ => UNKNOWN_LIBRARY.to_string(),
        }
    };

    let tallied: Vec<Record> = records
        .iter()
        .map(|record| Record {
            library: library_of(record),
            duplicate: record.flags & DUPLICATE != 0,
            secondary_or_supplementary: record.flags & (SECONDARY | SUPPLEMENTARY) != 0,
            unmapped: record.flags & READ_UNMAPPED != 0,
            paired: record.flags & READ_PAIRED != 0,
            mate_unmapped: record.flags & MATE_UNMAPPED != 0,
        })
        .collect();

    let mut file = MetricsFile::new();
    file.add_header("CollectDuplicateMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    for row in collect(&header_libraries, &tallied) {
        file.add_metric(&row);
    }
    std::fs::write(&metrics_file, file.write())?;
    Ok(())
}
