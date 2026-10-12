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
use picard_analysis::revert_sam::{revert_with, Options, RevertError, SortOrder};

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

/// A collection argument as Barclay builds it: the default, each value appended, `null` emptying it.
fn collection(args: &[String], keys: &[&str], default: Vec<String>, set: bool) -> Vec<String> {
    let mut out = default;
    for a in args {
        for key in keys {
            if let Some(value) = a.strip_prefix(key) {
                if value == "null" {
                    out.clear();
                } else if !set || !out.iter().any(|v| v == value) {
                    out.push(value.to_string());
                }
            }
        }
    }
    out
}

/// `customCommandLineValidation`: every message, in its order.
fn validation(
    options: &Options,
    by_read_group: bool,
    output: Option<&str>,
    output_map: Option<&str>,
) -> Vec<String> {
    let mut errors = Vec::new();
    if options.sanitize && options.sort_order != SortOrder::Queryname {
        errors.push(
            "SORT_ORDER must be queryname when sanitization is enabled with SANITIZE=true."
                .to_string(),
        );
    }
    if by_read_group {
        match (output, output_map) {
            (Some(output), _) => {
                if !std::path::Path::new(output).is_dir() {
                    errors.push(format!(
                        "When OUTPUT_BY_READGROUP=true and OUTPUT is provided, it must be a directory: {output}"
                    ));
                }
            }
            (None, None) => errors.push(
                "Must provide either OUTPUT or OUTPUT_MAP when OUTPUT_BY_READGROUP=true."
                    .to_string(),
            ),
            (None, Some(map)) => match std::fs::read_to_string(map) {
                Err(_) => errors.push(format!("Cannot read OUTPUT_MAP {map}")),
                Ok(text) => {
                    let labels: Vec<&str> = text.lines().next().unwrap_or("").split('\t').collect();
                    if labels.len() < 2 || labels[0] != "READ_GROUP_ID" || labels[1] != "OUTPUT" {
                        errors.push(format!("Invalid header: {map}. Must be a tab-separated file with READ_GROUP_ID as first column and OUTPUT as second column."));
                    }
                }
            },
        }
    } else {
        if output_map.is_some() {
            errors.push(
                "Cannot provide OUTPUT_MAP when OUTPUT_BY_READGROUP=false. Provide OUTPUT instead."
                    .to_string(),
            );
        }
        match output {
            None => errors.push("OUTPUT is required when OUTPUT_BY_READGROUP=false".to_string()),
            Some(output) if std::path::Path::new(output).is_dir() => errors.push(format!(
                "OUTPUT {output} should not be a directory when OUTPUT_BY_READGROUP=false"
            )),
            Some(_) => {}
        }
    }
    if !options.sanitize && options.keep_first_duplicate {
        errors.push("KEEP_FIRST_DUPLICATE cannot be used without SANITIZE".to_string());
    }
    errors
}

fn thrown(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
    std::process::exit(1);
}

/// The file the writer factory makes for a path: SAM text for `.sam`, a BAM otherwise.
fn write_records(path: &str, header: &SamHeader, records: &[BamRecord]) -> Result<(), String> {
    if path.ends_with(".sam") {
        let text = htsjdk_bam::sam_file::write_sam(header, records)
            .ok_or("records failed to re-encode as SAM")?;
        return std::fs::write(path, text).map_err(|e| format!("{e:?}"));
    }
    write_bam(path, header, records)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=").or_else(|| arg(&args, "O="));
    let output_map = arg(&args, "OUTPUT_MAP=");
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);
    let by_read_group = flag("OUTPUT_BY_READGROUP=", false) || flag("OBR=", false);
    let sort_order = match arg(&args, "SORT_ORDER=").as_deref() {
        None | Some("queryname") => SortOrder::Queryname,
        Some("coordinate") => SortOrder::Coordinate,
        Some(other) => return Err(format!("unknown SORT_ORDER: {other}").into()),
    };
    // Barclay's own checks come first: OUTPUT and OUTPUT_MAP exclude each other, and OUTPUT is
    // required unless OUTPUT_MAP is given.
    if output.is_some() && output_map.is_some() {
        eprintln!("Argument 'OUTPUT' cannot be used in conjunction with argument(s) OUTPUT_MAP");
        std::process::exit(1);
    }
    if output.is_none() && output_map.is_none() {
        eprintln!("Argument OUTPUT was missing: Argument 'OUTPUT' is required unless one of {{[OUTPUT_MAP]}} are provided");
        std::process::exit(1);
    }
    let (clear, reverse, complement) = picard_analysis::revert_sam::default_tags();
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
        attributes_to_clear: collection(&args, &["ATTRIBUTE_TO_CLEAR="], clear, false),
        attributes_to_reverse: collection(&args, &["ATTRIBUTE_TO_REVERSE=", "RV="], reverse, true),
        attributes_to_reverse_complement: collection(
            &args,
            &["ATTRIBUTE_TO_REVERSE_COMPLEMENT=", "RC="],
            complement,
            true,
        ),
    };

    // Barclay validates the command line before the tool is constructed, prints the usage, and
    // then every message it collected. The harness records that trailing block, so the messages
    // are printed here in the same order and nothing else is written.
    let errors = validation(
        &options,
        by_read_group,
        output.as_deref(),
        output_map.as_deref(),
    );
    if !errors.is_empty() {
        for message in &errors {
            eprintln!("{message}");
        }
        std::process::exit(1);
    }

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (mut header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
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

    // `validateHeaderOverrides`: an override needs every read group to agree on the field.
    let sample_alias = arg(&args, "SAMPLE_ALIAS=").or_else(|| arg(&args, "ALIAS="));
    let library_name = arg(&args, "LIBRARY_NAME=").or_else(|| arg(&args, "LIB="));
    let field = |group: &htsjdk_bam::header::ReadGroup, key: &str| {
        group.attributes.get(key).map(str::to_string)
    };
    if sample_alias.is_some() || library_name.is_some() {
        let groups = &header.read_groups;
        let same = |key: &str| {
            groups
                .iter()
                .skip(1)
                .all(|g| field(g, key) == field(&groups[0], key))
        };
        if sample_alias.is_some() && !same("SM") {
            thrown("Read groups have multiple values for sample.  A value for SAMPLE_ALIAS cannot be supplied.");
        }
        if library_name.is_some() && !same("LB") {
            thrown("Read groups have multiple values for library name.  A value for library name cannot be supplied.");
        }
    }
    for group in &mut header.read_groups {
        if let Some(sample) = &sample_alias {
            group.attributes.set("SM", sample);
        }
        if let Some(library) = &library_name {
            group.attributes.set("LB", library);
        }
    }

    // `OUTPUT_BY_READGROUP`: a writer per read group, at the path the map gives it.
    let mut outputs: Vec<(String, String)> = Vec::new();
    if by_read_group {
        if header.read_groups.is_empty() {
            thrown(&format!("{input} does not contain Read Groups"));
        }
        let extension = match arg(&args, "OUTPUT_BY_READGROUP_FILE_FORMAT=").as_deref() {
            None | Some("dynamic") => {
                if input.ends_with(".sam") {
                    ".sam".to_string()
                } else if input.ends_with(".cram") {
                    ".cram".to_string()
                } else {
                    ".bam".to_string()
                }
            }
            Some(format) => format!(".{format}"),
        };
        match &output_map {
            Some(map) => {
                let text = std::fs::read_to_string(map)?;
                for line in text.lines().skip(1).filter(|l| !l.is_empty()) {
                    let mut fields = line.split('\t');
                    let id = fields.next().unwrap_or("").to_string();
                    let path = fields.next().unwrap_or("");
                    outputs.push((id, picard_analysis::fingerprint::host_path(path)));
                }
            }
            None => {
                let directory = output.clone().unwrap_or_default();
                for group in &header.read_groups {
                    let path =
                        std::path::Path::new(&directory).join(format!("{}{extension}", group.id));
                    outputs.push((group.id.clone(), path.display().to_string()));
                }
            }
        }
    }

    // This one is a PicardException from doWork, thrown after the reader is opened and before any
    // record is written, so a refused row leaves no output behind.
    if options.restore_hardclips && !options.remove_alignment_information {
        thrown(&RevertError::HardclipsWithoutRemovingAlignment.message());
    }

    // The discard-fraction refusal happens after `out.close()`, so the output exists and then the
    // tool throws. Writing before the exit is not tidiness: a row that stopped short of the write
    // would differ from the reference on the file as well as on the message.
    let reverted = revert_with(&header, records, &options);
    if by_read_group {
        for (id, path) in &outputs {
            // `createHeaderMap`: the single header's fields, with this one read group.
            let mut one = reverted.header.clone();
            one.read_groups.retain(|group| &group.id == id);
            let mine: Vec<BamRecord> = reverted
                .records
                .iter()
                .filter(|record| {
                    matches!(record.tags.get(htsjdk_bam::tag::Tag::new(b"RG")),
                        Some(htsjdk_bam::tag::TagValue::Str(rg)) if rg == id)
                })
                .cloned()
                .collect();
            write_records(path, &one, &mine)?;
        }
    } else {
        write_records(
            output.as_deref().unwrap_or_default(),
            &reverted.header,
            &reverted.records,
        )?;
    }
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
