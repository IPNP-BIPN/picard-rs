//! `FilterSamReads` as a runnable binary: the covering array's port side.
//!
//! Eight filters, and each one takes a different argument that the other seven refuse. A covering
//! array carries every argument its domain holds on every row, so this tool cannot have all eight
//! measured at once: the repository declares a `READ_LIST_FILE` and no `INTERVAL_LIST`, `TAG` or
//! `JAVASCRIPT_FILE`, which makes the two read-list filters run and the other six refuse. Six of
//! the twenty-five rows produce a file, four of them distinct; nineteen are the tool saying no in
//! six different ways, and the messages name both the argument and the filter that clashed with
//! it.
//!
//! The refusals are Barclay's, collected and printed after the usage block rather than thrown.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::{coordinate, query_name};
use picard_analysis::filter_sam_reads::{
    keep_by_alignment_pairwise, keep_by_read_list, read_name_set, AlignedFilter, Filter,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `checkInputs`: an argument this filter does not take is a refusal, and so is one it needs and
/// has not been given. The wording is the reference's, list of filters included.
fn check_inputs(
    filters: &[&str],
    filter: &str,
    value: Option<&String>,
    name: &str,
) -> Option<String> {
    let listed = filters.contains(&filter);
    match (listed, value) {
        (true, None) => Some(format!(
            "{name} must be specified when using FILTER={filter}, but it was null."
        )),
        (false, Some(value)) => Some(format!(
            "{name} may only be specified when using FILTER from {}, FILTER value: {filter}, \
             {name} value: {value}",
            filters.join(", ")
        )),
        _ => None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let filter = arg(&args, "FILTER=").ok_or("FILTER= is required")?;
    let read_list = arg(&args, "READ_LIST_FILE=");
    let interval_list = arg(&args, "INTERVAL_LIST=");
    let javascript = arg(&args, "JAVASCRIPT_FILE=");
    let tag = arg(&args, "TAG=");

    let mut errors: Vec<String> = Vec::new();
    if input == output {
        errors.push("INPUT file and OUTPUT file must differ!".to_string());
    }
    for (filters, value, name) in [
        (
            vec!["includeReadList", "excludeReadList"],
            read_list.as_ref(),
            "READ_LIST_FILE",
        ),
        (
            vec!["includePairedIntervals"],
            interval_list.as_ref(),
            "INTERVAL_LIST",
        ),
        (
            vec!["includeJavascript"],
            javascript.as_ref(),
            "JAVASCRIPT_FILE",
        ),
        (
            vec!["includeTagValues", "excludeTagValues"],
            tag.as_ref(),
            "TAG",
        ),
    ] {
        if let Some(message) = check_inputs(&filters, &filter, value, name) {
            errors.push(message);
        }
    }
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

    let list_text = read_list
        .as_ref()
        .map(std::fs::read_to_string)
        .transpose()?;
    let mut kept = match filter.as_str() {
        "includeReadList" | "excludeReadList" => {
            let text = list_text.as_deref().unwrap_or("");
            let names = read_name_set(text);
            let which = if filter == "includeReadList" {
                Filter::IncludeReadList
            } else {
                Filter::ExcludeReadList
            };
            keep_by_read_list(&records, &names, which)
        }
        "includeAligned" | "excludeAligned" => match keep_by_alignment_pairwise(
            &records,
            if filter == "includeAligned" {
                AlignedFilter::IncludeAligned
            } else {
                AlignedFilter::ExcludeAligned
            },
        ) {
            Ok(kept) => kept,
            // `doWork` catches the assertion, logs it, deletes the output and returns 1.
            Err(_) => {
                let name = std::path::Path::new(&input)
                    .file_name()
                    .map(|name| name.to_string_lossy().into_owned())
                    .unwrap_or_default();
                eprintln!("ERROR\tFilterSamReads\tFailed to filter {name}");
                let _ = std::fs::remove_file(&output);
                std::process::exit(1);
            }
        },
        // `TagFilter`, record by record: a read is kept when whether its tag's value is one of
        // TAG_VALUE agrees with the filter's direction. The values are strings, so only a string
        // tag can equal one.
        "includeTagValues" | "excludeTagValues" => {
            let include = filter == "includeTagValues";
            let tag = arg(&args, "TAG=").unwrap_or_default();
            let values: Vec<String> = args
                .iter()
                .filter_map(|a| a.strip_prefix("TAG_VALUE="))
                .map(str::to_string)
                .collect();
            records
                .iter()
                .filter(|record| {
                    let listed = <&[u8; 2]>::try_from(tag.as_bytes())
                        .ok()
                        .is_some_and(|name| {
                            matches!(record.tags.get(htsjdk_bam::tag::Tag::new(name)),
                            Some(htsjdk_bam::tag::TagValue::Str(value)) if values.contains(value))
                        });
                    listed == include
                })
                .cloned()
                .collect()
        }
        // `JavascriptSamRecordFilter`, record by record: a read the script accepts is kept. Any
        // failure, the script's included, is caught by `doWork` and logged as a failed filter.
        "includeJavascript" => {
            let fail = || -> ! {
                let name = std::path::Path::new(&input)
                    .file_name()
                    .map(|name| name.to_string_lossy().into_owned())
                    .unwrap_or_default();
                eprintln!("ERROR\tFilterSamReads\tFailed to filter {name}");
                let _ = std::fs::remove_file(&output);
                std::process::exit(1);
            };
            let path = arg(&args, "JAVASCRIPT_FILE=").unwrap_or_default();
            let text = std::fs::read_to_string(&path).unwrap_or_else(|_| fail());
            let mut filter = picard_analysis::javascript::JavascriptFilter::new(&text)
                .unwrap_or_else(|_| fail());
            let mut kept = Vec::new();
            for record in &records {
                if filter
                    .accept_read(record, &header)
                    .unwrap_or_else(|_| fail())
                {
                    kept.push(record.clone());
                }
            }
            kept
        }
        other => return Err(format!("FILTER={other} is not ported").into()),
    };

    // `filterReads`: SORT_ORDER overwrites the header's order, and the writer is presorted only
    // when the input already carried the order asked for -- otherwise it sorts.
    let input_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    if let Some(order) = arg(&args, "SORT_ORDER=").or_else(|| arg(&args, "SO=")) {
        header.set_sort_order(&order);
        if order != input_order {
            match order.as_str() {
                "coordinate" => kept.sort_by(coordinate::compare),
                "queryname" => kept.sort_by(query_name::compare),
                other => return Err(format!("unknown SORT_ORDER: {other}").into()),
            }
        }
    }

    let mut writer =
        picard_analysis::output::bam_writer(Vec::new(), &header).map_err(|e| format!("{e:?}"))?;
    for record in &kept {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;
    Ok(())
}
