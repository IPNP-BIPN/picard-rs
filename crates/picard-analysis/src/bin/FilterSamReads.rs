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
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::{coordinate, query_name};
use picard_analysis::args::arg;
use picard_analysis::filter_sam_reads::{alignment_pairwise_indices, read_name_set, AlignedFilter};

/// An argument under its long name or its short one, as Barclay accepts either.
fn named(args: &[String], long: &str, short: &str) -> Option<String> {
    arg(args, &format!("{long}=")).or_else(|| arg(args, &format!("{short}=")))
}

/// `writeReadsFile`: `<OUTPUT's directory>/<basename of the file>.reads`, one `SAMRecord.toString()`
/// a line. The output's file is written with the output's header, which may carry a new `SO`.
fn write_reads_file(
    output: &str,
    file: &str,
    header: &SamHeader,
    records: &[BamRecord],
) -> std::io::Result<()> {
    let name = std::path::Path::new(file)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let base = match name.rfind('.') {
        Some(index) if index > 0 => &name[..index],
        _ => name.as_str(),
    };
    let path = std::path::Path::new(output)
        .parent()
        .unwrap_or(std::path::Path::new(""))
        .join(format!("{base}.reads"));
    std::fs::write(
        path,
        picard_analysis::filter_sam_reads::reads_file_text(header, records),
    )
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
    picard_analysis::argument_check::check("FilterSamReads");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let filter = arg(&args, "FILTER=").ok_or("FILTER= is required")?;
    let read_list = named(&args, "READ_LIST_FILE", "RLF");
    let interval_list = named(&args, "INTERVAL_LIST", "IL");
    let javascript = named(&args, "JAVASCRIPT_FILE", "JS");
    let tag = named(&args, "TAG", "T");
    let tag_values: Vec<String> = args
        .iter()
        .filter_map(|a| {
            a.strip_prefix("TAG_VALUE=")
                .or_else(|| a.strip_prefix("TV="))
        })
        .map(str::to_string)
        .collect();
    let write_reads_files = match named(&args, "WRITE_READS_FILES", "WRITE_READS_FILES") {
        None => false,
        Some(value) => match value.to_ascii_lowercase().as_str() {
            "true" | "t" | "yes" | "y" => true,
            "false" | "f" | "no" | "n" => false,
            _ => return Err(format!("unknown WRITE_READS_FILES: {value}").into()),
        },
    };

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

    // `doWork` catches every exception, logs it with the input's name, deletes the output and
    // returns 1.
    let input_name = std::path::Path::new(&input)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let fail = |message: &str| -> ! {
        eprintln!("ERROR\tFilterSamReads\tFailed to filter {input_name}");
        if !message.is_empty() {
            eprintln!("{message}");
        }
        let _ = std::fs::remove_file(&output);
        std::process::exit(1);
    };

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let raw = picard_analysis::cram_input::uncram(raw);
    let (mut header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = picard_analysis::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
        let reader =
            picard_analysis::record_validation::reader(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };
    if write_reads_files {
        write_reads_file(&output, &input, &header, &records)
            .unwrap_or_else(|e| fail(&e.to_string()));
    }
    let sort_order = named(&args, "SORT_ORDER", "SO");

    let list_text = read_list
        .as_ref()
        .map(std::fs::read_to_string)
        .transpose()?;
    // Each filter gives the indices of the records it keeps and, when it threw, the index of the
    // record it threw on. Nothing past that record is ever looked at.
    let (kept, error): (Vec<usize>, Option<(usize, String)>) = match filter.as_str() {
        "includeReadList" | "excludeReadList" => {
            let text = list_text.as_deref().unwrap_or("");
            let names = read_name_set(text);
            let include = filter == "includeReadList";
            let kept = (0..records.len())
                .filter(|&i| names.contains(records[i].read_name.as_str()) == include)
                .collect();
            (kept, None)
        }
        "includeAligned" | "excludeAligned" => alignment_pairwise_indices(
            &header,
            &records,
            if filter == "includeAligned" {
                AlignedFilter::IncludeAligned
            } else {
                AlignedFilter::ExcludeAligned
            },
        ),
        // `TagFilter`, record by record: a read is kept when whether its tag's value is one of
        // TAG_VALUE agrees with the filter's direction. The values are strings, so only a string
        // tag can equal one.
        "includeTagValues" | "excludeTagValues" => {
            let include = filter == "includeTagValues";
            let tag = tag.clone().unwrap_or_default();
            let values = &tag_values;
            let kept = (0..records.len())
                .filter(|&i| {
                    let listed = <&[u8; 2]>::try_from(tag.as_bytes())
                        .ok()
                        .is_some_and(|name| {
                            matches!(records[i].tags.get(htsjdk_bam::tag::Tag::new(name)),
                            Some(htsjdk_bam::tag::TagValue::Str(value)) if values.contains(value))
                        });
                    listed == include
                })
                .collect();
            (kept, None)
        }
        // `JavascriptSamRecordFilter`, record by record: a read the script accepts is kept. Any
        // failure, the script's included, is caught by `doWork` and logged as a failed filter.
        "includeJavascript" => {
            let path = javascript.clone().unwrap_or_default();
            let text = std::fs::read_to_string(&path).unwrap_or_else(|e| fail(&e.to_string()));
            let mut filter = picard_analysis::javascript::JavascriptFilter::new(&text)
                .unwrap_or_else(|_| fail(""));
            let mut kept = Vec::new();
            let mut error = None;
            for (i, record) in records.iter().enumerate() {
                match filter.accept_read(record, &header) {
                    Ok(true) => kept.push(i),
                    Ok(false) => {}
                    Err(_) => {
                        error = Some((i, String::new()));
                        break;
                    }
                }
            }
            (kept, error)
        }
        // `IntervalKeepPairFilter`, record by record. The list is read before the iterator is
        // built, so a bad list is the error a run reports first.
        "includePairedIntervals" => {
            let path = interval_list.clone().unwrap_or_default();
            if !std::path::Path::new(&path).is_file() {
                fail(&format!(
                    "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                    picard_analysis::metrics_cli::absolute(&path)
                ));
            }
            let text = std::fs::read_to_string(&path).unwrap_or_else(|e| fail(&e.to_string()));
            let intervals = picard_analysis::filter_sam_reads::PairedIntervals::parse(&text)
                .unwrap_or_else(|e| fail(&e));
            picard_analysis::filter_sam_reads::paired_interval_indices(
                &header, &records, &intervals,
            )
        }
        other => return Err(format!("unknown FILTER: {other}").into()),
    };

    // `FilteringSamIterator` finds its first record when it is built, so an error before the
    // first kept record is thrown before `filterReads` has looked at anything.
    if kept.is_empty() {
        if let Some((_, message)) = &error {
            fail(message);
        }
    }

    // `filterReads`: SORT_ORDER overwrites the header's order, and the writer is presorted only
    // when the input already carried the order asked for -- otherwise it sorts. The paired
    // interval filter checks the order after the overwrite, so an unsorted input asked to be
    // sorted by coordinate is accepted.
    let input_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    let order = sort_order.clone().unwrap_or_else(|| input_order.clone());
    if let Some(order) = &sort_order {
        header.set_sort_order(order);
    }
    if filter == "includePairedIntervals" && order != "coordinate" {
        fail(
            "java.lang.UnsupportedOperationException: Input must be coordinate sorted to use \
             includePairedIntervals",
        );
    }
    let presorted = order == input_order;
    let sam_output = output.ends_with(".sam");

    // The writer takes the records as the iterator hands them out, and the iterator runs one
    // kept record ahead: the record after the last kept one is where an error surfaces, so it
    // surfaces before the last kept record is written. A presorted writer checks each record
    // against the one before by the order's `fileOrderCompare`, and the message names the
    // rejected record twice, because `SAMSortOrderChecker.isSorted` has already made it the
    // previous one.
    let checked = presorted && matches!(order.as_str(), "coordinate" | "queryname");
    for (m, &index) in kept.iter().enumerate() {
        if m + 1 == kept.len() {
            if let Some((_, message)) = &error {
                fail(message);
            }
        }
        if checked && m > 0 {
            let (before, here) = (&records[kept[m - 1]], &records[index]);
            let out_of_order = if order == "coordinate" {
                let reference = |r: &BamRecord| {
                    if r.reference_index < 0 {
                        i64::MAX
                    } else {
                        r.reference_index as i64
                    }
                };
                (reference(before), before.alignment_start as i64)
                    .cmp(&(reference(here), here.alignment_start as i64))
                    .is_gt()
                    && !(before.reference_index < 0 && here.reference_index < 0)
            } else {
                query_name::compare_read_names(&before.read_name, &here.read_name).is_gt()
            };
            if out_of_order {
                let key = if order == "coordinate" {
                    let contig = usize::try_from(here.reference_index)
                        .ok()
                        .and_then(|i| header.sequences.get(i))
                        .map_or("*", |s| s.name.as_str());
                    format!("{contig}:{}", here.alignment_start)
                } else {
                    here.read_name.clone()
                };
                let file = if sam_output {
                    "null".to_string()
                } else {
                    format!("file://{}", picard_analysis::metrics_cli::absolute(&output))
                };
                fail(&format!(
                    "java.lang.IllegalArgumentException: Alignments added out of order in \
                     SAMFileWriterImpl.addAlignment for {file}. Sort order is {order}. Offending \
                     records are at [{key}] and [{key}]"
                ));
            }
        }
    }

    let mut kept: Vec<BamRecord> = kept.into_iter().map(|i| records[i].clone()).collect();
    if !presorted {
        match order.as_str() {
            "coordinate" => kept.sort_by(coordinate::compare),
            "queryname" => kept.sort_by(query_name::compare),
            "unsorted" => {}
            other => return Err(format!("unknown SORT_ORDER: {other}").into()),
        }
    }

    // The extension decides the format, and everything that is not `.sam` is a BAM.
    if sam_output {
        let text = write_sam(&header, &kept).ok_or("records failed to re-encode as SAM")?;
        std::fs::write(&output, text)?;
    } else {
        let mut writer = picard_analysis::output::bam_writer(Vec::new(), &header)
            .map_err(|e| format!("{e:?}"))?;
        for record in &kept {
            writer.write(record).map_err(|e| format!("{e:?}"))?;
        }
        // A `.cram` OUTPUT is written as CRAM, as Picard writes it.
        picard_analysis::cram_output::write(
            &output,
            writer.finish().map_err(|e| format!("{e:?}"))?,
        )?;
    }
    if write_reads_files {
        write_reads_file(&output, &output, &header, &kept).unwrap_or_else(|e| fail(&e.to_string()));
    }
    Ok(())
}
