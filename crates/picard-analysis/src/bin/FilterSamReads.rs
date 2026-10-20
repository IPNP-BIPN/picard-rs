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

use std::io::Write;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::{coordinate, query_name};
use picard_analysis::args::arg;
use picard_analysis::filter_sam_reads::{
    aligned_pairwise_step, read_name_set, AlignedFilter, Step,
};
use picard_analysis::stream_output::Output;

/// An argument under its long name or its short one, as Barclay accepts either.
fn named(args: &[String], long: &str, short: &str) -> Option<String> {
    arg(args, &format!("{long}=")).or_else(|| arg(args, &format!("{short}=")))
}

/// `writeReadsFile`'s file: `<OUTPUT's directory>/<basename of the file>.reads`.
fn reads_file_path(output: &str, file: &str) -> std::path::PathBuf {
    let name = std::path::Path::new(file)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let base = match name.rfind('.') {
        Some(index) if index > 0 => &name[..index],
        _ => name.as_str(),
    };
    std::path::Path::new(output)
        .parent()
        .unwrap_or(std::path::Path::new(""))
        .join(format!("{base}.reads"))
}

/// `writeReadsFile`'s lines, one `SAMRecord.toString()` each, written as the records go by. The
/// whole-file text was empty when any record could not be written as SAM, and so is the file.
struct ReadsFile {
    path: std::path::PathBuf,
    out: Option<std::io::BufWriter<std::fs::File>>,
}

impl ReadsFile {
    fn create(path: std::path::PathBuf) -> std::io::Result<ReadsFile> {
        let out = Some(std::io::BufWriter::new(std::fs::File::create(&path)?));
        Ok(ReadsFile { path, out })
    }

    fn write(&mut self, header: &SamHeader, record: &BamRecord) -> std::io::Result<()> {
        let Some(out) = &mut self.out else {
            return Ok(());
        };
        match picard_analysis::bam_stream::sam_line(header, record) {
            Some(line) => writeln!(out, "{line}"),
            None => {
                self.out = None;
                std::fs::File::create(&self.path).map(drop)
            }
        }
    }

    fn finish(self) -> std::io::Result<()> {
        self.out.map_or(Ok(()), |mut out| out.flush())
    }
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
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
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
    // The output's reads file is written beside it, and only kept when the run succeeds.
    let reads_path = reads_file_path(&output, &output);
    let reads_temporary = reads_path.with_extension("reads.partial");
    let fail = |message: &str| -> ! {
        let _ = std::fs::remove_file(&reads_temporary);
        eprintln!("ERROR\tFilterSamReads\tFailed to filter {input_name}");
        if !message.is_empty() {
            eprintln!("{message}");
        }
        let _ = std::fs::remove_file(&output);
        std::process::exit(1);
    };

    let opened = picard_analysis::bam_stream::Input::open(&input)?;
    let streamed = opened.streams();
    let (mut header, records) = opened.into_parts();
    // The input's reads file holds every record, written before the filter looks at any: a pass of
    // its own.
    let records = if write_reads_files {
        let mut reads = ReadsFile::create(reads_file_path(&output, &input))
            .unwrap_or_else(|e| fail(&e.to_string()));
        for record in records {
            reads
                .write(&header, &record)
                .unwrap_or_else(|e| fail(&e.to_string()));
        }
        reads.finish().unwrap_or_else(|e| fail(&e.to_string()));
        picard_analysis::bam_stream::Input::reopen(&input, streamed)?
            .into_parts()
            .1
    } else {
        records
    };
    let sort_order = named(&args, "SORT_ORDER", "SO");

    let list_text = read_list
        .as_ref()
        .map(std::fs::read_to_string)
        .transpose()?;
    let names = list_text.as_deref().map(read_name_set);
    // The filters read the input's header, before SORT_ORDER rewrites it.
    let input_header = header.clone();
    let mut javascript_filter = None;
    let mut paired_intervals = None;
    match filter.as_str() {
        "includeReadList" | "excludeReadList" | "includeAligned" | "excludeAligned"
        | "includeTagValues" | "excludeTagValues" => {}
        // `JavascriptSamRecordFilter`: any failure, the script's included, is caught by `doWork`
        // and logged as a failed filter.
        "includeJavascript" => {
            let path = javascript.clone().unwrap_or_default();
            let text = std::fs::read_to_string(&path).unwrap_or_else(|e| fail(&e.to_string()));
            javascript_filter = Some(
                picard_analysis::javascript::JavascriptFilter::new(&text)
                    .unwrap_or_else(|_| fail("")),
            );
        }
        // `IntervalKeepPairFilter`. The list is read before the iterator is built, so a bad list
        // is the error a run reports first.
        "includePairedIntervals" => {
            let path = interval_list.clone().unwrap_or_default();
            if !std::path::Path::new(&path).is_file() {
                fail(&format!(
                    "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                    picard_analysis::metrics_cli::absolute(&path)
                ));
            }
            let text = std::fs::read_to_string(&path).unwrap_or_else(|e| fail(&e.to_string()));
            paired_intervals = Some(
                picard_analysis::filter_sam_reads::PairedIntervals::parse(&text)
                    .unwrap_or_else(|e| fail(&e)),
            );
        }
        other => return Err(format!("unknown FILTER: {other}").into()),
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
    let unsupported = filter == "includePairedIntervals" && order != "coordinate";
    let presorted = order == input_order;
    let sam_output = output.ends_with(".sam");
    let checked = presorted && matches!(order.as_str(), "coordinate" | "queryname");
    // A presorted writer checks each record against the one before by the order's
    // `fileOrderCompare`, and the message names the rejected record twice, because
    // `SAMSortOrderChecker.isSorted` has already made it the previous one.
    let check_order = |before: &BamRecord, here: &BamRecord| {
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
    };

    // The extension decides the format, and everything that is not `.sam` is a BAM (a `.cram` is
    // written as CRAM, as Picard writes it). Any failure removes it, as `doWork` deletes it.
    let mut out = Output::create(&output, &header, false)?;
    let mut reads = if write_reads_files && presorted {
        Some(ReadsFile::create(reads_temporary.clone()).unwrap_or_else(|e| fail(&e.to_string())))
    } else {
        None
    };
    let mut to_sort: Vec<BamRecord> = Vec::new();
    let mut last_written: Option<BamRecord> = None;
    let mut write = |record: BamRecord| -> Result<(), Box<dyn std::error::Error>> {
        if !presorted {
            to_sort.push(record);
            return Ok(());
        }
        if checked {
            if let Some(before) = &last_written {
                check_order(before, &record);
            }
        }
        out.write(&record)?;
        if let Some(reads) = &mut reads {
            reads
                .write(&header, &record)
                .unwrap_or_else(|e| fail(&e.to_string()));
        }
        last_written = Some(record);
        Ok(())
    };

    // Each filter keeps records or throws on one, and nothing past that record is ever looked at.
    // The writer takes the records as the iterator hands them out, and the iterator runs one kept
    // record ahead: the record after the last kept one is where an error surfaces, so it surfaces
    // before the last kept record is written. `FilteringSamIterator` finds its first record when
    // it is built, so an error before the first kept record is thrown before `filterReads` has
    // looked at anything, and the paired interval filter's order check comes after that.
    let mut pending: Option<BamRecord> = None;
    let mut error: Option<String> = None;
    let mut records = records.peekable();
    while let Some(record) = records.next() {
        let keep = match filter.as_str() {
            "includeReadList" | "excludeReadList" => {
                let listed = names
                    .as_ref()
                    .is_some_and(|names| names.contains(record.read_name.as_str()));
                Ok(listed == (filter == "includeReadList"))
            }
            "includeAligned" | "excludeAligned" => {
                let which = if filter == "includeAligned" {
                    AlignedFilter::IncludeAligned
                } else {
                    AlignedFilter::ExcludeAligned
                };
                match aligned_pairwise_step(&input_header, &record, records.peek(), which) {
                    Ok(Step::Keep) => Ok(true),
                    Ok(Step::Drop) => Ok(false),
                    Ok(Step::DropWithMate) => {
                        records.next();
                        Ok(false)
                    }
                    Err(message) => Err(message),
                }
            }
            // `TagFilter`, record by record: a read is kept when whether its tag's value is one
            // of TAG_VALUE agrees with the filter's direction. The values are strings, so only a
            // string tag can equal one.
            "includeTagValues" | "excludeTagValues" => {
                let tag = tag.as_deref().unwrap_or_default();
                let listed = <&[u8; 2]>::try_from(tag.as_bytes())
                    .ok()
                    .is_some_and(|name| {
                        matches!(record.tags.get(htsjdk_bam::tag::Tag::new(name)),
                        Some(htsjdk_bam::tag::TagValue::Str(value)) if tag_values.contains(value))
                    });
                Ok(listed == (filter == "includeTagValues"))
            }
            "includeJavascript" => javascript_filter
                .as_mut()
                .expect("built above")
                .accept_read(&record, &input_header)
                .map_err(|_| String::new()),
            _ => paired_intervals
                .as_ref()
                .expect("parsed above")
                .keeps(&input_header, &record),
        };
        match keep {
            Ok(true) => {
                if pending.is_none() && unsupported {
                    break;
                }
                if let Some(previous) = pending.replace(record) {
                    write(previous)?;
                }
            }
            Ok(false) => {}
            Err(message) => {
                error = Some(message);
                break;
            }
        }
    }
    if let Some(message) = &error {
        fail(message);
    }
    if unsupported {
        fail(
            "java.lang.UnsupportedOperationException: Input must be coordinate sorted to use \
             includePairedIntervals",
        );
    }
    if let Some(last) = pending {
        write(last)?;
    }

    if !presorted {
        match order.as_str() {
            "coordinate" => to_sort.sort_by(coordinate::compare),
            "queryname" => to_sort.sort_by(query_name::compare),
            "unsorted" => {}
            other => return Err(format!("unknown SORT_ORDER: {other}").into()),
        }
        if write_reads_files {
            reads = Some(
                ReadsFile::create(reads_temporary.clone()).unwrap_or_else(|e| fail(&e.to_string())),
            );
        }
        for record in &to_sort {
            out.write(record)?;
            if let Some(reads) = &mut reads {
                reads
                    .write(&header, record)
                    .unwrap_or_else(|e| fail(&e.to_string()));
            }
        }
    }
    // `cram_output::write`: a CRAM takes `CREATE_INDEX` from the command line, a BAM no index.
    out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
    if let Some(reads) = reads {
        reads.finish().unwrap_or_else(|e| fail(&e.to_string()));
        std::fs::rename(&reads_temporary, &reads_path).unwrap_or_else(|e| fail(&e.to_string()));
    }
    Ok(())
}
