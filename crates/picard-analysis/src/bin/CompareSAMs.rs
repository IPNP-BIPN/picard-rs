//! `CompareSAMs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CompareSAMs.doWork` at tag 3.4.0. The comparison is
//! `picard_analysis::compare_sams`; this is the tool around it. The two files are positional, a
//! BAM is compared as its records decode, neither file held, `OUTPUT` is the one-row metrics file with the
//! mapping-quality histogram under `COMPARE_MQ`, and the exit code is the verdict: nought when the
//! files match, one when they differ, with the file written either way.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::args::arg;
use picard_analysis::compare_sams::SamComparisonMetric;
use picard_analysis::compare_sams::{compare_records, CompareError, CompareOptions};
use std::sync::mpsc;

/// How many records a decoding thread hands over at a time, and how many handfuls it may get
/// ahead of the comparison.
const BATCH: usize = 4096;
const AHEAD: usize = 4;

/// A file as `CompareSAMs` reads it: a BAM (or a CRAM, decoded as one) as its records decode, a
/// SAM as text for `compare_sams_with_options` to parse.
enum Input {
    Records(SamHeader, Vec<BamRecord>),
    Text(String),
}

fn read(path: &str) -> Result<Input, Box<dyn std::error::Error + Send + Sync>> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    let raw = picard_analysis::cram_input::uncram(raw);
    if !raw.starts_with(&[0x1f, 0x8b]) {
        return Ok(Input::Text(String::from_utf8(raw)?));
    }
    let plain = picard_analysis::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
    let reader =
        picard_analysis::record_validation::reader(&plain).map_err(|e| format!("{e:?}"))?;
    let header = reader.header.text.clone();
    let records: Vec<_> = reader
        .map(|r| r.map_err(|e| format!("{e:?}")))
        .collect::<Result<_, _>>()?;
    Ok(Input::Records(header, records))
}

/// Both BAMs compared as they decode, each on its own thread feeding the comparison, or `None`
/// when either is not a BAM or its single pass could not vouch for what the whole-file read
/// answers: the run then reads both whole, as it did.
fn streamed(
    files: [&str; 2],
    names: [&str; 2],
    options: &CompareOptions,
) -> Option<SamComparisonMetric> {
    for path in files {
        let mut magic = Vec::new();
        std::fs::File::open(path)
            .ok()?
            .take(2)
            .read_to_end(&mut magic)
            .ok()?;
        // shortcut: a SAM or CRAM is read whole, as the other tools read one; stream it if needed.
        if magic != [0x1f, 0x8b] {
            return None;
        }
    }
    std::thread::scope(|scope| {
        let (h1, left, left_clean) = decode(scope, files[0]);
        let (h2, right, right_clean) = decode(scope, files[1]);
        // The closure owns the record channels, and dropping them lets the decoders run on to
        // their answers whether or not the comparison read every record.
        let metric = (move || {
            Some(compare_records(
                &h1?, left, &h2?, right, names[0], names[1], options,
            ))
        })();
        let clean = left_clean.join().expect("the left decoder");
        let clean = right_clean.join().expect("the right decoder") && clean;
        metric.filter(|_| clean)
    })
}

/// One BAM decoded on its own thread: its header, its records as they decode, and whether the
/// single pass vouched for them. The records are dropped once the comparison stops reading them,
/// and the rest of the file is still decoded for its errors.
fn decode<'scope>(
    scope: &'scope std::thread::Scope<'scope, '_>,
    path: &'scope str,
) -> (
    Option<SamHeader>,
    impl Iterator<Item = BamRecord> + 'scope,
    std::thread::ScopedJoinHandle<'scope, bool>,
) {
    let (header_tx, header_rx) = mpsc::sync_channel(1);
    let (tx, rx) = mpsc::sync_channel::<Vec<BamRecord>>(AHEAD);
    let decoder = scope.spawn(move || {
        let mut batch = Vec::with_capacity(BATCH);
        let mut sending = true;
        let clean = picard_analysis::bam_stream::one_pass(
            path,
            |header| drop(header_tx.send(header.clone())),
            |record| {
                if !sending {
                    return;
                }
                batch.push(record);
                if batch.len() == BATCH {
                    let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH));
                    sending = tx.send(full).is_ok();
                }
            },
        );
        if sending && !batch.is_empty() {
            let _ = tx.send(batch);
        }
        clean
    });
    (header_rx.recv().ok(), rx.into_iter().flatten(), decoder)
}

/// Both files read whole, at once, one on its own thread; the left one's failure is still the one
/// reported when both fail.
fn whole(
    files: [&str; 2],
    names: [&str; 2],
    options: &CompareOptions,
) -> Result<SamComparisonMetric, Box<dyn std::error::Error>> {
    let (left, right) = std::thread::scope(|scope| {
        let right = scope.spawn(|| read(files[1]));
        (read(files[0]), right.join().expect("the reading thread"))
    });
    let unsend = |e: Box<dyn std::error::Error + Send + Sync>| -> Box<dyn std::error::Error> { e };
    let (left, right) = (left.map_err(unsend)?, right.map_err(unsend)?);
    let (left, right) = (parsed(left)?, parsed(right)?);
    Ok(compare_records(
        &left.0,
        left.1.into_iter(),
        &right.0,
        right.1.into_iter(),
        names[0],
        names[1],
        options,
    ))
}

/// The records of a SAM read as `compare_sams_with_options` reads them.
fn parsed(input: Input) -> Result<(SamHeader, Vec<BamRecord>), Box<dyn std::error::Error>> {
    match input {
        Input::Records(header, records) => Ok((header, records)),
        // The same refusal text the parse inside `compare_sams_with_options` gave.
        Input::Text(text) => Ok(read_sam_with(&text, ValidationStringency::Silent)
            .map_err(|e| format!("{:?}", CompareError::from(e)))?),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CompareSAMs");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process()
        .into_iter()
        .skip(1)
        .collect();
    let files: Vec<&String> = args.iter().filter(|a| !a.contains('=')).collect();
    if files.len() != 2 {
        return Err("exactly two SAM files are required".into());
    }
    let flag = |key: &str| {
        arg(&args, key)
            .map(|v| v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    };
    let options = CompareOptions {
        lenient_low_mq_alignment: flag("LENIENT_LOW_MQ_ALIGNMENT="),
        lenient_unknown_mq_alignment: flag("LENIENT_UNKNOWN_MQ_ALIGNMENT="),
        low_mq_threshold: arg(&args, "LOW_MQ_THRESHOLD=")
            .and_then(|v| v.parse().ok())
            .unwrap_or(3),
        compare_mq: flag("COMPARE_MQ="),
        lenient_header: flag("LENIENT_HEADER="),
        lenient_dup: flag("LENIENT_DUP="),
    };
    let absolute = |path: &str| {
        if std::path::Path::new(path).is_absolute() {
            path.to_string()
        } else {
            std::env::current_dir()
                .map(|cwd| cwd.join(path).display().to_string())
                .unwrap_or_else(|_| path.to_string())
        }
    };
    let (left_name, right_name) = (absolute(files[0]), absolute(files[1]));
    let metric = match streamed([files[0], files[1]], [&left_name, &right_name], &options) {
        Some(metric) => metric,
        None => whole([files[0], files[1]], [&left_name, &right_name], &options)?,
    };
    if let Some(exception) = &metric.thrown {
        for line in &metric.report {
            println!("{line}");
        }
        eprintln!("Exception in thread \"main\" {exception}");
        std::process::exit(1);
    }
    if let Some(output) = arg(&args, "OUTPUT=").or_else(|| arg(&args, "O=")) {
        let mut file = MetricsFile::new();
        for header in picard_analysis::command_line::default_headers("CompareSAMs") {
            file.add_header(&header);
        }
        file.add_metric(&metric);
        // `COMPARE_MQ`'s concordance histogram, keyed `mq1,mq2`, which the writer drops when empty.
        if !metric.mq_histogram.is_empty() {
            file.histograms.push(Histogram {
                bin_label: "BIN".to_string(),
                value_label: "VALUE".to_string(),
                key_class: "java.lang.String".to_string(),
                bins: metric
                    .mq_histogram
                    .iter()
                    .map(|(k, v)| (k.clone(), *v as f64))
                    .collect(),
            });
        }
        std::fs::write(output, file.write())?;
    }
    for line in &metric.report {
        println!("{line}");
    }
    println!("{}", picard_analysis::compare_sams::verdict(&metric));
    std::process::exit(if metric.are_equal { 0 } else { 1 });
}
