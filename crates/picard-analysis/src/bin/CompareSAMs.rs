//! `CompareSAMs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CompareSAMs.doWork` at tag 3.4.0. The comparison is
//! `picard_analysis::compare_sams`; this is the tool around it. The two files are positional, a
//! BAM is read as the records it holds, `OUTPUT` is the one-row metrics file with the
//! mapping-quality histogram under `COMPARE_MQ`, and the exit code is the verdict: nought when the
//! files match, one when they differ, with the file written either way.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::args::arg;
use picard_analysis::compare_sams::{compare_parsed, CompareError, CompareOptions};

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
    let args: Vec<String> = std::env::args().skip(1).collect();
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
    // The two files are read at once, one on its own thread; the left one's failure is still the
    // one reported when both fail.
    let (left, right) = std::thread::scope(|scope| {
        let right = scope.spawn(|| read(files[1]));
        (read(files[0]), right.join().expect("the reading thread"))
    });
    let unsend = |e: Box<dyn std::error::Error + Send + Sync>| -> Box<dyn std::error::Error> { e };
    let (left, right) = (left.map_err(unsend)?, right.map_err(unsend)?);
    let (left, right) = (parsed(left)?, parsed(right)?);
    let metric = compare_parsed(
        &left.0,
        &left.1,
        &right.0,
        &right.1,
        &absolute(files[0]),
        &absolute(files[1]),
        &options,
    );
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
