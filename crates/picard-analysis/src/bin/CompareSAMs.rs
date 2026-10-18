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
use htsjdk_bam::text_parse::{parse_line_with, ValidationStringency};
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::args::arg;
use picard_analysis::bam_stream::{Input, InputError};
use picard_analysis::compare_sams::{compare_records, CompareError, CompareOptions};

type Records = Box<dyn Iterator<Item = Result<BamRecord, String>>>;

/// The file's header and records as the comparison reads them, from its SAM text at SILENT: a
/// BAM's records are written as SAM lines and parsed back one at a time, as the whole text was.
/// A record that fails is the run's error, with the message the whole-file read gave.
fn as_sam(path: &str) -> Result<(SamHeader, Records), Box<dyn std::error::Error>> {
    let mut magic = Vec::new();
    std::fs::File::open(path)?.take(4).read_to_end(&mut magic)?;
    let bam = magic.starts_with(&[0x1f, 0x8b]) || picard_analysis::cram_input::is_cram(&magic);
    let (header, records) = match Input::open_with(path, ValidationStringency::Silent) {
        Err(InputError::Sam(e)) => return Err(format!("{:?}", CompareError::Parse(e)).into()),
        opened => opened?.into_parts(),
    };
    if !bam {
        return Ok((header, Box::new(records.map(Ok))));
    }
    let (text, _) = read_sam_with(&header.encode(), ValidationStringency::Silent)
        .map_err(|e| format!("{:?}", CompareError::Parse(e)))?;
    let parsed = text.clone();
    let records = records.map(move |record| {
        let line = picard_analysis::bam_stream::sam_line(&header, &record)
            .ok_or("records failed to re-encode as SAM")?;
        let index_of = |name: &str| {
            parsed
                .sequences
                .iter()
                .position(|s| s.name == name)
                .map(|i| i as i32)
        };
        parse_line_with(&line, index_of, ValidationStringency::Silent)
            .map_err(|e| format!("{:?}", CompareError::Parse(e)))
    });
    Ok((text, Box::new(records)))
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
    let (h1, left) = as_sam(files[0])?;
    let (h2, right) = as_sam(files[1])?;
    let (mut left_error, mut right_error) = (None, None);
    let metric = compare_records(
        &h1,
        left.map_while(|r| r.map_err(|e| left_error = Some(e)).ok()),
        &h2,
        right.map_while(|r| r.map_err(|e| right_error = Some(e)).ok()),
        &absolute(files[0]),
        &absolute(files[1]),
        &options,
    );
    if let Some(e) = left_error.or(right_error) {
        return Err(e.into());
    }
    let metric = metric.map_err(|e| format!("{e:?}"))?;
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
    println!("{}", picard_analysis::compare_sams::verdict(&metric));
    std::process::exit(if metric.are_equal { 0 } else { 1 });
}
