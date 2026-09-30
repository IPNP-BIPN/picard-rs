//! `CollectInsertSizeMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectInsertSizeMetrics` at tag 3.4.0. The collector is
//! `picard_analysis::insert_size`; this is `SinglePassSamProgram.makeItSo` around it:
//!
//! * the sort-order refusal, which `ASSUME_SORTED` (true by default) turns into a warning, and the
//!   reference walker that refuses to rewind whenever `REFERENCE_SEQUENCE` is given;
//! * the walk stops at the first unmapped record, because this program does not use them, or after
//!   `STOP_AFTER` records;
//! * the metrics file is written only when some orientation cleared `MINIMUM_PCT`; otherwise the
//!   tool warns and writes nothing, and the chart is not drawn.
//!
//! The chart itself is R's, and not ported: the port writes the metrics and no PDF. What it does
//! reproduce is R's failure when every histogram was trimmed away, which ends the run at one after
//! the metrics are written; the reference's log lines carry a timestamp the port's do not.
//! `METRIC_ACCUMULATION_LEVEL` is ported at `ALL_READS`, its default, alone.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::insert_size::{InsertSizeMetricsCollector, Options};
use picard_analysis::single_pass_rejections::{check_sort_order, walk_reference, SortOrder};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    if arg(&args, "Histogram_FILE=").is_none() && arg(&args, "H=").is_none() {
        eprintln!("Argument Histogram_FILE was missing: Argument 'Histogram_FILE' is required");
        std::process::exit(1);
    }
    let flag = |key: &str, default: bool| {
        arg(&args, key)
            .map(|v| v.eq_ignore_ascii_case("true"))
            .unwrap_or(default)
    };
    let optional_int = |key: &str| {
        arg(&args, key)
            .filter(|v| v != "null")
            .and_then(|v| v.parse::<i32>().ok())
    };
    let minimum_pct: f32 = arg(&args, "MINIMUM_PCT=")
        .or_else(|| arg(&args, "M="))
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.05);
    if !(0.0..=0.5).contains(&minimum_pct) {
        eprintln!(
            "MINIMUM_PCT was set to {minimum_pct:?}. It must be between 0 and 0.5 so all data categories don't get discarded."
        );
        std::process::exit(1);
    }
    for level in args.iter().filter_map(|a| {
        a.strip_prefix("METRIC_ACCUMULATION_LEVEL=")
            .or_else(|| a.strip_prefix("LEVEL="))
    }) {
        if level != "ALL_READS" {
            return Err(format!(
                "METRIC_ACCUMULATION_LEVEL={level} is not ported: only ALL_READS is"
            )
            .into());
        }
    }
    let options = Options {
        minimum_pct,
        deviations: arg(&args, "DEVIATIONS=")
            .and_then(|v| v.parse().ok())
            .unwrap_or(10.0),
        histogram_width: optional_int("HISTOGRAM_WIDTH=").or_else(|| optional_int("W=")),
        min_histogram_width: optional_int("MIN_HISTOGRAM_WIDTH=").or_else(|| optional_int("MW=")),
        include_duplicates: flag("INCLUDE_DUPLICATES=", false),
    };
    let stop_after: i64 = arg(&args, "STOP_AFTER=")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);

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
        read_sam(&String::from_utf8(raw)?).map_err(|e| format!("{e:?}"))?
    };

    let found = match header.attributes.get("SO") {
        Some("coordinate") => SortOrder::Coordinate,
        Some("queryname") => SortOrder::Queryname,
        Some("duplicate") => SortOrder::Duplicate,
        Some("unsorted") => SortOrder::Unsorted,
        _ => SortOrder::Unknown,
    };
    if let Err(rejection) = check_sort_order(&input, found, flag("ASSUME_SORTED=", true)) {
        eprintln!("Exception in thread \"main\" {}", rejection.thrown());
        std::process::exit(1);
    }
    let with_reference = arg(&args, "REFERENCE_SEQUENCE=").is_some() || arg(&args, "R=").is_some();
    let mut collector = InsertSizeMetricsCollector::new(options);
    let mut current: Option<i32> = None;
    let mut count = 0i64;
    for record in &records {
        if with_reference && record.reference_index != -1 {
            match walk_reference(current, record.reference_index) {
                Ok(next) => current = Some(next),
                Err(rejection) => {
                    eprintln!("Exception in thread \"main\" {}", rejection.thrown());
                    std::process::exit(1);
                }
            }
        }
        collector.accept(record);
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
        // `usesNoRefReads()` is false: the unmapped reads at the end are not walked.
        if record.reference_index == -1 {
            break;
        }
    }

    let results = collector.finish();
    if results.is_empty() {
        return Ok(());
    }
    let mut file = MetricsFile::new();
    file.add_header("CollectInsertSizeMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    for (metric, histogram) in results {
        file.add_metric(&metric);
        file.histograms.push(histogram);
    }
    std::fs::write(&output, file.write())?;
    // The chart reads the histogram section back, and a width that trimmed every histogram to
    // nothing leaves the file without one: R stops, and the tool fails after writing its metrics.
    if file.histograms.iter().all(|histogram| histogram.is_empty()) {
        eprintln!("ERROR\tProcessExecutor\tError in read.table(metricsFile, header = TRUE, sep = \"\\t\", skip = secondBlankLine,  : ");
        eprintln!("ERROR\tProcessExecutor\t  no lines available in input");
        eprintln!("ERROR\tProcessExecutor\tExecution halted");
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: R script picard/analysis/insertSizeHistogram.R failed with return code 1"
        );
        std::process::exit(1);
    }
    Ok(())
}
