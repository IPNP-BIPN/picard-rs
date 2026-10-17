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
//! `METRIC_ACCUMULATION_LEVEL` adds levels to its `ALL_READS` default: each level's units are the
//! header's read groups, in header order, and a read whose group lacks the field goes to an
//! `unknown` unit made when the first such read arrives.

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

/// `MultiLevelCollector`'s distributors, in output order: all reads, then one collector per
/// sample, library and platform unit the header names, in header order, plus an `unknown` one
/// made the first time a read has no read group or the group lacks the field.
struct Distributors {
    options: Options,
    levels: Vec<(&'static str, Vec<(String, InsertSizeMetricsCollector)>)>,
}

impl Distributors {
    fn new(levels: &[String], header: &SamHeader, options: Options) -> Self {
        let mut out = Vec::new();
        for level in ["ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"] {
            if !levels.iter().any(|l| l == level) {
                continue;
            }
            let mut collectors: Vec<(String, InsertSizeMetricsCollector)> = Vec::new();
            if level == "ALL_READS" {
                collectors.push((String::new(), InsertSizeMetricsCollector::new(options)));
            } else {
                for group in &header.read_groups {
                    let get = |tag: &str| group.attributes.get(tag).map(str::to_string);
                    let (key, unit) = match level {
                        "SAMPLE" => (get("SM"), (get("SM"), None, None)),
                        "LIBRARY" => (get("LB"), (get("SM"), get("LB"), None)),
                        _ => (get("PU"), (get("SM"), get("LB"), get("PU"))),
                    };
                    // A `LinkedHashMap` keyed by the field: a null key is a key too.
                    let key = key.unwrap_or_else(|| "\u{0}null".to_string());
                    if !collectors.iter().any(|(k, _)| *k == key) {
                        collectors.push((
                            key,
                            InsertSizeMetricsCollector::for_unit(options, unit.0, unit.1, unit.2),
                        ));
                    }
                }
            }
            out.push((level, collectors));
        }
        Distributors {
            options,
            levels: out,
        }
    }

    fn accept(&mut self, header: &SamHeader, record: &BamRecord) -> Result<(), String> {
        let group = match record.tags.get(htsjdk_bam::tag::Tag::new(b"RG")) {
            Some(htsjdk_bam::tag::TagValue::Str(id)) => {
                header.read_groups.iter().find(|g| g.id == *id)
            }
            _ => None,
        };
        for (level, collectors) in self.levels.iter_mut() {
            if *level == "ALL_READS" {
                collectors[0].1.accept(record);
                continue;
            }
            let tag = match *level {
                "SAMPLE" => "SM",
                "LIBRARY" => "LB",
                _ => "PU",
            };
            let key = group
                .and_then(|g| g.attributes.get(tag))
                .map_or("unknown".to_string(), str::to_string);
            let position = match collectors.iter().position(|(k, _)| *k == key) {
                Some(position) => position,
                None if key == "unknown" => {
                    let unknown = Some("unknown".to_string());
                    let unit = match *level {
                        "SAMPLE" => (unknown, None, None),
                        "LIBRARY" => (unknown.clone(), unknown, None),
                        _ => (unknown.clone(), unknown.clone(), unknown),
                    };
                    collectors.push((
                        key,
                        InsertSizeMetricsCollector::for_unit(self.options, unit.0, unit.1, unit.2),
                    ));
                    collectors.len() - 1
                }
                None => {
                    return Err(format!(
                        "picard.PicardException: Could not find collector for {key}"
                    ))
                }
            };
            collectors[position].1.accept(record);
        }
        Ok(())
    }

    fn finish(
        self,
    ) -> Vec<(
        picard_analysis::insert_size::InsertSizeMetrics,
        htsjdk_metrics::file::Histogram,
    )> {
        self.levels
            .into_iter()
            .flat_map(|(_, collectors)| collectors.into_iter().flat_map(|(_, c)| c.finish()))
            .collect()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    // Optional since Picard 3.5.0, which refuses it only where there is no R.
    picard_analysis::metrics_cli::refuse_chart_without_r(
        "CollectInsertSizeMetrics",
        arg(&args, "Histogram_FILE=").is_some() || arg(&args, "H=").is_some(),
    );
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
    // A collection argument is appended to its default, `ALL_READS`, and `null` empties it.
    let mut levels: Vec<String> = vec!["ALL_READS".to_string()];
    for level in args.iter().filter_map(|a| {
        a.strip_prefix("METRIC_ACCUMULATION_LEVEL=")
            .or_else(|| a.strip_prefix("LEVEL="))
    }) {
        if level == "null" {
            levels.clear();
        } else if !matches!(level, "ALL_READS" | "SAMPLE" | "LIBRARY" | "READ_GROUP") {
            return Err(format!("unknown METRIC_ACCUMULATION_LEVEL: {level}").into());
        } else {
            levels.push(level.to_string());
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
        let plain = picard_analysis::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
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
    let mut distributors = Distributors::new(&levels, &header, options);
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
        if let Err(thrown) = distributors.accept(&header, record) {
            eprintln!("Exception in thread \"main\" {thrown}");
            std::process::exit(1);
        }
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
        // `usesNoRefReads()` is false: the unmapped reads at the end are not walked.
        if record.reference_index == -1 {
            break;
        }
    }

    let results = distributors.finish();
    if results.is_empty() {
        return Ok(());
    }
    let mut file = MetricsFile::new();
    for header in picard_analysis::command_line::default_headers("CollectInsertSizeMetrics") {
        file.add_header(&header);
    }
    for (metric, histogram) in results {
        file.add_metric(&metric);
        file.histograms.push(histogram);
    }
    std::fs::write(&output, file.write())?;
    // With no `HISTOGRAM_WIDTH` but a `MIN_HISTOGRAM_WIDTH`, the chart's width is the largest bin
    // over every histogram, and `Histogram.getMax` on one with no bins left throws.
    let histogram_width = optional_int("HISTOGRAM_WIDTH=").or_else(|| optional_int("W="));
    let min_histogram_width = optional_int("MIN_HISTOGRAM_WIDTH=").or_else(|| optional_int("MW="));
    if histogram_width.is_none()
        && min_histogram_width.is_some()
        && file.histograms.iter().any(|histogram| histogram.is_empty())
    {
        eprintln!(
            "Exception in thread \"main\" java.lang.NullPointerException: Cannot invoke \"java.util.Map$Entry.getValue()\" because the return value of \"java.util.NavigableMap.lastEntry()\" is null"
        );
        std::process::exit(1);
    }
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
