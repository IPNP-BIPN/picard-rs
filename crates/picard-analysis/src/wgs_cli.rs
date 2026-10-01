//! The command line shared by `CollectWgsMetrics`, `CollectRawWgsMetrics` and
//! `CollectWgsMetricsWithNonZeroCoverage`: the three are one `doWork` with different defaults, a
//! different metric class, and (for the last) a second row and a chart.
//!
//! What runs is [`crate::wgs_walk`]. This is the order `doWork` does things in, because the order
//! decides which refusal a row that breaks two rules gets: the non-zero tool builds its collector
//! -- and so refuses a non-positive `COVERAGE_CAP` -- before it opens the locus iterator that
//! refuses a queryname-sorted input, and the other two do it the other way round.
//!
//! Ported from Picard 3.4.0 `CollectWgsMetrics.doWork`, `CollectRawWgsMetrics` and
//! `CollectWgsMetricsWithNonZeroCoverage.doWork`.

use std::io::Read;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::interval::IntervalList;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_metrics::format::format_double;

use crate::wgs_walk::{
    uniqued, walk, wgs_metrics, Source, WalkOptions, WgsMetricsRow, WGS_COLUMNS,
};

/// Which of the three tools is running.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    Wgs,
    Raw,
    NonZero,
}

impl Tool {
    fn name(self) -> &'static str {
        match self {
            Tool::Wgs => "CollectWgsMetrics",
            Tool::Raw => "CollectRawWgsMetrics",
            Tool::NonZero => "CollectWgsMetricsWithNonZeroCoverage",
        }
    }

    fn class_name(self) -> &'static str {
        match self {
            Tool::Wgs => "picard.analysis.WgsMetrics",
            Tool::Raw => "picard.analysis.CollectRawWgsMetrics$RawWgsMetrics",
            Tool::NonZero => {
                "picard.analysis.CollectWgsMetricsWithNonZeroCoverage$WgsMetricsWithNonZeroCoverage"
            }
        }
    }
}

struct Row<'a> {
    class_name: &'static str,
    columns: Vec<&'static str>,
    category: Option<&'a str>,
    metrics: &'a WgsMetricsRow,
}

impl MetricBean for Row<'_> {
    fn class_name(&self) -> &str {
        self.class_name
    }
    fn columns(&self) -> &[&'static str] {
        &self.columns
    }
    fn values(&self) -> Vec<Value> {
        let mut values = Vec::with_capacity(33);
        if let Some(category) = self.category {
            values.push(Value::Str(category.to_string()));
        }
        values.push(Value::Long(self.metrics.genome_territory));
        values.extend(self.metrics.doubles.iter().map(|&d| Value::Double(d)));
        values
    }
}

/// Every `NAME=value` the command line gave for one argument, under either of its names.
fn values_of<'a>(args: &'a [String], long: &str, short: Option<&str>) -> Vec<&'a str> {
    let mut out = Vec::new();
    for a in args {
        if let Some((name, value)) = a.split_once('=') {
            if name == long || Some(name) == short {
                out.push(value);
            }
        }
    }
    out
}

fn one<'a>(args: &'a [String], long: &str, short: Option<&str>) -> Option<&'a str> {
    values_of(args, long, short).last().copied()
}

fn fail(message: &str) -> ! {
    eprintln!("{message}");
    std::process::exit(1);
}

fn thrown(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" {message}");
    std::process::exit(1);
}

fn int_arg(args: &[String], long: &str, short: Option<&str>, default: i64) -> i64 {
    match one(args, long, short) {
        None => default,
        Some(v) => v.parse().unwrap_or_else(|_| {
            fail(&format!(
                "Argument '{long}' cannot be set to '{v}': it is not a number"
            ))
        }),
    }
}

fn bool_arg(args: &[String], long: &str, default: bool) -> bool {
    match one(args, long, None) {
        None => default,
        Some(v) if v.eq_ignore_ascii_case("true") => true,
        Some(v) if v.eq_ignore_ascii_case("false") => false,
        Some(v) => fail(&format!("Argument '{long}' cannot be set to '{v}'")),
    }
}

fn optional_path(args: &[String], long: &str, short: Option<&str>) -> Option<String> {
    one(args, long, short)
        .filter(|v| !v.eq_ignore_ascii_case("null"))
        .map(str::to_string)
}

/// `SamFiles.findIndex`: `<file>.bai`, then the `.bam` extension replaced by `.bai`.
fn has_index(input: &str) -> bool {
    if std::path::Path::new(&format!("{input}.bai")).exists() {
        return true;
    }
    if let Some(stem) = input.strip_suffix(".bam") {
        return std::path::Path::new(&format!("{stem}.bai")).exists();
    }
    false
}

fn read_input(input: &str) -> Result<(SamHeader, Vec<BamRecord>, bool), String> {
    let mut raw = Vec::new();
    std::fs::File::open(input)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(|e| format!("{e}"))?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((header, records, has_index(input)))
    } else {
        let text = String::from_utf8(raw).map_err(|e| format!("{e}"))?;
        let (header, records) = read_sam(&text).map_err(|e| format!("{e:?}"))?;
        Ok((header, records, false))
    }
}

/// The histogram section, written directly: the shared writer looks up each key by scanning,
/// which is quadratic, and `CollectRawWgsMetrics` writes a hundred thousand and one bins.
fn histogram_section(columns: &[(&str, Vec<(usize, i64)>)]) -> String {
    let mut out = String::new();
    out.push_str("## HISTOGRAM\tjava.lang.Integer\n");
    out.push_str("coverage");
    for (label, _) in columns {
        out.push('\t');
        out.push_str(label);
    }
    out.push('\n');
    let max = columns
        .iter()
        .flat_map(|(_, bins)| bins.iter().map(|(k, _)| *k))
        .max()
        .unwrap_or(0);
    let mut table: Vec<Vec<Option<i64>>> = vec![vec![None; max + 1]; columns.len()];
    for (c, (_, bins)) in columns.iter().enumerate() {
        for &(k, v) in bins {
            table[c][k] = Some(v);
        }
    }
    for key in 0..=max {
        if table.iter().all(|column| column[key].is_none()) {
            continue;
        }
        out.push_str(&key.to_string());
        for column in &table {
            out.push('\t');
            out.push_str(&format_double(column[key].unwrap_or(0) as f64));
        }
        out.push('\n');
    }
    out
}

fn bins(array: &[i64], from: usize) -> Vec<(usize, i64)> {
    array
        .iter()
        .enumerate()
        .skip(from)
        .map(|(i, &v)| (i, v))
        .collect()
}

/// Runs one of the three tools on `std::env::args()`, exiting the way Picard would.
pub fn main_for(tool: Tool) {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let input = one(&args, "INPUT", Some("I"))
        .map(str::to_string)
        .unwrap_or_else(|| fail("Argument 'INPUT' is required"));
    let output = one(&args, "OUTPUT", Some("O"))
        .map(str::to_string)
        .unwrap_or_else(|| fail("Argument 'OUTPUT' is required"));
    let reference = one(&args, "REFERENCE_SEQUENCE", Some("R"))
        .map(str::to_string)
        .unwrap_or_else(|| fail("Argument 'REFERENCE_SEQUENCE' is required"));
    if tool == Tool::NonZero && one(&args, "CHART_OUTPUT", Some("CHART")).is_none() {
        fail("Argument CHART_OUTPUT was missing: Argument 'CHART_OUTPUT' is required");
    }
    if optional_path(&args, "THEORETICAL_SENSITIVITY_OUTPUT", None).is_some() {
        fail("THEORETICAL_SENSITIVITY_OUTPUT is not ported");
    }
    if let Some(stringency) = one(&args, "VALIDATION_STRINGENCY", None) {
        if !matches!(stringency, "STRICT" | "LENIENT" | "SILENT") {
            fail(&format!("unknown VALIDATION_STRINGENCY: {stringency}"));
        }
    }
    let raw = tool == Tool::Raw;
    let minimum_mapping_quality = int_arg(
        &args,
        "MINIMUM_MAPPING_QUALITY",
        Some("MQ"),
        if raw { 0 } else { 20 },
    ) as i32;
    let minimum_base_quality = int_arg(
        &args,
        "MINIMUM_BASE_QUALITY",
        Some("Q"),
        if raw { 3 } else { 20 },
    ) as i32;
    let coverage_cap = int_arg(
        &args,
        "COVERAGE_CAP",
        Some("CAP"),
        if raw { 100_000 } else { 250 },
    ) as i32;
    let mut locus_accumulation_cap = int_arg(
        &args,
        "LOCUS_ACCUMULATION_CAP",
        None,
        if raw { 200_000 } else { 100_000 },
    ) as i32;
    let stop_after = int_arg(&args, "STOP_AFTER", None, -1);
    let include_bq_histogram = bool_arg(&args, "INCLUDE_BQ_HISTOGRAM", false);
    let count_unpaired = bool_arg(&args, "COUNT_UNPAIRED", false);
    let sample_size = int_arg(&args, "SAMPLE_SIZE", None, 10_000) as i32;
    let use_fast_algorithm = bool_arg(&args, "USE_FAST_ALGORITHM", false);
    let read_length = int_arg(&args, "READ_LENGTH", None, 150) as i32;
    let intervals_path = optional_path(&args, "INTERVALS", None);

    let (header, records, indexed) = read_input(&input).unwrap_or_else(|e| fail(&e));

    // CollectWgsMetricsWithNonZeroCoverage builds its collector up front, so the collector's own
    // refusal comes before anything the iterator would say.
    if tool == Tool::NonZero && coverage_cap <= 0 {
        thrown("java.lang.IllegalArgumentException: Coverage cap must be positive.");
    }
    if locus_accumulation_cap < coverage_cap {
        locus_accumulation_cap = coverage_cap;
    }

    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| fail(&format!("{e:?}")));

    // AbstractLocusIterator's constructor: an unsorted header is a warning, any order but
    // coordinate a refusal.
    match header.attributes.get("SO") {
        None | Some("unsorted") | Some("coordinate") => {}
        Some(_) => thrown(&format!(
            "htsjdk.samtools.SAMException: {} cannot operate on a SAM file that is not coordinate sorted.",
            if use_fast_algorithm {
                "EdgeReadIterator"
            } else {
                "SamLocusIterator"
            }
        )),
    }
    let names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let lengths: Vec<i32> = header.sequences.iter().map(|s| s.length).collect();
    let intervals = intervals_path.as_ref().map(|path| {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
        let list = IntervalList::parse_body(names.clone(), &text)
            .unwrap_or_else(|e| fail(&format!("{e:?}")));
        let triples = list
            .intervals
            .iter()
            .map(|iv| {
                let seq = names
                    .iter()
                    .position(|n| *n == iv.contig)
                    .unwrap_or_else(|| fail(&format!("unknown contig {}", iv.contig)));
                (seq as i32, iv.start, iv.end)
            })
            .collect::<Vec<_>>();
        uniqued(triples)
    });

    // SequenceDictionaryUtils.assertSequenceDictionariesEqual(input, reference).
    if !header.sequences.is_empty() {
        let same = contigs.len() == header.sequences.len()
            && contigs
                .iter()
                .zip(&header.sequences)
                .all(|(c, s)| c.name == s.name && c.bases.len() as i32 == s.length);
        if !same {
            fail("the input's sequence dictionary differs from the reference's");
        }
    }
    if coverage_cap <= 0 {
        thrown("java.lang.IllegalArgumentException: Coverage cap must be positive.");
    }

    let references: Vec<Vec<u8>> = contigs.into_iter().map(|c| c.bases).collect();
    let options = WalkOptions {
        minimum_mapping_quality,
        minimum_base_quality,
        coverage_cap,
        locus_accumulation_cap,
        stop_after,
        count_unpaired,
        use_fast_algorithm,
        fast_collector: use_fast_algorithm && tool != Tool::NonZero,
        read_length,
    };
    let source = match &intervals {
        None => Source::WholeFile,
        Some(list) => Source::Intervals {
            intervals: list,
            indexed,
        },
    };
    let result =
        walk(&records, &references, &lengths, source, &options).unwrap_or_else(|e| thrown(&e));

    let class_name = tool.class_name();
    let mut columns: Vec<&'static str> = Vec::new();
    if tool == Tool::NonZero {
        columns.push("CATEGORY");
    }
    columns.extend(WGS_COLUMNS);

    let mut file = MetricsFile::new();
    file.add_header(&format!("{} <command line>", tool.name()));
    file.add_header("Started on: <timestamp>");
    let mut histograms: Vec<(&str, Vec<(usize, i64)>)> = Vec::new();
    if tool == Tool::NonZero {
        let whole = wgs_metrics(
            &result.high_quality_depth,
            &result.unfiltered_depth,
            &result.unfiltered_baseq,
            &result,
            coverage_cap,
            sample_size,
        )
        .unwrap_or_else(|e| thrown(&e));
        let mut high_quality = result.high_quality_depth.clone();
        let mut unfiltered = result.unfiltered_depth.clone();
        high_quality[0] = 0;
        unfiltered[0] = 0;
        let non_zero = wgs_metrics(
            &high_quality,
            &unfiltered,
            &result.unfiltered_baseq,
            &result,
            coverage_cap,
            sample_size,
        )
        .unwrap_or_else(|e| thrown(&e));
        for (category, metrics) in [("WHOLE_GENOME", &whole), ("NON_ZERO_REGIONS", &non_zero)] {
            file.add_metric(&Row {
                class_name,
                columns: columns.clone(),
                category: Some(category),
                metrics,
            });
        }
        histograms.push(("count_WHOLE_GENOME", bins(&result.high_quality_depth, 0)));
        histograms.push((
            "count_NON_ZERO_REGIONS",
            bins(&result.high_quality_depth, 1),
        ));
    } else {
        let metrics = wgs_metrics(
            &result.high_quality_depth,
            &result.unfiltered_depth,
            &result.unfiltered_baseq,
            &result,
            coverage_cap,
            sample_size,
        )
        .unwrap_or_else(|e| thrown(&e));
        file.add_metric(&Row {
            class_name,
            columns,
            category: None,
            metrics: &metrics,
        });
        histograms.push((
            "high_quality_coverage_count",
            bins(&result.high_quality_depth, 0),
        ));
    }
    if include_bq_histogram {
        histograms.push(("unfiltered_baseq_count", bins(&result.unfiltered_baseq, 0)));
    }

    // The writer ends a file with no histogram in two newlines; the section goes between them.
    let mut text = file.write();
    text.pop();
    text.push_str(&histogram_section(&histograms));
    text.push('\n');
    if let Err(e) = std::fs::write(&output, text) {
        fail(&format!("{e}"));
    }
    // The non-zero tool then draws its chart through R. The chart is a rendering of the file
    // just written and is not ported; the metrics are the tool's answer.
}
