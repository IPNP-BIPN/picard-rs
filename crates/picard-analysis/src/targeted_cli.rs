//! The command line shared by `CollectHsMetrics` and `CollectTargetedPcrMetrics`, which are one
//! `CollectTargetedMetrics.doWork` with different probe arguments, defaults and metric classes.
//!
//! Ported from Picard 3.4.0 `CollectTargetedMetrics`, `CollectHsMetrics` and
//! `CollectTargetedPcrMetrics`. What runs is [`crate::targeted_metrics`].

use std::io::Read;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::interval::IntervalList;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_metrics::format::format_double;

use crate::targeted_metrics::{
    calculate_gc, collect, hs_penalty, library_size, Level, Options, Panel, ReadGroup, Span,
    TargetMetrics,
};
use crate::wgs_walk::uniqued;

/// Which of the two tools is running.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    Hs,
    TargetedPcr,
}

const PANEL_COLUMNS: [&str; 49] = [
    "TARGET_TERRITORY",
    "GENOME_SIZE",
    "TOTAL_READS",
    "PF_READS",
    "PF_BASES",
    "PF_UNIQUE_READS",
    "PF_UQ_READS_ALIGNED",
    "PF_BASES_ALIGNED",
    "PF_UQ_BASES_ALIGNED",
    "ON_TARGET_BASES",
    "PCT_PF_READS",
    "PCT_PF_UQ_READS",
    "PCT_PF_UQ_READS_ALIGNED",
    "MEAN_TARGET_COVERAGE",
    "MEDIAN_TARGET_COVERAGE",
    "MAX_TARGET_COVERAGE",
    "MIN_TARGET_COVERAGE",
    "ZERO_CVG_TARGETS_PCT",
    "PCT_EXC_DUPE",
    "PCT_EXC_ADAPTER",
    "PCT_EXC_MAPQ",
    "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP",
    "PCT_EXC_OFF_TARGET",
    "FOLD_80_BASE_PENALTY",
    "PCT_TARGET_BASES_1X",
    "PCT_TARGET_BASES_2X",
    "PCT_TARGET_BASES_10X",
    "PCT_TARGET_BASES_20X",
    "PCT_TARGET_BASES_30X",
    "PCT_TARGET_BASES_40X",
    "PCT_TARGET_BASES_50X",
    "PCT_TARGET_BASES_100X",
    "PCT_TARGET_BASES_250X",
    "PCT_TARGET_BASES_500X",
    "PCT_TARGET_BASES_1000X",
    "PCT_TARGET_BASES_2500X",
    "PCT_TARGET_BASES_5000X",
    "PCT_TARGET_BASES_10000X",
    "PCT_TARGET_BASES_25000X",
    "PCT_TARGET_BASES_50000X",
    "PCT_TARGET_BASES_100000X",
    "AT_DROPOUT",
    "GC_DROPOUT",
    "HET_SNP_SENSITIVITY",
    "HET_SNP_Q",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

const HS_COLUMNS: [&str; 20] = [
    "BAIT_SET",
    "BAIT_TERRITORY",
    "BAIT_DESIGN_EFFICIENCY",
    "ON_BAIT_BASES",
    "NEAR_BAIT_BASES",
    "OFF_BAIT_BASES",
    "PCT_SELECTED_BASES",
    "PCT_OFF_BAIT",
    "ON_BAIT_VS_SELECTED",
    "MEAN_BAIT_COVERAGE",
    "PCT_USABLE_BASES_ON_BAIT",
    "PCT_USABLE_BASES_ON_TARGET",
    "FOLD_ENRICHMENT",
    "HS_LIBRARY_SIZE",
    "HS_PENALTY_10X",
    "HS_PENALTY_20X",
    "HS_PENALTY_30X",
    "HS_PENALTY_40X",
    "HS_PENALTY_50X",
    "HS_PENALTY_100X",
];

const PCR_COLUMNS: [&str; 13] = [
    "CUSTOM_AMPLICON_SET",
    "AMPLICON_TERRITORY",
    "ON_AMPLICON_BASES",
    "NEAR_AMPLICON_BASES",
    "OFF_AMPLICON_BASES",
    "PCT_AMPLIFIED_BASES",
    "PCT_OFF_AMPLICON",
    "ON_AMPLICON_VS_SELECTED",
    "MEAN_AMPLICON_COVERAGE",
    "FOLD_ENRICHMENT",
    "PF_SELECTED_PAIRS",
    "PF_SELECTED_UNIQUE_PAIRS",
    "ON_TARGET_FROM_PAIR_BASES",
];

fn text(value: &Option<String>) -> Value {
    match value {
        Some(s) => Value::Str(s.clone()),
        None => Value::Null,
    }
}

fn panel_values(m: &TargetMetrics) -> Vec<Value> {
    let mut v = vec![
        Value::Long(m.target_territory),
        Value::Long(m.genome_size),
        Value::Long(m.total_reads),
        Value::Long(m.pf_reads),
        Value::Long(m.pf_bases),
        Value::Long(m.pf_unique_reads),
        Value::Long(m.pf_uq_reads_aligned),
        Value::Long(m.pf_bases_aligned),
        Value::Long(m.pf_uq_bases_aligned),
        Value::Long(m.on_target_bases),
        Value::Double(m.pct_pf_reads),
        Value::Double(m.pct_pf_uq_reads),
        Value::Double(m.pct_pf_uq_reads_aligned),
        Value::Double(m.mean_target_coverage),
        Value::Double(m.median_target_coverage),
        Value::Long(m.max_target_coverage),
        Value::Long(m.min_target_coverage),
        Value::Double(m.zero_cvg_targets_pct),
        Value::Double(m.pct_exc_dupe),
        Value::Double(m.pct_exc_adapter),
        Value::Double(m.pct_exc_mapq),
        Value::Double(m.pct_exc_baseq),
        Value::Double(m.pct_exc_overlap),
        Value::Double(m.pct_exc_off_target),
        Value::Double(m.fold_80_base_penalty),
    ];
    v.extend(m.pct_target_bases.iter().map(|&d| Value::Double(d)));
    v.extend([
        Value::Double(m.at_dropout),
        Value::Double(m.gc_dropout),
        Value::Double(m.het_snp_sensitivity),
        Value::Double(m.het_snp_q),
        text(&m.sample),
        text(&m.library),
        text(&m.read_group),
    ]);
    v
}

struct Row {
    tool: Tool,
    columns: Vec<&'static str>,
    values: Vec<Value>,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        match self.tool {
            Tool::Hs => "picard.analysis.directed.HsMetrics",
            Tool::TargetedPcr => "picard.analysis.directed.TargetedPcrMetrics",
        }
    }
    fn columns(&self) -> &[&'static str] {
        &self.columns
    }
    fn values(&self) -> Vec<Value> {
        self.values.clone()
    }
}

/// `HsMetricCollector.convertMetric` or `TargetedPcrMetricsCollector.convertMetric`.
fn convert(tool: Tool, m: &TargetMetrics) -> Row {
    let mut columns: Vec<&'static str> = Vec::new();
    let mut values: Vec<Value> = Vec::new();
    match tool {
        Tool::Hs => {
            let size = library_size(m);
            columns.extend(HS_COLUMNS);
            values.extend([
                text(&m.probe_set),
                Value::Long(m.probe_territory),
                Value::Double(m.target_territory as f64 / m.probe_territory as f64),
                Value::Long(m.on_probe_bases),
                Value::Long(m.near_probe_bases),
                Value::Long(m.off_probe_bases),
                Value::Double(m.pct_selected_bases),
                Value::Double(m.pct_off_probe),
                Value::Double(m.on_probe_vs_selected),
                Value::Double(m.mean_probe_coverage),
                Value::Double(m.on_probe_bases as f64 / m.pf_bases as f64),
                Value::Double(m.on_target_bases as f64 / m.pf_bases as f64),
                Value::Double(m.fold_enrichment),
                match size {
                    Some(s) => Value::Long(s),
                    None => Value::Null,
                },
            ]);
            for goal in [10, 20, 30, 40, 50, 100] {
                values.push(Value::Double(hs_penalty(size, m, goal)));
            }
        }
        Tool::TargetedPcr => {
            columns.extend(PCR_COLUMNS);
            values.extend([
                text(&m.probe_set),
                Value::Long(m.probe_territory),
                Value::Long(m.on_probe_bases),
                Value::Long(m.near_probe_bases),
                Value::Long(m.off_probe_bases),
                Value::Double(m.pct_selected_bases),
                Value::Double(m.pct_off_probe),
                Value::Double(m.on_probe_vs_selected),
                Value::Double(m.mean_probe_coverage),
                Value::Double(m.fold_enrichment),
                Value::Long(m.pf_selected_pairs),
                Value::Long(m.pf_selected_unique_pairs),
                Value::Long(m.on_target_from_pair_bases),
            ]);
        }
    }
    columns.extend(PANEL_COLUMNS);
    values.extend(panel_values(m));
    Row {
        tool,
        columns,
        values,
    }
}

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

fn int_arg(args: &[String], long: &str, short: Option<&str>, default: i32) -> i32 {
    match one(args, long, short) {
        None => default,
        Some(v) => v
            .parse()
            .unwrap_or_else(|_| fail(&format!("Argument '{long}' cannot be set to '{v}'"))),
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

fn not_null(v: &str) -> bool {
    !v.eq_ignore_ascii_case("null")
}

/// `CollectTargetedMetrics.renderProbeNameFromFile`: the file name up to its first period.
fn probe_name(path: &str) -> String {
    let name = std::path::Path::new(path)
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_default();
    match name.find('.') {
        Some(i) => name[..i].to_string(),
        None => name,
    }
}

fn read_input(input: &str) -> Result<(SamHeader, Vec<BamRecord>), String> {
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
        Ok((header, records))
    } else {
        let text = String::from_utf8(raw).map_err(|e| format!("{e}"))?;
        read_sam(&text).map_err(|e| format!("{e:?}"))
    }
}

/// `IntervalList.fromFiles`: every file's intervals, in file order, not uniqued.
fn read_intervals(paths: &[&str], names: &[String]) -> Vec<Span> {
    let mut out = Vec::new();
    for path in paths {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
        let list = IntervalList::parse_body(names.to_vec(), &text)
            .unwrap_or_else(|e| fail(&format!("{e:?}")));
        for iv in &list.intervals {
            let seq = names
                .iter()
                .position(|n| *n == iv.contig)
                .unwrap_or_else(|| fail(&format!("unknown contig {}", iv.contig)));
            out.push((seq as i32, iv.start, iv.end));
        }
    }
    out
}

fn histogram_section(columns: &[(&str, &[(i64, f64)])]) -> String {
    let mut out = String::from("## HISTOGRAM\tjava.lang.Integer\ncoverage_or_base_quality");
    for (label, _) in columns {
        out.push('\t');
        out.push_str(label);
    }
    out.push('\n');
    let mut keys: Vec<i64> = columns
        .iter()
        .flat_map(|(_, bins)| bins.iter().map(|(k, _)| *k))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let maps: Vec<std::collections::HashMap<i64, f64>> = columns
        .iter()
        .map(|(_, bins)| bins.iter().copied().collect())
        .collect();
    for key in keys {
        out.push_str(&key.to_string());
        for map in &maps {
            out.push('\t');
            out.push_str(&format_double(*map.get(&key).unwrap_or(&0.0)));
        }
        out.push('\n');
    }
    out
}

/// Runs one of the two tools on `std::env::args()`, exiting the way Picard would.
pub fn main_for(tool: Tool) {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let input = one(&args, "INPUT", Some("I"))
        .map(str::to_string)
        .unwrap_or_else(|| fail("Argument 'INPUT' is required"));
    let output = one(&args, "OUTPUT", Some("O"))
        .map(str::to_string)
        .unwrap_or_else(|| fail("Argument 'OUTPUT' is required"));
    let targets_paths: Vec<&str> = values_of(&args, "TARGET_INTERVALS", Some("TI"))
        .into_iter()
        .filter(|v| not_null(v))
        .collect();
    if targets_paths.is_empty() {
        fail("Argument 'TARGET_INTERVALS' is required");
    }
    let (probe_paths, probe_set_name): (Vec<&str>, String) = match tool {
        Tool::Hs => {
            let paths: Vec<&str> = values_of(&args, "BAIT_INTERVALS", Some("BI"))
                .into_iter()
                .filter(|v| not_null(v))
                .collect();
            if paths.is_empty() {
                fail("Argument 'BAIT_INTERVALS' is required");
            }
            let name = match one(&args, "BAIT_SET_NAME", Some("N")).filter(|v| not_null(v)) {
                Some(n) => n.to_string(),
                None => {
                    let mut names: Vec<String> = paths.iter().map(|p| probe_name(p)).collect();
                    names.sort();
                    names.dedup();
                    names.join(".")
                }
            };
            (paths, name)
        }
        Tool::TargetedPcr => {
            let path = one(&args, "AMPLICON_INTERVALS", Some("AI"))
                .filter(|v| not_null(v))
                .unwrap_or_else(|| fail("Argument 'AMPLICON_INTERVALS' is required"));
            let name =
                match one(&args, "CUSTOM_AMPLICON_SET_NAME", Some("N")).filter(|v| not_null(v)) {
                    Some(n) => n.to_string(),
                    None => probe_name(path),
                };
            (vec![path], name)
        }
    };
    for side in [
        "PER_TARGET_COVERAGE",
        "PER_BASE_COVERAGE",
        "THEORETICAL_SENSITIVITY_OUTPUT",
    ] {
        if one(&args, side, None).filter(|v| not_null(v)).is_some() {
            fail(&format!("{side} is not ported"));
        }
    }
    if let Some(stringency) = one(&args, "VALIDATION_STRINGENCY", None) {
        if !matches!(stringency, "STRICT" | "LENIENT" | "SILENT") {
            fail(&format!("unknown VALIDATION_STRINGENCY: {stringency}"));
        }
    }
    // A collection argument is appended to its default, so naming a level adds it to ALL_READS;
    // `null` clears the collection first.
    let mut levels = vec![Level::AllReads];
    for value in values_of(&args, "METRIC_ACCUMULATION_LEVEL", Some("LEVEL")) {
        let level = match value {
            "null" => {
                levels.clear();
                continue;
            }
            "ALL_READS" => Level::AllReads,
            "SAMPLE" => Level::Sample,
            "LIBRARY" => Level::Library,
            "READ_GROUP" => Level::ReadGroup,
            other => fail(&format!("unknown METRIC_ACCUMULATION_LEVEL: {other}")),
        };
        if !levels.contains(&level) {
            levels.push(level);
        }
    }
    let hs = tool == Tool::Hs;
    let options = Options {
        levels,
        near_distance: int_arg(&args, "NEAR_DISTANCE", None, 250),
        minimum_mapping_quality: int_arg(
            &args,
            "MINIMUM_MAPPING_QUALITY",
            Some("MQ"),
            if hs { 20 } else { 1 },
        ),
        minimum_base_quality: int_arg(
            &args,
            "MINIMUM_BASE_QUALITY",
            Some("Q"),
            if hs { 20 } else { 0 },
        ),
        clip_overlapping_reads: bool_arg(&args, "CLIP_OVERLAPPING_READS", hs),
        include_indels: bool_arg(&args, "INCLUDE_INDELS", false),
        coverage_cap: int_arg(&args, "COVERAGE_CAP", Some("covMax"), 200),
        sample_size: int_arg(&args, "SAMPLE_SIZE", None, 10_000),
        probe_set_name,
    };
    let reference = one(&args, "REFERENCE_SEQUENCE", Some("R"))
        .filter(|v| not_null(v))
        .map(str::to_string);

    let (header, records) = read_input(&input).unwrap_or_else(|e| fail(&e));
    let names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let raw_targets = read_intervals(&targets_paths, &names);
    let probes = uniqued(read_intervals(&probe_paths, &names));
    let targets = uniqued(raw_targets.clone());
    let genome_size: i64 = header.sequences.iter().map(|s| s.length as i64).sum();
    let target_gc = reference.as_ref().map(|path| {
        let contigs = read_fasta_file(path).unwrap_or_else(|e| fail(&format!("{e:?}")));
        targets
            .iter()
            .map(|&(seq, start, end)| {
                let bases = &contigs[seq as usize].bases;
                calculate_gc(&bases[(start - 1) as usize..end as usize])
            })
            .collect()
    });
    let panel = Panel {
        targets,
        probes,
        raw_target_count: raw_targets.len(),
        genome_size,
        target_gc,
    };
    let groups: Vec<ReadGroup> = header
        .read_groups
        .iter()
        .map(|g| ReadGroup {
            id: g.id.clone(),
            sample: g.attributes.get("SM").map(str::to_string),
            library: g.attributes.get("LB").map(str::to_string),
            platform_unit: g.attributes.get("PU").map(str::to_string),
        })
        .collect();
    let group_of = |record: &BamRecord| -> Option<usize> {
        let id = match record.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(s)) => s.to_string(),
            _ => return None,
        };
        groups.iter().position(|g| g.id == id)
    };
    let outcome =
        collect(&records, group_of, &groups, &panel, &options).unwrap_or_else(|e| thrown(&e));

    let mut file = MetricsFile::new();
    file.add_header(&format!(
        "{} <command line>",
        match tool {
            Tool::Hs => "CollectHsMetrics",
            Tool::TargetedPcr => "CollectTargetedPcrMetrics",
        }
    ));
    file.add_header("Started on: <timestamp>");
    let mut histograms: Vec<(&str, &[(i64, f64)])> = Vec::new();
    for row in &outcome.rows {
        file.add_metric(&convert(tool, row));
        histograms.push(("high_quality_coverage_count", &outcome.high_quality_depth));
        histograms.push(("unfiltered_baseq_count", &outcome.uncapped_base_q));
    }
    let mut text = file.write();
    text.pop();
    if !histograms.is_empty() {
        text.push_str(&histogram_section(&histograms));
    }
    text.push('\n');
    if let Err(e) = std::fs::write(&output, text) {
        fail(&format!("{e}"));
    }
}
