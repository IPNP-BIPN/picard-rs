//! `CollectIlluminaBasecallingMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.CollectIlluminaBasecallingMetrics` at tag 3.4.0 over
//! `picard_analysis::illumina_reader`. Each barcode (or the whole lane when no barcodes are
//! declared) is a row of cluster counts per tile: means and standard deviations rounded to whole
//! clusters, the PF percentage formatted by `DecimalFormat("#.##")` and read back, and the read and
//! base totals from the template reads. The rows are in barcode order (a `TreeMap`), then the lane.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::illumina_reader::{has_barcode_files, matched_barcodes, Run};
use picard_analysis::metrics_cli::{thrown, Args};
use std::collections::BTreeMap;

const COLUMNS: &[&str] = &[
    "LANE",
    "MOLECULAR_BARCODE_SEQUENCE_1",
    "MOLECULAR_BARCODE_NAME",
    "TOTAL_BASES",
    "PF_BASES",
    "TOTAL_READS",
    "PF_READS",
    "TOTAL_CLUSTERS",
    "PF_CLUSTERS",
    "MEAN_CLUSTERS_PER_TILE",
    "SD_CLUSTERS_PER_TILE",
    "MEAN_PCT_PF_CLUSTERS_PER_TILE",
    "SD_PCT_PF_CLUSTERS_PER_TILE",
    "MEAN_PF_CLUSTERS_PER_TILE",
    "SD_PF_CLUSTERS_PER_TILE",
];

struct Row(Vec<Value>);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.illumina.IlluminaBasecallingMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        self.0.clone()
    }
}

/// Two `Histogram<Integer>`s by tile: clusters and PF clusters.
#[derive(Default, Clone)]
struct Counts {
    name: Option<String>,
    clusters: BTreeMap<i32, f64>,
    pf: BTreeMap<i32, f64>,
}

fn mean(h: &BTreeMap<i32, f64>) -> f64 {
    let mut sum = 0.0;
    for v in h.values() {
        sum += v;
    }
    sum / h.len() as f64
}

fn sd(h: &BTreeMap<i32, f64>, mean: f64) -> f64 {
    let mut total = 0.0;
    for v in h.values() {
        total += (v - mean).powi(2);
    }
    (total / (h.len().saturating_sub(1).max(1)) as f64).sqrt()
}

/// `Math.round(double)`.
fn round(v: f64) -> f64 {
    if v.is_nan() {
        0.0
    } else {
        (v + 0.5).floor()
    }
}

/// `Double.valueOf(new DecimalFormat("#.##").format(v))`: two places, half-even.
fn two_places(v: f64) -> f64 {
    if v.is_nan() {
        return f64::NAN;
    }
    format!("{v:.2}").parse().unwrap_or(0.0)
}

fn row(lane: i32, barcode: Option<&str>, c: &Counts, templates: i64, template_bases: i64) -> Row {
    let m = mean(&c.clusters);
    let pm = mean(&c.pf);
    let pct: BTreeMap<i32, f64> = c.clusters.iter().map(|(k, v)| (*k, c.pf[k] / v)).collect();
    let pctm = mean(&pct);
    let total: f64 = c.clusters.values().sum();
    let pf: f64 = c.pf.values().sum();
    let (total, pf) = (total as i64, pf as i64);
    let text = |v: Option<&str>| v.map_or(Value::Null, |s| Value::Str(s.to_string()));
    Row(vec![
        Value::Str(lane.to_string()),
        text(barcode),
        text(c.name.as_deref()),
        Value::Long(total * template_bases),
        Value::Long(pf * template_bases),
        Value::Long(total * templates),
        Value::Long(pf * templates),
        Value::Long(total),
        Value::Long(pf),
        Value::Double(round(m)),
        Value::Double(round(sd(&c.clusters, m))),
        Value::Double(if pctm.is_nan() {
            0.0
        } else {
            two_places(pctm * 100.0)
        }),
        Value::Double(two_places(sd(&pct, pctm) * 100.0)),
        Value::Double(round(pm)),
        Value::Double(round(sd(&c.pf, pm))),
    ])
}

fn main() {
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("BARCODES_DIR", "BARCODES_DIR"),
        ("I", "INPUT"),
        ("L", "LANE"),
        ("RS", "READ_STRUCTURE"),
        ("O", "OUTPUT"),
    ]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let lane: i32 = args.required("LANE").parse().unwrap_or(0);
    let output = args.required("OUTPUT");
    let structure = parse_read_structure(&args.required("READ_STRUCTURE"))
        .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"));
    let sample_barcodes = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Barcode)
        .count();
    let templates: Vec<usize> = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Template)
        .map(|s| s.cycles)
        .collect();
    let barcodes_dir = args
        .get("BARCODES_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| basecalls.clone());

    let mut counts: BTreeMap<String, Counts> = BTreeMap::new();
    let mut barcode_length = 0usize;
    let declared = args.get("INPUT").map(str::to_string);
    if let Some(path) = &declared {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| thrown(&format!("{e}")));
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let col = |n: &str| header.iter().position(|h| *h == n);
        for line in lines.filter(|l| !l.trim().is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            let name = col("barcode_name")
                .and_then(|i| fields.get(i))
                .map(|s| s.to_string());
            let mut barcode = String::new();
            for i in 1..=sample_barcodes {
                if let Some(v) = col(&format!("barcode_sequence_{i}")).and_then(|c| fields.get(c)) {
                    barcode.push_str(v);
                }
                if barcode_length == 0 {
                    barcode_length = barcode.len();
                }
            }
            if !barcode.is_empty() {
                counts.insert(
                    barcode,
                    Counts {
                        name,
                        ..Counts::default()
                    },
                );
            }
        }
    }
    let with_barcodes = declared.is_some() && !counts.is_empty();
    if with_barcodes && !has_barcode_files(&barcodes_dir, lane) {
        thrown("picard.PicardException: Could not find a format with available files for the following data types: Barcodes");
    }
    let unmatched = "N".repeat(barcode_length);
    let run = Run::new(&basecalls, lane);
    let tiles = run.available_tiles().unwrap_or_else(|e| thrown(&e));
    for &tile in &tiles {
        let clusters = run.clusters(tile, &[], 2).unwrap_or_else(|e| thrown(&e));
        let matched = if with_barcodes {
            matched_barcodes(&barcodes_dir, lane, tile).unwrap_or_default()
        } else {
            Vec::new()
        };
        for (i, c) in clusters.iter().enumerate() {
            let barcode = matched
                .get(i)
                .cloned()
                .flatten()
                .unwrap_or_else(|| unmatched.clone());
            let entry = counts.entry(barcode).or_default();
            *entry.clusters.entry(tile).or_insert(0.0) += 1.0;
            *entry.pf.entry(tile).or_insert(0.0) += if c.pf { 1.0 } else { 0.0 };
        }
    }
    let template_bases: i64 = templates.iter().map(|c| *c as i64).sum();
    let n_templates = templates.len() as i64;
    let mut file = MetricsFile::new();
    file.add_header("CollectIlluminaBasecallingMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    let mut all = Counts::default();
    for (barcode, c) in &counts {
        file.add_metric(&row(lane, Some(barcode), c, n_templates, template_bases));
        for (k, v) in &c.clusters {
            *all.clusters.entry(*k).or_insert(0.0) += v;
        }
        for (k, v) in &c.pf {
            *all.pf.entry(*k).or_insert(0.0) += v;
        }
    }
    if !counts.contains_key("") {
        file.add_metric(&row(lane, None, &all, n_templates, template_bases));
    }
    if let Err(e) = std::fs::write(&output, file.write()) {
        thrown(&format!(
            "picard.PicardException: Error writing output file {output}: {e}"
        ));
    }
}
