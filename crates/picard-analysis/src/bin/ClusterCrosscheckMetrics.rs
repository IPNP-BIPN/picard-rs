//! `ClusterCrosscheckMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.ClusterCrosscheckMetrics.doWork` at tag 3.4.0, with the parts of
//! htsjdk 4.2.0 it reads and writes through. The graph and the rule that decides which rows come
//! back inside a cluster are `picard_analysis::cluster_crosscheck_metrics`; this is the table
//! around them, and the ORDER the rows are written in.
//!
//! # The order is a hash order, twice
//!
//! The reference does not sort what it writes:
//!
//! * the clusters come out of a `HashMap<Integer, Set<String>>`, so they are visited by bucket
//!   (the cluster number, modulo the table's width), and a number is the INDEX of the component's
//!   first node, not a count, so the numbers 0 and 3 are two clusters of three;
//! * within a cluster the rows are collected into a `HashSet<ClusteredCrosscheckMetric>`, which is
//!   visited by the bucket of `MetricBase.hashCode`: the hashes of every public field FORMATTED as
//!   it is written, folded with 31, in `getFields()` order (the subclass's two fields first).
//!
//! A `HashSet` also drops a row that equals an earlier one, which `MetricBase.equals` decides on
//! the same formatted values.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_metrics::format::{format_double, format_long};
use picard_analysis::cluster_crosscheck_metrics::Graph;
use picard_analysis::java_hash_map::{string_hash_code, JavaHashMap};
use picard_analysis::metrics_cli::{fail, thrown, Args};

const TOOL: &str = "ClusterCrosscheckMetrics";

/// The columns of `ClusteredCrosscheckMetric`, in `getFields()` order, with the Java type of each
/// as one letter: `s` string, `e` enum, `d` double, `i` integer.
const COLUMNS: &[(&str, char)] = &[
    ("CLUSTER", 'i'),
    ("CLUSTER_SIZE", 'i'),
    ("LEFT_GROUP_VALUE", 's'),
    ("RIGHT_GROUP_VALUE", 's'),
    ("RESULT", 'e'),
    ("DATA_TYPE", 'e'),
    ("LOD_SCORE", 'd'),
    ("LOD_SCORE_TUMOR_NORMAL", 'd'),
    ("LOD_SCORE_NORMAL_TUMOR", 'd'),
    ("LEFT_RUN_BARCODE", 's'),
    ("LEFT_LANE", 'i'),
    ("LEFT_MOLECULAR_BARCODE_SEQUENCE", 's'),
    ("LEFT_LIBRARY", 's'),
    ("LEFT_SAMPLE", 's'),
    ("LEFT_FILE", 's'),
    ("RIGHT_RUN_BARCODE", 's'),
    ("RIGHT_LANE", 'i'),
    ("RIGHT_MOLECULAR_BARCODE_SEQUENCE", 's'),
    ("RIGHT_LIBRARY", 's'),
    ("RIGHT_SAMPLE", 's'),
    ("RIGHT_FILE", 's'),
];

/// The names of every column, for the writer.
const NAMES: &[&str] = &[
    "CLUSTER",
    "CLUSTER_SIZE",
    "LEFT_GROUP_VALUE",
    "RIGHT_GROUP_VALUE",
    "RESULT",
    "DATA_TYPE",
    "LOD_SCORE",
    "LOD_SCORE_TUMOR_NORMAL",
    "LOD_SCORE_NORMAL_TUMOR",
    "LEFT_RUN_BARCODE",
    "LEFT_LANE",
    "LEFT_MOLECULAR_BARCODE_SEQUENCE",
    "LEFT_LIBRARY",
    "LEFT_SAMPLE",
    "LEFT_FILE",
    "RIGHT_RUN_BARCODE",
    "RIGHT_LANE",
    "RIGHT_MOLECULAR_BARCODE_SEQUENCE",
    "RIGHT_LIBRARY",
    "RIGHT_SAMPLE",
    "RIGHT_FILE",
];

/// One row of the input, every field already in the form it is written in.
#[derive(Clone)]
struct Metric {
    /// The 19 `CrosscheckMetric` fields, formatted, in the order of `COLUMNS[2..]`.
    fields: Vec<String>,
    lod: f64,
}

impl Metric {
    fn left(&self) -> &str {
        &self.fields[0]
    }
    fn right(&self) -> &str {
        &self.fields[1]
    }
}

/// A row of the output: the cluster and its size, then the input's fields.
struct Clustered {
    cluster: usize,
    size: usize,
    metric: Metric,
}

impl Clustered {
    fn formatted(&self) -> Vec<String> {
        let mut all = vec![
            format_long(self.cluster as i64),
            format_long(self.size as i64),
        ];
        all.extend(self.metric.fields.iter().cloned());
        all
    }

    /// `MetricBase.hashCode`.
    fn hash_code(&self) -> i32 {
        let mut result: i32 = 0;
        for text in self.formatted() {
            result = result
                .wrapping_mul(31)
                .wrapping_add(string_hash_code(&text));
        }
        result
    }
}

impl MetricBean for Clustered {
    fn class_name(&self) -> &str {
        "picard.fingerprint.ClusteredCrosscheckMetric"
    }
    fn columns(&self) -> &[&'static str] {
        NAMES
    }
    fn values(&self) -> Vec<Value> {
        self.formatted()
            .into_iter()
            .zip(COLUMNS)
            .map(|(text, (_, kind))| match kind {
                'i' if !text.is_empty() => Value::Long(text.parse().unwrap_or(0)),
                _ => Value::Str(text),
            })
            .collect()
    }
}

/// `HashMap.hash` followed by the table index: the bucket of an `int` hash code.
fn bucket(hash: i32, capacity: usize) -> usize {
    let h = hash as u32;
    ((h ^ (h >> 16)) as usize) & (capacity - 1)
}

/// The width a `HashMap` or `HashSet` has reached after `count` insertions: 16 to start, doubling
/// each time the size passes three quarters of it.
fn capacity_for(count: usize) -> usize {
    let mut capacity = 16;
    while count > capacity * 3 / 4 {
        capacity *= 2;
    }
    capacity
}

/// Items in the order a Java hash container would hand them out, given the order they went in:
/// by bucket, and in insertion order inside a bucket (a resize keeps each bucket's relative order).
fn in_hash_order<T>(items: Vec<T>, hash: impl Fn(&T) -> i32) -> Vec<T> {
    let capacity = capacity_for(items.len());
    let mut keyed: Vec<(usize, T)> = items
        .into_iter()
        .map(|item| (bucket(hash(&item), capacity), item))
        .collect();
    keyed.sort_by_key(|(index, _)| *index);
    keyed.into_iter().map(|(_, item)| item).collect()
}

/// `MetricsFile.read`, for a `CrosscheckMetric` file.
fn read_metrics(text: &str) -> Vec<Metric> {
    let mut lines = text.lines();
    let mut found = false;
    for line in lines.by_ref() {
        if line.starts_with("## METRICS CLASS") {
            found = true;
            break;
        }
    }
    if !found {
        return Vec::new();
    }
    let header: Vec<&str> = match lines.next() {
        Some(line) => line.split('\t').collect(),
        None => return Vec::new(),
    };
    let mut metrics = Vec::new();
    for line in lines {
        if line.is_empty() || line.starts_with("## HISTOGRAM") {
            break;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let mut fields = Vec::with_capacity(19);
        let mut lod = f64::NAN;
        for (name, kind) in &COLUMNS[2..] {
            let raw = header
                .iter()
                .position(|column| column == name)
                .and_then(|at| cells.get(at).copied())
                .unwrap_or("");
            let formatted = match kind {
                'd' if !raw.is_empty() => {
                    let value = if raw == "?" || raw == "-?" {
                        f64::NAN
                    } else {
                        raw.parse::<f64>().unwrap_or_else(|_| {
                            thrown(&format!(
                                "java.lang.NumberFormatException: For input string: \"{raw}\""
                            ))
                        })
                    };
                    if *name == "LOD_SCORE" {
                        lod = value;
                    }
                    format_double(value)
                }
                'i' if !raw.is_empty() => format_long(raw.parse::<i64>().unwrap_or_else(|_| {
                    thrown(&format!(
                        "java.lang.NumberFormatException: For input string: \"{raw}\""
                    ))
                })),
                _ => raw.to_string(),
            };
            fields.push(formatted);
        }
        metrics.push(Metric { fields, lod });
    }
    metrics
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("LOD", "LOD_THRESHOLD")]);
    let input = args.required("INPUT");
    let output = args.get("OUTPUT").map(str::to_string);
    let threshold = args.double("LOD_THRESHOLD", 0.0);
    if let Some(stringency) = args.get("VALIDATION_STRINGENCY") {
        if !matches!(stringency, "STRICT" | "LENIENT" | "SILENT") {
            fail(&format!(
                "Argument 'VALIDATION_STRINGENCY' cannot be set to '{stringency}'"
            ));
        }
    }

    let text = std::fs::read_to_string(&input).unwrap_or_else(|e| thrown(&format!("{e}")));
    let metrics = read_metrics(&text);

    let mut graph = Graph::new();
    for metric in &metrics {
        if metric.lod > threshold {
            graph.add_edge(metric.left(), metric.right());
        }
    }
    let clusters = graph.cluster();

    // `Collectors.toMap` into a HashMap<String, Integer>, in node order; then `groupingBy` over
    // its entries into a HashMap<Integer, ...>, whose keys are met in that map's order.
    let mut by_name: JavaHashMap<usize> = JavaHashMap::new();
    for node in graph.nodes() {
        by_name.put(node, clusters[node]);
    }
    let mut met: Vec<usize> = Vec::new();
    for (_, cluster) in by_name.iter() {
        if !met.contains(cluster) {
            met.push(*cluster);
        }
    }
    let ordered = in_hash_order(met, |cluster| *cluster as i32);

    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for cluster in ordered {
        let members: Vec<&String> = graph
            .nodes()
            .iter()
            .filter(|node| clusters[*node] == cluster)
            .collect();
        let size = members.len();
        // `collect(toSet())`: a HashSet, so a row equal to an earlier one is dropped, and the rest
        // come out in hash order.
        let mut rows: Vec<Clustered> = Vec::new();
        for metric in &metrics {
            if members.iter().any(|node| node.as_str() == metric.left())
                && members.iter().any(|node| node.as_str() == metric.right())
            {
                let row = Clustered {
                    cluster,
                    size,
                    metric: metric.clone(),
                };
                if !rows.iter().any(|held| held.formatted() == row.formatted()) {
                    rows.push(row);
                }
            }
        }
        for row in in_hash_order(rows, |row| row.hash_code()) {
            file.add_metric(&row);
        }
    }

    match output {
        Some(path) => {
            if let Err(e) = std::fs::write(&path, file.write()) {
                thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}"));
            }
        }
        None => print!("{}", file.write()),
    }
}
