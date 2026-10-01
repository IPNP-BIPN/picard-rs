//! `ClusterCrosscheckMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.ClusterCrosscheckMetrics.doWork` and `picard.util.GraphUtils` at tag
//! 3.4.0:
//!
//! * the input read as `MetricsFile<CrosscheckMetric>` (an empty cell is a null field, `?` a NaN);
//! * an edge for every row whose LOD is STRICTLY above `LOD_THRESHOLD` (a null LOD is the
//!   reference's `NullPointerException`), nodes numbered in first-seen order, union-find joining
//!   each neighbour's root under the node's own, and a cluster named by its root's index;
//! * the clusters visited in the order of a `HashMap<Integer, Set<String>>` built by
//!   `groupingBy` (which inserts at the head of a bucket) and then `toMap`, and within each cluster
//!   every input row whose two groups are both in it, collected into a `HashSet`: equal rows
//!   collapse and the rest come out in the order of `MetricBase.hashCode`, which hashes each
//!   field's FORMATTED value in declaration order (`CLUSTER` and `CLUSTER_SIZE` first).

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprinting::JavaMap;
use picard_analysis::java_hash_map::{string_hash_code, JavaHashMap};
use picard_analysis::metrics_cli::{fail, thrown, Args};

const TOOL: &str = "ClusterCrosscheckMetrics";

/// `CrosscheckMetric`'s columns, in declaration order.
const COLUMNS: [&str; 19] = [
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
const DOUBLES: [usize; 3] = [4, 5, 6];
const INTEGERS: [usize; 2] = [8, 14];
const RESULTS: [&str; 5] = [
    "EXPECTED_MATCH",
    "EXPECTED_MISMATCH",
    "UNEXPECTED_MATCH",
    "UNEXPECTED_MISMATCH",
    "INCONCLUSIVE",
];
const DATA_TYPES: [&str; 4] = ["FILE", "SAMPLE", "LIBRARY", "READGROUP"];

/// One row, each field as the value it was parsed into.
#[derive(Clone)]
struct Metric {
    fields: Vec<Value>,
}

impl Metric {
    fn string(&self, i: usize) -> Option<&str> {
        match &self.fields[i] {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }
}

#[derive(Clone)]
struct Clustered {
    cluster: i64,
    size: i64,
    metric: Metric,
}

impl Clustered {
    /// Every field formatted the way the file writes it, in `getFields()` order.
    fn formatted(&self) -> Vec<String> {
        let mut out = vec![
            Value::Long(self.cluster).format(),
            Value::Long(self.size).format(),
        ];
        out.extend(self.metric.fields.iter().map(Value::format));
        out
    }

    /// `MetricBase.hashCode`.
    fn hash(&self) -> i32 {
        let mut result: i32 = 0;
        for f in self.formatted() {
            result = result.wrapping_mul(31).wrapping_add(string_hash_code(&f));
        }
        result
    }
}

/// `MetricBase.equals`: every field's formatted value.
#[derive(Clone)]
struct Key(Vec<String>);

impl PartialEq for Key {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl MetricBean for Clustered {
    fn class_name(&self) -> &str {
        "picard.fingerprint.ClusteredCrosscheckMetric"
    }
    fn columns(&self) -> &[&'static str] {
        &[
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
        ]
    }
    fn values(&self) -> Vec<Value> {
        let mut out = vec![Value::Long(self.cluster), Value::Long(self.size)];
        out.extend(self.metric.fields.iter().cloned());
        out
    }
}

/// `FormatUtil.parseDouble`.
fn parse_double(text: &str) -> Option<f64> {
    if text == "?" || text == "-?" {
        return Some(f64::NAN);
    }
    htsjdk_vcf::genotype_likelihoods::parse_java_double(text)
}

/// `MetricsFile.read` for a `CrosscheckMetric` table.
fn read_metrics(text: &str) -> Vec<Metric> {
    let mut lines = text.lines();
    let mut columns: Option<Vec<String>> = None;
    for line in lines.by_ref() {
        if line.starts_with("## METRICS CLASS") {
            if let Some(header) = lines.next() {
                columns = Some(header.split('\t').map(str::to_string).collect());
            }
            break;
        }
    }
    let Some(columns) = columns else {
        return Vec::new();
    };
    let positions: Vec<usize> = columns
        .iter()
        .map(|c| {
            COLUMNS.iter().position(|k| k == c).unwrap_or_else(|| {
                thrown(&format!(
                    "htsjdk.samtools.SAMException: Could not find field with name {c} in metric class picard.fingerprint.CrosscheckMetric"
                ))
            })
        })
        .collect();
    let mut out = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            break;
        }
        let mut fields = vec![Value::Null; COLUMNS.len()];
        for (cell, &field) in line.split('\t').zip(positions.iter()) {
            if cell.is_empty() {
                continue;
            }
            fields[field] = if DOUBLES.contains(&field) {
                Value::Double(parse_double(cell).unwrap_or_else(|| {
                    thrown(&format!(
                        "java.lang.NumberFormatException: For input string: \"{cell}\""
                    ))
                }))
            } else if INTEGERS.contains(&field) {
                Value::Long(i64::from(cell.parse::<i32>().unwrap_or_else(|_| {
                    thrown(&format!(
                        "java.lang.NumberFormatException: For input string: \"{cell}\""
                    ))
                })))
            } else if field == 2 && !RESULTS.contains(&cell) {
                thrown(&format!(
                    "java.lang.IllegalArgumentException: No enum constant picard.fingerprint.CrosscheckMetric.FingerprintResult.{cell}"
                ))
            } else if field == 3 && !DATA_TYPES.contains(&cell) {
                thrown(&format!(
                    "java.lang.IllegalArgumentException: No enum constant picard.fingerprint.CrosscheckMetric.DataType.{cell}"
                ))
            } else {
                Value::Str(cell.to_string())
            };
        }
        out.push(Metric { fields });
    }
    out
}

/// `GraphUtils.Graph`.
#[derive(Default)]
struct Graph {
    nodes: Vec<String>,
    neighbors: Vec<Vec<usize>>,
}

impl Graph {
    fn add_node(&mut self, node: &str) -> usize {
        if let Some(i) = self.nodes.iter().position(|n| n == node) {
            return i;
        }
        self.nodes.push(node.to_string());
        self.neighbors.push(Vec::new());
        self.nodes.len() - 1
    }

    /// `addEdge`: its `left == right` guard compares references, and two strings parsed from a
    /// file are never the same object, so a self-comparison adds the node as its own neighbour.
    fn add_edge(&mut self, left: &str, right: &str) {
        let l = self.add_node(left);
        let r = self.add_node(right);
        if !self.neighbors[l].contains(&r) {
            self.neighbors[l].push(r);
        }
        if !self.neighbors[r].contains(&l) {
            self.neighbors[r].push(l);
        }
    }

    /// `cluster`: each node's representative index.
    fn cluster(&self) -> Vec<usize> {
        let mut grouping: Vec<usize> = (0..self.nodes.len()).collect();
        fn find(g: &mut [usize], mut node: usize) -> usize {
            let mut rep = node;
            while rep != g[rep] {
                rep = g[rep];
            }
            while node != rep {
                let next = g[node];
                g[node] = rep;
                node = next;
            }
            rep
        }
        for i in 0..self.neighbors.len() {
            for &j in &self.neighbors[i] {
                let a = find(&mut grouping, j);
                let b = find(&mut grouping, i);
                if a != b {
                    grouping[a] = b;
                }
            }
        }
        (0..self.nodes.len())
            .map(|n| find(&mut grouping, n))
            .collect()
    }
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("LOD", "LOD_THRESHOLD")]);
    let input = args.required("INPUT");
    let output = args.get("OUTPUT").map(str::to_string);
    let threshold = args.double("LOD_THRESHOLD", 0.0);

    let text = std::fs::read_to_string(&input).unwrap_or_else(|e| fail(&e.to_string()));
    let metrics = read_metrics(&text);

    let mut graph = Graph::default();
    for m in &metrics {
        let lod = match &m.fields[4] {
            Value::Double(d) => *d,
            _ => thrown(
                "java.lang.NullPointerException: Cannot invoke \"java.lang.Double.doubleValue()\" because \"metric.LOD_SCORE\" is null",
            ),
        };
        if lod > threshold {
            graph.add_edge(m.string(0).unwrap_or("null"), m.string(1).unwrap_or("null"));
        }
    }
    let reps = graph.cluster();

    // `clusters`: a HashMap<String, Integer> filled in node order.
    let mut clusters: JavaHashMap<usize> = JavaHashMap::new();
    for (node, rep) in graph.nodes.iter().zip(reps.iter()) {
        if clusters.get(node).is_none() {
            clusters.put(node, *rep);
        }
    }
    // `groupingBy(getValue)`: an Integer-keyed HashMap whose new keys go to their bucket's head.
    let mut grouped: Vec<(i32, Vec<String>)> = Vec::new();
    let mut order: JavaMap<i32, usize> = JavaMap::new();
    let mut heads: Vec<(i32, usize)> = Vec::new();
    for (node, rep) in clusters.iter() {
        let key = *rep as i32;
        if let Some(&slot) = order.get(key, &key) {
            grouped[slot].1.push(node.to_string());
        } else {
            grouped.push((key, vec![node.to_string()]));
            heads.push((key, grouped.len() - 1));
            order.put(key, key, grouped.len() - 1);
        }
    }
    // Rebuild the grouping map with head insertion to read its iteration order.
    let mut by_head: Vec<Vec<(i32, usize)>> = Vec::new();
    let mut capacity = 16usize;
    for (n, (key, slot)) in heads.iter().enumerate() {
        if by_head.is_empty() {
            by_head = vec![Vec::new(); capacity];
        } else if n > capacity * 3 / 4 {
            capacity *= 2;
            let mut grown = vec![Vec::new(); capacity];
            for bucket in by_head {
                for entry in bucket {
                    let idx = spread(entry.0) & (capacity - 1);
                    grown[idx].push(entry);
                }
            }
            by_head = grown;
        }
        let idx = spread(*key) & (capacity - 1);
        by_head[idx].insert(0, (*key, *slot));
    }
    // `toMap`: a second HashMap, filled in that order.
    let mut collection: JavaMap<i32, usize> = JavaMap::new();
    for (key, slot) in by_head.into_iter().flatten() {
        collection.put(key, key, slot);
    }

    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for (key, &slot) in collection.iter() {
        let members = &grouped[slot].1;
        let mut set: JavaMap<Key, Clustered> = JavaMap::new();
        for m in &metrics {
            let left = m.string(0).unwrap_or("null");
            let right = m.string(1).unwrap_or("null");
            if !(members.iter().any(|n| n == left) && members.iter().any(|n| n == right)) {
                continue;
            }
            let row = Clustered {
                cluster: i64::from(*key),
                size: members.len() as i64,
                metric: m.clone(),
            };
            let hash = row.hash();
            let k = Key(row.formatted());
            if !set.contains_key(hash, &k) {
                set.put(hash, k, row);
            }
        }
        for (_, row) in set.iter() {
            file.add_metric(row);
        }
    }
    let text = file.write();
    match &output {
        Some(path) => std::fs::write(path, text).unwrap_or_else(|e| fail(&e.to_string())),
        None => thrown("java.lang.NullPointerException"),
    }
}

/// `HashMap.hash` for an `Integer` key.
fn spread(h: i32) -> usize {
    let h = h as u32;
    (h ^ (h >> 16)) as usize
}
