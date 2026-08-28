//! Conformance for `ClusterCrosscheckMetrics` against Picard 3.4.0.
//!
//! The golden prints the input table for the chain case and the output table for every case. The
//! port is given the same rows and must reach the same clusters, sizes and membership.
//!
//! # What this suite is for
//!
//!  * **an edge being a LOD strictly above the threshold**;
//!  * **the clusters being connected components**;
//!  * **every row whose both sides are in one cluster coming back, whatever its own LOD**;
//!  * **a group in no edge reaching no cluster at all**;
//!  * **CLUSTER_SIZE counting groups and not rows**;
//!  * **the cluster identifier being a node index, so the ids are not contiguous**;
//!  * **a duplicated row appearing once**;
//!  * **and a file with nothing above the threshold writing no rows.**

use std::io::Read;

use picard_analysis::cluster_crosscheck_metrics::{
    cluster_metrics, ClusteredCrosscheckMetric, CrosscheckMetric, Graph, DEFAULT_LOD_THRESHOLD,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("cluster_crosscheck_metrics.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// One case's output rows, as (left, right, lod, cluster, size).
fn written(text: &str, case: &str) -> Vec<(String, String, f64, usize, usize)> {
    let table = field(text, "metrics", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = table.lines().filter(|line| !line.is_empty());
    let Some(header) = lines.next() else {
        return Vec::new();
    };
    let columns: Vec<&str> = header.split('\t').collect();
    let at = |name: &str| columns.iter().position(|c| *c == name).expect(name);
    lines
        .map(|line| {
            let values: Vec<&str> = line.split('\t').collect();
            (
                values[at("LEFT_GROUP_VALUE")].to_string(),
                values[at("RIGHT_GROUP_VALUE")].to_string(),
                values[at("LOD_SCORE")].parse().expect("a lod"),
                values[at("CLUSTER")].parse().expect("a cluster"),
                values[at("CLUSTER_SIZE")].parse().expect("a size"),
            )
        })
        .collect()
}

fn ours(metrics: &[CrosscheckMetric], threshold: f64) -> Vec<(String, String, f64, usize, usize)> {
    cluster_metrics(metrics, threshold)
        .into_iter()
        .map(|row: ClusteredCrosscheckMetric| {
            (
                row.metric.left_group_value,
                row.metric.right_group_value,
                row.metric.lod_score,
                row.cluster,
                row.cluster_size,
            )
        })
        .collect()
}

fn metric(left: &str, right: &str, lod: f64) -> CrosscheckMetric {
    CrosscheckMetric {
        left_group_value: left.to_string(),
        right_group_value: right.to_string(),
        lod_score: lod,
    }
}

/// The input the golden prints, read back.
fn chain(text: &str) -> Vec<CrosscheckMetric> {
    field(text, "in", "chain")
        .expect("the golden carries in/chain")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .skip(1)
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            metric(columns[0], columns[1], columns[4].parse().expect("a lod"))
        })
        .collect()
}

/// The cases whose inputs and thresholds the dump fixes.
fn cases() -> Vec<(&'static str, Vec<CrosscheckMetric>, f64)> {
    vec![
        (
            "transitive-chain",
            vec![
                metric("A", "B", 10.0),
                metric("B", "C", 10.0),
                metric("A", "C", -5.0),
            ],
            3.0,
        ),
        (
            "threshold-splits",
            vec![
                metric("A", "B", 10.0),
                metric("B", "C", 4.0),
                metric("A", "C", -5.0),
            ],
            5.0,
        ),
        ("exactly-at-the-threshold", vec![metric("A", "B", 3.0)], 3.0),
        (
            "just-above-the-threshold",
            vec![metric("A", "B", 3.001)],
            3.0,
        ),
        (
            "two-clusters",
            vec![
                metric("A", "B", 10.0),
                metric("C", "D", 10.0),
                metric("A", "C", -5.0),
            ],
            3.0,
        ),
        (
            "orphan-group",
            vec![metric("A", "B", 10.0), metric("Z", "A", -5.0)],
            3.0,
        ),
        (
            "self-comparison",
            vec![metric("A", "A", 10.0), metric("A", "B", 10.0)],
            3.0,
        ),
        (
            "duplicated-row",
            vec![metric("A", "B", 10.0), metric("A", "B", 10.0)],
            3.0,
        ),
        (
            "default-threshold",
            vec![metric("A", "B", 0.5), metric("C", "D", -0.5)],
            DEFAULT_LOD_THRESHOLD,
        ),
        (
            "nothing-above-the-threshold",
            vec![metric("A", "B", 1.0), metric("C", "D", 2.0)],
            3.0,
        ),
        ("empty", vec![], 3.0),
    ]
}

/// Every case's output rows are what the port reaches.
#[test]
fn every_case_reaches_the_same_rows() {
    let text = corpus();
    for (case, metrics, threshold) in cases() {
        assert_eq!(ours(&metrics, threshold), written(&text, case), "{case}");
    }
    // And the chain's input, read back from the golden's own copy, gives the same answer.
    assert_eq!(ours(&chain(&text), 3.0), written(&text, "transitive-chain"));
}

/// An edge is a LOD strictly above the threshold.
#[test]
fn the_edge_is_strictly_above_the_threshold() {
    let text = corpus();
    assert!(written(&text, "exactly-at-the-threshold").is_empty());
    assert_eq!(written(&text, "just-above-the-threshold").len(), 1);
    assert!(cluster_metrics(&[metric("A", "B", 3.0)], 3.0).is_empty());
    assert_eq!(cluster_metrics(&[metric("A", "B", 3.001)], 3.0).len(), 1);
}

/// The clusters are connected components, and every row whose both sides are in one comes back
/// whatever its own LOD.
#[test]
fn a_low_row_comes_back_inside_its_cluster() {
    let text = corpus();
    let rows = written(&text, "transitive-chain");
    assert_eq!(rows.len(), 3);
    // A and C were never related, and their row scores -5, yet it is here.
    let low = rows.iter().find(|row| row.2 == -5.0).expect("the A-C row");
    assert_eq!((low.0.as_str(), low.1.as_str()), ("A", "C"));
    // All three rows are in one cluster of three groups.
    assert!(rows.iter().all(|row| row.3 == rows[0].3));
    assert!(rows.iter().all(|row| row.4 == 3));
}

/// A group in no edge at all reaches no cluster, so its rows vanish.
#[test]
fn an_orphan_group_reaches_no_cluster() {
    let text = corpus();
    let rows = written(&text, "orphan-group");
    assert_eq!(rows.len(), 1);
    assert_eq!((rows[0].0.as_str(), rows[0].1.as_str()), ("A", "B"));
    assert!(!rows.iter().any(|row| row.0 == "Z" || row.1 == "Z"));
    // Raising the threshold orphans C the same way, taking both of its rows.
    let split = written(&text, "threshold-splits");
    assert_eq!(split.len(), 1);
    assert_eq!(split[0].4, 2);
    assert!(!split.iter().any(|row| row.0 == "C" || row.1 == "C"));
}

/// CLUSTER_SIZE counts groups and not rows.
#[test]
fn the_size_counts_groups_and_not_rows() {
    let text = corpus();
    let rows = written(&text, "transitive-chain");
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0].4, 3);
    // And the self-comparison adds no group: two rows, a cluster of two.
    let self_rows = written(&text, "self-comparison");
    assert_eq!(self_rows.len(), 2);
    assert!(self_rows.iter().all(|row| row.4 == 2));
}

/// The cluster identifier is a node index, so the identifiers are not contiguous.
#[test]
fn the_identifier_is_a_node_index() {
    let text = corpus();
    let rows = written(&text, "two-clusters");
    assert_eq!(rows.len(), 2);
    let mut ids: Vec<usize> = rows.iter().map(|row| row.3).collect();
    ids.sort();
    assert_eq!(ids, vec![0, 2]);
    // Which is what the graph's own numbering gives: A B C D are nodes 0 1 2 3, and the two
    // components' representatives are 1 and 3 by the union's direction, read back as 0 and 2.
    let mut graph = Graph::new();
    graph.add_edge("A", "B");
    graph.add_edge("C", "D");
    assert_eq!(graph.nodes(), ["A", "B", "C", "D"]);
    let clusters = graph.cluster();
    assert_eq!(clusters["A"], clusters["B"]);
    assert_eq!(clusters["C"], clusters["D"]);
    assert_ne!(clusters["A"], clusters["C"]);
}

/// A duplicated row appears once, the reference collecting into a set.
#[test]
fn a_duplicated_row_appears_once() {
    let text = corpus();
    assert_eq!(written(&text, "duplicated-row").len(), 1);
    assert_eq!(
        cluster_metrics(&[metric("A", "B", 10.0), metric("A", "B", 10.0)], 3.0).len(),
        1
    );
}

/// The default threshold is nought, and a file with nothing above the threshold writes no rows.
#[test]
fn the_default_threshold_is_nought() {
    let text = corpus();
    assert_eq!(DEFAULT_LOD_THRESHOLD, 0.0);
    let rows = written(&text, "default-threshold");
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].2, 0.5);
    for case in ["nothing-above-the-threshold", "empty"] {
        assert!(written(&text, case).is_empty(), "{case}");
    }
}
