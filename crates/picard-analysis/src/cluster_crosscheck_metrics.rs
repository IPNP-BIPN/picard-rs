//! `ClusterCrosscheckMetrics`: the connected components of a crosscheck table.
//!
//! Reading and writing the metrics file are not ported. The graph is, along with the rule that
//! decides which rows come back inside each cluster.
//!
//! Ported from `picard.fingerprint.ClusterCrosscheckMetrics` and
//! `picard.util.GraphUtils` in Picard 3.4.0.

use std::collections::BTreeMap;

/// `ClusterCrosscheckMetrics.LOD_THRESHOLD`.
pub const DEFAULT_LOD_THRESHOLD: f64 = 0.0;

/// One row of the input table, reduced to the three fields the clustering reads.
#[derive(Debug, Clone, PartialEq)]
pub struct CrosscheckMetric {
    pub left_group_value: String,
    pub right_group_value: String,
    pub lod_score: f64,
}

/// One row of the output, which is an input row with two fields added.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusteredCrosscheckMetric {
    pub metric: CrosscheckMetric,
    pub cluster: usize,
    /// The number of GROUPS in the cluster, not the number of rows.
    pub cluster_size: usize,
}

/// `GraphUtils.Graph`: nodes in first-seen order, and a union-find over their indices.
///
/// The cluster identifier a node ends up with is the INDEX of its component's representative, not
/// a counter, so the identifiers are not contiguous: two clusters of two are numbered 0 and 2.
#[derive(Debug, Clone, Default)]
pub struct Graph {
    nodes: Vec<String>,
    neighbours: Vec<Vec<usize>>,
}

impl Graph {
    pub fn new() -> Graph {
        Graph::default()
    }

    /// `addNode`, which answers the index and never adds a name twice.
    pub fn add_node(&mut self, node: &str) -> usize {
        if let Some(index) = self.nodes.iter().position(|held| held == node) {
            return index;
        }
        self.nodes.push(node.to_string());
        self.neighbours.push(Vec::new());
        self.nodes.len() - 1
    }

    /// `addEdge`, which is bidirectional.
    ///
    /// The reference guards against a self-edge with a REFERENCE comparison rather than an equality
    /// one, so two equal strings read separately from a file do not trip it. What that costs is
    /// nothing observable: a node joined to itself is already its own representative.
    pub fn add_edge(&mut self, left: &str, right: &str) {
        let left_index = self.add_node(left);
        let right_index = self.add_node(right);
        if left_index == right_index {
            return;
        }
        if !self.neighbours[left_index].contains(&right_index) {
            self.neighbours[left_index].push(right_index);
        }
        if !self.neighbours[right_index].contains(&left_index) {
            self.neighbours[right_index].push(left_index);
        }
    }

    pub fn nodes(&self) -> &[String] {
        &self.nodes
    }

    /// `cluster`: union-find over the node indices, answering each node's representative index.
    pub fn cluster(&self) -> BTreeMap<String, usize> {
        let mut representative: Vec<usize> = (0..self.nodes.len()).collect();
        for (i, neighbours) in self.neighbours.iter().enumerate() {
            for j in neighbours {
                join(&mut representative, *j, i);
            }
        }
        self.nodes
            .iter()
            .enumerate()
            .map(|(index, node)| (node.clone(), find(&representative, index)))
            .collect()
    }
}

fn find(representative: &[usize], mut node: usize) -> usize {
    while node != representative[node] {
        node = representative[node];
    }
    node
}

fn join(representative: &mut [usize], first: usize, second: usize) {
    let a = find(representative, first);
    let b = find(representative, second);
    if a == b {
        return;
    }
    representative[a] = b;
}

/// `clusterMetrics`: the graph of comparisons above the threshold, then every row whose BOTH
/// sides sit in one cluster.
///
/// An edge needs a LOD STRICTLY above the threshold, so a comparison exactly at it makes none. A
/// row's own LOD then plays no further part: a row well under the threshold comes back if the two
/// groups it names were joined by other rows. A group in no edge at all is not in the graph, so
/// its rows are dropped rather than forming a cluster of one.
///
/// The reference collects the rows into a HashSet, so a duplicated comparison appears once and the
/// order is a hash. This returns them sorted, by cluster and then by the two group names.
pub fn cluster_metrics(
    metrics: &[CrosscheckMetric],
    lod_threshold: f64,
) -> Vec<ClusteredCrosscheckMetric> {
    let mut graph = Graph::new();
    for metric in metrics {
        if metric.lod_score > lod_threshold {
            graph.add_edge(&metric.left_group_value, &metric.right_group_value);
        }
    }
    let clusters = graph.cluster();
    let mut sizes: BTreeMap<usize, usize> = BTreeMap::new();
    for cluster in clusters.values() {
        *sizes.entry(*cluster).or_default() += 1;
    }
    let mut rows: Vec<ClusteredCrosscheckMetric> = Vec::new();
    for metric in metrics {
        let (Some(left), Some(right)) = (
            clusters.get(&metric.left_group_value),
            clusters.get(&metric.right_group_value),
        ) else {
            continue;
        };
        if left != right {
            continue;
        }
        let row = ClusteredCrosscheckMetric {
            metric: metric.clone(),
            cluster: *left,
            cluster_size: sizes[left],
        };
        if !rows.contains(&row) {
            rows.push(row);
        }
    }
    rows.sort_by(|a, b| {
        a.cluster.cmp(&b.cluster).then(
            a.metric
                .left_group_value
                .cmp(&b.metric.left_group_value)
                .then(a.metric.right_group_value.cmp(&b.metric.right_group_value)),
        )
    });
    rows
}
