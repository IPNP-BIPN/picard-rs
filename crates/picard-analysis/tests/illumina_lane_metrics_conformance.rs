//! Conformance for `CollectIlluminaLaneMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/illumina-conformance/IlluminaLaneMetricsDump.java`: tile metrics files
//! written record by record, so every case is a statement about the metric CODES rather than about
//! a run.
//!
//! # What this suite is for
//!
//!  * **the codes: 100 and 101 the counts, 102 and 103 the densities, 200 and 201 the phasing**;
//!  * **the phasing code being offset by twice the read DESCRIPTOR's index, not the read's**;
//!  * **a lane being the MEAN over its tiles**;
//!  * **and the two refusals, each in the reference's own words.**

use std::io::Read;

use picard_analysis::illumina_files::TileMetric;
use picard_analysis::illumina_lane_metrics::{
    collect, phasing_code, Refusal, PHASING_BASE, PREPHASING_BASE,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/illumina_lane_metrics.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

/// The rows of a metrics table, by column name.
fn rows(text: &str, kind: &str, case: &str) -> Vec<Vec<String>> {
    match field(text, kind, case) {
        None => Vec::new(),
        Some(table) => table
            .lines()
            .skip(1)
            .filter(|line| !line.trim().is_empty())
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect(),
    }
}

/// The fixture's own tile: four counts and a phasing pair per template read.
fn tile(
    lane: u16,
    tile: u16,
    clusters: f32,
    passing: f32,
    descriptors: &[usize],
) -> Vec<TileMetric> {
    let mut metrics = vec![
        TileMetric {
            lane,
            tile,
            code: 100,
            value: clusters,
        },
        TileMetric {
            lane,
            tile,
            code: 101,
            value: passing,
        },
        TileMetric {
            lane,
            tile,
            code: 102,
            value: clusters * 10.0,
        },
        TileMetric {
            lane,
            tile,
            code: 103,
            value: passing * 10.0,
        },
    ];
    // The values are the DUMP's: a tenth and a fifth times the descriptor's index plus one, which
    // is what makes the second template read's numbers 0.3 and 0.6 rather than 0.2 and 0.4.
    for descriptor in descriptors {
        metrics.push(TileMetric {
            lane,
            tile,
            code: phasing_code(*descriptor, PHASING_BASE),
            value: 0.1 * (*descriptor + 1) as f32,
        });
        metrics.push(TileMetric {
            lane,
            tile,
            code: phasing_code(*descriptor, PREPHASING_BASE),
            value: 0.2 * (*descriptor + 1) as f32,
        });
    }
    metrics
}

/// A lane is the mean over its tiles, and two lanes are two rows.
#[test]
fn a_lane_is_the_mean_over_its_tiles() {
    let text = corpus();
    let (lanes, _) =
        collect(&tile(1, 1101, 1000.0, 800.0, &[0]), &[0]).expect("a well formed file");
    assert_eq!(lanes.len(), 1);
    assert_eq!(lanes[0].cluster_density, 1000.0);
    assert_eq!(
        rows(&text, "metrics", "one-tile.metrics.illumina_lane_metrics")[0][0],
        "1000"
    );

    let mut two = tile(1, 1101, 1000.0, 800.0, &[0]);
    two.extend(tile(1, 1102, 500.0, 400.0, &[0]));
    let (lanes, _) = collect(&two, &[0]).expect("a well formed file");
    assert_eq!(lanes[0].cluster_density, 750.0);
    assert_eq!(
        rows(&text, "metrics", "two-tiles.metrics.illumina_lane_metrics")[0][0],
        "750"
    );

    let mut three = two.clone();
    three.extend(tile(2, 1101, 200.0, 100.0, &[0]));
    let (lanes, _) = collect(&three, &[0]).expect("a well formed file");
    let recorded = rows(&text, "metrics", "two-lanes.metrics.illumina_lane_metrics");
    assert_eq!(lanes.len(), recorded.len());
    for (lane, row) in lanes.iter().zip(recorded) {
        assert_eq!(lane.cluster_density.to_string(), row[0]);
        assert_eq!(lane.lane.to_string(), row[1]);
    }
}

/// The phasing code is the DESCRIPTOR's index, which is why `4T8B4T` uses 204 and 205.
#[test]
fn the_phasing_code_counts_descriptors() {
    assert_eq!(phasing_code(0, PHASING_BASE), 200);
    assert_eq!(phasing_code(0, PREPHASING_BASE), 201);
    assert_eq!(phasing_code(2, PHASING_BASE), 204);
    assert_eq!(phasing_code(2, PREPHASING_BASE), 205);

    let text = corpus();
    // Two template reads at descriptors nought and two, which the golden's `4T8B4T` case wrote.
    let (_, phasing) =
        collect(&tile(1, 1101, 1000.0, 800.0, &[0, 2]), &[0, 2]).expect("a well formed file");
    let recorded = rows(
        &text,
        "metrics",
        "two-template-reads.metrics.illumina_phasing_metrics",
    );
    assert_eq!(phasing.len(), recorded.len());
    for (row, expected) in phasing.iter().zip(recorded) {
        assert_eq!(row.lane.to_string(), expected[0]);
        assert_eq!(
            match row.read {
                0 => "FIRST",
                _ => "SECOND",
            },
            expected[1]
        );
        let phasing_applied: f64 = expected[2].parse().expect("a number");
        assert!((phasing_applied - row.phasing_applied).abs() < 1e-4);
    }
    // A file whose second pair is at the OBVIOUS place rather than the descriptor's is refused.
    let wrong = tile(1, 1101, 1000.0, 800.0, &[0, 1]);
    let refusal = collect(&wrong, &[0, 2]).expect_err("the refusal");
    assert_eq!(
        refusal,
        Refusal::HalfAPhasingPair {
            which: "SECOND",
            cycle: 5,
            phasing: 204,
            prephasing: 205,
        }
    );
    // And the golden refused the same file, naming the same read.
    let recorded = field(&text, "error", "two-template-reads-with-one-phasing-pair")
        .expect("the golden's refusal");
    assert!(recorded.contains("Don't have both phasing and prephasing values for SECOND read"));
}

/// A file missing the counts or half a phasing pair is refused rather than reported with a gap.
#[test]
fn the_refusals_are_the_goldens() {
    let text = corpus();
    // Counts without densities, and densities without counts.
    let counts_only = vec![
        TileMetric {
            lane: 1,
            tile: 1101,
            code: 100,
            value: 1000.0,
        },
        TileMetric {
            lane: 1,
            tile: 1101,
            code: 101,
            value: 800.0,
        },
    ];
    assert_eq!(
        collect(&counts_only, &[0]).expect_err("the refusal"),
        Refusal::MissingCounts {
            lane: 1,
            tile: 1101
        }
    );
    assert!(field(&text, "error", "counts-without-densities")
        .expect("the golden")
        .contains("Expected to find cluster and density record codes (102 and 100)"));
    let densities_only = vec![
        TileMetric {
            lane: 1,
            tile: 1101,
            code: 102,
            value: 10000.0,
        },
        TileMetric {
            lane: 1,
            tile: 1101,
            code: 103,
            value: 8000.0,
        },
    ];
    assert!(matches!(
        collect(&densities_only, &[0]),
        Err(Refusal::MissingCounts { .. })
    ));

    // Half a phasing pair, for the first read.
    let mut half = counts_only.clone();
    half.push(TileMetric {
        lane: 1,
        tile: 1101,
        code: 102,
        value: 10000.0,
    });
    half.push(TileMetric {
        lane: 1,
        tile: 1101,
        code: 103,
        value: 8000.0,
    });
    half.push(TileMetric {
        lane: 1,
        tile: 1101,
        code: 200,
        value: 0.1,
    });
    assert_eq!(
        collect(&half, &[0]).expect_err("the refusal"),
        Refusal::HalfAPhasingPair {
            which: "FIRST",
            cycle: 1,
            phasing: 200,
            prephasing: 201,
        }
    );
    assert!(field(&text, "error", "half-a-phasing-pair")
        .expect("the golden")
        .contains("Don't have both phasing and prephasing values for FIRST read cycle 1"));
}
