//! Conformance for `CollectIlluminaLaneMetrics` against Picard, through the port's own binary.
//!
//! Every run of `tools/illumina-conformance/IlluminaLaneMetricsDump.java` is replayed with the same
//! tile metrics file, written record by record, and the same arguments; the rows the dump prints
//! are rebuilt from the metrics files the binary wrote (comments dropped, as the dump drops them)
//! or from its error. The rebuilt text must equal the golden line for line.

mod common;

use std::path::Path;

use common::{run, table, written, Dump};

/// One record of a tile metrics file: a lane, a tile, a metric code and a value.
type Metric = (u16, u16, u16, f32);

/// A tile metrics file of exactly the records given, version two with ten-byte records.
fn write_tile_metrics(file: &Path, metrics: &[Metric]) {
    let mut bytes = vec![2u8, 10];
    for (lane, tile, code, value) in metrics {
        bytes.extend_from_slice(&lane.to_le_bytes());
        bytes.extend_from_slice(&tile.to_le_bytes());
        bytes.extend_from_slice(&code.to_le_bytes());
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    std::fs::create_dir_all(file.parent().unwrap()).unwrap();
    std::fs::write(file, bytes).unwrap();
}

/// The dump's `tile`: the four counts, and a phasing pair per template read at descriptor indices
/// `{0}` for one template read and `{0, 2}` for two.
fn tile(lane: u16, tile: u16, clusters: f32, passing: f32, templates: usize) -> Vec<Metric> {
    let mut metrics = vec![
        (lane, tile, 100, clusters),
        (lane, tile, 101, passing),
        (lane, tile, 102, clusters * 10.0),
        (lane, tile, 103, passing * 10.0),
    ];
    let descriptors: &[u16] = if templates == 1 { &[0] } else { &[0, 2] };
    for &index in descriptors {
        metrics.push((lane, tile, 200 + index * 2, 0.1 * (index + 1) as f32));
        metrics.push((lane, tile, 201 + index * 2, 0.2 * (index + 1) as f32));
    }
    metrics
}

/// `IlluminaLaneMetricsDump.run`.
fn lane(dump: &mut Dump, case: &str, metrics: &[Metric], extra: &[&str]) {
    let dir = dump.dir();
    let run_dir = dir.join("run");
    write_tile_metrics(&run_dir.join("InterOp/TileMetricsOut.bin"), metrics);
    std::fs::create_dir_all(run_dir.join("Data/Intensities/BaseCalls")).unwrap();
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let mut argv = vec![
        format!("RUN_DIRECTORY={}", run_dir.display()),
        format!("OUTPUT_DIRECTORY={}", out.display()),
    ];
    argv.extend(extra.iter().map(|argument| argument.to_string()));
    if !extra.iter().any(|argument| argument.starts_with("O=")) {
        argv.push("O=metrics".to_string());
    }
    if !extra
        .iter()
        .any(|argument| argument.starts_with("READ_STRUCTURE="))
    {
        argv.push("READ_STRUCTURE=4T".to_string());
    }
    let ran = run(env!("CARGO_BIN_EXE_CollectIlluminaLaneMetrics"), &argv);
    if dump.failed(case, &dir, &ran, false) {
        return;
    }
    let files = written(&out);
    dump.emit("files", case, &files.join(" "));
    for file in files {
        dump.emit(
            "metrics",
            &format!("{case}.{file}"),
            &table(&out.join(&file)),
        );
    }
}

#[test]
fn lane_metrics_match_the_golden() {
    let mut dump = Dump::new("illumina_lane_metrics");
    lane(&mut dump, "one-tile", &tile(1, 1101, 1000.0, 800.0, 1), &[]);
    let mut two_tiles = tile(1, 1101, 1000.0, 800.0, 1);
    two_tiles.extend(tile(1, 1102, 500.0, 400.0, 1));
    lane(&mut dump, "two-tiles", &two_tiles, &[]);
    let mut two_lanes = two_tiles.clone();
    two_lanes.extend(tile(2, 1101, 200.0, 100.0, 1));
    lane(&mut dump, "two-lanes", &two_lanes, &[]);

    lane(
        &mut dump,
        "counts-without-densities",
        &[(1, 1101, 100, 1000.0), (1, 1101, 101, 800.0)],
        &[],
    );
    lane(
        &mut dump,
        "densities-without-counts",
        &[(1, 1101, 102, 10000.0), (1, 1101, 103, 8000.0)],
        &[],
    );
    lane(
        &mut dump,
        "a-code-nobody-reads",
        &[(1, 1101, 100, 1000.0), (1, 1101, 200, 42.0)],
        &[],
    );

    lane(
        &mut dump,
        "a-named-prefix",
        &tile(1, 1101, 1000.0, 800.0, 1),
        &["O=run17"],
    );
    lane(
        &mut dump,
        "no-read-structure",
        &tile(1, 1101, 1000.0, 800.0, 1),
        &["READ_STRUCTURE=null"],
    );
    lane(
        &mut dump,
        "two-template-reads",
        &tile(1, 1101, 1000.0, 800.0, 2),
        &["READ_STRUCTURE=4T8B4T"],
    );
    lane(
        &mut dump,
        "two-template-reads-with-one-phasing-pair",
        &tile(1, 1101, 1000.0, 800.0, 1),
        &["READ_STRUCTURE=4T8B4T"],
    );
    lane(
        &mut dump,
        "half-a-phasing-pair",
        &[
            (1, 1101, 100, 1000.0),
            (1, 1101, 101, 800.0),
            (1, 1101, 102, 10000.0),
            (1, 1101, 103, 8000.0),
            (1, 1101, 200, 0.1),
        ],
        &[],
    );
    dump.check("illumina_lane_metrics");
}
