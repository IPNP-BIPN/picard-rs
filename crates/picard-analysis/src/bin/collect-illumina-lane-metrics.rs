//! `CollectIlluminaLaneMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.CollectIlluminaLaneMetrics` and its `IlluminaLaneMetricsCollector` at tag
//! 3.4.0, with `TileMetricsUtil.parseTileMetrics` and `LanePhasingMetricsCollector`, for version-2
//! tile metrics.
//!
//! A lane's density is its clusters over its area, the area being each tile's clusters over its
//! density; a lane's phasing is the median over its tiles, times a hundred (version 2 stores a
//! fraction). The tiles reach both sums in the order of a `HashMap` keyed by `"lane:tile"`.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprint::{java_hash_order, reference_path};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::java_hash_map::string_hash_code;
use picard_analysis::metrics_cli::{thrown, Args};

struct Row {
    class: &'static str,
    columns: &'static [&'static str],
    values: Vec<Value>,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        self.class
    }
    fn columns(&self) -> &[&'static str] {
        self.columns
    }
    fn values(&self) -> Vec<Value> {
        self.values.clone()
    }
}

/// One tile: its density, its cluster count, and its phasing per template read (FIRST, SECOND).
struct Tile {
    lane: i32,
    density: f32,
    clusters: f32,
    phasing: Vec<(&'static str, f32, f32)>,
}

/// `MathUtil.median`.
fn median(values: &[f32]) -> f64 {
    let mut data: Vec<f64> = values.iter().map(|v| f64::from(*v)).collect();
    data.sort_by(|a, b| a.total_cmp(b));
    let middle = data.len() / 2;
    if data.len() % 2 == 1 {
        data[middle]
    } else {
        (data[middle - 1] + data[middle]) / 2.0
    }
}

/// The descriptors of `RunInfo.xml`: each `<Read>`'s cycle count, and whether it is an index.
fn run_info(path: &std::path::Path) -> Vec<(usize, bool)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
            reference_path(&path.display().to_string())
        ))
    });
    let attribute = |tag: &str, name: &str| -> String {
        let key = format!("{name}=\"");
        tag.find(&key)
            .map(|i| &tag[i + key.len()..])
            .and_then(|r| r.split('"').next())
            .unwrap_or("")
            .to_string()
    };
    let mut reads = Vec::new();
    for (i, tag) in text.split("<Read ").skip(1).enumerate() {
        let number: usize = attribute(tag, "Number").parse().unwrap_or(0);
        if number != i + 1 {
            thrown(&format!(
                "picard.PicardException: Read number in RunInfo.xml was out of order: {} != {number}",
                i + 1
            ));
        }
        let cycles = attribute(tag, "NumCycles").parse().unwrap_or(0);
        let indexed = attribute(tag, "IsIndexedRead").eq_ignore_ascii_case("Y");
        reads.push((cycles, indexed));
    }
    reads
}

fn main() {
    let args = Args::from_env(&[
        ("R", "RUN_DIRECTORY"),
        ("O", "OUTPUT_PREFIX"),
        ("RS", "READ_STRUCTURE"),
        ("EXT", "FILE_EXTENSION"),
    ]);
    let run = std::path::PathBuf::from(args.required("RUN_DIRECTORY"));
    let output_dir = std::path::PathBuf::from(args.required("OUTPUT_DIRECTORY"));
    let prefix = args.required("OUTPUT_PREFIX");
    let extension = args.get("FILE_EXTENSION").unwrap_or("").to_string();
    let strict = !matches!(
        args.get("VALIDATION_STRINGENCY"),
        Some("LENIENT") | Some("SILENT")
    );
    let lenient = args.get("VALIDATION_STRINGENCY") == Some("LENIENT");
    // The descriptors: template or not, in order.
    let templates: Vec<bool> = match args.get("READ_STRUCTURE") {
        Some(rs) => parse_read_structure(rs)
            .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"))
            .iter()
            .map(|s| s.kind == SegmentKind::Template)
            .collect(),
        None => run_info(&run.join("RunInfo.xml"))
            .iter()
            .map(|(_, indexed)| !indexed)
            .collect(),
    };

    let file = run.join("InterOp").join("TileMetricsOut.bin");
    let bytes = std::fs::read(&file).unwrap_or_else(|_| {
        thrown(&format!(
            "java.lang.IllegalStateException: No InterOp file found in {} or any of its cycle directories.\"",
            reference_path(&run.join("InterOp").display().to_string())
        ))
    });
    let records = picard_analysis::illumina_files::parse_tile_metrics(&bytes)
        .unwrap_or_else(|| thrown("picard.PicardException: Unsupported tile metrics version"));
    // `determineLastValueForLaneTileMetricsCode`, then by location.
    let mut last: Vec<((u16, u16, u16), f32)> = Vec::new();
    for r in &records {
        let key = (r.lane, r.tile, r.code);
        match last.iter_mut().find(|(k, _)| *k == key) {
            Some(slot) => slot.1 = r.value,
            None => last.push((key, r.value)),
        }
    }
    let mut locations: Vec<(String, Vec<(u16, f32)>)> = Vec::new();
    for ((lane, tile, code), value) in &last {
        let key = format!("{lane}:{tile}");
        match locations.iter_mut().find(|(k, _)| *k == key) {
            Some((_, codes)) => codes.push((*code, *value)),
            None => locations.push((key, vec![(*code, *value)])),
        }
    }
    let locations = java_hash_order(
        locations
            .into_iter()
            .map(|(k, v)| (string_hash_code(&k), (k, v)))
            .collect(),
    );
    let mut tiles: Vec<Tile> = Vec::new();
    for (key, codes) in &locations {
        let get = |c: u16| codes.iter().find(|(code, _)| *code == c).map(|(_, v)| *v);
        // `IlluminaMetricsCode`: 100 is DENSITY_ID and 102 is CLUSTER_ID.
        let (Some(density), Some(clusters)) = (get(100), get(102)) else {
            let mut observed: Vec<u16> = codes.iter().map(|(c, _)| *c).collect();
            observed = java_hash_order(observed.into_iter().map(|c| (i32::from(c), c)).collect());
            thrown(&format!(
                "picard.PicardException: Expected to find cluster and density record codes (102 and 100) in records read for tile location {key} (lane:tile), but found only [{}].",
                observed.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", ")
            ));
        };
        let mut phasing = Vec::new();
        let mut first = true;
        for (index, is_template) in templates.iter().enumerate() {
            if !is_template {
                continue;
            }
            let read = if first { "FIRST" } else { "SECOND" };
            let p = 200 + index as u16 * 2;
            let q = 201 + index as u16 * 2;
            let (a, b) = match (get(p), get(q)) {
                (Some(a), Some(b)) => (a, b),
                (pa, pb) => {
                    let message = format!(
                        "Don't have both phasing and prephasing values for {read} read cycle {}.  Phasing code was {p} and prephasing code was {q}.",
                        index + 1
                    );
                    if pa.is_none() && pb.is_none() && !strict {
                        if lenient {
                            eprintln!("WARN\t1970-01-01 00:00:00\tTileMetricsUtil\t{message}");
                        }
                        (0.0, 0.0)
                    } else {
                        thrown(&format!("picard.PicardException: {message}"));
                    }
                }
            };
            match phasing.iter_mut().find(|(r, _, _)| *r == read) {
                Some(slot) => *slot = (read, a, b),
                None => phasing.push((read, a, b)),
            }
            first = false;
        }
        let lane: i32 = key
            .split(':')
            .next()
            .and_then(|l| l.parse().ok())
            .unwrap_or(0);
        tiles.push(Tile {
            lane,
            density,
            clusters,
            phasing,
        });
    }
    tiles.retain(|t| t.lane > 0);
    let mut lanes: Vec<i32> = tiles.iter().map(|t| t.lane).collect();
    lanes.dedup();
    let lanes: Vec<i32> = java_hash_order({
        let mut seen: Vec<i32> = Vec::new();
        for l in lanes {
            if !seen.contains(&l) {
                seen.push(l);
            }
        }
        seen.into_iter().map(|l| (l, l)).collect()
    });

    let header = |file: &mut MetricsFile| {
        file.add_header("CollectIlluminaLaneMetrics <command line>");
        file.add_header("Started on: <timestamp>");
    };
    let mut lane_file = MetricsFile::new();
    header(&mut lane_file);
    let mut phasing_file = MetricsFile::new();
    header(&mut phasing_file);
    for &lane in &lanes {
        let of_lane: Vec<&Tile> = tiles.iter().filter(|t| t.lane == lane).collect();
        let (mut area, mut clusters) = (0.0f64, 0.0f64);
        for t in &of_lane {
            if t.density > 0.0 {
                area += f64::from(t.clusters / t.density);
            }
            clusters += f64::from(t.clusters);
        }
        lane_file.add_metric(&Row {
            class: "picard.illumina.IlluminaLaneMetrics",
            columns: &["CLUSTER_DENSITY", "LANE"],
            values: vec![
                Value::Double(if area > 0.0 { clusters / area } else { 0.0 }),
                Value::Long(i64::from(lane)),
            ],
        });
        for read in ["FIRST", "SECOND"] {
            let values: Vec<(f32, f32)> = of_lane
                .iter()
                .flat_map(|t| {
                    t.phasing
                        .iter()
                        .filter(|(r, _, _)| *r == read)
                        .map(|(_, a, b)| (*a, *b))
                })
                .collect();
            if values.is_empty() {
                continue;
            }
            let a: Vec<f32> = values.iter().map(|v| v.0).collect();
            let b: Vec<f32> = values.iter().map(|v| v.1).collect();
            phasing_file.add_metric(&Row {
                class: "picard.illumina.IlluminaPhasingMetrics",
                columns: &["LANE", "TYPE_NAME", "PHASING_APPLIED", "PREPHASING_APPLIED"],
                values: vec![
                    Value::Long(i64::from(lane)),
                    Value::Str(read.to_string()),
                    Value::Double(f64::from(median(&a) as f32 * 100.0)),
                    Value::Double(f64::from(median(&b) as f32 * 100.0)),
                ],
            });
        }
    }
    for (name, file) in [
        (
            format!("{prefix}.illumina_lane_metrics{extension}"),
            &lane_file,
        ),
        (
            format!("{prefix}.illumina_phasing_metrics{extension}"),
            &phasing_file,
        ),
    ] {
        if let Err(e) = std::fs::write(output_dir.join(name), file.write()) {
            thrown(&format!("htsjdk.samtools.SAMException: {e}"));
        }
    }
}
