//! `CollectIlluminaLaneMetrics` as a runnable binary: the covering array's port side.
//!
//! Two files from `InterOp/TileMetricsOut.bin`: a lane's cluster density, and the median phasing
//! and prephasing of each template read. What decides their bytes:
//!
//!  * the tile metrics are de-duplicated by (lane, tile, code) in a `HashMap` whose LAST record
//!    wins, then grouped by "lane:tile" in a `HashMap<String, List>` and by lane in a
//!    `HashMap<Integer, List>`, and those iteration orders are the order the lane density's float
//!    sums are taken in and the rows are written in;
//!  * code 100 is read as the DENSITY and code 102 as the CLUSTER COUNT (the reference's own
//!    `IlluminaMetricsCode` names), and the lane density is the cluster sum over the sum of each
//!    tile's `clusters / density`, a float division widened into a double sum;
//!  * a template read whose phasing codes are both missing is a refusal under STRICT and a zero
//!    otherwise; one of the two missing is a refusal whatever the stringency;
//!  * version 2 reports the phasing medians times one hundred (`usePercentage`), in float.
//!
//! The read structure comes from `READ_STRUCTURE` or from `RunInfo.xml`'s `Read` elements.
//! Version 3 tile metrics, which need `EmpiricalPhasingMetricsOut.bin`, are refused as
//! unsupported.
//!
//! Ported from `picard.illumina.CollectIlluminaLaneMetrics`, `IlluminaLaneMetrics`,
//! `IlluminaPhasingMetrics`, `LanePhasingMetricsCollector`, `parser.TileMetricsUtil` and
//! `parser.Tile` in Picard 3.4.0.

use std::path::{Path, PathBuf};

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::illumina_dir::{
    absolute, assert_file_is_readable, read_tile_metrics, unsupported, ReadDescriptor,
    ReadStructure, ReadType, Thrown, TileMetricRecord,
};
use picard_analysis::java_hash_map::{JavaHashMap, JavaHashMapBy};
use picard_analysis::metrics_cli::{fail, thrown, Args};

struct LaneRow {
    cluster_density: f64,
    lane: i64,
}

impl MetricBean for LaneRow {
    fn class_name(&self) -> &str {
        "picard.illumina.IlluminaLaneMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &["CLUSTER_DENSITY", "LANE"]
    }
    fn values(&self) -> Vec<Value> {
        vec![Value::Double(self.cluster_density), Value::Long(self.lane)]
    }
}

struct PhasingRow {
    lane: i64,
    type_name: &'static str,
    phasing: f64,
    prephasing: f64,
}

impl MetricBean for PhasingRow {
    fn class_name(&self) -> &str {
        "picard.illumina.IlluminaPhasingMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &["LANE", "TYPE_NAME", "PHASING_APPLIED", "PREPHASING_APPLIED"]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Long(self.lane),
            Value::Str(self.type_name.to_string()),
            Value::Double(self.phasing),
            Value::Double(self.prephasing),
        ]
    }
}

/// `Tile`: the density (code 100), the cluster count (code 102) and the per-read phasing.
struct Tile {
    lane: i32,
    density: f32,
    clusters: f32,
    /// FIRST then SECOND, as an `EnumMap` holds them.
    phasing: [Option<(f32, f32)>; 2],
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Stringency {
    Strict,
    Lenient,
    Silent,
}

fn main() {
    let args = Args::from_env(&[
        ("O", "OUTPUT_PREFIX"),
        ("RS", "READ_STRUCTURE"),
        ("EXT", "FILE_EXTENSION"),
    ]);
    let run_directory = PathBuf::from(args.required("RUN_DIRECTORY"));
    let output_directory = PathBuf::from(args.required("OUTPUT_DIRECTORY"));
    let output_prefix = args.required("OUTPUT_PREFIX");
    let extension = args.get("FILE_EXTENSION").unwrap_or("").to_string();
    let stringency = match args.get("VALIDATION_STRINGENCY").unwrap_or("STRICT") {
        "STRICT" => Stringency::Strict,
        "LENIENT" => Stringency::Lenient,
        "SILENT" => Stringency::Silent,
        other => fail(&format!(
            "Argument 'VALIDATION_STRINGENCY' cannot be set to '{other}'"
        )),
    };
    // Barclay builds the argument through `new ReadStructure(String)`.
    let given = args.get("READ_STRUCTURE").map(|text| {
        ReadStructure::parse(text).unwrap_or_else(|e| {
            fail(&format!(
                "Failed to parse value '{text}' for argument READ_STRUCTURE: {}",
                e.message
            ))
        })
    });

    let structure = match given {
        Some(structure) => structure,
        None => match read_structure_from_run_info(&run_directory) {
            Ok(structure) => structure,
            Err(e) => thrown(&e.render()),
        },
    };
    let (lanes, version) = match collect(&run_directory, &structure, stringency) {
        Ok(result) => result,
        Err(e) => thrown(&e.render()),
    };

    let mut lane_file = metrics_file();
    for (lane, tiles) in &lanes {
        let mut area = 0f64;
        let mut clusters = 0f64;
        for tile in tiles {
            if tile.density > 0.0 {
                area += (tile.clusters / tile.density) as f64;
            }
            clusters += tile.clusters as f64;
        }
        lane_file.add_metric(&LaneRow {
            cluster_density: if area > 0.0 { clusters / area } else { 0.0 },
            lane: *lane as i64,
        });
    }
    write(
        &output_directory,
        &format!("{output_prefix}.illumina_lane_metrics{extension}"),
        &lane_file,
    );

    let mut phasing_file = metrics_file();
    for (lane, tiles) in &lanes {
        for (read, name) in [(0usize, "FIRST"), (1, "SECOND")] {
            let values: Vec<(f32, f32)> = tiles.iter().filter_map(|t| t.phasing[read]).collect();
            if values.is_empty() {
                continue;
            }
            let median = |pick: &dyn Fn(&(f32, f32)) -> f32| -> f32 {
                let mut data: Vec<f64> = values.iter().map(|v| pick(v) as f64).collect();
                data.sort_by(|a, b| a.total_cmp(b));
                let middle = data.len() / 2;
                let m = if data.len() % 2 == 1 {
                    data[middle]
                } else {
                    (data[middle - 1] + data[middle]) / 2.0
                } as f32;
                if version == 2 {
                    m * 100.0
                } else {
                    m
                }
            };
            phasing_file.add_metric(&PhasingRow {
                lane: *lane as i64,
                type_name: name,
                phasing: median(&|v| v.0) as f64,
                prephasing: median(&|v| v.1) as f64,
            });
        }
    }
    write(
        &output_directory,
        &format!("{output_prefix}.illumina_phasing_metrics{extension}"),
        &phasing_file,
    );
}

fn metrics_file() -> MetricsFile {
    let mut file = MetricsFile::new();
    file.add_header("CollectIlluminaLaneMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    file
}

fn write(directory: &Path, name: &str, file: &MetricsFile) {
    if let Err(e) = std::fs::write(directory.join(name), file.write()) {
        fail(&format!("{e}"));
    }
}

/// `RunInfo.xml`'s `Read` elements, in document order.
fn read_structure_from_run_info(run_directory: &Path) -> Result<ReadStructure, Thrown> {
    let run_info = PathBuf::from(format!("{}/RunInfo.xml", run_directory.display()));
    assert_file_is_readable(&run_info)?;
    let text = std::fs::read_to_string(&run_info).map_err(|e| Thrown::picard(e.to_string()))?;
    let wrap = |e: Thrown| Thrown::picard(e.message);
    let mut descriptors = Vec::new();
    for (i, attributes) in read_elements(&text).into_iter().enumerate() {
        let attribute = |name: &str| -> Result<String, Thrown> {
            attributes
                .iter()
                .find(|(k, _)| k == name)
                .map(|(_, v)| v.clone())
                .ok_or_else(|| {
                    Thrown::picard(
                        "Cannot invoke \"org.w3c.dom.Node.getNodeValue()\" because the return \
                         value of \"org.w3c.dom.NamedNodeMap.getNamedItem(String)\" is null",
                    )
                })
        };
        let parse = |value: String| -> Result<i32, Thrown> {
            value
                .parse()
                .map_err(|_| Thrown::picard(format!("For input string: \"{value}\"")))
        };
        let number = parse(attribute("Number")?)?;
        let cycles = parse(attribute("NumCycles")?)?;
        let indexed = attribute("IsIndexedRead")?.to_uppercase() == "Y";
        if number != i as i32 + 1 {
            return Err(Thrown::picard(format!(
                "Read number in RunInfo.xml was out of order: {} != {number}",
                i + 1
            )));
        }
        descriptors.push(ReadDescriptor {
            length: cycles,
            kind: if indexed { ReadType::B } else { ReadType::T },
        });
    }
    ReadStructure::new(descriptors).map_err(wrap)
}

/// The attributes of every `<Read ...>` start tag, in order.
fn read_elements(text: &str) -> Vec<Vec<(String, String)>> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(at) = rest.find('<') {
        rest = &rest[at + 1..];
        let end = rest.find('>').unwrap_or(rest.len());
        let tag = &rest[..end];
        rest = &rest[end.min(rest.len())..];
        let name_end = tag
            .find(|c: char| c.is_whitespace() || c == '/')
            .unwrap_or(tag.len());
        if &tag[..name_end] != "Read" {
            continue;
        }
        let mut attributes = Vec::new();
        let mut body = &tag[name_end..];
        while let Some(eq) = body.find('=') {
            let name = body[..eq].trim().to_string();
            let after = body[eq + 1..].trim_start();
            let Some(quote) = after.chars().next() else {
                break;
            };
            let value_rest = &after[1..];
            let close = value_rest.find(quote).unwrap_or(value_rest.len());
            let value = value_rest[..close]
                .replace("&lt;", "<")
                .replace("&gt;", ">")
                .replace("&quot;", "\"")
                .replace("&apos;", "'")
                .replace("&amp;", "&");
            attributes.push((name, value));
            body = &value_rest[(close + 1).min(value_rest.len())..];
        }
        out.push(attributes);
    }
    out
}

/// `TileMetricsUtil.findTileMetricsFiles`.
fn find_tile_metrics_files(run_directory: &Path, cycles: i32) -> Result<Vec<PathBuf>, Thrown> {
    let interop = run_directory.join("InterOp");
    let mut paths = vec![interop.join("TileMetricsOut.bin")];
    for cycle in (1..=cycles).rev() {
        paths.push(interop.join(format!("C{cycle}.1/TileMetricsOut.bin")));
    }
    let found: Vec<PathBuf> = paths.into_iter().filter(|p| p.exists()).collect();
    if found.is_empty() {
        return Err(Thrown::new(
            "java.lang.IllegalStateException",
            format!(
                "No InterOp file found in {} or any of its cycle directories.\"",
                absolute(&interop)
            ),
        ));
    }
    Ok(found)
}

/// `IlluminaLaneTileCode.hashCode`.
fn lane_tile_code_hash(record: &TileMetricRecord) -> i32 {
    let mut result = record.lane;
    result = result.wrapping_mul(31).wrapping_add(record.tile);
    result.wrapping_mul(31).wrapping_add(record.code)
}

type Lanes = Vec<(i32, Vec<Tile>)>;

/// `collectLaneMetrics` up to the writing: the lanes in `HashMap<Integer>` order, each with its
/// tiles in the order they were grouped, and the tile metrics version.
fn collect(
    run_directory: &Path,
    structure: &ReadStructure,
    stringency: Stringency,
) -> Result<(Lanes, i32), Thrown> {
    let files = find_tile_metrics_files(run_directory, structure.total_cycles)?;
    let version = read_tile_metrics(&files[0])?.version;
    for file in &files {
        if read_tile_metrics(file)?.version != version {
            return Err(Thrown::picard(format!(
                "Not all tile metrics files match expected version: {version}"
            )));
        }
    }
    if version == 3 {
        return Err(unsupported(
            "Version 3 tile metrics (EmpiricalPhasingMetricsOut.bin)",
        ));
    }
    let metrics = read_tile_metrics(&files[0])?;

    // `determineLastValueForLaneTileMetricsCode`: the last record of each (lane, tile, code).
    let mut last: JavaHashMapBy<(i32, i32, i32), TileMetricRecord> = JavaHashMapBy::new();
    for record in &metrics.records {
        last.put(
            (record.lane, record.tile, record.code),
            lane_tile_code_hash(record),
            *record,
        );
    }
    // `partitionTileMetricsByLocation`.
    let mut by_location: JavaHashMap<Vec<TileMetricRecord>> = JavaHashMap::new();
    for (_, record) in last.iter() {
        let key = format!("{}:{}", record.lane, record.tile);
        match by_location.get_mut(&key) {
            Some(list) => list.push(*record),
            None => by_location.put(&key, vec![*record]),
        }
    }

    let mut tiles: Vec<Tile> = Vec::new();
    for (location, records) in by_location.iter() {
        let mut by_code: JavaHashMapBy<i32, Vec<TileMetricRecord>> = JavaHashMapBy::new();
        for record in records {
            match by_code.get_mut(&record.code, record.code) {
                Some(list) => list.push(*record),
                None => by_code.put(record.code, record.code, vec![*record]),
            }
        }
        let codes: Vec<i32> = by_code.iter().map(|(code, _)| *code).collect();
        let value = |code: i32| -> Option<TileMetricRecord> {
            by_code
                .iter()
                .find(|(c, _)| **c == code)
                .map(|(_, list)| list[0])
        };
        let (Some(density), Some(cluster)) = (value(100), value(102)) else {
            let observed: Vec<String> = codes.iter().map(|c| c.to_string()).collect();
            return Err(Thrown::picard(format!(
                "Expected to find cluster and density record codes (102 and 100) in records read for tile location {location} (lane:tile), but found only [{}].",
                observed.join(", ")
            )));
        };
        let mut phasing = [None, None];
        let mut first = true;
        for (index, descriptor) in structure.descriptors.iter().enumerate() {
            if descriptor.kind != ReadType::T {
                continue;
            }
            let read = if first { 0 } else { 1 };
            let name = if first { "FIRST" } else { "SECOND" };
            let phasing_code = 200 + index as i32 * 2;
            let prephasing_code = 201 + index as i32 * 2;
            let pair = match (value(phasing_code), value(prephasing_code)) {
                (Some(p), Some(q)) => (p.value, q.value),
                (p, q) => {
                    let message = format!(
                        "Don't have both phasing and prephasing values for {name} read cycle {}.  Phasing code was {phasing_code} and prephasing code was {prephasing_code}.",
                        index + 1
                    );
                    if p.is_none() && q.is_none() && stringency != Stringency::Strict {
                        (0.0, 0.0)
                    } else {
                        return Err(Thrown::picard(message));
                    }
                }
            };
            // `Tile` keeps one value per read; a structure with three template reads hands it
            // two SECONDs, which `getSoleElement` refuses.
            if phasing[read].is_some() {
                return Err(Thrown::illegal_argument(
                    "Expected a single element in the phasing values of a template read",
                ));
            }
            phasing[read] = Some(pair);
            first = false;
        }
        tiles.push(Tile {
            lane: density.lane,
            density: density.value,
            clusters: cluster.value,
            phasing,
        });
    }

    // `groupingBy(Tile::getLaneNumber)`, lanes above zero.
    let mut lanes: JavaHashMapBy<i32, Vec<Tile>> = JavaHashMapBy::new();
    for tile in tiles.into_iter().filter(|t| t.lane > 0) {
        let lane = tile.lane;
        match lanes.get_mut(&lane, lane) {
            Some(list) => list.push(tile),
            None => lanes.put(lane, lane, vec![tile]),
        }
    }
    let mut out = Vec::new();
    let keys: Vec<i32> = lanes.iter().map(|(k, _)| *k).collect();
    for key in keys {
        let list = std::mem::take(lanes.get_mut(&key, key).unwrap());
        out.push((key, list));
    }
    Ok((out, version))
}
