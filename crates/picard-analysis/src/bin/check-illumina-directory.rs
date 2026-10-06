//! `CheckIlluminaDirectory` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.CheckIlluminaDirectory.doWork` and `verifyLane` at tag 3.4.0 for
//! per-tile runs, over `picard_analysis::illumina_dir`. The tool's answer is its log and its exit
//! code: each data type is matched to the first format whose files are there, each format checks
//! its files for the tiles `TileMetricsOut.bin` lists and the cycles the read structure keeps, and
//! every failure is one log line.
//!
//! The reference also writes the failure count to `./errors.count` in its working directory, which
//! is outside what a run is compared on; the port leaves that file unwritten. `FAKE_FILES` and
//! `LINK_LOCS` write into the basecalls directory and are not ported.

use picard_analysis::illumina_dir::{
    expected_tiles, log_info, long_lane, Layout, PerTile, PerTilePerCycle,
};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::metrics_cli::{thrown, Args};

const CLASS: &str = "CheckIlluminaDirectory";

/// `IlluminaDataType`, in enum order, which is the order a `TreeSet` of them walks.
const DATA_TYPES: [&str; 5] = ["Position", "BaseCalls", "QualityScores", "PF", "Barcodes"];

fn main() {
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("RS", "READ_STRUCTURE"),
        ("L", "LANES"),
        ("T", "TILE_NUMBERS"),
        ("DT", "DATA_TYPES"),
    ]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let structure = parse_read_structure(&args.required("READ_STRUCTURE"))
        .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"));
    let lanes: Vec<i32> = args
        .collection("LANES", &[])
        .iter()
        .filter_map(|v| v.parse().ok())
        .collect();
    let tile_numbers: Vec<i32> = args
        .collection("TILE_NUMBERS", &[])
        .iter()
        .filter_map(|v| v.parse().ok())
        .collect();
    let mut types: Vec<String> = args.collection("DATA_TYPES", &[]);
    if types.is_empty() {
        types = ["BaseCalls", "QualityScores", "Position", "PF"]
            .iter()
            .map(|s| s.to_string())
            .collect();
    }
    let types: Vec<&str> = DATA_TYPES
        .iter()
        .copied()
        .filter(|t| types.iter().any(|w| w == t))
        .collect();

    // `OutputMapping.getOutputCycles`: every cycle not in a skip.
    let mut cycles: Vec<i32> = Vec::new();
    let mut cycle = 1;
    for segment in &structure {
        for _ in 0..segment.cycles {
            if segment.kind != SegmentKind::Skip {
                cycles.push(cycle);
            }
            cycle += 1;
        }
    }
    let join = |v: &[i32], sep: &str| {
        v.iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(sep)
    };
    let absolute = |p: &std::path::Path| {
        picard_analysis::fingerprint::reference_path(&p.display().to_string())
    };
    log_info(
        CLASS,
        &format!(
            "Checking lanes({} in basecalls directory ({})\n",
            join(&lanes, ","),
            absolute(&basecalls)
        ),
    );
    log_info(CLASS, &format!("Expected cycles: {}", join(&cycles, ", ")));
    let layout = Layout::new(&basecalls);
    let mut failing: Vec<i32> = Vec::new();
    let mut total = 0usize;
    for &lane in &lanes {
        let lane_dir = basecalls.join(long_lane(lane));
        if !lane_dir.is_dir() {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Directory does not exist: {}",
                absolute(&lane_dir)
            ));
        }
        let mut tiles = expected_tiles(&layout, lane).unwrap_or_else(|e| thrown(&e));
        if !tile_numbers.is_empty() {
            tiles.retain(|t| tile_numbers.contains(t));
        }
        log_info(CLASS, &format!("Checking lane {lane}"));
        log_info(CLASS, &format!("Expected tiles: {}", join(&tiles, ", ")));
        if tiles.is_empty() {
            thrown("picard.PicardException: 0 input tiles were specified!  Check to make sure this lane is in the InterOp file!");
        }
        if cycles.is_empty() {
            thrown("picard.PicardException: 0 output cycles were specified!");
        }
        // The formats, in `SupportedIlluminaFormat` order; the multi-tile and CBCL ones are
        // never available without their tile index or compressed cycle files.
        let bcl = PerTilePerCycle::new(".bcl", &lane_dir, lane);
        let intensity_lane = layout.intensities.join(long_lane(lane));
        let mut locs = PerTile::or_per_run(".locs", &intensity_lane, lane);
        let mut clocs = PerTile::new(".clocs", &intensity_lane, lane);
        let mut pos = PerTile::new("_pos.txt", &layout.intensities, lane);
        let mut filter = PerTile::new(".filter", &lane_dir, lane);
        let mut barcode = PerTile::new("_barcode.txt", &basecalls, lane);
        let format_of = |t: &str| -> Option<&'static str> {
            match t {
                "BaseCalls" | "QualityScores" => bcl.files_available().then_some("Bcl"),
                "PF" => filter.files_available().then_some("Filter"),
                "Position" => {
                    if locs.files_available() {
                        Some("Locs")
                    } else if clocs.files_available() {
                        Some("Clocs")
                    } else if pos.files_available() {
                        Some("Pos")
                    } else {
                        None
                    }
                }
                _ => barcode.files_available().then_some("Barcode"),
            }
        };
        let chosen: Vec<(&str, Option<&str>)> = types.iter().map(|t| (*t, format_of(t))).collect();
        let mut failures = 0usize;
        let unmatched: Vec<&str> = chosen
            .iter()
            .filter(|(_, f)| f.is_none())
            .map(|(t, _)| *t)
            .collect();
        if !unmatched.is_empty() {
            log_info(
                CLASS,
                &format!(
                    "Could not find a format with available files for the following data types: {}",
                    unmatched.join(", ")
                ),
            );
            failures += unmatched.len();
        }
        for format in ["Bcl", "Locs", "Clocs", "Pos", "Filter", "Barcode"] {
            if !chosen.iter().any(|(_, f)| *f == Some(format)) {
                continue;
            }
            let found = match format {
                "Bcl" => bcl.verify(&tiles, &cycles),
                "Locs" => locs.verify(&tiles),
                "Clocs" => clocs.verify(&tiles),
                "Pos" => pos.verify(&tiles),
                "Filter" => filter.verify(&tiles),
                _ => barcode.verify(&tiles),
            };
            failures += found.len();
            for f in &found {
                log_info(CLASS, f);
            }
        }
        if failures > 0 {
            log_info(
                CLASS,
                &format!("Lane {lane} FAILED  Total Errors: {failures}"),
            );
            failing.push(lane);
            total += failures;
        } else {
            log_info(CLASS, &format!("Lane {lane} SUCCEEDED "));
        }
    }
    if total == 0 {
        log_info(
            CLASS,
            "SUCCEEDED!  All required files are present and non-empty.",
        );
    } else {
        log_info(
            CLASS,
            &format!(
                "FAILED! There were {total} in the following lanes: {}",
                join(&failing, ", ")
            ),
        );
        std::process::exit(1);
    }
}
