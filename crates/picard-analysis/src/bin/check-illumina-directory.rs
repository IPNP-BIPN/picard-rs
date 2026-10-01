//! `CheckIlluminaDirectory` as a runnable binary: the covering array's port side.
//!
//! The tool's answer is its exit status and, when anything failed, the NUMBER of failures written
//! to `./errors.count` in the working directory; every failure message is logged and nothing else
//! is written. The count is what [`picard_analysis::illumina_dir`] decides: which formats the
//! requested data types resolve to, and what each format's `verify` finds missing. A tile list
//! that ends up empty, a lane directory that is not there and an unreadable tile metrics file are
//! exceptions rather than counts.
//!
//! `FAKE_FILES` and `LINK_LOCS` write into the run directory and are refused here rather than
//! ported; the CBCL branch is refused too (see `illumina_dir`).
//!
//! Ported from `picard.illumina.CheckIlluminaDirectory` in Picard 3.4.0.

use std::collections::BTreeSet;
use std::path::PathBuf;

use picard_analysis::illumina_dir::{
    assert_directory_is_readable, determine_formats, has_cbcls, long_lane_str, unmatched_types,
    unsupported, DataType, IlluminaFileUtil, OutputMapping, ReadStructure, Thrown, Util,
};
use picard_analysis::metrics_cli::{fail, refuse_validation, thrown, Args};

fn main() {
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("DT", "DATA_TYPES"),
        ("RS", "READ_STRUCTURE"),
        ("L", "LANES"),
        ("T", "TILE_NUMBERS"),
        ("F", "FAKE_FILES"),
        ("X", "LINK_LOCS"),
    ]);
    let basecalls = PathBuf::from(args.required("BASECALLS_DIR"));
    let read_structure = args.required("READ_STRUCTURE");
    let lanes: Vec<i32> = args
        .all("LANES")
        .iter()
        .map(|v| {
            v.parse()
                .unwrap_or_else(|_| fail(&format!("Argument 'LANES' cannot be set to '{v}'")))
        })
        .collect();
    let tile_numbers: Vec<i32> = args
        .all("TILE_NUMBERS")
        .iter()
        .map(|v| {
            v.parse().unwrap_or_else(|_| {
                fail(&format!("Argument 'TILE_NUMBERS' cannot be set to '{v}'"))
            })
        })
        .collect();
    let mut data_types: BTreeSet<DataType> = BTreeSet::new();
    for value in args.all("DATA_TYPES") {
        match DataType::parse(&value) {
            Some(t) => {
                data_types.insert(t);
            }
            None => fail(&format!("Argument 'DATA_TYPES' cannot be set to '{value}'")),
        }
    }
    if args.bool("FAKE_FILES", false) || args.bool("LINK_LOCS", false) {
        thrown(&unsupported("FAKE_FILES and LINK_LOCS (writing into the run directory)").render());
    }

    // `customCommandLineValidation`.
    if let Err(e) = assert_directory_is_readable(&basecalls) {
        thrown(&e.render());
    }
    if lanes.iter().any(|lane| *lane < 1) {
        let joined: Vec<String> = lanes.iter().map(|l| l.to_string()).collect();
        refuse_validation(
            "CheckIlluminaDirectory",
            &[format!(
                "LANES must be greater than or equal to 1.  LANES passed in {}",
                joined.join(", ")
            )],
        );
    }

    match run(
        &basecalls,
        &read_structure,
        &lanes,
        &tile_numbers,
        data_types,
    ) {
        Ok(0) => {}
        Ok(failures) => {
            if let Err(e) = std::fs::write("./errors.count", failures.to_string()) {
                eprintln!("Unable to write number of errors to file: {e}");
            }
            std::process::exit(1);
        }
        Err(e) => thrown(&e.render()),
    }
}

/// `doWork`: the total number of failures over the lanes.
fn run(
    basecalls: &std::path::Path,
    read_structure: &str,
    lanes: &[i32],
    tile_numbers: &[i32],
    mut data_types: BTreeSet<DataType>,
) -> Result<usize, Thrown> {
    let structure = ReadStructure::parse(read_structure)?;
    if data_types.is_empty() {
        data_types = [
            DataType::BaseCalls,
            DataType::QualityScores,
            DataType::Position,
            DataType::PF,
        ]
        .into_iter()
        .collect();
    }
    let mapping = OutputMapping::new(&structure)?;
    let cycles = mapping.output_cycles.clone();
    let mut total = 0;
    for lane in lanes {
        assert_directory_is_readable(&basecalls.join(long_lane_str(*lane)))?;
        if has_cbcls(basecalls, *lane) {
            return Err(unsupported("A CBCL run directory"));
        }
        let mut util = IlluminaFileUtil::new(basecalls, None, *lane);
        let mut expected = util.expected_tiles()?;
        if !tile_numbers.is_empty() {
            expected.retain(|tile| tile_numbers.contains(tile));
        }
        total += verify_lane(&mut util, &expected, &cycles, &data_types)?;
    }
    Ok(total)
}

/// `verifyLane`.
fn verify_lane(
    util: &mut IlluminaFileUtil,
    expected_tiles: &[i32],
    cycles: &[i32],
    data_types: &BTreeSet<DataType>,
) -> Result<usize, Thrown> {
    if expected_tiles.is_empty() {
        return Err(Thrown::picard(
            "0 input tiles were specified!  Check to make sure this lane is in the InterOp file!",
        ));
    }
    if cycles.is_empty() {
        return Err(Thrown::picard("0 output cycles were specified!"));
    }
    let formats = determine_formats(data_types, util)?;
    let mut failures = unmatched_types(data_types, &formats).len();
    for format in formats.keys() {
        let verified = match util.util(*format)? {
            Util::PerTile(per_tile) => {
                per_tile.set_tiles_for_per_run_file(expected_tiles);
                per_tile.verify(expected_tiles)
            }
            Util::PerTilePerCycle(per_cycle) => per_cycle.verify(expected_tiles, cycles),
            Util::Cbcl(_) => {
                return Err(Thrown::new(
                    "java.lang.UnsupportedOperationException",
                    "`verify()` is not implemented for CBCLs",
                ))
            }
            Util::Absent => return Err(unsupported(&format!("{format:?}"))),
        };
        failures += verified.len();
    }
    Ok(failures)
}
