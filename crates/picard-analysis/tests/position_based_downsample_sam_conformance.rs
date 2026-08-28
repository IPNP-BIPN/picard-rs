//! Conformance for `PositionBasedDownsampleSam` against Picard 3.4.0.
//!
//! Each case carries the file the tool read, as SAM without its header, and the names of the
//! reads it kept. The read names are Illumina's, so the port reads the tile and the two
//! coordinates straight out of them and must keep the same reads.
//!
//! # What this suite is for
//!
//!  * **the selection carrying no randomness, so two runs keep the same reads**;
//!  * **a fraction over a half inverting the mask rather than growing the circle**;
//!  * **each tile being masked by its own extent**;
//!  * **the extent starting at zero and then being widened by its span over the read count**;
//!  * **`STOP_AFTER` cutting both passes, so it moves the extent as well as the count**;
//!  * **the duplicate flag being cleared on the reads that are kept, unless told otherwise**;
//!  * **a second run being refused unless it is allowed**;
//!  * **and a fraction outside the unit interval being an exit code of one.**

use std::io::Read;

use picard_analysis::position_based_downsample_sam::{
    keep, misses_the_fraction, tile_extents, CircleSelector, Coord, PhysicalLocation,
    ACCEPTABLE_FUDGE_FACTOR, PG_PROGRAM_NAME,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("position_based_downsample_sam.txt.gz");
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

/// One case's read names, in file order.
fn names(text: &str, case: &str) -> Vec<String> {
    field(text, "sam", case)
        .unwrap_or_else(|| panic!("{case} has an input"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').next().expect("a name").to_string())
        .collect()
}

/// The names the tool kept, in the order it wrote them.
fn kept_names(text: &str, case: &str) -> Vec<String> {
    field(text, "kept", case)
        .unwrap_or_else(|| panic!("{case} kept something"))
        .split(',')
        .filter(|name| !name.is_empty())
        .map(str::to_string)
        .collect()
}

fn kept_flags(text: &str, case: &str) -> Vec<i32> {
    field(text, "flags", case)
        .unwrap_or_else(|| panic!("{case} has flags"))
        .split(',')
        .filter(|flag| !flag.is_empty())
        .map(|flag| flag.parse().expect("a flag"))
        .collect()
}

/// `RUN:1:<tile>:<x>:<y>`, which the default read-name regex reads the last three fields of.
fn location(name: &str) -> PhysicalLocation {
    let fields: Vec<&str> = name.split(':').collect();
    if fields.len() < 5 {
        // The reference's parser gives an unparseable name a tile of -1 and no coordinates.
        return PhysicalLocation {
            tile: -1,
            x: 0,
            y: 0,
        };
    }
    PhysicalLocation {
        tile: fields[2].parse().expect("a tile"),
        x: fields[3].parse().expect("an x"),
        y: fields[4].parse().expect("a y"),
    }
}

fn ours(text: &str, case: &str, fraction: f64, stop_after: Option<usize>) -> Vec<String> {
    let names = names(text, case);
    let locations: Vec<PhysicalLocation> = names.iter().map(|name| location(name)).collect();
    let kept = keep(&locations, fraction, stop_after);
    names
        .into_iter()
        .zip(kept)
        .filter(|(_, keep)| *keep)
        .map(|(name, _)| name)
        .collect()
}

/// The cases whose fraction and STOP_AFTER the dump names on the command line.
const CASES: &[(&str, f64, Option<usize>)] = &[
    ("fraction-one-tenth", 0.1, None),
    ("fraction-one-half", 0.5, None),
    ("fraction-nine-tenths", 0.9, None),
    ("fraction-one", 1.0, None),
    ("fraction-zero", 0.0, None),
    ("repeatable", 0.3, None),
    ("repeatable-again", 0.3, None),
    ("two-tiles", 0.3, None),
    ("duplicates-cleared", 0.3, None),
    ("duplicates-kept", 0.3, None),
    ("stop-after-ten", 0.3, Some(10)),
    ("unparseable-name", 0.3, None),
];

/// Every case keeps exactly the reads the port selects.
#[test]
fn every_case_keeps_the_same_reads() {
    let text = corpus();
    for (case, fraction, stop_after) in CASES {
        assert_eq!(
            ours(&text, case, *fraction, *stop_after),
            kept_names(&text, case),
            "{case}"
        );
    }
}

/// The selection carries no randomness: the same input at the same fraction keeps the same reads.
#[test]
fn the_selection_carries_no_randomness() {
    let text = corpus();
    let once = kept_names(&text, "repeatable");
    let again = kept_names(&text, "repeatable-again");
    assert_eq!(once, again);
    assert!(!once.is_empty());
}

/// A fraction over a half inverts the mask rather than growing the circle: the selector is built
/// for `1 - fraction` and its sense is reversed.
#[test]
fn a_fraction_over_a_half_inverts_the_mask() {
    let text = corpus();
    let low = CircleSelector::new(0.1);
    let high = CircleSelector::new(0.9);
    assert!(low.positive_selection);
    assert!(!high.positive_selection);
    // 1 - 0.9 is not 0.1 to the last bit, so the two radii agree to within a rounding.
    assert!((low.radius_squared - high.radius_squared).abs() < 1e-15);
    assert!((low.offset - high.offset).abs() < 1e-15);
    // Which is why the two keep complementary sets of the same grid.
    let tenth: std::collections::HashSet<String> = kept_names(&text, "fraction-one-tenth")
        .into_iter()
        .collect();
    let ninth: std::collections::HashSet<String> = kept_names(&text, "fraction-nine-tenths")
        .into_iter()
        .collect();
    let all: std::collections::HashSet<String> =
        names(&text, "fraction-one-tenth").into_iter().collect();
    assert_eq!(tenth.union(&ninth).count(), all.len());
    assert!(tenth.is_disjoint(&ninth));
}

/// A fraction of one keeps everything and a fraction of nought keeps nothing.
#[test]
fn the_two_extremes_keep_all_and_none() {
    let text = corpus();
    assert_eq!(
        kept_names(&text, "fraction-one").len(),
        names(&text, "fraction-one").len()
    );
    assert!(kept_names(&text, "fraction-zero").is_empty());
}

/// Each tile is masked by its own extent, so the same coordinates are kept on one and dropped on
/// another.
#[test]
fn each_tile_is_masked_by_its_own_extent() {
    let text = corpus();
    let kept = kept_names(&text, "two-tiles");
    let on = |tile: &str| -> Vec<String> {
        kept.iter()
            .filter(|name| name.starts_with(&format!("RUN:1:{tile}:")))
            .map(|name| name[format!("RUN:1:{tile}:").len()..].to_string())
            .collect()
    };
    let first = on("1101");
    let second = on("1102");
    assert_ne!(first, second);
    // The second tile has an extra read far out, which widens it and keeps more of the grid.
    assert!(second.len() > first.len(), "{first:?} {second:?}");
}

/// The extent starts at zero and is then widened by its own span over the read count.
#[test]
fn the_extent_starts_at_zero_and_is_widened() {
    let locations: Vec<PhysicalLocation> = (1..=5)
        .map(|i| PhysicalLocation {
            tile: 1,
            x: i * 1000,
            y: i * 1000,
        })
        .collect();
    let tiles = tile_extents(&locations);
    let coord = tiles[&1];
    // Nothing sits below 1000, and the minimum is still nought before the widening.
    assert_eq!(coord.count, 5);
    // The span was 5000 and the count 5, so each side moved by 1000.
    assert_eq!(coord.min_x, -1000);
    assert_eq!(coord.max_x, 6000);
    assert_eq!(coord.min_y, -1000);
    assert_eq!(coord.max_y, 6000);
    // A tile of ONE read still has a span, because the minimum starts at zero: its extent runs
    // from 0 to 1000 before the widening, and the widening then doubles it in both directions.
    let single = tile_extents(&locations[..1]);
    assert_eq!(
        single[&1],
        Coord {
            min_x: -1000,
            min_y: -1000,
            max_x: 2000,
            max_y: 2000,
            count: 1
        }
    );
}

/// STOP_AFTER cuts both passes, so it moves the tile's extent as well as the read count: the
/// first ten reads of the grid keep two, where the whole grid keeps four of those same ten.
#[test]
fn stop_after_cuts_both_passes() {
    let text = corpus();
    let stopped = kept_names(&text, "stop-after-ten");
    assert_eq!(stopped.len(), 2);
    let whole = kept_names(&text, "repeatable");
    let first_ten: Vec<String> = names(&text, "repeatable").into_iter().take(10).collect();
    let whole_of_first_ten: Vec<&String> = whole
        .iter()
        .filter(|name| first_ten.contains(name))
        .collect();
    assert_eq!(whole_of_first_ten.len(), 4);
}

/// The duplicate flag is cleared on the reads that are kept, unless the argument says otherwise.
#[test]
fn the_duplicate_flag_is_cleared_unless_kept() {
    let text = corpus();
    assert_eq!(
        kept_names(&text, "duplicates-cleared"),
        kept_names(&text, "duplicates-kept")
    );
    assert!(kept_flags(&text, "duplicates-cleared")
        .iter()
        .all(|flag| flag & 0x400 == 0));
    assert!(kept_flags(&text, "duplicates-kept")
        .iter()
        .all(|flag| flag & 0x400 != 0));
}

/// A read name the regex cannot parse is not refused: it is masked with the other unparseable
/// ones, and here it survives.
#[test]
fn an_unparseable_name_is_not_refused() {
    let text = corpus();
    let kept = kept_names(&text, "unparseable-name");
    assert!(kept.contains(&"no-colons-at-all".to_string()), "{kept:?}");
    assert_eq!(location("no-colons-at-all").tile, -1);
}

/// A second run over an already-downsampled file is refused by a message naming the program, and
/// allowed when the argument says so.
#[test]
fn a_second_run_is_refused_unless_allowed() {
    let text = corpus();
    let error = field(&text, "error", "downsampled-twice").expect("the refusal");
    assert!(error.contains("has been downsampled already"), "{error}");
    assert!(!kept_names(&text, "downsampled-twice-allowed").is_empty());
    // The guard looks for the tool's own name in the @PG records.
    let programs = field(&text, "pg", "downsampled-twice-allowed").expect("the programs");
    assert!(programs.contains(PG_PROGRAM_NAME), "{programs}");
}

/// A fraction outside the unit interval is an exit code of one rather than an exception.
#[test]
fn a_fraction_out_of_range_is_an_exit_code() {
    let text = corpus();
    for case in ["fraction-too-big", "fraction-negative"] {
        assert_eq!(
            field(&text, "error", case).as_deref(),
            Some("exit 1"),
            "{case}"
        );
    }
}

/// The closing check warns when the rate misses the fraction by more than a fifth of the smaller
/// of the two, which a positional mask does easily on a small file.
#[test]
fn the_closing_check_measures_the_miss() {
    assert_eq!(ACCEPTABLE_FUDGE_FACTOR, 0.2);
    assert!(misses_the_fraction(1, 25, 0.1));
    assert!(!misses_the_fraction(13, 25, 0.5));
    // Nothing read at all is not a miss.
    assert!(!misses_the_fraction(0, 0, 0.5));
}

/// The repeating part of the mask, which is what makes it tile the flowcell rather than sit once
/// in the middle of it.
#[test]
fn the_rounded_part_is_signed() {
    assert_eq!(CircleSelector::rounded_part(0.2), 0.2);
    assert_eq!(CircleSelector::rounded_part(0.0), 0.0);
    assert!((CircleSelector::rounded_part(0.8) + 0.2).abs() < 1e-12);
    assert!((CircleSelector::rounded_part(1.2) - 0.2).abs() < 1e-12);
    // Exactly halfway rounds UP, which is what Math.round does.
    assert_eq!(CircleSelector::rounded_part(0.5), -0.5);
}
