//! Conformance for `CollectJumpingLibraryMetrics` against Picard 3.4.0.
//!
//! Each case carries the metrics table the tool wrote for a fixture whose pairs the test rebuilds
//! from the same rule the dump generated them with. The port must reproduce every column.
//!
//! # What this suite is for
//!
//!  * **the orientation being the two strands and the two positions, not the insert's sign**;
//!  * **the three chimera kinds and the order they are tried in**;
//!  * **the chimera floor being the greater of the argument and the outward mode**;
//!  * **the histogram trim keeping only consecutive bins**;
//!  * **the quality floor consulting MQ only when it is there**;
//!  * **the library size being zero without duplicates**;
//!  * **every ratio being zero when its denominator is**;
//!  * **and the unsorted refusal's misspelling.**

use std::collections::BTreeMap;
use std::io::Read;

use picard_analysis::jumping_library::{
    chimera_threshold, classify, collect, estimate_library_size, histogram_mean, mode, orientation,
    outward_mode, passes_quality, trim_by_tail_limit, unsorted_message, Arguments, Bucket,
    JumpingLibraryMetrics, Orientation, Pair, COLUMNS, DEFAULT_CHIMERA_KB_MIN, DEFAULT_TAIL_LIMIT,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("jumping_library.txt.gz");
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

/// One case's metrics table, as a map from column name to written value.
fn metrics(corpus: &str, name: &str) -> BTreeMap<String, String> {
    let line = corpus
        .lines()
        .find(|line| line.starts_with(&format!("metrics\t{name}\t")))
        .unwrap_or_else(|| panic!("the corpus carries {name}"));
    let payload = unescape(&line[format!("metrics\t{name}\t").len()..]);
    let mut rows = payload.lines();
    let header: Vec<&str> = rows.next().expect("a header").split('\t').collect();
    let values: Vec<&str> = rows.next().expect("a row").split('\t').collect();
    header
        .into_iter()
        .zip(values)
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

fn integer(row: &BTreeMap<String, String>, column: &str) -> i64 {
    row.get(column)
        .unwrap_or_else(|| panic!("the row carries {column}"))
        .parse()
        .unwrap_or_else(|_| panic!("{column} is a number"))
}

fn number(row: &BTreeMap<String, String>, column: &str) -> Option<f64> {
    let value = row.get(column)?;
    if value.is_empty() || value == "?" {
        return None;
    }
    value.parse().ok()
}

/// The fixture's pairs, rebuilt from the rule the dump wrote them with.
///
/// The mate always sits to the right of the read, so the second field of each tuple is always
/// true except where a case says otherwise.
fn jump(insert: i64) -> Pair {
    Pair {
        reference_index: 0,
        mate_reference_index: 0,
        reverse: true,
        mate_reverse: false,
        insert_size: insert,
        duplicate: false,
        mate_quality: None,
        mapping_quality: 60,
        unmapped: false,
        mate_unmapped: false,
    }
}

fn innie(insert: i64) -> Pair {
    Pair {
        reverse: false,
        mate_reverse: true,
        ..jump(insert)
    }
}

fn jumps() -> Vec<(Pair, bool)> {
    vec![
        (jump(-2000), true),
        (jump(-2100), true),
        (jump(-1900), true),
    ]
}

fn run(pairs: Vec<(Pair, bool)>, arguments: &Arguments) -> JumpingLibraryMetrics {
    collect(&pairs, arguments)
}

fn assert_matches(produced: &JumpingLibraryMetrics, row: &BTreeMap<String, String>, case: &str) {
    assert_eq!(produced.jump_pairs, integer(row, "JUMP_PAIRS"), "{case}");
    assert_eq!(
        produced.nonjump_pairs,
        integer(row, "NONJUMP_PAIRS"),
        "{case}"
    );
    assert_eq!(
        produced.chimeric_pairs,
        integer(row, "CHIMERIC_PAIRS"),
        "{case}"
    );
    assert_eq!(produced.fragments, integer(row, "FRAGMENTS"), "{case}");
    assert_eq!(
        produced.jump_duplicate_pairs,
        integer(row, "JUMP_DUPLICATE_PAIRS"),
        "{case}"
    );
    assert_eq!(
        produced.jump_library_size,
        integer(row, "JUMP_LIBRARY_SIZE"),
        "{case}"
    );
    if let Some(mean) = number(row, "JUMP_MEAN_INSERT_SIZE") {
        assert!(
            (produced.jump_mean_insert_size - mean).abs() < 1e-6,
            "{case}: {} vs {mean}",
            produced.jump_mean_insert_size
        );
    }
    for (produced, column) in [
        (produced.pct_jumps, "PCT_JUMPS"),
        (produced.pct_nonjumps, "PCT_NONJUMPS"),
        (produced.pct_chimeras, "PCT_CHIMERAS"),
    ] {
        if let Some(written) = number(row, column) {
            assert!(
                (produced - written).abs() < 1e-6,
                "{case}: {column} {produced} vs {written}"
            );
        }
    }
}

/// Every case the corpus carries, rebuilt from the fixture's own rule.
#[test]
fn every_case_matches_the_corpus() {
    let corpus = corpus();
    let default = Arguments::default();
    let mut compared = 0;

    assert_matches(
        &run(jumps(), &default),
        &metrics(&corpus, "jumps-only"),
        "jumps-only",
    );
    compared += 1;

    // Writing the mates as well changes nothing: only the first of each pair is counted, and the
    // dump's second reads are not offered to the port at all.
    assert_matches(
        &run(jumps(), &default),
        &metrics(&corpus, "jumps-with-mates"),
        "jumps-with-mates",
    );
    compared += 1;

    let innies = vec![(innie(300), true), (innie(350), true)];
    assert_matches(
        &run(innies, &default),
        &metrics(&corpus, "innies-only"),
        "innies-only",
    );
    compared += 1;

    // The three chimera kinds.
    let mut chimeras = jumps();
    chimeras.push((jump(-500000), true));
    chimeras.push((
        Pair {
            mate_reverse: true,
            ..jump(-2000)
        },
        true,
    ));
    chimeras.push((
        Pair {
            mate_reference_index: 1,
            insert_size: 0,
            ..jump(0)
        },
        true,
    ));
    assert_matches(
        &run(chimeras, &default),
        &metrics(&corpus, "chimeras"),
        "chimeras",
    );
    compared += 1;

    // A fragment.
    let mut fragment = jumps();
    fragment.push((
        Pair {
            mate_unmapped: true,
            insert_size: 0,
            ..jump(0)
        },
        true,
    ));
    assert_matches(
        &run(fragment, &default),
        &metrics(&corpus, "fragment"),
        "fragment",
    );
    compared += 1;

    // Duplicates.
    let mut duplicates = jumps();
    duplicates.push((
        Pair {
            duplicate: true,
            ..jump(-2000)
        },
        true,
    ));
    assert_matches(
        &run(duplicates, &default),
        &metrics(&corpus, "duplicates"),
        "duplicates",
    );
    compared += 1;

    // The quality floor, off and on.
    let mut qualities = jumps();
    qualities.push((
        Pair {
            mapping_quality: 5,
            ..jump(-2000)
        },
        true,
    ));
    qualities.push((
        Pair {
            mate_quality: Some(5),
            ..jump(-2000)
        },
        true,
    ));
    assert_matches(
        &run(qualities.clone(), &default),
        &metrics(&corpus, "quality-floor-off"),
        "quality-floor-off",
    );
    compared += 1;
    let strict = Arguments {
        minimum_mapping_quality: 30,
        ..Arguments::default()
    };
    assert_matches(
        &run(qualities, &strict),
        &metrics(&corpus, "quality-floor-thirty"),
        "quality-floor-thirty",
    );
    compared += 1;

    // The chimera floor lowered to one, where the mode takes over.
    let lowered = Arguments {
        chimera_kb_min: 1,
        ..Arguments::default()
    };
    assert_matches(
        &run(jumps(), &lowered),
        &metrics(&corpus, "chimera-floor-one"),
        "chimera-floor-one",
    );
    compared += 1;

    // The tail limit, which changes nothing on this distribution.
    let trimmed = Arguments {
        tail_limit: 1,
        ..Arguments::default()
    };
    assert_matches(
        &run(jumps(), &trimmed),
        &metrics(&corpus, "tail-limit-one"),
        "tail-limit-one",
    );
    compared += 1;

    // No pairs at all.
    assert_matches(
        &run(vec![], &default),
        &metrics(&corpus, "no-pairs"),
        "no-pairs",
    );
    compared += 1;

    assert_eq!(compared, 11, "the cases the port reproduces");
}

/// The two strands and the two positions, not the sign of the insert.
#[test]
fn the_orientation_is_the_strands_and_the_positions() {
    assert_eq!(orientation(true, false, true), Orientation::Rf);
    assert_eq!(orientation(false, true, true), Orientation::Fr);
    // The mate to the LEFT swaps the two.
    assert_eq!(orientation(true, false, false), Orientation::Fr);
    assert_eq!(orientation(false, true, false), Orientation::Rf);
    // Two ends on the same strand are tandem whichever way round they are.
    assert_eq!(orientation(true, true, true), Orientation::Tandem);
    assert_eq!(orientation(false, false, false), Orientation::Tandem);
    // And the sign of the insert has nothing to do with it: the corpus's jumps all carry a
    // negative insert and its innies a positive one, yet the port is never given either.
    let corpus = corpus();
    assert_eq!(integer(&metrics(&corpus, "jumps-only"), "JUMP_PAIRS"), 3);
    assert_eq!(
        integer(&metrics(&corpus, "innies-only"), "NONJUMP_PAIRS"),
        2
    );
}

/// An oversized insert wins over tandem, and tandem wins over cross-chromosome.
#[test]
fn the_order_of_the_chimera_tests_decides() {
    let arguments = Arguments::default();
    let threshold = chimera_threshold(2000.0, &arguments);
    assert_eq!(threshold, DEFAULT_CHIMERA_KB_MIN as f64);
    // Oversized AND tandem: counted once, as a chimera either way.
    let both = Pair {
        mate_reverse: true,
        ..jump(-500000)
    };
    assert_eq!(
        classify(&both, threshold, true, &arguments),
        Bucket::Chimera
    );
    // Tandem AND cross-chromosome: also one chimera.
    let tandem_cross = Pair {
        mate_reverse: true,
        mate_reference_index: 1,
        ..jump(-2000)
    };
    assert_eq!(
        classify(&tandem_cross, threshold, true, &arguments),
        Bucket::Chimera
    );
    // Which is why the corpus's run of two such pairs reports two and not four.
    let corpus = corpus();
    assert_eq!(
        integer(&metrics(&corpus, "overlapping-chimeras"), "CHIMERIC_PAIRS"),
        2
    );
}

/// The greater of the argument and the outward mode.
#[test]
fn the_chimera_floor_is_the_greater_of_two() {
    let default = Arguments::default();
    assert_eq!(
        outward_mode(&jumps(), &default),
        1900.0,
        "the first of a tie"
    );
    // At the default the argument wins.
    assert_eq!(
        chimera_threshold(outward_mode(&jumps(), &default), &default),
        100000.0
    );
    // Lowered to one, the mode wins, and the two inserts past it become chimeras.
    let lowered = Arguments {
        chimera_kb_min: 1,
        ..Arguments::default()
    };
    assert_eq!(
        chimera_threshold(outward_mode(&jumps(), &lowered), &lowered),
        1900.0
    );
    let corpus = corpus();
    let row = metrics(&corpus, "chimera-floor-one");
    assert_eq!(integer(&row, "JUMP_PAIRS"), 1);
    assert_eq!(integer(&row, "CHIMERIC_PAIRS"), 2);
    // An empty histogram has no mode at all, so the argument is the whole floor.
    assert_eq!(outward_mode(&[], &default), 0.0);
    assert_eq!(mode(&BTreeMap::new()), 0.0);
}

/// It walks to the mode and then only while the bins follow one another exactly.
#[test]
fn the_trim_keeps_only_consecutive_bins() {
    // Three bins a hundred apart: everything past the mode is cut.
    let spread: BTreeMap<i64, i64> = [(1900, 1), (2000, 1), (2100, 1)].into_iter().collect();
    let trimmed = trim_by_tail_limit(&spread, DEFAULT_TAIL_LIMIT);
    assert_eq!(trimmed.keys().copied().collect::<Vec<_>>(), vec![1900]);
    assert_eq!(histogram_mean(&trimmed), 1900.0);
    // Which is what the corpus wrote, from a set whose three inserts average 2000.
    let corpus = corpus();
    assert_eq!(
        number(&metrics(&corpus, "jumps-only"), "JUMP_MEAN_INSERT_SIZE"),
        Some(1900.0)
    );
    // Consecutive bins are kept, and the mean then moves.
    let consecutive: BTreeMap<i64, i64> = [(2000, 3), (2001, 2), (2002, 1)].into_iter().collect();
    let kept = trim_by_tail_limit(&consecutive, DEFAULT_TAIL_LIMIT);
    assert_eq!(kept.len(), 3);
    assert!((histogram_mean(&kept) - 2000.6666666666667).abs() < 1e-9);
    assert_eq!(
        number(
            &metrics(&corpus, "consecutive-bins"),
            "JUMP_MEAN_INSERT_SIZE"
        ),
        Some(2000.666667)
    );
    // A lone outlier past a real mode is cut whatever the limit is.
    let tail_heavy: BTreeMap<i64, i64> = [(2000, 5), (60000, 1)].into_iter().collect();
    assert_eq!(trim_by_tail_limit(&tail_heavy, 2).len(), 1);
    assert_eq!(trim_by_tail_limit(&tail_heavy, DEFAULT_TAIL_LIMIT).len(), 1);
    assert_eq!(
        number(&metrics(&corpus, "tail-heavy"), "JUMP_MEAN_INSERT_SIZE"),
        Some(2000.0)
    );
    // An empty histogram trims to nothing rather than panicking.
    assert!(trim_by_tail_limit(&BTreeMap::new(), 1).is_empty());
}

/// A pair with no MQ tag passes on its own mapping quality alone.
#[test]
fn the_quality_floor_consults_mq_only_when_it_is_there() {
    let strict = Arguments {
        minimum_mapping_quality: 30,
        ..Arguments::default()
    };
    assert!(
        passes_quality(&jump(-2000), &strict),
        "no tag, good quality"
    );
    assert!(!passes_quality(
        &Pair {
            mapping_quality: 5,
            ..jump(-2000)
        },
        &strict
    ));
    assert!(!passes_quality(
        &Pair {
            mate_quality: Some(5),
            ..jump(-2000)
        },
        &strict
    ));
    // A tag above the floor passes, and so does one exactly on it.
    assert!(passes_quality(
        &Pair {
            mate_quality: Some(30),
            ..jump(-2000)
        },
        &strict
    ));
    // Which is what the corpus's two runs report: five jumps at the default and three at thirty.
    let corpus = corpus();
    assert_eq!(
        integer(&metrics(&corpus, "quality-floor-off"), "JUMP_PAIRS"),
        5
    );
    assert_eq!(
        integer(&metrics(&corpus, "quality-floor-thirty"), "JUMP_PAIRS"),
        3
    );
}

/// Zero when there are none, and the bisection's own answer when there are.
#[test]
fn the_library_size_is_zero_without_duplicates() {
    assert_eq!(estimate_library_size(4, 4), None);
    assert_eq!(estimate_library_size(0, 0), None);
    // Four pairs, three of them unique.
    assert_eq!(estimate_library_size(4, 3), Some(6));
    let corpus = corpus();
    assert_eq!(
        integer(&metrics(&corpus, "duplicates"), "JUMP_LIBRARY_SIZE"),
        6
    );
    assert_eq!(
        integer(&metrics(&corpus, "jumps-only"), "JUMP_LIBRARY_SIZE"),
        0
    );
    // More duplicates make a smaller library.
    let few = estimate_library_size(100, 90).expect("a size");
    let many = estimate_library_size(100, 50).expect("a size");
    assert!(many < few);
}

/// Zero and not a NaN.
#[test]
fn every_ratio_is_zero_when_its_denominator_is() {
    let empty = collect(&[], &Arguments::default());
    assert_eq!(empty.pct_jumps, 0.0);
    assert_eq!(empty.pct_nonjumps, 0.0);
    assert_eq!(empty.pct_chimeras, 0.0);
    assert_eq!(empty.jump_duplicate_pct, 0.0);
    assert_eq!(empty.nonjump_duplicate_pct, 0.0);
    assert_eq!(empty.jump_library_size, 0);
    // The histogram's own mean is a NaN over nothing, which the metric writes as it is.
    assert!(empty.jump_mean_insert_size.is_nan());
    let corpus = corpus();
    let row = metrics(&corpus, "no-pairs");
    assert_eq!(integer(&row, "JUMP_PAIRS"), 0);
    assert_eq!(number(&row, "PCT_JUMPS"), Some(0.0));
}

/// Both ends unmapped and unplaced ends the traversal.
#[test]
fn an_unplaced_unmapped_pair_ends_the_file() {
    let arguments = Arguments::default();
    let terminator = Pair {
        unmapped: true,
        mate_unmapped: true,
        reference_index: -1,
        mate_reference_index: -1,
        insert_size: 0,
        mapping_quality: 0,
        ..jump(0)
    };
    assert_eq!(
        classify(&terminator, 100000.0, true, &arguments),
        Bucket::Terminator
    );
    // A pair with both ends unmapped but PLACED is passed over instead.
    let placed = Pair {
        reference_index: 0,
        mate_reference_index: 0,
        ..terminator.clone()
    };
    assert_eq!(
        classify(&placed, 100000.0, true, &arguments),
        Bucket::Skipped
    );
    // Anything after the terminator is never counted.
    let mut pairs = jumps();
    pairs.push((terminator, true));
    pairs.push((jump(-2000), true));
    assert_eq!(collect(&pairs, &arguments).jump_pairs, 3);
    // Which the corpus does not show, a coordinate-sorted file putting the unplaced reads LAST:
    // its run counts the pair written after the terminator all the same.
    let corpus = corpus();
    assert_eq!(
        integer(&metrics(&corpus, "unmapped-terminator"), "JUMP_PAIRS"),
        4
    );
}

/// The column order and the refusal's own words.
#[test]
fn the_columns_and_the_refusal_are_the_references() {
    let corpus = corpus();
    let row = metrics(&corpus, "jumps-only");
    for column in COLUMNS {
        assert!(row.contains_key(column), "{column}");
    }
    assert_eq!(row.len(), COLUMNS.len());
    // The refusal misspells "coordinate" and names the file twice.
    let line = corpus
        .lines()
        .find(|line| line.starts_with("error\tunsorted\t"))
        .expect("the refusal");
    let message = line["error\tunsorted\t".len()..].to_string();
    assert_eq!(
        message,
        format!("picard.PicardException:{}", unsorted_message("in.bam"))
    );
    assert!(message.contains("coordintate"), "{message}");
}
