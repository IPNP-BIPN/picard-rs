//! Conformance for `CollectRawWgsMetrics` against Picard 3.4.0.
//!
//! The tool is `CollectWgsMetrics` with four defaults changed and nothing else, so the same
//! accounting is ported once and this suite asks it the same questions under the other defaults.
//! The golden's fixture is the same reference and the same shapes of read.
//!
//! # What this suite is for
//!
//!  * **the mapping-quality floor being zero, so a read the other tool excludes whole is counted**;
//!  * **the base-quality floor being three, which quality two and quality five straddle**;
//!  * **an N base still being excluded, quality zero being under three as under twenty**;
//!  * **the coverage cap being past anything the fixture reaches**;
//!  * **and the duplicate, unpaired and overlap rules being untouched.**

use std::io::Read;

use picard_analysis::collect_wgs_metrics::{
    base_fate, genome_territory, mean_coverage, raw_arguments, read_fate, Arguments, Counts, Fate,
    Read as WgsRead, RAW_COVERAGE_CAP, RAW_LOCUS_ACCUMULATION_CAP, RAW_MINIMUM_BASE_QUALITY,
    RAW_MINIMUM_MAPPING_QUALITY,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_raw_wgs_metrics.txt.gz");
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

fn metrics(text: &str, case: &str) -> std::collections::HashMap<String, String> {
    let table = field(text, "metrics", case).unwrap_or_else(|| panic!("{case}"));
    let mut lines = table.lines().filter(|line| !line.is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let values: Vec<&str> = lines.next().expect("a value line").split('\t').collect();
    header
        .into_iter()
        .zip(values)
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

fn value(text: &str, case: &str, name: &str) -> f64 {
    metrics(text, case)
        .get(name)
        .unwrap_or_else(|| panic!("{case}/{name}"))
        .parse()
        .expect("a number")
}

fn reference() -> Vec<u8> {
    let mut bases: Vec<u8> = (0..190).map(|i| b"ACGT"[i % 4]).collect();
    bases.extend(std::iter::repeat_n(b'N', 10));
    bases
}

struct Record {
    read: WgsRead,
    start: usize,
    bases: Vec<u8>,
    qualities: Vec<i32>,
    pair: String,
}

fn records(text: &str, case: &str) -> Vec<Record> {
    field(text, "sam", case)
        .unwrap_or_else(|| panic!("{case} has an input"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("a flag word");
            Record {
                read: WgsRead {
                    adapter: false,
                    mapping_quality: columns[4].parse().expect("a mapping quality"),
                    duplicate: flags & 0x400 != 0,
                    paired: flags & 0x1 != 0,
                    mate_unmapped: flags & 0x8 != 0,
                },
                start: columns[3].parse::<usize>().expect("a start") - 1,
                bases: columns[9].bytes().collect(),
                qualities: columns[10].bytes().map(|b| i32::from(b) - 33).collect(),
                pair: columns[0].to_string(),
            }
        })
        .collect()
}

fn walk(text: &str, case: &str, arguments: &Arguments) -> Counts {
    let reference = reference();
    let mut depths = vec![0i32; reference.len()];
    let mut covered_by: Vec<Option<String>> = vec![None; reference.len()];
    let mut counts = Counts::default();
    for record in records(text, case) {
        if let Some(fate) = read_fate(&record.read, arguments) {
            for _ in 0..record.bases.len() {
                counts.add(fate);
            }
            continue;
        }
        for (offset, base) in record.bases.iter().enumerate() {
            let at = record.start + offset;
            let already = covered_by[at].as_deref() == Some(record.pair.as_str());
            let fate = base_fate(
                *base,
                record.qualities[offset],
                already,
                depths[at],
                arguments,
            );
            counts.add(fate);
            if fate == Fate::Counted {
                depths[at] += 1;
                covered_by[at] = Some(record.pair.clone());
            }
        }
    }
    counts
}

fn arguments(case: &str) -> Arguments {
    let mut arguments = raw_arguments();
    if matches!(
        case,
        "low-mapping-quality"
            | "quality-two"
            | "quality-five"
            | "n-bases"
            | "one-unpaired-read-counted"
            | "duplicate"
            | "deep"
    ) {
        arguments.count_unpaired = true;
    }
    arguments
}

const CASES: &[&str] = &[
    "low-mapping-quality",
    "quality-two",
    "quality-five",
    "n-bases",
    "one-unpaired-read",
    "one-unpaired-read-counted",
    "duplicate",
    "pair-overlapping",
    "pair-disjoint",
    "deep",
    "empty",
];

/// Every case's metrics are what the port reaches under the raw defaults.
#[test]
fn every_case_reaches_the_same_metrics() {
    let text = corpus();
    let territory = genome_territory(&reference());
    for case in CASES {
        let counts = walk(&text, case, &arguments(case));
        let round = |x: f64| (x * 1e6).round() / 1e6;
        assert_eq!(
            round(mean_coverage(counts.counted, territory)),
            value(&text, case, "MEAN_COVERAGE"),
            "{case} mean"
        );
        for (name, count) in [
            ("PCT_EXC_MAPQ", counts.mapping_quality),
            ("PCT_EXC_DUPE", counts.duplicate),
            ("PCT_EXC_UNPAIRED", counts.unpaired),
            ("PCT_EXC_BASEQ", counts.base_quality),
            ("PCT_EXC_OVERLAP", counts.overlap),
            ("PCT_EXC_CAPPED", counts.capped),
        ] {
            assert_eq!(
                round(counts.fraction(count)),
                value(&text, case, name),
                "{case} {name}"
            );
        }
    }
}

/// The four defaults are the whole of the difference.
#[test]
fn the_four_defaults_are_the_whole_difference() {
    let raw = raw_arguments();
    let plain = Arguments::default();
    assert_eq!(raw.minimum_mapping_quality, RAW_MINIMUM_MAPPING_QUALITY);
    assert_eq!(raw.minimum_base_quality, RAW_MINIMUM_BASE_QUALITY);
    assert_eq!(raw.coverage_cap, RAW_COVERAGE_CAP);
    assert_eq!(RAW_LOCUS_ACCUMULATION_CAP, 200_000);
    assert_eq!(raw.count_unpaired, plain.count_unpaired);
    assert_eq!(RAW_MINIMUM_MAPPING_QUALITY, 0);
    assert_eq!(RAW_MINIMUM_BASE_QUALITY, 3);
}

/// A read at mapping quality five is counted here, where the other tool excludes it whole.
#[test]
fn the_mapping_quality_floor_is_nothing() {
    let text = corpus();
    assert_eq!(value(&text, "low-mapping-quality", "PCT_EXC_MAPQ"), 0.0);
    assert_eq!(
        value(&text, "low-mapping-quality", "MEAN_COVERAGE"),
        0.105263
    );
    let low = WgsRead {
        adapter: false,
        mapping_quality: 5,
        duplicate: false,
        paired: false,
        mate_unmapped: false,
    };
    let counted = Arguments {
        count_unpaired: true,
        ..raw_arguments()
    };
    assert_eq!(read_fate(&low, &counted), None);
    // The other tool's defaults exclude the same read whole.
    assert_eq!(
        read_fate(
            &low,
            &Arguments {
                count_unpaired: true,
                ..Arguments::default()
            }
        ),
        Some(Fate::MappingQuality)
    );
}

/// The base-quality floor is three, which quality two and quality five straddle.
#[test]
fn the_base_quality_floor_is_three() {
    let text = corpus();
    assert_eq!(value(&text, "quality-two", "PCT_EXC_BASEQ"), 0.5);
    assert_eq!(value(&text, "quality-five", "PCT_EXC_BASEQ"), 0.0);
    assert_eq!(value(&text, "quality-five", "MEAN_COVERAGE"), 0.105263);
    let raw = raw_arguments();
    assert_eq!(base_fate(b'A', 2, false, 0, &raw), Fate::BaseQuality);
    assert_eq!(base_fate(b'A', 3, false, 0, &raw), Fate::Counted);
    assert_eq!(base_fate(b'A', 5, false, 0, &raw), Fate::Counted);
    // Where the other tool's floor excludes both.
    let plain = Arguments::default();
    assert_eq!(base_fate(b'A', 5, false, 0, &plain), Fate::BaseQuality);
}

/// An N base is still excluded, quality zero being under three as it is under twenty.
#[test]
fn an_n_base_is_still_excluded() {
    let text = corpus();
    assert_eq!(value(&text, "n-bases", "PCT_EXC_BASEQ"), 0.25);
    assert_eq!(
        base_fate(b'N', 40, false, 0, &raw_arguments()),
        Fate::BaseQuality
    );
}

/// The cap is past anything the fixture reaches, so ten reads deep is not capped.
#[test]
fn the_cap_is_out_of_reach() {
    let text = corpus();
    assert_eq!(value(&text, "deep", "PCT_EXC_CAPPED"), 0.0);
    assert_eq!(value(&text, "deep", "MEAN_COVERAGE"), 1.052632);
    assert_eq!(
        base_fate(b'A', 40, false, 10, &raw_arguments()),
        Fate::Counted
    );
    // The default cap is a hundred thousand, which the one case that leaves it alone shows in the
    // trailer counting the bins past the deepest one anything reached: 99,990 of them against the
    // 240 the same reads leave under a cap of two hundred and fifty, which every other case names
    // to keep the reference from writing a hundred thousand histogram lines apiece.
    let histogram =
        |case: &str| field(&text, "histogram", case).unwrap_or_else(|| panic!("histogram/{case}"));
    assert!(histogram("default-coverage-cap").contains("# 99990 further bins"));
    assert!(histogram("deep").contains("# 240 further bins"));
    // The two runs are the same reads, so their counted bins agree up to where the cap cuts.
    assert_eq!(
        value(&text, "default-coverage-cap", "MEAN_COVERAGE"),
        value(&text, "deep", "MEAN_COVERAGE")
    );
}

/// The duplicate, unpaired and overlap rules are untouched.
#[test]
fn the_other_three_rules_are_untouched() {
    let text = corpus();
    assert_eq!(value(&text, "duplicate", "PCT_EXC_DUPE"), 1.0);
    assert_eq!(value(&text, "one-unpaired-read", "PCT_EXC_UNPAIRED"), 1.0);
    assert_eq!(
        value(&text, "one-unpaired-read-counted", "PCT_EXC_UNPAIRED"),
        0.0
    );
    assert_eq!(value(&text, "pair-overlapping", "PCT_EXC_OVERLAP"), 0.5);
    assert_eq!(value(&text, "pair-disjoint", "PCT_EXC_OVERLAP"), 0.0);
    assert_eq!(value(&text, "empty", "GENOME_TERRITORY"), 190.0);
}
