//! Conformance for `CollectWgsMetrics` against Picard 3.4.0.
//!
//! Each case carries the file the tool read, as SAM without its header, and the metrics table it
//! wrote. The fixture's reference is two hundred bases whose last ten are Ns, so the territory is
//! known by arithmetic, and every read is twenty bases long, so the port can replay the chain
//! read by read and reach the same eight counters.
//!
//! # What this suite is for
//!
//!  * **the territory being the non-N bases and not the covered ones**;
//!  * **the seven exclusions partitioning the bases that did not count**;
//!  * **the mapping quality and the duplicate flag taking whole reads**;
//!  * **the base quality taking single bases, and an N being one of them**;
//!  * **an unpaired read and a half-mapped pair reaching the same counter**;
//!  * **a pair's overlap counting once**;
//!  * **the cap truncating the depth and counting the remainder**;
//!  * **the mean being over the territory**;
//!  * **and the histogram argument adding a column rather than a table.**

use std::io::Read;

use picard_analysis::collect_wgs_metrics::{
    base_fate, genome_territory, mean_coverage, read_fate, Arguments, Counts, Fate, Read as WgsRead,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_wgs_metrics.txt.gz");
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

/// One case's metrics, as a name-to-value map.
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

/// The reference the dump wrote: 190 bases of `ACGT` then ten Ns.
fn reference() -> Vec<u8> {
    let mut bases: Vec<u8> = (0..190).map(|i| b"ACGT"[i % 4]).collect();
    bases.extend(std::iter::repeat_n(b'N', 10));
    bases
}

/// One record of a case, as the chain reads it.
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

/// The whole chain over one case, which is what the tool's two passes amount to here.
fn walk(text: &str, case: &str, arguments: &Arguments) -> (Counts, Vec<i32>) {
    let reference = reference();
    let mut depths = vec![0i32; reference.len()];
    // Which pair has already covered each base, so the second end's overlap is seen.
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
    (counts, depths)
}

/// The arguments each case ran under, which the dump names on the command line.
fn arguments(case: &str) -> Arguments {
    let mut arguments = Arguments::default();
    if matches!(
        case,
        "one-unpaired-read-counted"
            | "low-mapping-quality"
            | "low-base-quality"
            | "n-bases"
            | "duplicate"
            | "deep-uncapped"
            | "deep-capped"
    ) {
        arguments.count_unpaired = true;
    }
    if case == "deep-capped" {
        arguments.coverage_cap = 2;
    }
    arguments
}

const CASES: &[&str] = &[
    "one-unpaired-read",
    "one-unpaired-read-counted",
    "pair-disjoint",
    "pair-overlapping",
    "low-mapping-quality",
    "low-base-quality",
    "n-bases",
    "duplicate",
    "mate-unmapped",
    "deep-uncapped",
    "deep-capped",
    "with-histogram",
    "empty",
];

/// Every case's exclusion fractions and mean coverage are what the port reaches.
#[test]
fn every_case_reaches_the_same_metrics() {
    let text = corpus();
    let territory = genome_territory(&reference());
    for case in CASES {
        let (counts, _) = walk(&text, case, &arguments(case));
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
        assert_eq!(
            round(counts.excluded_fraction()),
            value(&text, case, "PCT_EXC_TOTAL"),
            "{case} total"
        );
    }
}

/// The territory is the non-N bases of the reference and not the covered ones.
#[test]
fn the_territory_is_the_non_n_bases() {
    let text = corpus();
    assert_eq!(reference().len(), 200);
    assert_eq!(genome_territory(&reference()), 190);
    for case in CASES {
        assert_eq!(value(&text, case, "GENOME_TERRITORY"), 190.0, "{case}");
    }
}

/// The seven exclusions partition the bases that did not count, so their sum is PCT_EXC_TOTAL.
#[test]
fn the_exclusions_are_a_partition() {
    let text = corpus();
    for case in CASES {
        let m = metrics(&text, case);
        let sum: f64 = [
            "PCT_EXC_ADAPTER",
            "PCT_EXC_MAPQ",
            "PCT_EXC_DUPE",
            "PCT_EXC_UNPAIRED",
            "PCT_EXC_BASEQ",
            "PCT_EXC_OVERLAP",
            "PCT_EXC_CAPPED",
        ]
        .iter()
        .map(|name| m[*name].parse::<f64>().expect("a number"))
        .sum();
        let total: f64 = m["PCT_EXC_TOTAL"].parse().expect("a number");
        assert!((sum - total).abs() < 1e-6, "{case}: {sum} against {total}");
    }
}

/// The mapping quality and the duplicate flag take whole reads, and reach different counters.
#[test]
fn a_whole_read_goes_to_one_counter() {
    let text = corpus();
    assert_eq!(value(&text, "low-mapping-quality", "PCT_EXC_MAPQ"), 1.0);
    assert_eq!(value(&text, "low-mapping-quality", "PCT_EXC_DUPE"), 0.0);
    assert_eq!(value(&text, "duplicate", "PCT_EXC_DUPE"), 1.0);
    assert_eq!(value(&text, "duplicate", "PCT_EXC_MAPQ"), 0.0);
    let counted = Arguments {
        count_unpaired: true,
        ..Arguments::default()
    };
    let low = WgsRead {
        adapter: false,
        mapping_quality: 5,
        duplicate: false,
        paired: false,
        mate_unmapped: false,
    };
    assert_eq!(read_fate(&low, &counted), Some(Fate::MappingQuality));
    // The order is a chain: a duplicate that is also under the floor goes to the floor.
    let both = WgsRead {
        duplicate: true,
        ..low
    };
    assert_eq!(read_fate(&both, &counted), Some(Fate::MappingQuality));
}

/// The base quality takes single bases, so a read contributes some of them and not others, and an
/// N is one of those bases whatever its quality says.
#[test]
fn the_base_quality_takes_single_bases() {
    let text = corpus();
    assert_eq!(value(&text, "low-base-quality", "PCT_EXC_BASEQ"), 0.5);
    assert_eq!(value(&text, "low-base-quality", "MEAN_COVERAGE"), 0.052632);
    // Five Ns of twenty, all at quality 'I' which is forty.
    assert_eq!(value(&text, "n-bases", "PCT_EXC_BASEQ"), 0.25);
    let arguments = Arguments::default();
    assert_eq!(base_fate(b'N', 40, false, 0, &arguments), Fate::BaseQuality);
    assert_eq!(base_fate(b'A', 40, false, 0, &arguments), Fate::Counted);
    assert_eq!(base_fate(b'A', 2, false, 0, &arguments), Fate::BaseQuality);
}

/// An unpaired read and a pair with one end unmapped reach the same counter.
#[test]
fn an_unmapped_mate_is_unpaired() {
    let text = corpus();
    assert_eq!(value(&text, "one-unpaired-read", "PCT_EXC_UNPAIRED"), 1.0);
    assert_eq!(value(&text, "mate-unmapped", "PCT_EXC_UNPAIRED"), 1.0);
    // And COUNT_UNPAIRED lets the first through.
    assert_eq!(
        value(&text, "one-unpaired-read-counted", "PCT_EXC_UNPAIRED"),
        0.0
    );
    assert_eq!(
        value(&text, "one-unpaired-read-counted", "MEAN_COVERAGE"),
        0.105263
    );
    let half = WgsRead {
        adapter: false,
        mapping_quality: 60,
        duplicate: false,
        paired: true,
        mate_unmapped: true,
    };
    assert_eq!(
        read_fate(&half, &Arguments::default()),
        Some(Fate::Unpaired)
    );
    assert_eq!(
        read_fate(
            &half,
            &Arguments {
                count_unpaired: true,
                ..Arguments::default()
            }
        ),
        None
    );
}

/// A pair's overlap counts once: two twenty-base ends on the same span cover it at depth one.
#[test]
fn the_overlap_counts_once() {
    let text = corpus();
    assert_eq!(value(&text, "pair-overlapping", "PCT_EXC_OVERLAP"), 0.5);
    assert_eq!(value(&text, "pair-overlapping", "MEAN_COVERAGE"), 0.105263);
    // A disjoint pair of the same two reads covers twice as much and overlaps nowhere.
    assert_eq!(value(&text, "pair-disjoint", "PCT_EXC_OVERLAP"), 0.0);
    assert_eq!(value(&text, "pair-disjoint", "MEAN_COVERAGE"), 0.210526);
    let (_, depths) = walk(&text, "pair-overlapping", &arguments("pair-overlapping"));
    assert!(depths[..20].iter().all(|d| *d == 1));
}

/// The cap truncates the depth and counts the remainder.
#[test]
fn the_cap_truncates_and_counts_the_remainder() {
    let text = corpus();
    assert_eq!(value(&text, "deep-uncapped", "PCT_EXC_CAPPED"), 0.0);
    assert_eq!(value(&text, "deep-uncapped", "MEAN_COVERAGE"), 1.052632);
    // Ten reads of twenty bases under a cap of two: forty counted, a hundred and sixty capped.
    assert_eq!(value(&text, "deep-capped", "PCT_EXC_CAPPED"), 0.8);
    assert_eq!(value(&text, "deep-capped", "MEAN_COVERAGE"), 0.210526);
    let (counts, depths) = walk(&text, "deep-capped", &arguments("deep-capped"));
    assert_eq!(counts.counted, 40);
    assert_eq!(counts.capped, 160);
    assert!(depths[..20].iter().all(|d| *d == 2));
}

/// The histogram argument adds a column rather than a table, and changes no metric.
#[test]
fn the_histogram_argument_adds_a_column() {
    let text = corpus();
    let plain = field(&text, "histogram", "pair-disjoint").expect("a histogram");
    let with = field(&text, "histogram", "with-histogram").expect("a histogram");
    assert_eq!(
        plain.lines().next(),
        Some("coverage\thigh_quality_coverage_count")
    );
    assert_eq!(
        with.lines().next(),
        Some("coverage\thigh_quality_coverage_count\tunfiltered_baseq_count")
    );
    assert_eq!(
        metrics(&text, "with-histogram"),
        metrics(&text, "pair-disjoint")
    );
}

/// A file with no reads still reports its territory, every depth being zero.
#[test]
fn a_file_with_no_reads_still_has_a_territory() {
    let text = corpus();
    assert_eq!(value(&text, "empty", "GENOME_TERRITORY"), 190.0);
    assert_eq!(value(&text, "empty", "MEAN_COVERAGE"), 0.0);
    assert_eq!(value(&text, "empty", "PCT_EXC_TOTAL"), 0.0);
    let (counts, depths) = walk(&text, "empty", &Arguments::default());
    assert_eq!(counts, Counts::default());
    assert_eq!(counts.total(), 0);
    // A total of nothing is a fraction of nought and not a division by zero.
    assert_eq!(counts.excluded_fraction(), 0.0);
    assert!(depths.iter().all(|d| *d == 0));
    // The histogram's first bin is the whole territory.
    let histogram = field(&text, "histogram", "empty").expect("a histogram");
    assert_eq!(histogram.lines().nth(1), Some("0\t190"));
}
