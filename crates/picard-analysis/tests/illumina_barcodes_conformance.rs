//! Conformance for the Illumina binary formats and `ExtractIlluminaBarcodes` against Picard 3.4.0.
//!
//! Goldens from `tools/illumina-conformance/`: the directory check, the lane metrics and the
//! barcodes, all three over a basecalls directory written byte by byte by the same fixture writer.
//!
//! # What this suite is for
//!
//!  * **the basecall byte: two bits of base, six of quality, and zero meaning a no-call**;
//!  * **the read structure cutting the same cycles into different segments**;
//!  * **the barcode match being two tests, so equidistant barcodes match neither**;
//!  * **a quality floor turning a base into a no-call rather than a mismatch**;
//!  * **and the metrics counting PF apart from the rest.**

use std::io::Read;

use picard_analysis::extract_illumina_barcodes::{
    best_match, metrics, mismatches, Observed, Options,
};
use picard_analysis::illumina_files::{
    decode_basecall, parse_read_structure, segment_cycles, total_cycles, SegmentKind,
};

fn corpus(name: &str) -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.txt.gz"));
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

/// The fixture's own bases, which the writer put there: four clusters over four cycles.
fn cluster(index: usize) -> Vec<u8> {
    // Cycle one and two are `ACGT`, cycle three is `AACC` and cycle four is `GGTT`.
    let cycles = ["ACGT", "ACGT", "AACC", "GGTT"];
    cycles.iter().map(|cycle| cycle.as_bytes()[index]).collect()
}

/// A basecall byte is the base and the quality, and zero is neither.
#[test]
fn a_basecall_byte_is_two_bits_and_six() {
    // The fixture writes quality 30, so `A` is `30 << 2` and `T` is that plus three.
    assert_eq!(decode_basecall(30 << 2).base, b'A');
    assert_eq!(decode_basecall(30 << 2).quality, 30);
    assert_eq!(decode_basecall((30 << 2) | 3).base, b'T');
    assert_eq!(decode_basecall((30 << 2) | 3).quality, 30);
    // A byte of zero is a no-call, and its quality is not zero-the-quality but no quality at all.
    assert_eq!(decode_basecall(0).base, b'.');
    assert_eq!(decode_basecall(0).quality, 0);
}

/// The read structure cuts the cycles, and the golden's cases are its arithmetic.
#[test]
fn the_read_structure_cuts_the_cycles() {
    let four = parse_read_structure("4T").expect("a structure");
    assert_eq!(total_cycles(&four), 4);
    assert_eq!(segment_cycles(&four), vec![vec![1, 2, 3, 4]]);

    let two_and_two = parse_read_structure("2T2T").expect("a structure");
    assert_eq!(segment_cycles(&two_and_two), vec![vec![1, 2], vec![3, 4]]);

    let barcoded = parse_read_structure("2T2B").expect("a structure");
    assert_eq!(barcoded[1].kind, SegmentKind::Barcode);
    let skipped = parse_read_structure("2T2S").expect("a structure");
    assert_eq!(skipped[1].kind, SegmentKind::Skip);

    // A longer structure asks for more cycles than the fixture has, which is what the directory
    // check refuses: `more-cycles-than-there-are` is a failure and `3T` over the same files is not.
    let text = corpus("check_illumina_directory");
    assert_eq!(
        field(&text, "code", "a-complete-directory").as_deref(),
        Some("0")
    );
    assert_eq!(
        field(&text, "code", "more-cycles-than-there-are").as_deref(),
        Some("1")
    );
    assert_eq!(
        total_cycles(&parse_read_structure("6T").expect("a structure")),
        6
    );
    assert_eq!(
        field(&text, "code", "a-missing-cycle-not-asked-for").as_deref(),
        Some("0")
    );
    assert_eq!(
        total_cycles(&parse_read_structure("3T").expect("a structure")),
        3
    );
    // And a letter the reference does not know is not a structure at all.
    assert!(parse_read_structure("4X").is_none());
    assert!(parse_read_structure("T").is_none());
}

/// The barcodes the fixture's clusters carry, read as `2T2B`.
fn observed(index: usize) -> Observed {
    let bases = cluster(index);
    Observed {
        bases: bases[2..4].to_vec(),
        qualities: vec![30, 30],
    }
}

/// Every barcode case of the golden, decided by the port.
#[test]
fn the_barcodes_are_the_goldens() {
    let text = corpus("extract_illumina_barcodes");
    // The fixture's four clusters carry `AG`, `AG`, `CT` and `CT`.
    let carried: Vec<String> = (0..4)
        .map(|index| String::from_utf8(observed(index).bases).expect("bases"))
        .collect();
    assert_eq!(carried, ["AG", "AG", "CT", "CT"]);

    // Two declared barcodes: every cluster matches one of them exactly.
    let two = vec![b"AG".to_vec(), b"CT".to_vec()];
    for index in 0..4 {
        let decision = best_match(&observed(index), &two, &Options::default());
        assert!(decision.matched, "{index}");
        assert_eq!(decision.mismatches, 0);
    }
    let rows = metrics(
        &(0..4)
            .map(|index| (observed(index), index != 3))
            .collect::<Vec<_>>(),
        &two,
        &Options::default(),
    );
    let recorded = field(&text, "metrics", "two-barcodes").expect("the golden");
    for row in &rows {
        let line = recorded
            .lines()
            .find(|line| line.starts_with(&format!("{}\t", row.barcode)))
            .unwrap_or_else(|| panic!("{}", row.barcode));
        let columns: Vec<&str> = line.split('\t').collect();
        assert_eq!(columns[4].parse::<i64>().expect("reads"), row.reads);
        assert_eq!(columns[5].parse::<i64>().expect("pf reads"), row.pf_reads);
        assert_eq!(
            columns[6].parse::<i64>().expect("perfect"),
            row.perfect_matches
        );
    }

    // One declared barcode: two clusters match it and two match nothing.
    let one = vec![b"AG".to_vec()];
    let matched = (0..4)
        .filter(|index| best_match(&observed(*index), &one, &Options::default()).matched)
        .count();
    assert_eq!(matched, 2);
}

/// The two thresholds, each side of their cut.
#[test]
fn the_match_is_two_tests() {
    // A barcode one base from every cluster takes them all at the default allowance.
    let near = vec![b"AT".to_vec()];
    for index in 0..4 {
        assert!(best_match(&observed(index), &near, &Options::default()).matched);
    }
    // And none at zero.
    let strict = Options {
        max_mismatches: 0,
        ..Options::default()
    };
    for index in 0..4 {
        assert!(!best_match(&observed(index), &near, &strict).matched);
    }

    // Two barcodes EQUIDISTANT from `AG` match neither, because the better wins by nothing.
    let equidistant = vec![b"AA".to_vec(), b"GG".to_vec()];
    let decision = best_match(&observed(0), &equidistant, &Options::default());
    assert_eq!(decision.mismatches, 1);
    assert_eq!(decision.mismatch_delta, 0);
    assert!(!decision.matched);
    // With no delta required, the FIRST of the two takes it.
    let no_delta = Options {
        min_mismatch_delta: 0,
        ..Options::default()
    };
    let decision = best_match(&observed(0), &equidistant, &no_delta);
    assert!(decision.matched);
    assert_eq!(decision.barcode, Some(0));

    // And the golden agrees about both, in its per-tile file: `N` with a lower-case barcode where
    // the delta refused it, `Y` with an upper-case one where it did not.
    let text = corpus("extract_illumina_barcodes");
    let refused = field(
        &text,
        "barcodes",
        "two-equidistant-barcodes.s_1_1101_barcode.txt",
    )
    .expect("the golden");
    assert!(refused
        .lines()
        .next()
        .expect("a line")
        .contains("\tN\taa\t"));
    let taken = field(
        &text,
        "barcodes",
        "two-equidistant-barcodes-with-no-delta.s_1_1101_barcode.txt",
    )
    .expect("the golden");
    assert!(taken.lines().next().expect("a line").contains("\tY\tAA\t"));
}

/// A base below the quality floor is a no-call rather than a mismatch.
#[test]
fn the_quality_floor_makes_no_calls() {
    let two = vec![b"AG".to_vec(), b"CT".to_vec()];
    let floor = Options {
        minimum_base_quality: 40,
        ..Options::default()
    };
    // The fixture's qualities are thirty, so a floor of forty rejects every base of every cluster.
    for index in 0..4 {
        assert_eq!(mismatches(&observed(index), &two[0], 40), 2);
        assert!(!best_match(&observed(index), &two, &floor).matched);
    }
    // Which is what the golden says: every read goes to the unmatched row.
    let text = corpus("extract_illumina_barcodes");
    let recorded = field(&text, "metrics", "a-quality-floor-above-the-bases").expect("the golden");
    let unmatched = recorded
        .lines()
        .find(|line| line.starts_with("NN\t"))
        .expect("the unmatched row");
    assert_eq!(unmatched.split('\t').nth(4), Some("4"));
    // And a floor below them changes nothing.
    let low = Options {
        minimum_base_quality: 10,
        ..Options::default()
    };
    assert!(best_match(&observed(0), &two, &low).matched);
}
