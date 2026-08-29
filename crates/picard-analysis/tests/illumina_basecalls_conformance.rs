//! Conformance for `IlluminaBasecallsToFastq` and `IlluminaBasecallsToSam` against Picard 3.4.0.
//!
//! Goldens from `tools/illumina-conformance/`, both over the same four-cluster basecalls
//! directory, which is what lets the two be compared with each other as well as with the
//! reference.
//!
//! # What this suite is for
//!
//!  * **the read structure deciding what is written and what is only read**;
//!  * **the read name carrying the position and the filter's verdict, inverted**;
//!  * **two template reads becoming a PAIR in a BAM and two files in a FASTQ**;
//!  * **the failing cluster being flagged in one and absent from neither**;
//!  * **and the qualities being the BCL's six bits, plus thirty-three in a FASTQ.**

use std::io::Read;

use picard_analysis::illumina_basecalls::{
    phred33, position_in_name, read_group, read_name, sam_flags, segment, written_segments,
    Cluster, ReadNameFormat, Run,
};
use picard_analysis::illumina_files::{decode_basecall, parse_read_structure, BaseCall};

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

/// The fixture's four clusters, as the writer wrote them.
fn clusters() -> Vec<Cluster> {
    let cycles = ["ACGT", "ACGT", "AACC", "GGTT"];
    (0..4)
        .map(|index| Cluster {
            calls: cycles
                .iter()
                .map(|cycle| decode_basecall(30 << 2 | base_code(cycle.as_bytes()[index])))
                .collect::<Vec<BaseCall>>(),
            // Three of the four passed; the last did not.
            passed_filter: index != 3,
            // The `.locs` file carries `100.0 * (index + 1)` and `200.0 * (index + 1)`; the NAME
            // carries ten times that plus a thousand.
            x: position_in_name(100.0 * (index as f32 + 1.0), 0.0).0,
            y: position_in_name(0.0, 200.0 * (index as f32 + 1.0)).1,
        })
        .collect()
}

fn base_code(base: u8) -> u8 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}

fn run() -> Run {
    Run {
        machine: "machine".to_string(),
        run_barcode: "run17".to_string(),
        flowcell: "flowcell".to_string(),
        lane: 1,
        tile: 1101,
    }
}

/// The FASTQ a run of `4T` writes, record for record.
#[test]
fn the_fastq_is_the_goldens() {
    let text = corpus("basecalls_to_fastq");
    let structure = parse_read_structure("4T").expect("a structure");
    let recorded = field(&text, "fastq", "one-read-of-four.reads.1.fastq").expect("the golden");
    let mut lines = recorded.lines();
    for cluster in clusters() {
        let name = lines.next().expect("a name");
        let bases = lines.next().expect("bases");
        let _plus = lines.next().expect("a plus");
        let qualities = lines.next().expect("qualities");

        assert_eq!(
            name,
            format!(
                "@{}",
                read_name(&run(), &cluster, ReadNameFormat::Casava18, None, None)
            )
        );
        let (produced, produced_qualities) = segment(&cluster, &structure, 0);
        assert_eq!(String::from_utf8(produced).expect("bases"), bases);
        assert_eq!(
            String::from_utf8(phred33(&produced_qualities)).expect("qualities"),
            qualities
        );
    }
    // The cluster that failed the filter says so in its name, with `Y` rather than `N`.
    assert!(recorded.lines().any(|line| line.contains(" :Y:0:")));
}

/// A barcode segment is read and not written, and a skip is neither.
#[test]
fn the_structure_decides_what_is_written() {
    let text = corpus("basecalls_to_fastq");
    for (case, structure, files) in [
        ("one-read-of-four", "4T", 1),
        ("two-reads-of-two", "2T2T", 2),
        ("a-skipped-segment", "2T2S", 1),
    ] {
        let parsed = parse_read_structure(structure).expect("a structure");
        assert_eq!(written_segments(&parsed).len(), files, "{case}");
        let names = field(&text, "files", case).expect("the golden");
        assert_eq!(names.split_whitespace().count(), files, "{case}: {names}");
    }
    // A barcode is written too, but to a file of its own rather than as a read.
    let names = field(&text, "files", "a-barcode-segment").expect("the golden");
    assert_eq!(names, "reads.1.fastq reads.barcode_1.fastq");
    let parsed = parse_read_structure("2T2B").expect("a structure");
    assert_eq!(written_segments(&parsed), vec![0]);
    // And the two reads of `2T2T` carry their ordinal in the name, where one read carries none.
    let two = field(&text, "fastq", "two-reads-of-two.reads.1.fastq").expect("the golden");
    assert!(two.lines().next().expect("a name").contains(" 1:N:0:"));
    let one = field(&text, "fastq", "one-read-of-four.reads.1.fastq").expect("the golden");
    assert!(one.lines().next().expect("a name").contains(" :N:0:"));
}

/// The coordinates in a name are ten times the file's, plus a thousand.
#[test]
fn the_position_in_a_name_is_not_the_position_in_the_file() {
    assert_eq!(position_in_name(100.0, 200.0), (2000, 3000));
    assert_eq!(position_in_name(500.0, 900.0), (6000, 10000));
    let text = corpus("basecalls_to_fastq");
    let recorded = field(&text, "fastq", "one-read-of-four.reads.1.fastq").expect("the golden");
    assert!(recorded
        .lines()
        .next()
        .expect("a name")
        .contains(":2000:3000 "));
}

/// The other read name format, which drops everything but the run and the position.
#[test]
fn the_other_name_format_is_shorter() {
    let text = corpus("basecalls_to_fastq");
    let recorded =
        field(&text, "fastq", "the-other-read-name-format.reads.1.fastq").expect("the golden");
    let name = recorded.lines().next().expect("a name");
    assert_eq!(
        name,
        format!(
            "@{}",
            read_name(&run(), &clusters()[0], ReadNameFormat::Illumina, None, None)
        )
    );
}

/// The BAM says what the FASTQ cannot: the pairing, and the vendor check.
#[test]
fn the_bam_flags_are_the_goldens() {
    let text = corpus("basecalls_to_sam");
    // One template read: unmapped and nothing else.
    let recorded = field(&text, "sam", "one-read-of-four.reads.bam").expect("the golden");
    for (line, cluster) in recorded.lines().zip(clusters()) {
        let flags: u16 = line
            .split('\t')
            .nth(1)
            .expect("flags")
            .parse()
            .expect("a number");
        assert_eq!(flags, sam_flags(1, 0, cluster.passed_filter));
    }
    // The failing cluster is 516, which is unmapped plus the vendor check.
    assert_eq!(sam_flags(1, 0, false), 0x4 | 0x200);
    assert!(recorded
        .lines()
        .any(|line| line.split('\t').nth(1) == Some("516")));

    // Two template reads: a pair, flagged 77 and 141.
    let recorded = field(&text, "sam", "two-reads-of-two.reads.bam").expect("the golden");
    let flags: Vec<u16> = recorded
        .lines()
        .take(2)
        .map(|line| {
            line.split('\t')
                .nth(1)
                .expect("flags")
                .parse()
                .expect("a number")
        })
        .collect();
    assert_eq!(flags, vec![sam_flags(2, 0, true), sam_flags(2, 1, true)]);
    assert_eq!(flags, vec![77, 141]);

    // And the run's identity is written once, in the read group.
    let group = field(&text, "header", "one-read-of-four.reads.bam").expect("the golden");
    let (id, unit) = read_group(&run(), None, false);
    assert!(group.contains(&format!("ID:{id}")), "{group}");
    assert!(group.contains(&format!("PU:{unit}")), "{group}");
    // With a barcode, the platform unit carries it.
    let group = field(&text, "header", "split-by-barcode.first.bam").expect("the golden");
    let (_, unit) = read_group(&run(), Some("AG"), true);
    assert!(group.contains(&format!("PU:{unit}")), "{group}");
    assert_eq!(unit, "run17.1.AG");
}

/// The cluster that failed the filter is dropped only where the argument says so.
#[test]
fn the_failing_cluster_is_flagged_or_dropped() {
    let fastq = corpus("basecalls_to_fastq");
    let with = field(&fastq, "fastq", "one-read-of-four.reads.1.fastq").expect("the golden");
    let without =
        field(&fastq, "fastq", "without-the-non-pf-reads.reads.1.fastq").expect("the golden");
    assert_eq!(with.lines().count(), 16);
    assert_eq!(without.lines().count(), 12);

    let sam = corpus("basecalls_to_sam");
    let with = field(&sam, "sam", "one-read-of-four.reads.bam").expect("the golden");
    let without = field(&sam, "sam", "without-the-non-pf-reads.reads.bam").expect("the golden");
    assert_eq!(with.lines().count(), 4);
    assert_eq!(without.lines().count(), 3);
    // And what is dropped is exactly the record the flag would have marked.
    assert!(with
        .lines()
        .any(|line| line.split('\t').nth(1) == Some("516")));
    assert!(!without
        .lines()
        .any(|line| line.split('\t').nth(1) == Some("516")));
}
