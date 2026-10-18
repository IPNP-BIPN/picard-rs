//! Conformance for `UmiAwareMarkDuplicatesWithMateCigar` on reads with no read group, against
//! Picard 3.5.0, run through the binary.
//!
//! Golden from `tools/matecigardup-conformance/UmiNoReadGroupDump.java`. Each input is rebuilt
//! here as SAM text from the dump's own case list:
//!
//!  * **a duplicate set whose representative has no RG tag throws a NullPointerException**,
//!    single-end, unmapped and with no `@RG` line alike;
//!  * **what the output holds at the throw**: a SAM output its header and no record, a BAM output
//!    nothing;
//!  * **and the queryname refusal coming first, with no output at all.**

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn corpus() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/umi_no_read_group.txt.gz");
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
                .replace("\\\\", "\\")
        })
}

fn header(groups: bool, order: &str) -> String {
    let mut text = format!("@HD\tVN:1.6\tSO:{order}\n@SQ\tSN:chr1\tLN:2000\n");
    if groups {
        text.push_str("@RG\tID:rg1\tSM:sample1\tLB:lib1\n");
    }
    text
}

/// name, flags, start (0 when unplaced), mate start (0 when the mate is unmapped), RG.
type Spec = (&'static str, u16, i32, i32, Option<&'static str>);

fn line((name, flags, start, mate_start, group): Spec) -> String {
    let mapped = flags & 0x4 == 0;
    let (reference, position) = if start > 0 { ("chr1", start) } else { ("*", 0) };
    let mut tags = String::new();
    let (mut mate_reference, mut mate_position) = ("*", 0);
    if flags & 0x1 != 0 {
        mate_reference = "=";
        mate_position = if mate_start == 0 { start } else { mate_start };
        if flags & 0x8 == 0 {
            tags.push_str("\tMC:Z:20M");
        }
    }
    tags.push_str("\tRX:Z:ACGT");
    if let Some(group) = group {
        tags.push_str(&format!("\tRG:Z:{group}"));
    }
    format!(
        "{name}\t{flags}\t{reference}\t{position}\t{}\t{}\t{mate_reference}\t{mate_position}\t0\t{}\t{}{tags}\n",
        if mapped { 60 } else { 0 },
        if mapped { "20M" } else { "*" },
        "A".repeat(20),
        "I".repeat(20),
    )
}

/// The dump's two pairs at one position, in coordinate order.
fn pairs(first: Option<&'static str>, second: Option<&'static str>) -> Vec<Spec> {
    vec![
        ("HWI:1:FC:1:1:1:1", 99, 100, 300, first),
        ("HWI:1:FC:1:1:9000:9000", 99, 100, 300, second),
        ("HWI:1:FC:1:1:1:1", 147, 300, 100, first),
        ("HWI:1:FC:1:1:9000:9000", 147, 300, 100, second),
    ]
}

/// Every case of the dump but `undeclared-rg`, which htsjdk's strict reader refuses before the
/// tool runs: that refusal is the reader's, not this tool's.
fn cases() -> Vec<(&'static str, String, Vec<Spec>, &'static str)> {
    let coordinate = |groups| header(groups, "coordinate");
    vec![
        (
            "pairs-with-rg",
            coordinate(true),
            pairs(Some("rg1"), Some("rg1")),
            "out.sam",
        ),
        (
            "pairs-without-rg",
            coordinate(true),
            pairs(None, None),
            "out.sam",
        ),
        (
            "pairs-without-rg-and-no-rg-lines",
            coordinate(false),
            pairs(None, None),
            "out.sam",
        ),
        (
            "one-pair-with-rg-one-without",
            coordinate(true),
            pairs(Some("rg1"), None),
            "out.sam",
        ),
        (
            "one-pair-without-rg-one-with",
            coordinate(true),
            pairs(None, Some("rg1")),
            "out.sam",
        ),
        (
            "singles-without-rg",
            coordinate(true),
            vec![
                ("HWI:1:FC:1:1:1:1", 0, 100, 0, None),
                ("HWI:1:FC:1:1:9000:9000", 0, 100, 0, None),
            ],
            "out.sam",
        ),
        (
            "rg-reads-then-one-without",
            coordinate(true),
            vec![
                ("HWI:1:FC:1:1:1:1", 0, 100, 0, Some("rg1")),
                ("HWI:1:FC:1:1:2:2", 0, 200, 0, Some("rg1")),
                ("HWI:1:FC:1:1:3:3", 0, 300, 0, Some("rg1")),
                ("HWI:1:FC:1:1:4:4", 0, 900, 0, None),
            ],
            "out.sam",
        ),
        (
            "mate-unmapped-without-rg",
            coordinate(true),
            vec![
                ("HWI:1:FC:1:1:1:1", 0x1 | 0x8 | 0x40, 100, 0, None),
                ("HWI:1:FC:1:1:1:1", 0x1 | 0x4 | 0x80, 100, 0, None),
            ],
            "out.sam",
        ),
        (
            "unplaced-unmapped-without-rg",
            coordinate(true),
            vec![
                ("HWI:1:FC:1:1:1:1", 0, 100, 0, Some("rg1")),
                ("HWI:1:FC:1:1:2:2", 0x4, 0, 0, None),
            ],
            "out.sam",
        ),
        (
            "pairs-without-rg-to-bam",
            coordinate(true),
            pairs(None, None),
            "out.bam",
        ),
        (
            "queryname-pairs-without-rg",
            header(true, "queryname"),
            vec![
                ("HWI:1:FC:1:1:1:1", 99, 100, 300, None),
                ("HWI:1:FC:1:1:1:1", 147, 300, 100, None),
                ("HWI:1:FC:1:1:9000:9000", 99, 100, 300, None),
                ("HWI:1:FC:1:1:9000:9000", 147, 300, 100, None),
            ],
            "out.sam",
        ),
    ]
}

/// A directory of its own per case.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("umi-no-read-group-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Every case refuses, or marks, as the reference did, and leaves the file it left.
#[test]
fn every_case_ends_as_the_reference_ended() {
    let text = corpus();
    for (case, header, reads, output) in cases() {
        let dir = scratch(case);
        let input: String = header + &reads.into_iter().map(line).collect::<String>();
        std::fs::write(dir.join("in.sam"), input).unwrap();
        let run = Command::new(env!(
            "CARGO_BIN_EXE_umi-aware-mark-duplicates-with-mate-cigar"
        ))
        .current_dir(&dir)
        .args([
            "I=in.sam",
            &format!("O={output}"),
            "M=metrics.txt",
            "UMI_METRICS=umi.txt",
        ])
        .output()
        .unwrap();
        let out = dir.join(output);
        match field(&text, "error", case) {
            Some(error) => {
                assert!(!run.status.success(), "{case}");
                // The dump writes `class:message`, the JVM `class: message`.
                let (class, message) = error.split_once(':').expect("class:message");
                assert_eq!(
                    String::from_utf8_lossy(&run.stderr),
                    format!("Exception in thread \"main\" {class}: {message}\n"),
                    "{case}"
                );
                let partial = field(&text, "partial", case).expect("the partial output");
                match partial.as_str() {
                    "absent" => assert!(!out.exists(), "{case}"),
                    "0 bytes" => assert_eq!(std::fs::read(&out).unwrap().len(), 0, "{case}"),
                    sam => assert_eq!(std::fs::read_to_string(&out).unwrap(), sam, "{case}"),
                }
            }
            None => {
                assert!(run.status.success(), "{case}");
                let (_, records) = picard_analysis::metrics_cli::read_input(out.to_str().unwrap());
                let produced: Vec<(String, u16)> = records
                    .iter()
                    .map(|record| (record.read_name.clone(), record.flags))
                    .collect();
                let expected: Vec<(String, u16)> = field(&text, "marked", case)
                    .expect("the golden")
                    .lines()
                    .map(|line| {
                        let columns: Vec<&str> = line.split('\t').collect();
                        (columns[0].to_string(), columns[1].parse().unwrap())
                    })
                    .collect();
                assert_eq!(produced, expected, "{case}");
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
