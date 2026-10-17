//! Conformance for `SortSam` spilling to `TMP_DIR` against Picard 3.5.0.
//!
//! The corpus carries an unsorted SAM of 600 records, most of them equal under the comparators,
//! and what Picard wrote sorting it with `MAX_RECORDS_IN_RAM=50`: twelve runs merged. The dump
//! only emits once the same sorts in memory wrote the same bytes. The port's binary sorts the same
//! input with the same threshold, so it spills as often, and must write every file byte for byte.

use std::io::Read;
use std::process::Command;

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sort_sam_spill.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\n", "\n")
        .replace("\\t", "\t")
        .replace("\\\\", "\\")
}

fn from_hex(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}

fn payload(kind: &str, case: &str) -> String {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .find_map(|l| {
            let mut it = l.splitn(3, '\t');
            (it.next()? == kind && it.next()? == case).then(|| it.next().unwrap_or("").to_string())
        })
        .unwrap_or_else(|| panic!("no {kind} {case} row"))
}

/// Runs the port's `SortSam` on the corpus input and returns the output and, if it wrote one,
/// the index.
fn sort(order: &str, ext: &str) -> (Vec<u8>, Option<Vec<u8>>) {
    let dir = std::env::temp_dir().join(format!(
        "sort-sam-spill-{}-{order}-{ext}",
        std::process::id()
    ));
    let tmp = dir.join("tmp");
    std::fs::create_dir_all(&tmp).unwrap();
    let input = dir.join("in.sam");
    std::fs::write(&input, unescape(&payload("sam", "input"))).unwrap();
    let output = dir.join(format!("out.{ext}"));
    let status = Command::new(env!("CARGO_BIN_EXE_SortSam"))
        .arg(format!("INPUT={}", input.display()))
        .arg(format!("OUTPUT={}", output.display()))
        .arg(format!("SORT_ORDER={order}"))
        .arg("MAX_RECORDS_IN_RAM=50")
        .arg(format!("TMP_DIR={}", tmp.display()))
        .arg(format!("CREATE_INDEX={}", order == "coordinate"))
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(
        std::fs::read_dir(&tmp).unwrap().count(),
        0,
        "a run was left in TMP_DIR"
    );
    let out = std::fs::read(&output).unwrap();
    let index = std::fs::read(dir.join("out.bai")).ok();
    std::fs::remove_dir_all(&dir).unwrap();
    (out, index)
}

#[test]
fn the_spilled_coordinate_sort_and_its_index_are_byte_identical() {
    let (bam, bai) = sort("coordinate", "bam");
    assert!(
        bam == from_hex(&payload("bam", "coordinate-bam")),
        "BAM differs"
    );
    assert!(
        bai.unwrap() == from_hex(&payload("bai", "coordinate-bam")),
        "BAI differs"
    );
}

#[test]
fn the_spilled_queryname_sort_is_byte_identical_as_bam() {
    let (bam, _) = sort("queryname", "bam");
    assert!(
        bam == from_hex(&payload("bam", "queryname-bam")),
        "BAM differs"
    );
}

#[test]
fn the_spilled_queryname_sort_is_byte_identical_as_sam() {
    let (sam, _) = sort("queryname", "sam");
    assert_eq!(
        String::from_utf8(sam).unwrap(),
        unescape(&payload("sam", "queryname-sam"))
    );
}
