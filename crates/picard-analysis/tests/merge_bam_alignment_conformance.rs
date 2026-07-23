//! Conformance for `MergeBamAlignment`'s merged record (unpaired, single primary hit, coordinate
//! output) against Picard 3.4.0.
//!
//! Each case carries the unmapped SAM, the aligned SAM, and the merged record data line(s) Picard
//! produced. The port reads the unmapped and aligned records, runs `merge_aligned_fragment` (which does
//! the transfer, the `PG` linkage, and the reference-based `NM`/`MD`/`UQ`), and must reproduce each
//! record line byte-for-byte. The output header (`@SQ`/`@RG`/`@PG`) is a later slice and is not
//! compared. The reference is a fixed `chr1` of forty bases.

use std::collections::HashMap;
use std::io::Read;

use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use picard_analysis::merge_bam_alignment::{merge_aligned_fragment, merge_bam_alignment_records};

const REF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/merge_bam_alignment.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

#[derive(Default)]
struct Case {
    unmapped: String,
    aligned: String,
    records: Vec<String>,
}

fn cases() -> Vec<(String, Case)> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, Case> = HashMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let name = it.next().unwrap().to_string();
        let payload = unescape(it.next().unwrap_or(""));
        let case = map.entry(name.clone()).or_insert_with(|| {
            order.push(name.clone());
            Case::default()
        });
        match kind {
            "unmapped" => case.unmapped = payload,
            "aligned" => case.aligned = payload,
            "record" => case.records.push(payload),
            "rc" => {}
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|n| (n.clone(), map.remove(&n).unwrap()))
        .collect()
}

/// The single-contig output header (`chr1`), used both to resolve the aligned record's reference name
/// and to render the merged record.
fn out_sequences() -> Vec<SequenceRecord> {
    vec![SequenceRecord::new("chr1", 40)]
}

fn record_line(rec: &BamRecord) -> String {
    let mut header = SamHeader::new();
    header.sequences = out_sequences();
    write_sam(&header, std::slice::from_ref(rec))
        .unwrap()
        .lines()
        .find(|l| !l.starts_with('@'))
        .unwrap()
        .to_string()
}

#[test]
fn every_merged_record_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 4, "case count");
    for (name, case) in &cases {
        let (_, unmapped) = read_sam(&case.unmapped).unwrap();
        let (aligned_header, aligned) = read_sam(&case.aligned).unwrap();
        assert_eq!(unmapped.len(), aligned.len());

        for (u, a) in unmapped.iter().zip(&aligned) {
            // The aligned record's reference name, resolved through the aligned file's dictionary.
            let ref_name = if a.reference_index < 0 {
                "*".to_string()
            } else {
                aligned_header.sequences[a.reference_index as usize]
                    .name
                    .clone()
            };
            let merged =
                merge_aligned_fragment(u, a, &ref_name, &out_sequences(), REF, Some("bwa"))
                    .unwrap();
            let got = record_line(&merged);
            assert!(
                case.records.contains(&got),
                "{name}: produced record not in golden set:\n  got:    {got}\n  golden: {:?}",
                case.records
            );
        }
        assert_eq!(unmapped.len(), case.records.len(), "{name} record count");
    }
}

/// The whole-file producer must match reads by name and reproduce the coordinate-sorted **order**.
#[test]
fn the_whole_file_producer_matches_and_coordinate_sorts() {
    let case = cases()
        .into_iter()
        .find(|(n, _)| n == "multi_read")
        .unwrap()
        .1;
    let (_, unmapped) = read_sam(&case.unmapped).unwrap();
    let (aligned_header, aligned) = read_sam(&case.aligned).unwrap();

    let reference_bases = HashMap::from([("chr1".to_string(), REF.to_vec())]);
    let merged = merge_bam_alignment_records(
        &unmapped,
        &aligned,
        &aligned_header.sequences,
        &out_sequences(),
        &reference_bases,
        Some("bwa"),
    )
    .unwrap();

    let got: Vec<String> = merged.iter().map(record_line).collect();
    assert_eq!(got, case.records, "coordinate-sorted record order");
}
