//! Conformance for `SamToFastqWithTags`' per-tag-group FASTQ output against Picard 3.4.0.
//!
//! The corpus carries the input SAM and, per case, each tag FASTQ Picard produced (`out\t<case>:<file>`
//! rows). The port runs `sam_to_fastq_with_tags_unpaired` on the same records and tag groups and must
//! reproduce every file byte-for-byte. FASTQ has no header, so each file is compared raw. The base
//! read FASTQ is the already-ported `SamToFastq` path and is not covered here.

use std::collections::HashMap;
use std::io::Read;

use htsjdk_bam::sam_file::read_sam;
use picard_analysis::sam_to_fastq_with_tags::{
    sam_to_fastq_with_tags_paired, sam_to_fastq_with_tags_unpaired, TagGroup,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/sam_to_fastq_with_tags.txt.gz");
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

/// Parse the corpus into `input[case]` and `out["case:file"]` maps.
fn parse() -> (HashMap<String, String>, HashMap<String, String>) {
    let text = corpus();
    let mut inputs = HashMap::new();
    let mut outs = HashMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let key = it.next().unwrap().to_string();
        let payload = unescape(it.next().unwrap_or(""));
        match kind {
            "input" => {
                inputs.insert(key, payload);
            }
            "out" => {
                outs.insert(key, payload);
            }
            "rc" => {} // every case returns 0
            other => panic!("unexpected row kind {other}"),
        }
    }
    (inputs, outs)
}

#[test]
fn every_tag_fastq_is_byte_identical() {
    let (inputs, outs) = parse();

    // Each case's tag groups, matching the harness's SEQUENCE_TAG_GROUP / QUALITY_TAG_GROUP args.
    let cases: Vec<(&str, Vec<TagGroup>)> = vec![
        (
            "two_groups",
            vec![
                TagGroup::new("CR").with_quality("CY"),
                TagGroup::new("CB,UR").with_quality("CY,UY"),
            ],
        ),
        ("no_quality", vec![TagGroup::new("CB,UR")]),
    ];

    for (case, groups) in &cases {
        let sam = inputs.get(*case).expect("case input");
        let records = read_sam(sam).unwrap().1;
        let produced = sam_to_fastq_with_tags_unpaired(&records, groups).expect("tag fastq");
        assert_files(case, &produced, &outs);
    }

    // The paired case: the same tags but reads split across first/second of pair, into one file.
    let sam = inputs.get("paired").expect("paired input");
    let records = read_sam(sam).unwrap().1;
    let groups = [TagGroup::new("CR").with_quality("CY")];
    let produced = sam_to_fastq_with_tags_paired(&records, &groups).expect("paired tag fastq");
    assert_files("paired", &produced, &outs);
}

fn assert_files(case: &str, produced: &[(String, String)], outs: &HashMap<String, String>) {
    for (file_name, text) in produced {
        let key = format!("{case}:{file_name}");
        let golden = outs
            .get(&key)
            .unwrap_or_else(|| panic!("no golden for {key}"));
        assert_eq!(text, golden, "{key}");
    }
}
