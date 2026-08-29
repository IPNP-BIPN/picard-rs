//! Conformance for `MarkIlluminaAdapters` against Picard 3.4.0.
//!
//! Golden from `tools/markadapters-conformance/MarkIlluminaAdaptersDump.java`, fifteen runs whose
//! marked records are in the golden as SAM text.
//!
//! # What this suite is for
//!
//!  * **the tag being one-based**;
//!  * **the search taking the LAST start that matches, not the first**;
//!  * **the error allowance being truncated from the overlap's length**;
//!  * **the minimum below which a read is never marked**;
//!  * **and a tag the input carried not surviving a run that finds nothing.**

use std::io::Read;

use picard_analysis::mark_illumina_adapters::{
    clipped_bases, find_index_of_clip_sequence, first_matching_adapter, three_prime, xt_tag,
    AN_EXISTING_TAG_SURVIVES, DEFAULT_ADAPTERS, MAX_ERROR_RATE, MIN_MATCH_BASES, NO_MATCH,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mark_illumina_adapters.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> String {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
        .unwrap_or_else(|| panic!("{kind}/{case}"))
}

/// The records of one case: the bases, and the XT tag if it has one.
fn records(text: &str, kind: &str, case: &str) -> Vec<(String, Option<i32>)> {
    field(text, kind, case)
        .split('\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let tag = columns
                .iter()
                .find_map(|column| column.strip_prefix("XT:i:"))
                .map(|value| value.parse::<i32>().expect("a tag"));
            (columns[9].to_string(), tag)
        })
        .collect()
}

fn tag(text: &str, case: &str) -> Option<i32> {
    records(text, "marked", case)[0].1
}

fn bases(text: &str, case: &str) -> String {
    records(text, "sam", case)[0].0.clone()
}

/// The tag is the index plus one, and the search takes the last start that matches.
#[test]
fn the_tag_is_one_based_and_the_last_match() {
    let text = corpus();
    let adapter = three_prime("INDEXED").expect("the indexed adapter");
    // An adapter at the very first base is `XT:i:1` and not `XT:i:0`.
    let read = bases(&text, "adapter-at-the-first-base");
    let index = find_index_of_clip_sequence(
        read.as_bytes(),
        adapter.as_bytes(),
        MIN_MATCH_BASES,
        MAX_ERROR_RATE,
    );
    assert_eq!(index, 0);
    assert_eq!(xt_tag(index), tag(&text, "adapter-at-the-first-base"));
    assert_eq!(tag(&text, "adapter-at-the-first-base"), Some(1));
    // And one halfway is where it was put.
    let read = bases(&text, "adapter-halfway");
    let index = find_index_of_clip_sequence(
        read.as_bytes(),
        adapter.as_bytes(),
        MIN_MATCH_BASES,
        MAX_ERROR_RATE,
    );
    assert_eq!(xt_tag(index), tag(&text, "adapter-halfway"));
    assert_eq!(tag(&text, "adapter-halfway"), Some(31));
    // A read with no adapter carries no tag at all.
    let read = bases(&text, "no-adapter");
    assert_eq!(
        find_index_of_clip_sequence(
            read.as_bytes(),
            adapter.as_bytes(),
            MIN_MATCH_BASES,
            MAX_ERROR_RATE
        ),
        NO_MATCH
    );
    assert_eq!(tag(&text, "no-adapter"), None);
    // The custom adapter is the clearest case of the search's direction: a poly-T run placed at
    // offset thirty is found at thirty-three, because every later start matches too and the loop
    // returns the last one it can.
    let read = bases(&text, "a-custom-adapter");
    let custom = "TTTTTTTTTTTTTTTTTTTT";
    let index = find_index_of_clip_sequence(
        read.as_bytes(),
        custom.as_bytes(),
        MIN_MATCH_BASES,
        MAX_ERROR_RATE,
    );
    assert_eq!(xt_tag(index), tag(&text, "a-custom-adapter"));
    assert_eq!(tag(&text, "a-custom-adapter"), Some(34));
}

/// The minimum, below which a read is never marked.
#[test]
fn a_read_shorter_than_the_minimum_is_never_marked() {
    let text = corpus();
    let adapter = three_prime("INDEXED").expect("the indexed adapter");
    // Eleven bases of adapter at the end of the read: the loop starts at `len - 12`, which is one
    // base before the adapter, and no start from there on matches.
    let read = bases(&text, "eleven-bases-of-adapter");
    assert_eq!(
        find_index_of_clip_sequence(read.as_bytes(), adapter.as_bytes(), 12, MAX_ERROR_RATE),
        NO_MATCH
    );
    assert_eq!(tag(&text, "eleven-bases-of-adapter"), None);
    // Twelve of them are found, and so are eleven once the minimum is lowered.
    let read = bases(&text, "twelve-bases-of-adapter");
    let index =
        find_index_of_clip_sequence(read.as_bytes(), adapter.as_bytes(), 12, MAX_ERROR_RATE);
    assert_eq!(xt_tag(index), tag(&text, "twelve-bases-of-adapter"));
    let read = bases(&text, "eleven-bases-with-a-lower-minimum");
    let index =
        find_index_of_clip_sequence(read.as_bytes(), adapter.as_bytes(), 11, MAX_ERROR_RATE);
    assert_eq!(
        xt_tag(index),
        tag(&text, "eleven-bases-with-a-lower-minimum")
    );
    // A read shorter than the minimum entirely is refused before the loop runs.
    assert_eq!(
        find_index_of_clip_sequence(b"ACGT", adapter.as_bytes(), 12, MAX_ERROR_RATE),
        NO_MATCH
    );
}

/// The error allowance is truncated, and it is computed from the overlap.
#[test]
fn the_allowance_is_truncated_from_the_overlap() {
    let text = corpus();
    let adapter = three_prime("INDEXED").expect("the indexed adapter");
    // Twelve bases at a tenth allow one mismatch.
    let read = bases(&text, "one-mismatch-in-twelve");
    let index =
        find_index_of_clip_sequence(read.as_bytes(), adapter.as_bytes(), 12, MAX_ERROR_RATE);
    assert_eq!(xt_tag(index), tag(&text, "one-mismatch-in-twelve"));
    assert!(tag(&text, "one-mismatch-in-twelve").is_some());
    // Two do not, until the rate is widened.
    let read = bases(&text, "two-mismatches-in-twelve");
    assert_eq!(
        find_index_of_clip_sequence(read.as_bytes(), adapter.as_bytes(), 12, MAX_ERROR_RATE),
        NO_MATCH
    );
    assert_eq!(tag(&text, "two-mismatches-in-twelve"), None);
    let read = bases(&text, "two-mismatches-with-a-wider-rate");
    let index = find_index_of_clip_sequence(read.as_bytes(), adapter.as_bytes(), 12, 0.2);
    assert_eq!(
        xt_tag(index),
        tag(&text, "two-mismatches-with-a-wider-rate")
    );
    // The truncation itself: twelve at a tenth is one, nine is none.
    assert_eq!((12.0 * MAX_ERROR_RATE) as usize, 1);
    assert_eq!((9.0 * MAX_ERROR_RATE) as usize, 0);
}

/// The adapter list's order decides which pair is found, and the histogram counts bases.
#[test]
fn the_list_is_tried_in_order() {
    let text = corpus();
    let read = bases(&text, "adapter-halfway");
    let (position, index) = first_matching_adapter(
        read.as_bytes(),
        &DEFAULT_ADAPTERS,
        MIN_MATCH_BASES,
        MAX_ERROR_RATE,
    )
    .expect("a match");
    assert_eq!(DEFAULT_ADAPTERS[position], "INDEXED");
    assert_eq!(xt_tag(index), tag(&text, "adapter-halfway"));
    // A list of one that does not contain it finds nothing, which is what says the list matters.
    assert_eq!(
        first_matching_adapter(
            read.as_bytes(),
            &["PAIRED_END"],
            MIN_MATCH_BASES,
            MAX_ERROR_RATE
        ),
        None
    );
    assert_eq!(tag(&text, "one-adapter-named"), None);
    // The histogram counts the bases a read would lose, which is its length minus the index.
    let histogram = field(&text, "metrics", "adapter-halfway");
    assert!(histogram.contains("30\t1"), "{histogram}");
    assert_eq!(clipped_bases(read.len(), 31), 30);
}

/// A tag the input carried does not survive a run that finds nothing.
#[test]
fn an_existing_tag_does_not_survive() {
    let text = corpus();
    // The input carried `XT:i:7` and the output carries no tag at all.
    assert_eq!(records(&text, "sam", "an-existing-tag")[0].1, Some(7));
    assert_eq!(tag(&text, "an-existing-tag"), None);
    // The constant says what the golden says, rather than being asserted for its own sake.
    assert_eq!(
        AN_EXISTING_TAG_SURVIVES,
        tag(&text, "an-existing-tag") == records(&text, "sam", "an-existing-tag")[0].1
    );
    // Where a run that DOES find an adapter overwrites it with its own answer.
    assert_eq!(
        records(&text, "sam", "an-existing-tag-and-an-adapter")[0].1,
        Some(7)
    );
    assert_eq!(tag(&text, "an-existing-tag-and-an-adapter"), Some(31));
    // A pair is marked as a PAIR and not read by read: one read carrying an adapter and one not
    // gives TWO tags, both at the position the match was found at, because the paired path finds
    // one index and writes it onto both reads.
    let pair = records(&text, "marked", "a-pair-with-one-adapter");
    assert_eq!(pair.len(), 2);
    assert_eq!(pair.iter().filter(|(_, tag)| tag.is_some()).count(), 2);
    assert_eq!(pair[0].1, pair[1].1);
    // The read that carries no adapter is the second one, and its bases say so.
    assert!(!pair[1].0.contains("AGATCGGAAGAGC"));
    // The pair whose reads both carry one is marked the same way.
    let pair = records(&text, "marked", "a-pair");
    assert_eq!(pair.iter().filter(|(_, tag)| tag.is_some()).count(), 2);
    assert_eq!(pair[0].1, pair[1].1);
}
