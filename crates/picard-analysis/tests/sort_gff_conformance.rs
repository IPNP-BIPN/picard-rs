//! Conformance for `SortGff` against Picard 3.4.0.
//!
//! Each case carries the input GFF and the sorted file the tool wrote for it. The port sorts the
//! same features and must reach the same order.
//!
//! # What this suite is for
//!
//!  * **contigs sorting lexicographically without a dictionary**;
//!  * **and by the dictionary's own order with one**;
//!  * **a contig the dictionary does not name sorting first**;
//!  * **the order within a contig being the start alone, stably**;
//!  * **the version directive being the codec's own**;
//!  * **the record count changing nothing about the output**;
//!  * **and a file with no feature being refused like one that is not GFF.**

use std::io::Read;

use picard_analysis::sort_gff::{
    cannot_decode_message, compare, sequence_index, sort, Feature, GFF_VERSION_DIRECTIVE,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("sort_gff.txt.gz");
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

fn field(corpus: &str, kind: &str, name: &str) -> Option<String> {
    corpus
        .lines()
        .find(|line| line.starts_with(&format!("{kind}\t{name}\t")))
        .map(|line| unescape(&line[format!("{kind}\t{name}\t").len()..]))
}

/// The features of one of the corpus's files, in the order they appear.
fn features(text: &str) -> Vec<Feature> {
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .enumerate()
        .map(|(index, line)| {
            let columns: Vec<&str> = line.split('\t').collect();
            Feature {
                contig: columns[0].to_string(),
                start: columns[3].parse().expect("a start"),
                end: columns[4].parse().expect("an end"),
                index,
            }
        })
        .collect()
}

/// The IDs of a file's features, in order.
fn ids(text: &str) -> Vec<String> {
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            line.split('\t')
                .nth(8)
                .expect("attributes")
                .split(';')
                .next()
                .expect("an ID")
                .to_string()
        })
        .collect()
}

/// The order of a feature is its contig and its start, and nothing else.
fn placed(features: &[Feature]) -> Vec<String> {
    features
        .iter()
        .map(|feature| format!("{}:{}", feature.contig, feature.start))
        .collect()
}

fn output_placed(corpus: &str, name: &str) -> Vec<String> {
    placed(&features(&field(corpus, "out", name).expect("its output")))
}

/// `chr10` before `chr2`.
#[test]
fn the_contigs_sort_lexicographically_without_a_dictionary() {
    let corpus = corpus();
    let input = features(&field(&corpus, "gff", "mixed").expect("the input"));
    let sorted = sort(&input, None);
    assert_eq!(placed(&sorted), output_placed(&corpus, "lexicographic"));
    assert_eq!(
        placed(&sorted),
        vec!["chr1:100", "chr1:100", "chr1:900", "chr10:100", "chr2:500"]
    );
    // The record count changes nothing about the output.
    assert_eq!(
        output_placed(&corpus, "spill-to-disk"),
        output_placed(&corpus, "lexicographic")
    );
}

/// And by its order rather than by the name.
#[test]
fn a_dictionary_puts_the_contigs_in_its_own_order() {
    let corpus = corpus();
    let input = features(&field(&corpus, "gff", "mixed").expect("the input"));
    let dictionary: Vec<String> = ["chr1", "chr2", "chr10"]
        .iter()
        .map(|c| c.to_string())
        .collect();
    let sorted = sort(&input, Some(&dictionary));
    assert_eq!(placed(&sorted), output_placed(&corpus, "dictionary-order"));
    assert_eq!(
        placed(&sorted),
        vec!["chr1:100", "chr1:100", "chr1:900", "chr2:500", "chr10:100"]
    );
    assert_eq!(sequence_index(&dictionary, "chr1"), 0);
    assert_eq!(sequence_index(&dictionary, "chr10"), 2);
}

/// It gets index -1, which is before every contig the dictionary does name.
#[test]
fn a_contig_the_dictionary_lacks_sorts_first() {
    let corpus = corpus();
    let input = features(&field(&corpus, "gff", "mixed").expect("the input"));
    let partial: Vec<String> = ["chr2", "chr10"].iter().map(|c| c.to_string()).collect();
    assert_eq!(sequence_index(&partial, "chr1"), -1);
    let sorted = sort(&input, Some(&partial));
    assert_eq!(
        placed(&sorted),
        output_placed(&corpus, "dictionary-partial")
    );
    // chr1 is not in the dictionary and comes first all the same.
    assert!(placed(&sorted)[0].starts_with("chr1:"));
    assert_eq!(
        placed(&sorted),
        vec!["chr1:100", "chr1:100", "chr1:900", "chr2:500", "chr10:100"]
    );
}

/// The start alone, and stably.
#[test]
fn the_order_within_a_contig_is_the_start_alone() {
    let corpus = corpus();
    // The fixture has two features starting at 100 on chr1, one ending at 300 and one at 200.
    let written = ids(&field(&corpus, "out", "lexicographic").expect("its output"));
    assert_eq!(written[0], "ID=c2");
    assert_eq!(written[1], "ID=c3");
    // The one that ends LATER comes first, because it was read first: the ends never meet.
    let a = Feature {
        contig: "chr1".to_string(),
        start: 100,
        end: 300,
        index: 0,
    };
    let b = Feature {
        end: 200,
        index: 1,
        ..a.clone()
    };
    assert_eq!(compare(&a, &b, None), std::cmp::Ordering::Equal);
    assert_eq!(
        placed(&sort(&[a.clone(), b.clone()], None)),
        placed(&[a, b])
    );
    // A parent written after its child is put wherever its own coordinates say.
    let child_first = ids(&field(&corpus, "out", "child-before-parent").expect("its output"));
    assert_eq!(child_first, vec!["ID=g1", "ID=e1"]);
}

/// The codec writes its own version, not the input's.
#[test]
fn the_version_directive_is_the_codecs_own() {
    let corpus = corpus();
    let input = field(&corpus, "gff", "mixed").expect("the input");
    assert!(input.starts_with("##gff-version 3.1.26"), "{input}");
    let output = field(&corpus, "out", "lexicographic").expect("its output");
    assert!(output.starts_with(GFF_VERSION_DIRECTIVE), "{output}");
    assert_eq!(GFF_VERSION_DIRECTIVE, "##gff-version 3.1.25");
    // The comment line is carried over.
    assert!(
        output.contains("#a comment the sorter carries over"),
        "{output}"
    );
}

/// Both get the same refusal, naming the input.
#[test]
fn a_file_with_no_feature_is_refused_like_one_that_is_not_gff() {
    let corpus = corpus();
    for name in ["no-features", "not-a-gff"] {
        let line = corpus
            .lines()
            .find(|line| line.starts_with(&format!("error\t{name}\t")))
            .unwrap_or_else(|| panic!("the corpus carries error/{name}"));
        let message = &line[format!("error\t{name}\t").len()..];
        assert_eq!(
            message,
            format!(
                "java.lang.IllegalArgumentException:{}",
                cannot_decode_message("<dir>/in.gff3")
            ),
            "{name}"
        );
    }
    // The empty file did carry a directive and a comment, so it is not blank.
    let empty = field(&corpus, "gff", "empty").expect("the input");
    assert!(empty.contains("##gff-version"), "{empty}");
    assert!(features(&empty).is_empty());
}

/// The attributes survive the round trip through the codec.
#[test]
fn the_attributes_are_rewritten_by_the_codec() {
    let corpus = corpus();
    let output = field(&corpus, "out", "escaped-attributes").expect("its output");
    let line = output
        .lines()
        .find(|line| !line.starts_with('#') && !line.is_empty())
        .expect("its feature");
    let attributes = line.split('\t').nth(8).expect("attributes");
    // The escaped comma stays escaped and the space stays a space.
    assert!(attributes.contains("Note=a%2Cb"), "{attributes}");
    assert!(attributes.contains("with space"), "{attributes}");
    assert!(attributes.starts_with("ID=g1"), "{attributes}");
}
