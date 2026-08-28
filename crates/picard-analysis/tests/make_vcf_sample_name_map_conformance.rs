//! Conformance for `MakeVcfSampleNameMap` against Picard 3.4.0.
//!
//! Each case carries the map file the tool wrote for a list of inputs. The inputs are named
//! relative to the working directory, which is what makes the line order reproducible: the map is
//! keyed by the path string and its order is that key's hash.
//!
//! # What this suite is for
//!
//!  * **the line being the path first and the name second**;
//!  * **the order being the map's and not the arguments', so three inputs named forwards and
//!    backwards produce the same file**;
//!  * **the key being the path string, so `a.vcf` and `./a.vcf` are two entries**;
//!  * **the same string twice being one entry**;
//!  * **two paths naming one sample both being kept**;
//!  * **an input with no sample or with two being refused by a message naming the count**;
//!  * **a bad input stopping the run, so nothing at all is written**;
//!  * **and the file always ending on a newline.**

use std::io::Read;

use picard_analysis::make_vcf_sample_name_map::{
    build, hash_map_order, java_string_hash, line, render, table_size_for,
    wrong_sample_count_message, Entry,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("make_vcf_sample_name_map.txt.gz");
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

/// The sample each fixture VCF names, which the dump's own file names encode.
fn sample_of(path: &str) -> &'static str {
    match path.trim_start_matches("./") {
        "a.vcf" | "d.vcf" => "SAMPLE_ONE",
        "b.vcf" => "SAMPLE_TWO",
        "c.vcf" => "SAMPLE_THREE",
        other => panic!("{other} is not a single-sample fixture"),
    }
}

fn entries(paths: &[&str]) -> Vec<Entry> {
    paths
        .iter()
        .map(|path| Entry {
            path: (*path).to_string(),
            sample: sample_of(path).to_string(),
        })
        .collect()
}

/// The cases that write a file, with the inputs each was given.
const WRITTEN: &[(&str, &[&str])] = &[
    ("one-input", &["a.vcf"]),
    ("three-inputs", &["a.vcf", "b.vcf", "c.vcf"]),
    ("three-inputs-reversed", &["c.vcf", "b.vcf", "a.vcf"]),
    ("same-path-twice", &["a.vcf", "a.vcf"]),
    ("same-sample-two-paths", &["a.vcf", "d.vcf"]),
    ("unnormalised-path", &["a.vcf", "./a.vcf"]),
];

/// Every case that writes a file writes exactly what the port renders.
#[test]
fn every_case_writes_the_same_file() {
    let text = corpus();
    for (case, paths) in WRITTEN {
        let expected = field(&text, "out", case).unwrap_or_else(|| panic!("{case}"));
        assert_eq!(render(&build(&entries(paths))), expected, "{case}");
    }
}

/// The line is the path first and the name second.
#[test]
fn the_line_is_the_path_then_the_name() {
    let text = corpus();
    assert_eq!(line("a.vcf", "SAMPLE_ONE"), "a.vcf\tSAMPLE_ONE");
    let one = field(&text, "out", "one-input").expect("one-input");
    assert_eq!(one, "a.vcf\tSAMPLE_ONE\n");
    assert!(one.starts_with("a.vcf"));
    assert!(!one.starts_with("SAMPLE_ONE"));
}

/// The order is the map's and not the arguments': the same three inputs named backwards produce
/// the same file, byte for byte.
#[test]
fn the_order_does_not_follow_the_arguments() {
    let text = corpus();
    let forwards = field(&text, "out", "three-inputs").expect("three-inputs");
    let backwards = field(&text, "out", "three-inputs-reversed").expect("three-inputs-reversed");
    assert_eq!(forwards, backwards);
    // And it is not the arguments' order in the other case either: d.vcf comes before a.vcf.
    let two = field(&text, "out", "same-sample-two-paths").expect("same-sample-two-paths");
    assert_eq!(two, "d.vcf\tSAMPLE_ONE\na.vcf\tSAMPLE_ONE\n");
    assert_eq!(
        hash_map_order(&["a.vcf".to_string(), "d.vcf".to_string()], 2),
        vec!["d.vcf".to_string(), "a.vcf".to_string()]
    );
}

/// The key is the path STRING, so one file named two ways is two entries while one string twice
/// is one.
#[test]
fn the_key_is_the_path_string() {
    let text = corpus();
    let twice = field(&text, "out", "same-path-twice").expect("same-path-twice");
    assert_eq!(twice.lines().count(), 1);
    let dotted = field(&text, "out", "unnormalised-path").expect("unnormalised-path");
    assert_eq!(dotted.lines().count(), 2);
    assert_eq!(dotted, "a.vcf\tSAMPLE_ONE\n./a.vcf\tSAMPLE_ONE\n");
}

/// Two paths naming one sample are both kept: the duplicate is only warned about.
#[test]
fn two_paths_naming_one_sample_are_both_kept() {
    let text = corpus();
    let two = field(&text, "out", "same-sample-two-paths").expect("same-sample-two-paths");
    assert_eq!(two.lines().count(), 2);
    assert!(two.lines().all(|line| line.ends_with("SAMPLE_ONE")));
}

/// An input that does not name exactly one sample is refused by a message carrying the count, and
/// a bad input stops the run so nothing at all is written.
#[test]
fn an_input_without_one_sample_is_refused() {
    let text = corpus();
    for (case, path, count) in [
        ("no-sample", "empty.vcf", 0usize),
        ("two-samples", "pair.vcf", 2),
        ("good-then-bad", "pair.vcf", 2),
    ] {
        let error = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        assert_eq!(
            error,
            format!(
                "picard.PicardException:{}",
                wrong_sample_count_message(path, count)
            ),
            "{case}"
        );
    }
    // The good input that came first left no file behind.
    assert_eq!(field(&text, "out", "good-then-bad"), None);
}

/// The file always ends on a newline, `Files.write` writing one after every line.
#[test]
fn the_file_ends_on_a_newline() {
    let text = corpus();
    for (case, _) in WRITTEN {
        let written = field(&text, "out", case).unwrap_or_else(|| panic!("{case}"));
        assert!(written.ends_with('\n'), "{case}");
    }
    assert_eq!(render(&[]), "");
}

/// The order the map iterates in is the spread hash modulo a table whose size comes from the
/// INPUT count, which is why two keys and three keys are not laid out the same way.
#[test]
fn the_order_is_the_hash_modulo_the_table() {
    assert_eq!(java_string_hash(""), 0);
    assert_eq!(java_string_hash("a"), 97);
    assert_eq!(java_string_hash("ab"), 97 * 31 + 98);
    assert_eq!(table_size_for(1), 1);
    assert_eq!(table_size_for(2), 2);
    assert_eq!(table_size_for(3), 4);
    assert_eq!(table_size_for(4), 4);
    assert_eq!(table_size_for(5), 8);
    // Three inputs start in a table of four and are never resized.
    let three = [
        "a.vcf".to_string(),
        "b.vcf".to_string(),
        "c.vcf".to_string(),
    ];
    assert_eq!(hash_map_order(&three, 3), three.to_vec());
    // Two start in a table of two and the second insertion doubles it.
    assert_eq!(
        hash_map_order(&["a.vcf".to_string(), "./a.vcf".to_string()], 2),
        vec!["a.vcf".to_string(), "./a.vcf".to_string()]
    );
}
