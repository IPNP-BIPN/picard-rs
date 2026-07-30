//! Conformance for `IntervalListTools`'s `PADDING` and `BREAK_BANDS_AT_MULTIPLES_OF`, against
//! Picard 3.4.0.
//!
//! Golden from `tools/intervallisttools-conformance/IntervalListToolsPadBreakDump.java`, produced
//! by the pinned container on real x86-64. The one this slice first produced came from a laptop and
//! was refused, which is what decision 0008 exists for.
//!
//! Both options are clamps and renames rather than arithmetic, and the golden pins both ends:
//!
//! ```text
//! chr1  5  20   +  A     padded by 10 becomes  chr1  1   30   +  A
//! chr1  95 100  +  C     padded by 10 becomes  chr1  85  100  +  C
//! chr1  5  25   +  A     broken at 10 becomes  chr1  5   9    +  A.1
//!                                              chr1  10  19   +  A.2
//!                                              chr1  20  25   +  A.3
//! ```
//!
//! Padding is clamped to `[1, contig length]`, so the interval starting at 5 does not reach -5 and
//! the one ending at 100 does not reach 110. Breaking renames every piece `A.1`, `A.2`, `A.3`, so
//! the first piece is renamed even though it starts where the original did, and the pieces are cut
//! at the band multiples rather than into equal parts: the first is five bases, the second ten, the
//! third six.
//!
//! `@PG` and `UR:` are stripped, as everywhere: the first is the command line and the second the
//! reference's path, both carrying the run's temp directory.

use std::collections::HashMap;
use std::io::Read;

use picard_analysis::interval_list_tools::{interval_list_tools, Options};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/interval_list_tools_padbreak.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars();
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

/// Drop `@PG` lines and every `UR:` field: the command line and the reference's path, both of
/// which carry the temp directory the run happened in.
fn canonicalize(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("@PG"))
        .map(|line| {
            let kept: Vec<&str> = line.split('\t').filter(|f| !f.starts_with("UR:")).collect();
            format!("{}\n", kept.join("\t"))
        })
        .collect()
}

#[derive(Default)]
struct Case {
    input1: String,
    padding: i32,
    break_bands: i32,
    output: String,
}

fn cases() -> Vec<(String, Case)> {
    let text = corpus();
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, Case> = HashMap::new();
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.splitn(3, '\t');
        let kind = fields.next().expect("a kind");
        let name = fields.next().expect("a case").to_string();
        let payload = unescape(fields.next().unwrap_or(""));
        let case = map.entry(name.clone()).or_insert_with(|| {
            order.push(name.clone());
            Case::default()
        });
        match kind {
            "input1" => case.input1 = payload,
            "padding" => case.padding = payload.parse().expect("an integer"),
            "break_bands" => case.break_bands = payload.parse().expect("an integer"),
            "output" => case.output = payload,
            other => panic!("unexpected row kind {other}"),
        }
    }
    order
        .into_iter()
        .map(|name| {
            let case = map.remove(&name).expect("a case");
            (name, case)
        })
        .collect()
}

#[test]
fn padding_and_break_bands_are_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 2, "case count");
    for (name, case) in &cases {
        let options = Options {
            padding: case.padding,
            break_bands_at_multiples_of: case.break_bands,
            ..Options::default()
        };
        let got = interval_list_tools(&[&case.input1], &[], &options).expect("tool");
        assert_eq!(canonicalize(&got), canonicalize(&case.output), "{name}");
    }
}

/// The two clamps and the rename, asserted on the golden itself so a change in the fixture cannot
/// quietly make the test above vacuous.
#[test]
fn the_golden_pins_both_clamps_and_the_rename() {
    let text = corpus();

    // Padding does not run past the start of the contig, nor past its length.
    assert!(
        text.contains(r"chr1\t1\t30\t+\tA"),
        "clamped at the contig start"
    );
    assert!(
        text.contains(r"chr1\t85\t100\t+\tC"),
        "clamped at the contig length"
    );

    // Breaking cuts at the band multiples, not into equal parts, and renames every piece.
    assert!(text.contains(r"chr1\t5\t9\t+\tA.1"));
    assert!(text.contains(r"chr1\t10\t19\t+\tA.2"));
    assert!(text.contains(r"chr1\t20\t25\t+\tA.3"));
}
