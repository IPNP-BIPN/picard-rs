//! Conformance for `CrosscheckReadGroupFingerprints` against Picard 3.4.0.
//!
//! Golden from `tools/checkfingerprint-conformance/`: the deprecated wrapper's roll-ups and the
//! two command lines it refuses.
//!
//! # What this suite is for
//!
//!  * **the roll-up being two booleans of its own** rather than the parent's `--CROSSCHECK_BY`;
//!  * **a roll-up moving the output**, because what it writes is a matrix and not a table;
//!  * **`--CROSSCHECK_BY` being refused outright**, in the tool's own sentence;
//!  * **and the two expectations being mutually exclusive**, which is the parser's refusal and
//!    not the tool's.

use std::io::Read;

use picard_analysis::crosscheck_fingerprints::DataType;
use picard_analysis::crosscheck_read_group_fingerprints::{
    destination, expect_all_groups_to_match, mutex_refusal, refusals, Options,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/fingerprint_metrics.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

/// The roll-up is asked for by a boolean, and it decides both what is compared and where it goes.
#[test]
fn the_roll_up_moves_the_output() {
    let by_read_group = destination("metrics.txt", &Options::default());
    assert_eq!(by_read_group.crosscheck_by, DataType::ReadGroup);
    assert_eq!(by_read_group.output, "metrics.txt");
    assert_eq!(by_read_group.matrix_output, None);

    for (options, expected) in [
        (
            Options {
                crosscheck_samples: true,
                ..Options::default()
            },
            DataType::Sample,
        ),
        (
            Options {
                crosscheck_libraries: true,
                ..Options::default()
            },
            DataType::Library,
        ),
    ] {
        let rolled = destination("metrics.txt", &options);
        assert_eq!(rolled.crosscheck_by, expected);
        // The named file becomes the matrix, and the table is thrown away.
        assert_eq!(rolled.matrix_output.as_deref(), Some("metrics.txt"));
        assert_eq!(rolled.output, "/dev/null");
    }

    // Both at once is the library one: it is checked first.
    let both = destination(
        "metrics.txt",
        &Options {
            crosscheck_samples: true,
            crosscheck_libraries: true,
            ..Options::default()
        },
    );
    assert_eq!(both.crosscheck_by, DataType::Library);

    // And what each roll-up writes is a matrix of the level it rolled up to, which is what the
    // golden holds.
    let text = corpus();
    let samples = field(&text, "metrics", "crosscheck-rolled-up-to-samples").expect("the golden");
    assert_eq!(samples.lines().next(), Some("SAMPLE\tsample1"));
    let libraries =
        field(&text, "metrics", "crosscheck-rolled-up-to-libraries").expect("the golden");
    assert_eq!(
        libraries.lines().next(),
        Some("LIBRARY\tsample1::lib1\tsample1::lib2")
    );
}

/// The wrapper refuses the argument it sets itself, in its own words.
#[test]
fn crosscheck_by_is_refused() {
    let text = corpus();
    let refused = refusals(&Options {
        crosscheck_by: Some(DataType::Sample),
        ..Options::default()
    });
    assert_eq!(refused.len(), 1);
    assert_eq!(
        refused[0],
        field(&text, "refusal", "crosscheck-by-is-refused").expect("the golden")
    );
    assert_eq!(
        field(&text, "code", "crosscheck-by-is-refused").as_deref(),
        Some("1")
    );
    // A command line that names none of the three is not refused at all.
    assert!(refusals(&Options::default()).is_empty());
}

/// The two expectations are mutually exclusive, which the parser refuses before the tool runs.
#[test]
fn the_two_expectations_are_exclusive() {
    let text = corpus();
    assert_eq!(
        mutex_refusal(
            "EXPECT_ALL_READ_GROUPS_TO_MATCH",
            &["EXPECT_ALL_GROUPS_TO_MATCH"]
        ),
        field(&text, "refusal", "the-two-expectations-are-exclusive").expect("the golden")
    );
    assert_eq!(
        field(&text, "code", "the-two-expectations-are-exclusive").as_deref(),
        Some("1")
    );

    // Set on its own, the wrapper's flag becomes the parent's.
    assert!(expect_all_groups_to_match(
        false,
        &Options {
            expect_all_read_groups_to_match: true,
            ..Options::default()
        }
    ));
    assert!(!expect_all_groups_to_match(false, &Options::default()));
}
