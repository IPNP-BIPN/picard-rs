//! Conformance for `CollectMultipleMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/multiplemetrics-conformance/CollectMultipleMetricsDump.java`, which ran the
//! tool eighteen times and recorded the files each run landed on, the tables it wrote and the
//! refusals it made.
//!
//! # What this suite is for
//!
//!  * **the default set being five programs of the nine the enum declares**;
//!  * **the files each program lands on, and the extension landing on half of them**;
//!  * **what a program needs before a record is read**;
//!  * **one EXTRA_ARGUMENT reaching one program and no other**;
//!  * **and the charts being R's, which is why their bytes are not a claim.**

use std::io::Read;

use picard_analysis::collect_multiple_metrics::{
    extra_argument, file_names, plan, Program, Refusal, CHART_BYTES_ARE_REPRODUCIBLE,
    DEFAULT_PROGRAMS, NO_PROGRAMS_MESSAGE, PROGRAMS,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/collect_multiple_metrics.txt.gz");
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

fn files(text: &str, case: &str) -> Vec<String> {
    field(text, "files", case)
        .unwrap_or_else(|| panic!("files/{case}"))
        .split(' ')
        .map(str::to_string)
        .collect()
}

/// What the port would write for a set of programs, sorted the way the dump sorted a directory.
fn written(programs: &[Program], extension: Option<&str>) -> Vec<String> {
    let mut names: Vec<String> = programs
        .iter()
        .flat_map(|program| file_names("m", *program, extension))
        .collect();
    names.sort();
    names
}

/// The default set is five programs, and the files it lands on are the golden's.
#[test]
fn the_default_set_is_five_of_the_nine() {
    let text = corpus();
    assert_eq!(PROGRAMS.len(), 9);
    assert_eq!(DEFAULT_PROGRAMS.len(), 5);
    assert_eq!(
        files(&text, "default-programs"),
        written(&DEFAULT_PROGRAMS, None)
    );
    // None of the five needs a reference, which is why the same run without one writes the same
    // files rather than being refused.
    assert!(!DEFAULT_PROGRAMS
        .iter()
        .any(|program| program.needs_reference_sequence()));
    assert_eq!(
        files(&text, "default-programs-without-a-reference"),
        written(&DEFAULT_PROGRAMS, None)
    );
    // A run names its own set by emptying the default one first, and an empty set is refused.
    assert_eq!(
        files(&text, "one-program"),
        written(&[Program::QualityScoreDistribution], None)
    );
    assert_eq!(
        files(&text, "two-programs"),
        written(
            &[
                Program::QualityScoreDistribution,
                Program::MeanQualityByCycle
            ],
            None
        )
    );
    let refusal = field(&text, "error", "no-programs").expect("the refusal");
    assert!(refusal.contains(NO_PROGRAMS_MESSAGE), "{refusal}");
    assert_eq!(plan(&[], true, true, &[]), Err(Refusal::NoPrograms));
}

/// The extension lands on the metrics files and not on the charts.
#[test]
fn the_extension_lands_on_half_the_files() {
    let text = corpus();
    let with = files(&text, "file-extension");
    assert_eq!(with, written(&DEFAULT_PROGRAMS, Some(".txt")));
    // Five renamed and five left alone, which is what makes this worth an assertion.
    assert_eq!(with.iter().filter(|name| name.ends_with(".txt")).count(), 5);
    assert_eq!(with.iter().filter(|name| name.ends_with(".pdf")).count(), 5);
    // And the chart's own name is unchanged between the two runs.
    for name in files(&text, "default-programs") {
        if name.ends_with(".pdf") {
            assert!(with.contains(&name), "{name}");
        }
    }
}

/// A program that needs more than the reads is refused before a record is read.
#[test]
fn a_program_that_needs_more_is_refused_by_name() {
    let text = corpus();
    let cases = [
        (
            "gc-bias-without-a-reference",
            Refusal::NeedsReferenceSequence(Program::CollectGcBiasMetrics),
        ),
        (
            "artifacts-without-a-reference",
            Refusal::NeedsReferenceSequence(Program::CollectSequencingArtifactMetrics),
        ),
        (
            "rna-seq-without-a-refflat",
            Refusal::NeedsRefflatFile(Program::RnaSeqMetrics),
        ),
    ];
    for (case, refusal) in cases {
        let written = field(&text, "error", case).unwrap_or_else(|| panic!("{case}"));
        assert!(written.ends_with(&refusal.message()), "{written}");
    }
    assert_eq!(
        plan(&[Program::CollectGcBiasMetrics], false, true, &[]),
        Err(Refusal::NeedsReferenceSequence(
            Program::CollectGcBiasMetrics
        ))
    );
    // With the reference it runs, and writes three files rather than a pair.
    assert_eq!(
        files(&text, "gc-bias"),
        written(&[Program::CollectGcBiasMetrics], None)
    );
    assert_eq!(Program::CollectGcBiasMetrics.extensions().len(), 3);
}

/// One extra argument reaches one program, and one for a program that is not running is an error.
#[test]
fn an_extra_argument_reaches_one_program() {
    let text = corpus();
    let parsed = extra_argument("CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200").expect("a value");
    assert_eq!(parsed.program, Program::CollectInsertSizeMetrics);
    assert_eq!(parsed.values, vec!["HISTOGRAM_WIDTH=200".to_string()]);
    // The new-parser spelling is two entries, the reluctant group taking as little as it can.
    let split = extra_argument("CollectInsertSizeMetrics::--HISTOGRAM_WIDTH 200").expect("a value");
    assert_eq!(
        split.values,
        vec!["--HISTOGRAM_WIDTH".to_string(), "200".to_string()]
    );
    assert_eq!(
        files(&text, "extra-argument"),
        written(&[Program::CollectInsertSizeMetrics], None)
    );
    // The leftover check runs after the loop, so it is the program that is named and not the
    // argument that the message carries.
    let leftover = field(&text, "error", "extra-argument-for-another-program").expect("a refusal");
    assert!(
        leftover.contains(
            &Refusal::ExtraArgumentNotRequested(Program::CollectInsertSizeMetrics).message()
        ),
        "{leftover}"
    );
    assert_eq!(
        plan(
            &[Program::QualityScoreDistribution],
            true,
            true,
            &["CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200"]
        ),
        Err(Refusal::ExtraArgumentNotRequested(
            Program::CollectInsertSizeMetrics
        ))
    );
    // A value with no `::` and one whose program does not resolve are two different refusals.
    let malformed = field(&text, "error", "extra-argument-malformed").expect("a refusal");
    assert!(
        malformed
            .contains(&Refusal::ExtraArgumentMalformed("HISTOGRAM_WIDTH=200".into()).message()),
        "{malformed}"
    );
    let unknown = field(&text, "error", "extra-argument-unknown-program").expect("a refusal");
    assert!(unknown.contains("NoSuchProgram"), "{unknown}");
    assert_eq!(
        extra_argument("NoSuchProgram::X=1"),
        Err(Refusal::ExtraArgumentUnknownProgram("NoSuchProgram".into()))
    );
    // And an argument the pass owns is accepted and ignored: the run writes its files.
    assert_eq!(
        files(&text, "extra-argument-the-pass-owns"),
        written(&[Program::QualityScoreDistribution], None)
    );
}

/// The charts are R's output, which two runs of one fixture do not agree on.
#[test]
fn the_charts_are_not_bytes_a_golden_can_hold() {
    let text = corpus();
    let stability: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("chart-stability\t"))
        .collect();
    assert_eq!(stability.len(), 5);
    for line in &stability {
        assert!(line.ends_with("=differs"), "{line}");
    }
    // The constant says what the golden says, rather than being asserted for its own sake.
    assert_eq!(
        CHART_BYTES_ARE_REPRODUCIBLE,
        stability.iter().all(|line| line.ends_with("=same"))
    );
    // The same fixture run twice lands on the same file NAMES, which is what the port claims.
    assert_eq!(
        files(&text, "default-programs"),
        files(&text, "default-programs-again")
    );
}

/// The dispatcher changes no number: one program alone writes what the pass writes for it.
#[test]
fn one_pass_is_the_standalone_tools_numbers() {
    let text = corpus();
    let alone = field(&text, "metrics", "standalone/.alignment_summary_metrics").expect("alone");
    let together = field(
        &text,
        "metrics",
        "default-programs/.alignment_summary_metrics",
    )
    .expect("together");
    assert_eq!(alone, together);
    // And truncating the pass truncates it for every program at once, which is what makes it one
    // pass: with four records the insert-size program writes nothing at all.
    let truncated = files(&text, "stop-after");
    assert!(!truncated.iter().any(|name| name.contains("insert_size")));
    assert!(truncated
        .iter()
        .any(|name| name.contains("quality_by_cycle")));
}
