//! Conformance for `CollectSequencingArtifactMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/artifactmetrics-conformance`. Each case carries the input as SAM, the five
//! file names the prefix produced, the row counts of each file, the summary tables in full and the
//! detail rows that counted anything.
//!
//! # What this suite is for
//!
//!  * **the output argument being a prefix for five files**;
//!  * **`--FILE_EXTENSION` appending to all five**;
//!  * **the detail files holding a row per substitution per context**;
//!  * **`--CONTEXTS_TO_PRINT` cutting the detail files and not the counting**;
//!  * **a pre-adapter artifact following the read and a bait-bias one the strand**;
//!  * **`--TANDEM_READS` swapping which end counts as which**;
//!  * **the rates being floored so the Q is finite**;
//!  * **the base-quality floor reading the `OQ` tag unless it is turned off**;
//!  * **and each whole-read filter dropping its read, three of them with an argument that puts it
//!    back.**

use std::io::Read;

use picard_analysis::collect_sequencing_artifact_metrics::{
    bait_bias, bait_bias_error_rates, detail_rows, file_names, phred_from_error_probability,
    pre_adapter, pre_adapter_error_rate, transitions, Alignment, MIN_ERROR,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_sequencing_artifact_metrics.txt.gz");
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

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// One `<kind>=<payload>` row, of which a case has one per file.
fn keyed(text: &str, kind: &str, case: &str, file: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t{file}=");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn rows(text: &str, case: &str, file: &str) -> usize {
    keyed(text, "rows", case, file)
        .unwrap_or_else(|| panic!("rows/{case}/{file}"))
        .parse()
        .expect("a count")
}

/// The detail rows that counted anything, as maps.
fn detail(text: &str, case: &str, file: &str) -> Vec<std::collections::HashMap<String, String>> {
    let body = keyed(text, "detail", case, file).unwrap_or_else(|| panic!("detail/{case}/{file}"));
    let mut lines = body.lines().filter(|line| !line.is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

/// The one pre-adapter row a case's substitution lands in.
///
/// The context is part of the key: the fixture covers many G sites with the reference base, so
/// several `G>T` rows carry counts and only the one whose context is `CGT` carries the
/// substitution.
fn artifact(
    text: &str,
    case: &str,
    reference: &str,
    alternate: &str,
) -> Option<std::collections::HashMap<String, String>> {
    let context = if reference == "G" { "CGT" } else { "ACG" };
    detail(text, case, "pre_adapter_detail_metrics")
        .into_iter()
        .find(|row| {
            row["REF_BASE"] == reference
                && row["ALT_BASE"] == alternate
                && row["CONTEXT"] == context
        })
}

/// The output argument is a prefix for five files, and the extension appends to all five.
#[test]
fn the_output_is_a_prefix_for_five_files() {
    let text = corpus();
    let written = field(&text, "files", "plain").expect("the file names");
    let names: std::collections::BTreeSet<String> =
        written.split(',').map(str::to_string).collect();
    assert_eq!(names.len(), 5);
    assert_eq!(
        file_names("out", None)
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>(),
        names
    );
    // The extension goes on the end of each, rather than replacing anything.
    let extended = field(&text, "files", "file-extension").expect("the file names");
    assert!(extended.split(',').all(|name| name.ends_with(".txt")));
    assert_eq!(
        file_names("out", Some(".txt"))
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>(),
        extended
            .split(',')
            .map(str::to_string)
            .collect::<std::collections::BTreeSet<_>>()
    );
}

/// The detail files hold a row per substitution per context; the summaries hold twelve.
#[test]
fn the_detail_files_are_a_row_per_substitution_per_context() {
    let text = corpus();
    assert_eq!(transitions().len(), 12);
    assert_eq!(detail_rows(1), 192);
    assert_eq!(detail_rows(0), 12);
    for file in ["pre_adapter_detail_metrics", "bait_bias_detail_metrics"] {
        assert_eq!(
            rows(&text, "alt-read-one-forward", file),
            detail_rows(1),
            "{file}"
        );
        assert_eq!(
            rows(&text, "context-size-zero", file),
            detail_rows(0),
            "{file}"
        );
    }
    for file in ["pre_adapter_summary_metrics", "bait_bias_summary_metrics"] {
        assert_eq!(rows(&text, "alt-read-one-forward", file), 12, "{file}");
        assert_eq!(rows(&text, "context-size-zero", file), 12, "{file}");
    }
}

/// `--CONTEXTS_TO_PRINT` cuts the detail files and the error summary, and not the counting.
#[test]
fn naming_a_context_cuts_the_files_and_not_the_counting() {
    let text = corpus();
    assert_eq!(
        rows(&text, "contexts-to-print", "pre_adapter_detail_metrics"),
        3
    );
    assert_eq!(
        rows(&text, "contexts-to-print", "bait_bias_detail_metrics"),
        3
    );
    assert_eq!(rows(&text, "contexts-to-print", "error_summary_metrics"), 3);
    assert_eq!(
        rows(&text, "alt-read-one-forward", "error_summary_metrics"),
        6
    );
    // The summary rows are the same table, which is what says the counting was untouched.
    for file in ["pre_adapter_summary_metrics", "bait_bias_summary_metrics"] {
        assert_eq!(
            keyed(&text, "summary", "contexts-to-print", file),
            keyed(&text, "summary", "alt-read-one-forward", file),
            "{file}"
        );
    }
}

/// A pre-adapter artifact follows the read; the end and the strand decide together.
#[test]
fn the_end_and_the_strand_decide_the_side() {
    let text = corpus();
    let forward = artifact(&text, "alt-read-one-forward", "G", "T").expect("a G>T row");
    assert_eq!(forward["PRO_ALT_BASES"], "1");
    assert_eq!(forward["CON_ALT_BASES"], "0");
    let read_two = artifact(&text, "alt-read-two-forward", "G", "T").expect("a G>T row");
    assert_eq!(read_two["PRO_ALT_BASES"], "0");
    assert_eq!(read_two["CON_ALT_BASES"], "1");
    let reverse = artifact(&text, "alt-read-one-reverse", "G", "T").expect("a G>T row");
    assert_eq!(reverse["PRO_ALT_BASES"], "0");
    assert_eq!(reverse["CON_ALT_BASES"], "1");
    // Which is what the port's own fold gives for the same alignments.
    let one_forward = Alignment {
        r1_pos: 1,
        ..Alignment::default()
    };
    let counts = pre_adapter(
        Alignment::default(),
        one_forward,
        Alignment::default(),
        Alignment::default(),
        false,
    );
    assert_eq!((counts.pro_alt, counts.con_alt), (1, 0));
    let two_forward = Alignment {
        r2_pos: 1,
        ..Alignment::default()
    };
    let counts = pre_adapter(
        Alignment::default(),
        two_forward,
        Alignment::default(),
        Alignment::default(),
        false,
    );
    assert_eq!((counts.pro_alt, counts.con_alt), (0, 1));
}

/// `--TANDEM_READS` swaps read two's half of every sum.
#[test]
fn tandem_reads_swap_read_twos_half() {
    let text = corpus();
    // A pair whose ends are on opposite strands lands on the same side under either convention,
    // so the swap is invisible there.
    let opposite = artifact(&text, "alt-read-one-forward", "G", "T").expect("a G>T row");
    let opposite_tandem = artifact(&text, "tandem-reads", "G", "T").expect("a G>T row");
    assert_eq!(opposite, opposite_tandem);
    // It is plain on a substitution carried by read two on the forward strand: that end changes
    // side, and the Q with it.
    let read_two = artifact(&text, "alt-read-two-forward", "G", "T").expect("a G>T row");
    let read_two_tandem = artifact(&text, "tandem-read-two-forward", "G", "T").expect("a G>T row");
    assert_eq!(
        (&read_two["PRO_ALT_BASES"], &read_two["CON_ALT_BASES"]),
        (&"0".to_string(), &"1".to_string())
    );
    assert_eq!(
        (
            &read_two_tandem["PRO_ALT_BASES"],
            &read_two_tandem["CON_ALT_BASES"]
        ),
        (&"1".to_string(), &"0".to_string())
    );
    assert_eq!(read_two["QSCORE"], "100");
    assert_eq!(read_two_tandem["QSCORE"], "11");
    // Which is the port's own fold of the same alignment.
    let two_forward = Alignment {
        r2_pos: 1,
        ..Alignment::default()
    };
    let plain = pre_adapter(
        Alignment::default(),
        two_forward,
        Alignment::default(),
        Alignment::default(),
        false,
    );
    let tandem = pre_adapter(
        Alignment::default(),
        two_forward,
        Alignment::default(),
        Alignment::default(),
        true,
    );
    assert_eq!((plain.pro_alt, plain.con_alt), (0, 1));
    assert_eq!((tandem.pro_alt, tandem.con_alt), (1, 0));
}

/// The bait-bias counters sum both ends and both strands away.
#[test]
fn the_bait_bias_counters_sum_the_pair_away() {
    let text = corpus();
    for case in [
        "alt-read-one-forward",
        "alt-read-two-forward",
        "alt-read-one-reverse",
    ] {
        let row = detail(&text, case, "bait_bias_detail_metrics")
            .into_iter()
            .find(|row| row["REF_BASE"] == "G" && row["ALT_BASE"] == "T" && row["CONTEXT"] == "CGT")
            .unwrap_or_else(|| panic!("a G>T row in {case}"));
        let forward: i64 = row["FWD_CXT_ALT_BASES"].parse().expect("a count");
        let reverse: i64 = row["REV_CXT_ALT_BASES"].parse().expect("a count");
        assert_eq!(forward + reverse, 1, "{case}");
    }
    let counts = bait_bias(
        Alignment::default(),
        Alignment {
            r2_pos: 1,
            ..Alignment::default()
        },
        Alignment::default(),
        Alignment::default(),
    );
    assert_eq!((counts.fwd_alt, counts.rev_alt), (1, 0));
}

/// The rates are floored, so a row nothing was seen for still has a finite Q of a hundred.
#[test]
fn the_rates_are_floored_and_the_q_is_an_integer() {
    let text = corpus();
    let row = artifact(&text, "alt-read-one-forward", "C", "A").expect("a C>A row");
    assert_eq!(row["QSCORE"], "100");
    assert_eq!(phred_from_error_probability(MIN_ERROR), 100);
    let counts = pre_adapter(
        Alignment::default(),
        Alignment::default(),
        Alignment::default(),
        Alignment::default(),
        false,
    );
    assert_eq!(pre_adapter_error_rate(&counts), MIN_ERROR);
    let (forward, reverse, difference) = bait_bias_error_rates(&bait_bias(
        Alignment::default(),
        Alignment::default(),
        Alignment::default(),
        Alignment::default(),
    ));
    assert_eq!(
        (forward, reverse, difference),
        (MIN_ERROR, MIN_ERROR, MIN_ERROR)
    );
    // And the Q column is rounded to an integer, which is why it reads 11 and not 11.4.
    let seen = artifact(&text, "alt-read-one-forward", "G", "T").expect("a G>T row");
    assert!(!seen["QSCORE"].contains('.'));
    let rate: f64 = seen["ERROR_RATE"].parse().expect("a rate");
    assert_eq!(
        seen["QSCORE"].parse::<i32>().expect("a q"),
        phred_from_error_probability(rate)
    );
}

/// The base-quality floor reads the `OQ` tag unless it is turned off.
#[test]
fn the_quality_floor_reads_the_oq_tag() {
    let text = corpus();
    // The row is still there, with the reference counts the clean end left: what the floor takes
    // is the alternate observation.
    for case in ["low-base-quality", "original-qualities"] {
        let row = artifact(&text, case, "G", "T").unwrap_or_else(|| panic!("{case}"));
        assert_eq!(row["PRO_ALT_BASES"], "0", "{case}");
        assert_eq!(row["CON_ALT_BASES"], "0", "{case}");
    }
    let ignored = artifact(&text, "original-qualities-ignored", "G", "T").expect("a G>T row");
    assert_eq!(ignored["PRO_ALT_BASES"], "1");
}

/// Each whole-read filter drops its read; three of them have an argument that puts it back.
#[test]
fn each_filter_drops_its_read() {
    let text = corpus();
    for case in [
        "low-mapping-quality",
        "insert-too-small",
        "insert-too-large",
        "duplicate",
        "secondary",
        "fails-vendor",
        "unpaired",
    ] {
        let row = artifact(&text, case, "G", "T").unwrap_or_else(|| panic!("{case}"));
        assert_eq!(row["PRO_ALT_BASES"], "0", "{case}");
        assert_eq!(row["CON_ALT_BASES"], "0", "{case}");
    }
    for case in [
        "duplicate-included",
        "fails-vendor-included",
        "unpaired-included",
    ] {
        let row = artifact(&text, case, "G", "T").unwrap_or_else(|| panic!("{case}"));
        assert_eq!(row["PRO_ALT_BASES"], "1", "{case}");
    }
    // An unpaired read counts as read one, which is what makes it a propitious observation.
    let mut counts = Alignment::default();
    counts.count(false, false, false);
    assert_eq!(counts.r1_pos, 1);
}
