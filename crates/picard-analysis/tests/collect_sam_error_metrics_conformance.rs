//! Conformance for `CollectSamErrorMetrics` against Picard 3.4.0.
//!
//! Golden from `tools/samerrormetrics-conformance/`: thirteen runs over reads on a reference of
//! repeating `ACGT`, plus the refusal a list of metrics is refused with.
//!
//! # What this suite is for
//!
//!  * **the rate being Bayesian rather than a ratio**, so no errors still reports a quality;
//!  * **the two thresholds dropping observations**, one by base and one by read;
//!  * **the VCF being subtractive**, taking bases out of the denominator;
//!  * **a stratifier splitting the rows and nothing else**;
//!  * **a deletion counted once for the read rather than once per locus**;
//!  * **the indel table's `Q_SCORE` staying zero**, because that metric derives its own fields
//!    and not the one it inherits;
//!  * **and a run of poorly mapped reads writing a table with no header at all.**

use std::io::Read as _;

use picard_analysis::collect_sam_error_metrics::{
    aggregation_suffix, collect, gc_content, pileup, prior_error, processed_loci, q_score, render,
    suffixes, Calculator, Options, Read, Refusal, Stratifier, Table, DEFAULT_ERROR_METRICS,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/sam_error_metrics.txt.gz");
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

/// The contig: `ACGT` repeating, so a mismatch is any other base.
fn reference() -> Vec<u8> {
    (0..600).map(|index| b"ACGT"[index % 4]).collect()
}

/// The reference's own bases over a window, which a read copies before it is edited.
fn window(start: usize, length: usize) -> Vec<u8> {
    (0..length)
        .map(|offset| b"ACGT"[(start - 1 + offset) % 4])
        .collect()
}

/// Eight reads over one window, the given number of them carrying a mismatch at one position.
fn reads(mismatches: usize, mismatch_at: usize, quality: u8, mapping_quality: u8) -> Vec<Read> {
    (0..8)
        .map(|index| {
            let mut bases = window(101, 20);
            let mut qualities = vec![40u8; 20];
            if index < mismatches {
                let offset = mismatch_at - 101;
                bases[offset] = if bases[offset] == b'A' { b'C' } else { b'A' };
                qualities[offset] = quality;
            }
            Read {
                name: format!("r{index}"),
                start: 101,
                bases,
                qualities,
                flags: 0,
                mate_start: 0,
                cigar: vec![(20, 'M')],
                mapping_quality,
            }
        })
        .collect()
}

/// Run one case the way the tool runs it, and render the one table it writes.
fn run(
    reads: &[Read],
    calculator: Calculator,
    stratifier: Stratifier,
    options: &Options,
) -> String {
    let loci = processed_loci(pileup(reads, options), options);
    render(&collect(
        reads,
        &reference(),
        1,
        &loci,
        calculator,
        stratifier,
        options,
    ))
}

/// The prior is a pseudo-count, so a file with no errors reports a finite quality that moves when
/// the prior does.
#[test]
fn the_rate_is_bayesian() {
    let text = corpus();
    let plain = reads(2, 105, 40, 60);
    assert_eq!(
        run(
            &plain,
            Calculator::Error,
            Stratifier::All,
            &Options::default()
        ),
        field(&text, "metrics", "two-mismatches.errors.error_by_all").expect("the golden")
    );

    let clean = reads(0, 105, 40, 60);
    assert_eq!(
        run(
            &clean,
            Calculator::Error,
            Stratifier::All,
            &Options::default()
        ),
        field(&text, "metrics", "no-mismatches.errors.error_by_all").expect("the golden")
    );
    let another_prior = Options {
        prior_q: 10,
        ..Options::default()
    };
    assert_eq!(
        run(&clean, Calculator::Error, Stratifier::All, &another_prior),
        field(
            &text,
            "metrics",
            "no-mismatches-with-another-prior.errors.error_by_all"
        )
        .expect("the golden")
    );

    // The same bases, three priors, three qualities: none of them infinite.
    assert_eq!(q_score(0, 160, prior_error(30)), 52);
    assert_eq!(q_score(0, 160, prior_error(10)), 32);
    assert_eq!(q_score(2, 160, prior_error(30)), 19);
}

/// A stratifier splits the rows and nothing else.
#[test]
fn a_stratifier_splits_the_rows() {
    let text = corpus();
    let plain = reads(2, 105, 40, 60);
    let options = Options::default();

    for (stratifier, case) in [
        (
            Stratifier::BaseQuality,
            "stratified-by-base-quality.errors.error_by_base_quality",
        ),
        (
            Stratifier::Cycle,
            "stratified-by-cycle.errors.error_by_cycle",
        ),
        (
            Stratifier::GcContent,
            "two-metrics-two-files.errors.error_by_gc",
        ),
    ] {
        assert_eq!(
            run(&plain, Calculator::Error, stratifier, &options),
            field(&text, "metrics", case).expect("the golden"),
            "{case}"
        );
    }

    // The totals survive the split: twenty cycles of eight bases is the whole hundred and sixty.
    let Table::Base(rows) = collect(
        &plain,
        &reference(),
        1,
        &processed_loci(pileup(&plain, &options), &options),
        Calculator::Error,
        Stratifier::Cycle,
        &options,
    ) else {
        panic!("a simple error table")
    };
    assert_eq!(rows.len(), 20);
    assert_eq!(rows.iter().map(|row| row.total_bases).sum::<u64>(), 160);
    assert_eq!(rows.iter().map(|row| row.error_bases).sum::<u64>(), 2);

    // The edited reads carry a different GC, which is what splits that table in two.
    assert_eq!(gc_content(&window(101, 20)), 0.5);
    let mut edited = window(101, 20);
    edited[4] = b'C';
    assert_eq!(gc_content(&edited), 0.55);
}

/// The VCF takes bases out of the denominator rather than out of the rows.
#[test]
fn a_known_site_is_not_an_error() {
    let text = corpus();
    let plain = reads(2, 105, 40, 60);
    let known = Options {
        known_sites: vec![105],
        ..Options::default()
    };
    assert_eq!(
        run(&plain, Calculator::Error, Stratifier::All, &known),
        field(&text, "metrics", "a-known-site.errors.error_by_all").expect("the golden")
    );
}

/// One threshold per base and one per read, and a run that keeps no read at all.
#[test]
fn the_thresholds_drop_observations() {
    let text = corpus();

    let low_quality = reads(2, 105, 2, 60);
    assert_eq!(
        run(
            &low_quality,
            Calculator::Error,
            Stratifier::All,
            &Options::default()
        ),
        field(
            &text,
            "metrics",
            "a-low-quality-mismatch.errors.error_by_all"
        )
        .expect("the golden")
    );
    let lower = Options {
        min_base_q: 2,
        ..Options::default()
    };
    assert_eq!(
        run(&low_quality, Calculator::Error, Stratifier::All, &lower),
        field(
            &text,
            "metrics",
            "a-low-quality-mismatch-with-a-lower-threshold.errors.error_by_all"
        )
        .expect("the golden")
    );

    let poorly_mapped = reads(2, 105, 40, 5);
    // Every read is dropped, so the table has no rows, and a table with no rows is a file with no
    // header rather than a header over nothing.
    assert_eq!(
        run(
            &poorly_mapped,
            Calculator::Error,
            Stratifier::All,
            &Options::default()
        ),
        field(&text, "metrics", "poorly-mapped-reads.errors.error_by_all").expect("the golden")
    );
    let lower = Options {
        min_mapping_q: 1,
        ..Options::default()
    };
    assert_eq!(
        run(&poorly_mapped, Calculator::Error, Stratifier::All, &lower),
        field(
            &text,
            "metrics",
            "poorly-mapped-reads-with-a-lower-threshold.errors.error_by_all"
        )
        .expect("the golden")
    );
}

/// A pair that overlaps itself, counted in the metric that exists for it.
#[test]
fn an_overlapping_pair_has_its_own_metric() {
    let text = corpus();
    let mut first = window(101, 20);
    first[4] = if first[4] == b'A' { b'C' } else { b'A' };
    let pair = vec![
        Read {
            name: "p1".to_string(),
            start: 101,
            bases: first,
            qualities: vec![40; 20],
            flags: 0x1 | 0x2 | 0x40 | 0x20,
            mate_start: 106,
            cigar: vec![(20, 'M')],
            mapping_quality: 60,
        },
        Read {
            name: "p1".to_string(),
            start: 106,
            bases: window(106, 20),
            qualities: vec![40; 20],
            flags: 0x1 | 0x2 | 0x80 | 0x10,
            mate_start: 101,
            cigar: vec![(20, 'M')],
            mapping_quality: 60,
        },
    ];
    let options = Options::default();
    assert_eq!(
        run(&pair, Calculator::Error, Stratifier::All, &options),
        field(&text, "metrics", "an-overlapping-pair.errors.error_by_all").expect("the golden")
    );
    assert_eq!(
        run(
            &pair,
            Calculator::OverlappingError,
            Stratifier::All,
            &options
        ),
        field(
            &text,
            "metrics",
            "an-overlapping-pair.errors.overlapping_error_by_all"
        )
        .expect("the golden")
    );

    // The overlap is fifteen positions read twice, and the mismatch is outside it, which is why
    // the simple table counts an error the overlapping one does not.
    let Table::Overlapping(rows) = collect(
        &pair,
        &reference(),
        1,
        &processed_loci(pileup(&pair, &options), &options),
        Calculator::OverlappingError,
        Stratifier::All,
        &options,
    ) else {
        panic!("an overlapping table")
    };
    assert_eq!(rows[0].bases_with_overlapping_reads, 30);
    assert_eq!(rows[0].disagrees_with_reference_only, 0);
}

/// An indel has a metric, a denominator, and a column it never fills in.
#[test]
fn a_deletion_is_counted_once_for_the_read() {
    let text = corpus();
    let indels: Vec<Read> = (0..4)
        .map(|index| {
            let mut bases = window(101, 10);
            bases.extend(window(113, 8));
            Read {
                name: format!("d{index}"),
                start: 101,
                bases,
                qualities: vec![40; 18],
                flags: 0,
                mate_start: 0,
                cigar: vec![(10, 'M'), (2, 'D'), (8, 'M')],
                mapping_quality: 60,
            }
        })
        .collect();
    let options = Options::default();
    assert_eq!(
        run(&indels, Calculator::IndelError, Stratifier::All, &options),
        field(&text, "metrics", "a-deletion.errors.indel_error_by_all").expect("the golden")
    );

    let Table::Indel(rows) = collect(
        &indels,
        &reference(),
        1,
        &processed_loci(pileup(&indels, &options), &options),
        Calculator::IndelError,
        Stratifier::All,
        &options,
    ) else {
        panic!("an indel table")
    };
    // Four reads, one deletion of two bases each: the deletion spans two loci and is still one
    // deletion.
    assert_eq!(rows[0].deletions, 4);
    assert_eq!(rows[0].deleted_bases, 8);
    // The deleted bases are the errors, and the deleted loci are not part of the denominator.
    assert_eq!(rows[0].error_bases, 8);
    assert_eq!(rows[0].total_bases, 72);
    // And the inherited column stays at zero: this metric derives its own fields and does not
    // derive that one.
    assert_eq!(rows[0].q_score, 0);
}

/// The cap counts what is left after the known sites have been taken out.
#[test]
fn the_cap_stops_early() {
    let text = corpus();
    let plain = reads(2, 105, 40, 60);
    let capped = Options {
        max_loci: 5,
        ..Options::default()
    };
    assert_eq!(
        run(&plain, Calculator::Error, Stratifier::All, &capped),
        field(&text, "metrics", "a-cap-on-the-loci.errors.error_by_all").expect("the golden")
    );
}

/// The list of metrics is appended to, so naming one it already carries is a refusal.
#[test]
fn the_default_list_is_appended_to() {
    let text = corpus();
    let mut appended: Vec<String> = DEFAULT_ERROR_METRICS
        .iter()
        .map(|s| s.to_string())
        .collect();
    appended.push("ERROR".to_string());
    let refusal = suffixes(&appended).expect_err("a duplicated suffix");
    assert_eq!(
        refusal,
        Refusal::DuplicatedSuffix {
            suffix: "error_by_all".to_string(),
            class: "class picard.sam.SamErrorMetric.BaseErrorAggregation".to_string(),
        }
    );
    let recorded = field(&text, "error", "the-default-list-is-appended-to").expect("the golden");
    assert_eq!(
        recorded,
        format!("java.lang.IllegalArgumentException:{}", refusal.message())
    );

    // The default list on its own writes twenty-seven files, all differently named.
    let default: Vec<String> = DEFAULT_ERROR_METRICS
        .iter()
        .map(|s| s.to_string())
        .collect();
    let written = suffixes(&default).expect("no duplicate");
    assert_eq!(written.len(), 27);
    assert_eq!(written[0], "error_by_all");
    assert_eq!(written[26], "indel_error_by_all");
    // Several stratifiers fold into one suffix from the left.
    assert_eq!(
        aggregation_suffix("ERROR:READ_ORDINALITY:CYCLE").as_deref(),
        Some("error_by_read_ordinality_and_cycle")
    );
    // And the suffix is not the name lower-cased.
    assert_eq!(
        aggregation_suffix("ERROR:GC_CONTENT").as_deref(),
        Some("error_by_gc")
    );
}

/// Two metrics in one run write two files, and the file is named for the directive.
#[test]
fn one_file_per_metric() {
    let text = corpus();
    let recorded = field(&text, "files", "two-metrics-two-files").expect("the golden");
    let written =
        suffixes(&["ERROR".to_string(), "ERROR:GC_CONTENT".to_string()]).expect("no duplicate");
    assert_eq!(
        recorded,
        written
            .iter()
            .map(|suffix| format!("errors.{suffix}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
}
