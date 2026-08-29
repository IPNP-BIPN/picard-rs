//! Conformance for which alignment `MergeBamAlignment` calls primary, against Picard 3.4.0.
//!
//! Golden from `tools/mergebamalignment-conformance/`: ten runs over three alignments of one read
//! and two pairings of one pair.
//!
//! # What this suite is for
//!
//!  * **three strategies taking the mapping quality and one taking the earliest read base**;
//!  * **the aligner's own primary standing where it named exactly one**;
//!  * **the losers being written as secondary rather than dropped**;
//!  * **a pairing being chosen as a pairing**, by quality or by distance;
//!  * **and `EarliestFragment` refusing a paired read outright.**

use std::io::Read;

use picard_analysis::merge_bam_alignment_hits::{
    aligners_choice_stands, combine_mapqs, earliest_fragment_refusal, index_of_first_aligned_base,
    pair_distance, primary_candidates, tally_primary_alignments, written, Alignment, Hit,
    NumPrimaryAlignmentState, Strategy,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/merge_bam_alignment_hits.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn records(text: &str, case: &str) -> Vec<Vec<String>> {
    let prefix = format!("record\t{case}\t");
    text.lines()
        .filter(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .split('\t')
                .map(str::to_string)
                .collect()
        })
        .collect()
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].replace("\\t", "\t"))
}

/// Which of a case's records is the primary one, by its flag.
fn primary_row(text: &str, case: &str) -> usize {
    records(text, case)
        .iter()
        .position(|record| {
            let flags: u16 = record[1].parse().expect("a flag");
            flags & 0x100 == 0
        })
        .expect("a primary")
}

fn alignment(
    start: i32,
    mapping_quality: i32,
    cigar: &[(usize, char)],
    secondary: bool,
) -> Alignment {
    Alignment {
        reference: "chr1".to_string(),
        start,
        mapping_quality,
        cigar: cigar.to_vec(),
        negative_strand: false,
        secondary,
    }
}

/// The three alignments the fixture gives one read: the best quality is NOT on the alignment that
/// maps the read's first base.
fn three_hits(named_primary: Option<usize>) -> Vec<Hit> {
    [
        (121, 60, vec![(5, 'S'), (15, 'M')]),
        (61, 30, vec![(20, 'M')]),
        (41, 10, vec![(10, 'S'), (10, 'M')]),
    ]
    .iter()
    .enumerate()
    .map(|(index, (start, quality, cigar))| Hit {
        first: Some(alignment(
            *start,
            *quality,
            cigar,
            named_primary != Some(index),
        )),
        second: None,
    })
    .collect()
}

/// Three strategies take the quality and one takes the earliest base.
#[test]
fn the_strategies_disagree() {
    let text = corpus();
    let hits = three_hits(None);

    // Nobody named a primary, so every strategy chooses.
    assert!(!aligners_choice_stands(&hits));

    for (strategy, case) in [
        (Strategy::BestMapq, "three-hits-best-mapq"),
        (Strategy::BestEndMapq, "three-hits-best-end-mapq"),
        (Strategy::MostDistant, "three-hits-most-distant"),
    ] {
        let candidates = primary_candidates(&hits, strategy);
        assert_eq!(candidates, vec![0], "{case}");
        assert_eq!(primary_row(&text, case), 0, "{case}");
    }

    // The earliest-base strategy takes the second alignment, whose quality is half the first's.
    assert_eq!(
        primary_candidates(&hits, Strategy::EarliestFragment),
        vec![1]
    );
    assert_eq!(primary_row(&text, "three-hits-earliest-fragment"), 1);
    // Which base each alignment starts at is what decides it.
    assert_eq!(
        index_of_first_aligned_base(&alignment(121, 60, &[(5, 'S'), (15, 'M')], true)),
        5
    );
    assert_eq!(
        index_of_first_aligned_base(&alignment(61, 30, &[(20, 'M')], true)),
        0
    );
}

/// The aligner's own primary stands where it named exactly one.
#[test]
fn the_aligners_choice_stands() {
    let text = corpus();
    let hits = three_hits(Some(1));
    assert_eq!(
        tally_primary_alignments(&hits, true),
        NumPrimaryAlignmentState::One
    );
    assert!(aligners_choice_stands(&hits));
    assert_eq!(primary_candidates(&hits, Strategy::BestMapq), vec![1]);
    assert_eq!(primary_row(&text, "the-aligners-own-primary"), 1);

    // Two named primaries are as good as none: the strategy chooses again.
    let two = vec![
        Hit {
            first: Some(alignment(121, 60, &[(20, 'M')], false)),
            second: None,
        },
        Hit {
            first: Some(alignment(61, 30, &[(20, 'M')], false)),
            second: None,
        },
    ];
    assert_eq!(
        tally_primary_alignments(&two, true),
        NumPrimaryAlignmentState::MoreThanOne
    );
    assert!(!aligners_choice_stands(&two));
}

/// The losers are written as secondary rather than dropped.
#[test]
fn the_losers_are_written_as_secondary() {
    let text = corpus();
    let hits = three_hits(None);
    assert_eq!(written(&hits, 0, true), vec![0, 1, 2]);
    assert_eq!(records(&text, "three-hits-best-mapq").len(), 3);
    let flags: Vec<u16> = records(&text, "three-hits-best-mapq")
        .iter()
        .map(|record| record[1].parse().expect("a flag"))
        .collect();
    assert_eq!(flags, vec![0, 256, 256]);

    assert_eq!(written(&hits, 0, false), vec![0]);
    assert_eq!(
        records(&text, "three-hits-without-the-secondaries").len(),
        1
    );
}

/// A pairing is chosen as a pairing, by quality or by distance.
#[test]
fn a_pairing_is_chosen_as_a_pairing() {
    let text = corpus();
    let pairing = |first_start: i32, second_start: i32, quality: i32| Hit {
        first: Some(alignment(first_start, quality, &[(20, 'M')], true)),
        second: Some(alignment(second_start, quality, &[(20, 'M')], true)),
    };
    // Two pairings: the near one is better mapped, the far one reaches further.
    let hits = vec![pairing(41, 81, 60), pairing(41, 161, 30)];

    assert_eq!(primary_candidates(&hits, Strategy::BestMapq), vec![0]);
    assert_eq!(primary_row(&text, "a-pair-with-two-hits"), 0);

    assert_eq!(primary_candidates(&hits, Strategy::MostDistant), vec![1]);
    assert_eq!(primary_row(&text, "a-pair-with-two-hits-most-distant"), 0);
    // The chosen pairing is written first under that strategy, which is why its row is the first
    // one: the reference moves it to the head of the list rather than marking it in place.
    let chosen = &records(&text, "a-pair-with-two-hits-most-distant")[0];
    assert_eq!(chosen[7], "161");

    // The distance is the whole span of the pairing, and the quality of a pair combines both ends
    // a hundred to one.
    assert_eq!(
        pair_distance(
            &alignment(41, 60, &[(20, 'M')], true),
            &alignment(161, 30, &[(20, 'M')], true)
        ),
        140
    );
    assert_eq!(combine_mapqs(60, 30), 9000);
    // An unknown quality of 255 counts as one rather than as the largest.
    assert_eq!(combine_mapqs(255, 0), 1);
}

/// The fragment strategy refuses a paired read outright.
#[test]
fn earliest_fragment_refuses_a_pair() {
    let text = corpus();
    let recorded = field(&text, "error", "a-pair-under-earliest-fragment").expect("the golden");
    assert_eq!(
        recorded,
        format!(
            "java.lang.UnsupportedOperationException:{}",
            earliest_fragment_refusal("p", "1/2 20b aligned to chr1:41-60.")
        )
    );
}
