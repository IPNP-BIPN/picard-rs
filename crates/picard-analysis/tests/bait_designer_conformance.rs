//! Conformance for `BaitDesigner` against Picard 3.4.0.
//!
//! Golden from `tools/baitdesigner-conformance/`: eleven designs over a contig of repeating
//! `ACGT`, so a bait's sequence is known by arithmetic.
//!
//! # What this suite is for
//!
//!  * **the baits being laid at a fixed offset**, overlapping where the target is long;
//!  * **a target shorter than a bait still getting the minimum number of them**, hanging off both
//!    ends;
//!  * **the three strategies laying the same short target differently**;
//!  * **the primers being on the ordered sequence and on neither interval list**;
//!  * **two targets a bait could cover together being merged into one design**;
//!  * **and the pool being the design repeated until the plate is full**, every second copy
//!    reverse complemented and the numbering restarting each time.

use std::io::Read;

use picard_analysis::bait_designer::{
    copies, design, estimate_baits, make_bait_name, pool_rows, prepare_targets, render_bait,
    render_pool, reverse_complement, written_files, Options, Strategy, Target,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bait_designer.txt.gz");
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

const LENGTH: usize = 2000;

/// The contig: `ACGT` repeating.
fn reference() -> Vec<u8> {
    (0..LENGTH).map(|index| b"ACGT"[index % 4]).collect()
}

fn target(start: i32, end: i32, index: usize) -> Target {
    Target {
        contig: "chr1".to_string(),
        start,
        end,
        negative_strand: false,
        name: format!("target{index}"),
    }
}

/// One case: the bait interval list it writes.
fn baits_of(targets: &[Target], options: &Options) -> String {
    design(targets, &reference(), options)
        .iter()
        .map(render_bait)
        .collect::<Vec<_>>()
        .join("\n")
}

/// The offset is fixed, so a target longer than one bait is covered by several that overlap.
#[test]
fn the_baits_are_laid_at_a_fixed_offset() {
    let text = corpus();
    let long = [target(201, 500, 1)];
    assert_eq!(
        baits_of(&long, &Options::default()),
        field(&text, "baits", "a-target-longer-than-a-bait").expect("the golden")
    );

    // Every bait is the requested length, and each begins one offset after the last.
    let baits = design(&long, &reference(), &Options::default());
    assert!(baits.iter().all(|bait| bait.length() == 120));
    let starts: Vec<i32> = baits.iter().map(|bait| bait.start).collect();
    assert_eq!(starts, vec![171, 251, 331, 411]);

    // The bait's shape is what moves them.
    let smaller = Options {
        bait_size: 60,
        ..Options::default()
    };
    assert_eq!(
        baits_of(&long, &smaller),
        field(&text, "baits", "a-smaller-bait").expect("the golden")
    );
    let wider = Options {
        bait_offset: 120,
        ..Options::default()
    };
    assert_eq!(
        baits_of(&long, &wider),
        field(&text, "baits", "a-wider-offset").expect("the golden")
    );
    let padded = Options {
        padding: 50,
        ..Options::default()
    };
    assert_eq!(
        baits_of(&long, &padded),
        field(&text, "baits", "with-padding").expect("the golden")
    );
}

/// A target shorter than a bait still gets the minimum number of them.
#[test]
fn a_short_target_gets_the_minimum() {
    let text = corpus();
    let short = [target(201, 260, 1)];
    assert_eq!(
        baits_of(&short, &Options::default()),
        field(&text, "baits", "a-target-shorter-than-a-bait").expect("the golden")
    );
    // Two baits over sixty bases: the target is widened to what two baits at the offset would
    // tile before anything is laid, so the first starts well before it and the last ends well
    // after it.
    let baits = design(&short, &reference(), &Options::default());
    assert_eq!(baits.len(), 2);
    assert!(baits[0].start < 201);
    assert!(baits[1].end > 260);

    let one = Options {
        minimum_baits_per_target: 1,
        ..Options::default()
    };
    assert_eq!(
        baits_of(&short, &one),
        field(&text, "baits", "a-short-target-with-one-bait").expect("the golden")
    );
    assert_eq!(design(&short, &reference(), &one).len(), 1);
}

/// The three strategies lay the same short target differently.
#[test]
fn the_strategies_differ() {
    let text = corpus();
    let short = [target(201, 260, 1)];
    for (strategy, case) in [
        (
            Strategy::CenteredConstrained,
            "strategy-centeredconstrained",
        ),
        (Strategy::FixedOffset, "strategy-fixedoffset"),
        (Strategy::Simple, "strategy-simple"),
    ] {
        let options = Options {
            strategy,
            ..Options::default()
        };
        assert_eq!(
            baits_of(&short, &options),
            field(&text, "baits", case).expect("the golden"),
            "{case}"
        );
    }

    // The constrained strategy centres one bait on a target it cannot tile; the simple one starts
    // at the target's own start and stops as soon as a bait would begin past its end.
    let centred = design(
        &short,
        &reference(),
        &Options {
            strategy: Strategy::CenteredConstrained,
            ..Options::default()
        },
    );
    assert_eq!(centred.len(), 1);
    assert_eq!(centred[0].start, 171);
    let simple = design(
        &short,
        &reference(),
        &Options {
            strategy: Strategy::Simple,
            ..Options::default()
        },
    );
    assert_eq!(simple.len(), 1);
    assert_eq!(simple[0].start, 201);
}

/// Two targets a bait could cover together become one design, unless the merging is turned off.
#[test]
fn nearby_targets_are_merged() {
    let text = corpus();
    let two = [target(201, 260, 1), target(301, 360, 2)];
    assert_eq!(
        baits_of(&two, &Options::default()),
        field(&text, "baits", "two-nearby-targets").expect("the golden")
    );
    let unmerged = Options {
        merge_nearby_targets: false,
        ..Options::default()
    };
    assert_eq!(
        baits_of(&two, &unmerged),
        field(&text, "baits", "two-nearby-targets-unmerged").expect("the golden")
    );

    // Merged, the two intervals are one target, and both baits are named for the first of them.
    let merged = prepare_targets(&two, LENGTH as i32, &Options::default());
    assert_eq!(merged.len(), 1);
    assert_eq!((merged[0].start, merged[0].end), (201, 360));
    assert!(baits_of(&two, &Options::default()).contains("target1_bait#2"));
    // Unmerged, each target is designed on its own and the second one says so.
    assert!(baits_of(&two, &unmerged).contains("target2_bait#"));

    // The estimate that decides it reads like a ceiling and is not one: the rounding is applied
    // before the division and the division itself is truncated.
    assert_eq!(estimate_baits(201, 260, &Options::default()), 2);
    assert_eq!(estimate_baits(201, 400, &Options::default()), 2);
    // Two hundred and twenty bases is one and a quarter offsets past the first bait, and the
    // estimate is two: a ceiling would have said three.
    assert_eq!(estimate_baits(201, 420, &Options::default()), 2);
    assert_eq!(estimate_baits(201, 481, &Options::default()), 3);
}

/// The pool is the design repeated until the plate is full.
#[test]
fn the_pool_is_the_design_repeated() {
    let text = corpus();
    let options = Options::default();
    let baits = design(&[target(201, 500, 1)], &reference(), &options);
    let rows = pool_rows(&baits, &options);
    assert_eq!(
        render_pool(&rows),
        field(&text, "pool", "a-target-longer-than-a-bait").expect("the golden")
    );

    // Four baits fill the plate thirteen thousand seven hundred and fifty times over.
    assert_eq!(copies(baits.len(), &options), 13_750);
    assert_eq!(rows.len(), 55_000);
    // The numbering restarts with each copy rather than running on.
    assert_eq!(rows[0].0, "design_000001");
    assert_eq!(rows[4].0, "design_000001");
    // And every second copy is reverse complemented, which is the same sequence read the other
    // way.
    assert_eq!(
        rows[4].1,
        String::from_utf8(reverse_complement(rows[0].1.as_bytes())).expect("bases")
    );

    // A design that does not divide the plate evenly stops short of it.
    let wider = Options {
        bait_offset: 120,
        ..Options::default()
    };
    let three = design(&[target(201, 500, 1)], &reference(), &wider);
    assert_eq!(three.len(), 3);
    assert_eq!(pool_rows(&three, &wider).len(), 54_999);
}

/// The primers are on the ordered sequence and on neither interval list.
#[test]
fn the_primers_are_on_the_order_only() {
    let text = corpus();
    let options = Options::default();
    let baits = design(&[target(201, 500, 1)], &reference(), &options);
    let rows = pool_rows(&baits, &options);

    assert!(rows[0].1.starts_with(&options.left_primer));
    assert!(rows[0].1.ends_with(&options.right_primer));
    assert_eq!(
        rows[0].1.len(),
        options.left_primer.len() + 120 + options.right_primer.len()
    );
    let listed = field(&text, "baits", "a-target-longer-than-a-bait").expect("the golden");
    assert!(!listed.contains(&options.left_primer));
}

/// A bait's name is its target's, with the index padded to the width of the count.
#[test]
fn the_names_are_padded_to_the_count() {
    assert_eq!(make_bait_name("target1", 1, 1), "target1_bait#1");
    assert_eq!(make_bait_name("target1", 2, 4), "target1_bait#2");
    assert_eq!(make_bait_name("target1", 9, 12), "target1_bait#09");
}

/// The design writes a directory of files named after it.
#[test]
fn the_output_is_a_directory() {
    let text = corpus();
    assert_eq!(
        written_files(&Options::default()).join(" "),
        field(&text, "files", "a-target-longer-than-a-bait").expect("the golden")
    );
}
