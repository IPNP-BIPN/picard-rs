//! Conformance for what `MergeBamAlignment` does to a pair, against Picard 3.4.0.
//!
//! Golden from `tools/mergebamalignment-conformance/`: eleven runs whose aligned records carry
//! deliberately wrong mate fields, so that what the merger decides is visible rather than
//! inherited.
//!
//! # What this suite is for
//!
//!  * **the mate fields being written from the other end**, whatever the aligner put there;
//!  * **the insert size being computed and signed by which end starts first**;
//!  * **the proper-pair flag being the merger's**, unless the aligner's are kept;
//!  * **`MC` and `MQ` being added, and dropped where asked**;
//!  * **an overlapping pair being clipped, which MOVES the end that reaches back**;
//!  * **hard clipping keeping what it removed, in `XB` and `XQ`**;
//!  * **and an end with no alignment being placed at its mate's coordinates.**

use std::io::Read;

use picard_analysis::merge_bam_alignment_pair::{
    clip_for_overlapping_reads, fix_up_pair, insert_size, is_proper_pair, pair_orientation,
    set_mate_info, End, Options, PairOrientation,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/merge_bam_alignment_pair.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

/// The record lines of one case, in the order the file wrote them.
fn records(text: &str, case: &str) -> Vec<String> {
    let prefix = format!("record\t{case}\t");
    text.lines()
        .filter(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].replace("\\t", "\t"))
        .collect()
}

/// One aligned end, as the fixture hands it to the merger: the mate fields are the aligner's, and
/// they are wrong.
fn end(first_of_pair: bool, start: i32, cigar: &[(usize, char)], negative: bool) -> End {
    End {
        reference: Some("chr1".to_string()),
        start,
        mapping_quality: 60,
        negative_strand: negative,
        unmapped: false,
        cigar: cigar.to_vec(),
        bases: b"ACGTACGT".to_vec(),
        qualities: vec![40; 8],
        proper_pair: false,
        mate_reference: Some("chr2".to_string()),
        mate_start: 39,
        mate_negative_strand: false,
        mate_unmapped: false,
        insert_size: if first_of_pair { 999 } else { -999 },
        mate_mapping_quality: None,
        mate_cigar: None,
        clipped_bases: None,
        clipped_qualities: None,
    }
}

/// The record as the file writes it, for the fields this suite is about.
fn line(name: &str, flags: u16, record: &End) -> String {
    let mut tags = Vec::new();
    if let Some(bases) = &record.clipped_bases {
        tags.push(format!("XB:Z:{bases}"));
    }
    if let Some(cigar) = &record.mate_cigar {
        tags.push(format!("MC:Z:{cigar}"));
    }
    tags.push("PG:Z:bwa".to_string());
    tags.push("RG:Z:rg1".to_string());
    if let Some(quality) = record.mate_mapping_quality {
        tags.push(format!("MQ:i:{quality}"));
    }
    if let Some(qualities) = &record.clipped_qualities {
        tags.push(format!("XQ:Z:{qualities}"));
    }
    let bases = String::from_utf8(record.bases.clone()).expect("bases are ASCII");
    let qualities: String = record
        .qualities
        .iter()
        .map(|quality| (quality + 33) as char)
        .collect();
    format!(
        "{name}\t{flags}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{bases}\t{qualities}\t{}",
        record.reference.clone().unwrap_or_else(|| "*".to_string()),
        record.start,
        record.mapping_quality,
        record.cigar_string(),
        match (&record.mate_reference, &record.reference) {
            (Some(mate), Some(own)) if mate == own => "=".to_string(),
            (Some(mate), _) => mate.clone(),
            (None, _) => "*".to_string(),
        },
        record.mate_start,
        record.insert_size,
        tags.join("\t")
    )
}

/// The mate fields the aligner wrote are replaced, and the insert size is computed.
#[test]
fn the_mate_fields_come_from_the_other_end() {
    let text = corpus();
    let mut first = end(true, 1, &[(8, 'M')], false);
    let mut second = end(false, 25, &[(8, 'M')], false);
    fix_up_pair(&mut first, &mut second, &Options::default());

    assert_eq!(
        vec![line("p", 65, &first), line("p", 129, &second)],
        records(&text, "mate-fields-from-the-other-end")
    );
    // The mate contig the aligner named is gone, the position is the other end's, and the mate's
    // mapping quality and cigar are on the record.
    assert_eq!(first.mate_reference.as_deref(), Some("chr1"));
    assert_eq!(first.mate_start, 25);
    assert_eq!(first.mate_mapping_quality, Some(60));
    assert_eq!(first.mate_cigar.as_deref(), Some("8M"));
    assert_eq!((first.insert_size, second.insert_size), (25, -25));

    // The same pair the other way round signs it the other way.
    let mut first = end(true, 25, &[(8, 'M')], false);
    let mut second = end(false, 1, &[(8, 'M')], false);
    fix_up_pair(&mut first, &mut second, &Options::default());
    assert_eq!(
        vec![line("p", 65, &first), line("p", 129, &second)],
        records(&text, "the-second-end-first")
    );
    assert_eq!((first.insert_size, second.insert_size), (-25, 25));
}

/// The proper-pair flag is the merger's own, unless the aligner's are kept.
#[test]
fn the_proper_pair_flag_is_the_mergers() {
    let text = corpus();
    let across = |aligner_flags: bool| {
        let mut first = End {
            reference: Some("chr1".to_string()),
            proper_pair: true,
            ..end(true, 1, &[(8, 'M')], false)
        };
        let mut second = End {
            reference: Some("chr2".to_string()),
            proper_pair: true,
            ..end(false, 1, &[(8, 'M')], false)
        };
        let options = Options {
            aligner_proper_pair_flags: aligner_flags,
            ..Options::default()
        };
        fix_up_pair(&mut first, &mut second, &options);
        (first, second)
    };

    let (first, second) = across(false);
    assert!(!first.proper_pair);
    assert_eq!(
        vec![line("p", 65, &first), line("p", 129, &second)],
        records(&text, "across-two-contigs")
    );
    // Two contigs is no insert size either.
    assert_eq!(first.insert_size, 0);

    let (first, second) = across(true);
    assert!(first.proper_pair);
    assert_eq!(
        vec![line("p", 67, &first), line("p", 131, &second)],
        records(&text, "across-two-contigs-with-the-aligners-flags")
    );

    // A pair on one contig pointing at each other is proper; the same pair pointing away is not.
    let mut first = end(true, 1, &[(8, 'M')], false);
    let mut second = end(false, 25, &[(8, 'M')], true);
    set_mate_info(&mut first, &mut second, true);
    assert_eq!(pair_orientation(&first), PairOrientation::FR);
    assert!(is_proper_pair(&first, &second, &[PairOrientation::FR]));
    assert!(!is_proper_pair(&first, &second, &[PairOrientation::RF]));
}

/// The mate's cigar is added unless the run says not to.
#[test]
fn the_mate_cigar_is_added_unless_it_is_not() {
    let text = corpus();
    let mut first = end(true, 1, &[(8, 'M')], false);
    let mut second = end(false, 25, &[(8, 'M')], false);
    let options = Options {
        add_mate_cigar: false,
        ..Options::default()
    };
    fix_up_pair(&mut first, &mut second, &options);
    assert_eq!(first.mate_cigar, None);
    assert_eq!(
        vec![line("p", 65, &first), line("p", 129, &second)],
        records(&text, "without-the-mate-cigar")
    );
    // The mate's mapping quality is not part of that bargain: MQ stays.
    assert_eq!(first.mate_mapping_quality, Some(60));
}

/// Clipping an overlapping pair moves the end that reaches back past its mate.
#[test]
fn an_overlapping_pair_is_clipped() {
    let text = corpus();
    let overlapping = || {
        (
            end(true, 5, &[(8, 'M')], false),
            End {
                negative_strand: true,
                ..end(false, 3, &[(8, 'M')], true)
            },
        )
    };

    let (mut first, mut second) = overlapping();
    fix_up_pair(&mut first, &mut second, &Options::default());
    assert_eq!(
        vec![line("p", 99, &first), line("p", 147, &second)],
        records(&text, "overlapping-ends")
    );
    // The forward end keeps its start and loses two bases from its far end; the reverse end's
    // START MOVES from 3 to 5, the clip being at its own five-prime end.
    assert_eq!((first.start, first.cigar_string()), (5, "6M2S".to_string()));
    assert_eq!(
        (second.start, second.cigar_string()),
        (5, "2S6M".to_string())
    );
    assert_eq!(first.insert_size, 6);

    // Turned off, the ends stay where the aligner put them.
    let (mut first, mut second) = overlapping();
    let unclipped = Options {
        clip_overlapping_reads: false,
        ..Options::default()
    };
    fix_up_pair(&mut first, &mut second, &unclipped);
    assert_eq!(
        vec![line("p", 99, &first), line("p", 147, &second)],
        records(&text, "overlapping-ends-unclipped")
    );
    assert_eq!((second.start, second.cigar_string()), (3, "8M".to_string()));

    // A pair that does not overlap is not touched.
    let mut apart = end(true, 1, &[(8, 'M')], false);
    let mut away = End {
        negative_strand: true,
        ..end(false, 25, &[(8, 'M')], true)
    };
    clip_for_overlapping_reads(&mut apart, &mut away, false);
    assert_eq!(apart.cigar_string(), "8M");
    assert_eq!(away.cigar_string(), "8M");
}

/// Hard clipping keeps what it removed, so the bases are recoverable from the record.
#[test]
fn hard_clipping_keeps_what_it_removed() {
    let text = corpus();
    let mut first = end(true, 5, &[(8, 'M')], false);
    let mut second = End {
        negative_strand: true,
        ..end(false, 3, &[(8, 'M')], true)
    };
    let options = Options {
        hard_clip_overlapping_reads: true,
        ..Options::default()
    };
    fix_up_pair(&mut first, &mut second, &options);
    assert_eq!(
        vec![line("p", 99, &first), line("p", 147, &second)],
        records(&text, "overlapping-ends-hard-clipped")
    );
    // The record holds six bases now, and the two it lost are in XB with their qualities in XQ,
    // in the order the sequencer produced them.
    assert_eq!(first.bases.len(), 6);
    assert_eq!(first.clipped_bases.as_deref(), Some("GT"));
    assert_eq!(first.clipped_qualities.as_deref(), Some("II"));
    assert_eq!(first.cigar_string(), "6M2H");
    assert_eq!(second.cigar_string(), "2H6M");
}

/// An end with no alignment is placed at its mate's coordinates.
#[test]
fn an_unaligned_end_is_placed_at_its_mates() {
    let text = corpus();
    let mut mapped = end(true, 1, &[(8, 'M')], false);
    let mut unaligned = End {
        reference: None,
        start: 0,
        mapping_quality: 0,
        unmapped: true,
        cigar: Vec::new(),
        mate_reference: None,
        mate_start: 0,
        insert_size: 0,
        ..end(false, 0, &[], false)
    };
    fix_up_pair(&mut mapped, &mut unaligned, &Options::default());
    assert_eq!(
        vec![line("p", 73, &mapped), line("p", 133, &unaligned)],
        records(&text, "one-end-unaligned")
    );
    // The unmapped end sits where its mate does, and it is the one that carries MC and MQ: the
    // mapped end has no mate cigar, its mate having no cigar to give.
    assert_eq!(unaligned.reference.as_deref(), Some("chr1"));
    assert_eq!(unaligned.start, 1);
    assert_eq!(unaligned.mate_cigar.as_deref(), Some("8M"));
    assert_eq!(mapped.mate_cigar, None);
    assert_eq!(mapped.mate_mapping_quality, None);
    assert_eq!(insert_size(&mapped, &unaligned), 0);
}

/// A soft-clipped end is measured from where it aligns, not from where it would.
#[test]
fn the_insert_size_is_measured_from_the_alignment() {
    let text = corpus();
    let mut first = end(true, 1, &[(2, 'S'), (6, 'M')], false);
    let mut second = end(false, 25, &[(8, 'M')], false);
    fix_up_pair(&mut first, &mut second, &Options::default());
    assert_eq!(
        vec![line("p", 65, &first), line("p", 129, &second)],
        records(&text, "a-soft-clipped-end")
    );
    // The clip does not move the start, so the distance is the same twenty-five as an unclipped
    // pair's, and the mate cigar carries the clip.
    assert_eq!(first.insert_size, 25);
    assert_eq!(second.mate_cigar.as_deref(), Some("2S6M"));
}

/// The deprecated argument changes nothing.
#[test]
fn the_deprecated_argument_is_inert() {
    let text = corpus();
    assert_eq!(
        records(&text, "the-deprecated-paired-run"),
        records(&text, "mate-fields-from-the-other-end")
    );
}
