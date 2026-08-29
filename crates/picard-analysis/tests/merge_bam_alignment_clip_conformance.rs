//! Conformance for what `MergeBamAlignment` clips and unmaps, against Picard 3.4.0.
//!
//! Golden from `tools/mergebamalignment-conformance/`: eleven runs over a two-hundred-base contig.
//!
//! # What this suite is for
//!
//!  * **an alignment off the end of its contig being clipped to fit**, and not clipped twice;
//!  * **the adapter `XT` marked being clipped from that base on**, unless the run says otherwise;
//!  * **contamination taking two things**: too few aligned bases AND clips at both ends;
//!  * **and the four unmapping strategies keeping four different things.**

use std::io::Read;

use picard_analysis::merge_bam_alignment_clip::{
    cigar_string, clip_adapter, clip_off_the_end, encode_mapping_information, is_contaminant,
    unmap_contaminant, Record, UnmappingReadStrategy, CONTAMINATION_COMMENT,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/merge_bam_alignment_clip.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

/// The record lines of one case, split into fields.
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

/// One record's flag, position, mapping quality and cigar.
fn shape(record: &[String]) -> (u16, i32, i32, String) {
    (
        record[1].parse().expect("a flag"),
        record[3].parse().expect("a position"),
        record[4].parse().expect("a mapping quality"),
        record[5].clone(),
    )
}

/// An alignment that runs off the end of its contig is clipped to fit, once.
#[test]
fn an_alignment_off_the_end_is_clipped() {
    let text = corpus();
    // Forty bases at 171 of a two-hundred-base contig: ten of them have nowhere to go.
    let clipped = clip_off_the_end(171, &[(40, 'M')], 200).expect("a clip");
    assert_eq!(cigar_string(&clipped), "30M10S");
    assert_eq!(
        shape(&records(&text, "off-the-end-of-the-contig")[0]).3,
        cigar_string(&clipped)
    );

    // The same read already ending in a five-base soft clip is clipped to the same cigar, the
    // existing clip being taken off what is added rather than added to it.
    let already = clip_off_the_end(171, &[(35, 'M'), (5, 'S')], 200).expect("a clip");
    assert_eq!(cigar_string(&already), "30M10S");
    assert_eq!(
        shape(&records(&text, "off-the-end-with-a-soft-clip")[0]).3,
        cigar_string(&already)
    );

    // An alignment that fits is not touched at all.
    assert_eq!(clip_off_the_end(41, &[(40, 'M')], 200), None);
    assert_eq!(shape(&records(&text, "inside-the-contig")[0]).3, "40M");
}

/// The adapter is clipped from the base `XT` names to the three-prime end.
#[test]
fn the_adapter_is_clipped_from_where_it_starts() {
    let text = corpus();
    let clipped = clip_adapter(&[(40, 'M')], false, 31);
    assert_eq!(cigar_string(&clipped), "30M10S");
    assert_eq!(
        shape(&records(&text, "an-adapter-marked-in-the-unmapped-bam")[0]).3,
        cigar_string(&clipped)
    );
    // Turned off, the same record keeps its whole alignment, and the tag stays on it either way.
    assert_eq!(shape(&records(&text, "an-adapter-left-alone")[0]).3, "40M");
    assert!(records(&text, "an-adapter-left-alone")[0]
        .iter()
        .any(|field| field == "XT:i:31"));
}

/// Contamination takes two things: too few aligned bases, and clips at both ends.
#[test]
fn contamination_takes_two_things() {
    let text = corpus();
    // Ten aligned bases between two clips, under a threshold of thirty-two.
    let both_ends = [(15, 'S'), (10, 'M'), (15, 'S')];
    assert!(is_contaminant(&both_ends, 32));
    assert_eq!(shape(&records(&text, "a-short-alignment-unmapped")[0]).0, 4);

    // The same ten bases clipped at one end only are not a contaminant, however short they are.
    let one_end = [(10, 'M'), (30, 'S')];
    assert!(!is_contaminant(&one_end, 32));
    assert_eq!(shape(&records(&text, "clipped-at-one-end-only")[0]).0, 0);

    // And the threshold is the whole of the first test: lowered under the aligned length, the same
    // alignment is kept.
    assert!(!is_contaminant(&both_ends, 8));
    assert_eq!(
        shape(&records(&text, "a-short-alignment-under-a-lower-threshold")[0]).3,
        "15S10M15S"
    );
    // Without the flag at all, the contaminant is kept as it was.
    assert_eq!(
        shape(&records(&text, "a-short-alignment-kept")[0]),
        (0, 41, 60, "15S10M15S".to_string())
    );
}

/// The four strategies keep four different things.
#[test]
fn the_strategies_keep_different_things() {
    let text = corpus();
    let aligned = || Record {
        reference: Some("chr1".to_string()),
        start: 41,
        mapping_quality: 60,
        cigar: vec![(15, 'S'), (10, 'M'), (15, 'S')],
        unmapped: false,
        edit_distance: None,
        original_alignment: None,
        comment: None,
    };

    // The default keeps the coordinates and loses the cigar and the mapping quality.
    let mut record = aligned();
    unmap_contaminant(&mut record, UnmappingReadStrategy::DoNotChange);
    assert_eq!(
        (
            record.unmapped,
            record.reference.clone(),
            record.start,
            record.mapping_quality,
            cigar_string(&record.cigar)
        ),
        (true, Some("chr1".to_string()), 41, 0, "*".to_string())
    );
    assert_eq!(record.original_alignment, None);
    let golden = &records(&text, "a-short-alignment-unmapped")[0];
    assert_eq!(shape(golden), (4, 41, 0, "*".to_string()));
    assert!(golden
        .iter()
        .any(|field| field == &format!("CO:Z:{CONTAMINATION_COMMENT}")));

    // Copying to the tag keeps the coordinates and records the alignment beside them.
    let mut record = aligned();
    unmap_contaminant(&mut record, UnmappingReadStrategy::CopyToTag);
    assert_eq!(
        record.original_alignment.as_deref(),
        Some("chr1,41,15S10M15S,60,;")
    );
    let golden = &records(&text, "a-contaminant-copy_to_tag")[0];
    assert_eq!(shape(golden), (4, 41, 0, "*".to_string()));
    assert!(golden
        .iter()
        .any(|field| field == "OA:Z:chr1,41,15S10M15S,60,;"));

    // Moving to the tag takes the coordinates away as well.
    let mut record = aligned();
    unmap_contaminant(&mut record, UnmappingReadStrategy::MoveToTag);
    assert_eq!(record.reference, None);
    assert_eq!(record.start, 0);
    assert_eq!(
        record.original_alignment.as_deref(),
        Some("chr1,41,15S10M15S,60,;")
    );
    assert_eq!(
        shape(&records(&text, "a-contaminant-move_to_tag")[0]),
        (4, 0, 0, "*".to_string())
    );

    // And the invalid strategy leaves the cigar and the mapping quality on a record flagged
    // unmapped, which is what its name is about.
    let mut record = aligned();
    unmap_contaminant(&mut record, UnmappingReadStrategy::DoNotChangeInvalid);
    assert_eq!(record.mapping_quality, 60);
    assert_eq!(cigar_string(&record.cigar), "15S10M15S");
    assert_eq!(
        shape(&records(&text, "a-contaminant-do_not_change_invalid")[0]),
        (4, 41, 60, "15S10M15S".to_string())
    );

    // A record that already carries a comment keeps it, with the new one after a pipe.
    let mut record = Record {
        comment: Some("something else".to_string()),
        ..aligned()
    };
    unmap_contaminant(&mut record, UnmappingReadStrategy::DoNotChange);
    assert_eq!(
        record.comment.as_deref(),
        Some("something else | Cross-species contamination")
    );

    // The tag's format is the SA tag's, and a record with no edit distance leaves that field
    // empty rather than filling it in.
    assert_eq!(
        encode_mapping_information(&aligned()),
        "chr1,41,15S10M15S,60,;"
    );
}
