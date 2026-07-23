//! `CompareSAMs`.
//!
//! Ports `picard.sam.CompareSAMs` + `picard.sam.util.SamComparison` at tag 3.4.0, for the
//! **default strict mode over any single shared sort order**. Compares two SAM files and reports, in
//! a [`SamComparisonMetric`] row, how many primary alignments match, differ, are unmapped on one or
//! both sides, are missing on one side, or disagree on duplicate marking; the tool prints
//! `SAM files match.` or `SAM files differ.` accordingly.
//!
//! Scope of this slice: inputs that share a sort order, dispatched on it to the coordinate,
//! queryname, or (any other value) unsorted comparison path; strict comparison (`LENIENT_*` all
//! false); `COMPARE_MQ=false` (so no mapping-quality histogram). Each path matches reads by
//! `PrimaryAlignmentKey` rather than by position, so order within a coordinate (and, in the unsorted
//! path, order at all) does not matter; the per-type counts are commutative, so the port reproduces
//! htsjdk's totals without reproducing its exact `LinkedHashMap` iteration order. Header comparison
//! is a structural equality of the two headers, which is correct for equal headers; htsjdk's finer
//! field-by-field `compareHeaders` (and thus header-difference reporting), mismatched sort orders,
//! the lenient modes, and the mapping-quality histogram are deferred.
//!
//! Two Picard behaviours are reproduced deliberately. `alignmentsMatch` compares
//! `s1.getReadNegativeStrandFlag() == s1.getReadNegativeStrandFlag()` (the left record against
//! itself), so strand is effectively ignored: two reads at the same reference and start match even
//! on opposite strands. And when the left file has trailing records with no right counterpart,
//! htsjdk throws (it reads `it1.getCurrent()` after exhausting the iterator); the port instead
//! counts them as `MISSING_RIGHT`, a divergence on that malformed edge which the conformance corpus
//! does not cover.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// Why `CompareSAMs` could not run this slice.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompareError {
    Parse(ParseError),
    /// A sort order this slice does not yet handle (only `queryname` is supported here).
    UnsupportedSortOrder(String),
}

impl From<ParseError> for CompareError {
    fn from(e: ParseError) -> Self {
        CompareError::Parse(e)
    }
}

/// `picard.sam.SamComparisonMetric`, the one-row comparison result.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SamComparisonMetric {
    pub left_file: String,
    pub right_file: String,
    pub mappings_match: i64,
    pub mappings_differ: i64,
    pub unmapped_both: i64,
    pub unmapped_left: i64,
    pub unmapped_right: i64,
    pub missing_left: i64,
    pub missing_right: i64,
    pub duplicate_markings_differ: i64,
    pub are_equal: bool,
}

impl SamComparisonMetric {
    /// `allVisitedAlignmentsEqual`: none missing, differing, or unmapped on a single side.
    fn all_visited_equal(&self) -> bool {
        self.missing_left == 0
            && self.missing_right == 0
            && self.mappings_differ == 0
            && self.unmapped_left == 0
            && self.unmapped_right == 0
    }
}

impl MetricBean for SamComparisonMetric {
    fn class_name(&self) -> &str {
        "picard.sam.SamComparisonMetric"
    }

    fn columns(&self) -> &[&'static str] {
        &[
            "LEFT_FILE",
            "RIGHT_FILE",
            "MAPPINGS_MATCH",
            "MAPPINGS_DIFFER",
            "UNMAPPED_BOTH",
            "UNMAPPED_LEFT",
            "UNMAPPED_RIGHT",
            "MISSING_LEFT",
            "MISSING_RIGHT",
            "DUPLICATE_MARKINGS_DIFFER",
            "ARE_EQUAL",
        ]
    }

    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.left_file.clone()),
            Value::Str(self.right_file.clone()),
            Value::Long(self.mappings_match),
            Value::Long(self.mappings_differ),
            Value::Long(self.unmapped_both),
            Value::Long(self.unmapped_left),
            Value::Long(self.unmapped_right),
            Value::Long(self.missing_left),
            Value::Long(self.missing_right),
            Value::Long(self.duplicate_markings_differ),
            Value::Bool(self.are_equal),
        ]
    }
}

/// `PrimaryAlignmentKey`: a read's name plus its pairing status, ordered `UNPAIRED < FIRST < SECOND`.
fn primary_alignment_key(rec: &BamRecord) -> (String, u8) {
    let pair_status = if rec.flags & READ_PAIRED != 0 {
        if rec.flags & SECOND_OF_PAIR != 0 {
            2
        } else {
            1
        }
    } else {
        0
    };
    (rec.read_name.clone(), pair_status)
}

fn is_secondary_or_supplementary(rec: &BamRecord) -> bool {
    rec.flags & (SECONDARY | SUPPLEMENTARY) != 0
}

fn reference_name<'a>(header: &'a SamHeader, rec: &BamRecord) -> Option<&'a str> {
    if rec.reference_index < 0 {
        None
    } else {
        Some(header.sequences[rec.reference_index as usize].name.as_str())
    }
}

/// `alignmentsMatch` (strict): same reference and start. Strand is not compared, reproducing
/// htsjdk's `s1 == s1` self-comparison.
fn alignments_match(h1: &SamHeader, s1: &BamRecord, h2: &SamHeader, s2: &BamRecord) -> bool {
    reference_name(h1, s1) == reference_name(h2, s2) && s1.alignment_start == s2.alignment_start
}

/// `compareAlignmentRecords` followed by `updateMetric`, plus `catalogDuplicateDifferences`.
fn tally(
    h1: &SamHeader,
    s1: &BamRecord,
    h2: &SamHeader,
    s2: &BamRecord,
    m: &mut SamComparisonMetric,
) {
    if (s1.flags & DUPLICATE) != (s2.flags & DUPLICATE) {
        m.duplicate_markings_differ += 1;
    }
    let s1_unmapped = s1.flags & READ_UNMAPPED != 0;
    let s2_unmapped = s2.flags & READ_UNMAPPED != 0;
    if s1_unmapped && s2_unmapped {
        m.unmapped_both += 1;
    } else if s1_unmapped {
        m.unmapped_left += 1;
    } else if s2_unmapped {
        m.unmapped_right += 1;
    } else if alignments_match(h1, s1, h2, s2) {
        m.mappings_match += 1;
    } else {
        m.mappings_differ += 1;
    }
}

/// `CompareSAMs`/`SamComparison` for queryname-sorted, strict input. `left_name`/`right_name` fill
/// the `LEFT_FILE`/`RIGHT_FILE` columns.
pub fn compare_sams(
    left_sam: &str,
    right_sam: &str,
    left_name: &str,
    right_name: &str,
) -> Result<SamComparisonMetric, CompareError> {
    let (h1, recs1) = read_sam_with(left_sam, ValidationStringency::Silent)?;
    let (h2, recs2) = read_sam_with(right_sam, ValidationStringency::Silent)?;

    // htsjdk requires both files to share a sort order and dispatches on it: coordinate, queryname,
    // or (any other value) the unsorted path. Mismatched sort orders are deferred.
    let so1 = h1.attributes.get("SO").unwrap_or("unsorted");
    let so2 = h2.attributes.get("SO").unwrap_or("unsorted");
    if so1 != so2 {
        return Err(CompareError::UnsupportedSortOrder(format!("{so1}/{so2}")));
    }

    let mut m = SamComparisonMetric {
        left_file: left_name.to_string(),
        right_file: right_name.to_string(),
        ..Default::default()
    };

    let left: Vec<&BamRecord> = recs1
        .iter()
        .filter(|r| !is_secondary_or_supplementary(r))
        .collect();
    let right: Vec<&BamRecord> = recs2
        .iter()
        .filter(|r| !is_secondary_or_supplementary(r))
        .collect();

    match so1 {
        "queryname" => compare_queryname(&left, &right, &h1, &h2, &mut m),
        "coordinate" => compare_coordinate(&left, &right, &h1, &h2, &mut m),
        _ => compare_unsorted(&left, &right, &h1, &h2, &mut m),
    }

    let headers_equal = h1 == h2;
    m.are_equal = headers_equal && m.all_visited_equal() && m.duplicate_markings_differ == 0;
    Ok(m)
}

/// `compareQueryNameSortedAlignments`: a merge of the two name-ordered streams by `PrimaryAlignmentKey`.
fn compare_queryname(
    left: &[&BamRecord],
    right: &[&BamRecord],
    h1: &SamHeader,
    h2: &SamHeader,
    m: &mut SamComparisonMetric,
) {
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.len() {
        if j >= right.len() {
            // htsjdk throws here; the port counts the trailing lefts (see the module note).
            m.missing_right += (left.len() - i) as i64;
            break;
        }
        let lk = primary_alignment_key(left[i]);
        let rk = primary_alignment_key(right[j]);
        match lk.cmp(&rk) {
            std::cmp::Ordering::Less => {
                m.missing_right += 1;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                m.missing_left += 1;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                tally(h1, left[i], h2, right[j], m);
                i += 1;
                j += 1;
            }
        }
    }
    if j < right.len() {
        m.missing_left += (right.len() - j) as i64;
    }
}

/// `compareAlignmentCoordinates`: unmapped reads (no reference) sort last and equal to each other,
/// else order by reference index then alignment start.
fn compare_alignment_coordinates(left: &BamRecord, right: &BamRecord) -> std::cmp::Ordering {
    let lu = left.reference_index < 0;
    let ru = right.reference_index < 0;
    match (lu, ru) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) => left
            .reference_index
            .cmp(&right.reference_index)
            .then(left.alignment_start.cmp(&right.alignment_start)),
    }
}

/// `compareCoordinateSortedAlignments`: within each coordinate, reads are matched by
/// `PrimaryAlignmentKey` rather than by position, so order within a coordinate does not matter. The
/// per-type counts are commutative, so the port matches htsjdk's totals without reproducing its
/// exact `LinkedHashMap` iteration order.
fn compare_coordinate(
    left: &[&BamRecord],
    right: &[&BamRecord],
    h1: &SamHeader,
    h2: &SamHeader,
    m: &mut SamComparisonMetric,
) {
    use std::collections::HashMap;

    let (mut i, mut j) = (0usize, 0usize);
    let mut left_unmatched: HashMap<(String, u8), &BamRecord> = HashMap::new();
    let mut right_unmatched: HashMap<(String, u8), &BamRecord> = HashMap::new();

    while i < left.len() {
        if j >= right.len() {
            // Right exhausted: match remaining lefts against saved rights, else MISSING_RIGHT.
            for l in &left[i..] {
                match right_unmatched.remove(&primary_alignment_key(l)) {
                    Some(r) => tally(h1, l, h2, r, m),
                    None => m.missing_right += 1,
                }
            }
            break;
        }
        // Grab all lefts sharing this coordinate.
        let anchor = left[i];
        let mut left_here: HashMap<(String, u8), &BamRecord> = HashMap::new();
        while i < left.len()
            && compare_alignment_coordinates(anchor, left[i]) == std::cmp::Ordering::Equal
        {
            left_here.insert(primary_alignment_key(left[i]), left[i]);
            i += 1;
        }
        // Advance right past everything ordered before this coordinate, saving those rights.
        while j < right.len()
            && compare_alignment_coordinates(anchor, right[j]) == std::cmp::Ordering::Greater
        {
            right_unmatched.insert(primary_alignment_key(right[j]), right[j]);
            j += 1;
        }
        // Rights at this coordinate: match against left_here by key, else save.
        while j < right.len()
            && compare_alignment_coordinates(anchor, right[j]) == std::cmp::Ordering::Equal
        {
            let rk = primary_alignment_key(right[j]);
            match left_here.remove(&rk) {
                Some(l) => tally(h1, l, h2, right[j], m),
                None => {
                    right_unmatched.insert(rk, right[j]);
                }
            }
            j += 1;
        }
        // Unmatched lefts at this coordinate carry over.
        for (k, l) in left_here {
            left_unmatched.insert(k, l);
        }
    }

    // Remaining rights: match against saved lefts, else MISSING_LEFT.
    for r in &right[j..] {
        match left_unmatched.remove(&primary_alignment_key(r)) {
            Some(l) => tally(h1, l, h2, r, m),
            None => m.missing_left += 1,
        }
    }
    // Saved lefts: match against saved rights, else MISSING_RIGHT.
    let keys: Vec<(String, u8)> = left_unmatched.keys().cloned().collect();
    for k in keys {
        let l = left_unmatched[&k];
        match right_unmatched.remove(&k) {
            Some(r) => tally(h1, l, h2, r, m),
            None => m.missing_right += 1,
        }
    }
    m.missing_left += right_unmatched.len() as i64;
}

/// `compareUnsortedAlignments`: with no order assumed, index every left by `PrimaryAlignmentKey`,
/// match each right against it, and treat whatever is left over as missing. As above, the counts are
/// commutative, so a plain map matches htsjdk's totals.
fn compare_unsorted(
    left: &[&BamRecord],
    right: &[&BamRecord],
    h1: &SamHeader,
    h2: &SamHeader,
    m: &mut SamComparisonMetric,
) {
    use std::collections::HashMap;

    let mut left_unmatched: HashMap<(String, u8), &BamRecord> = HashMap::new();
    for l in left {
        left_unmatched.insert(primary_alignment_key(l), l);
    }
    for r in right {
        match left_unmatched.remove(&primary_alignment_key(r)) {
            Some(l) => tally(h1, l, h2, r, m),
            None => m.missing_left += 1,
        }
    }
    m.missing_right += left_unmatched.len() as i64;
}

/// The stdout line `CompareSAMs.doWork` prints.
pub fn verdict(metric: &SamComparisonMetric) -> &'static str {
    if metric.are_equal {
        "SAM files match."
    } else {
        "SAM files differ."
    }
}

/// `SamComparison.writeReport`: the metrics file (without the command-line/timestamp banner, which
/// the caller supplies). One `SamComparisonMetric` row, no histogram (`COMPARE_MQ=false`).
pub fn write_report(metric: &SamComparisonMetric) -> String {
    let mut mf = MetricsFile::new();
    mf.add_metric(metric);
    mf.write()
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: &str = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n";

    #[test]
    fn identical_files_match() {
        let sam = format!("{H}a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let m = compare_sams(&sam, &sam, "L", "R").unwrap();
        assert_eq!(m.mappings_match, 1);
        assert!(m.are_equal);
        assert_eq!(verdict(&m), "SAM files match.");
    }

    #[test]
    fn opposite_strand_still_matches() {
        // Reproduces the s1==s1 strand no-op: same ref+start, opposite strand, still a match.
        let l = format!("{H}a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let r = format!("{H}a\t16\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let m = compare_sams(&l, &r, "L", "R").unwrap();
        assert_eq!(m.mappings_match, 1);
        assert!(m.are_equal);
    }

    #[test]
    fn coordinate_matching_is_order_independent_within_a_coordinate() {
        // Two reads at the same coordinate, in opposite order across the files, still match.
        let h = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n";
        let p1 = "p\t99\tchr1\t10\t60\t4M\t=\t10\t0\tACGT\tIIII\n";
        let p2 = "p\t147\tchr1\t10\t60\t4M\t=\t10\t0\tACGT\tIIII\n";
        let l = format!("{h}{p1}{p2}");
        let r = format!("{h}{p2}{p1}");
        let m = compare_sams(&l, &r, "L", "R").unwrap();
        assert_eq!(m.mappings_match, 2);
        assert!(m.are_equal);
    }

    #[test]
    fn unsorted_matching_ignores_order() {
        let h = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n";
        let a = "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let b = "b\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let m = compare_sams(&format!("{h}{a}{b}"), &format!("{h}{b}{a}"), "L", "R").unwrap();
        assert_eq!(m.mappings_match, 2);
        assert!(m.are_equal);
    }

    #[test]
    fn differing_duplicate_marks_are_counted() {
        let l = format!("{H}a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let r = format!("{H}a\t1024\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let m = compare_sams(&l, &r, "L", "R").unwrap();
        assert_eq!(m.duplicate_markings_differ, 1);
        assert!(!m.are_equal);
        assert_eq!(verdict(&m), "SAM files differ.");
    }
}
