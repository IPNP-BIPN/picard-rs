//! `CompareSAMs`.
//!
//! Ports `picard.sam.CompareSAMs` + `picard.sam.util.SamComparison` at tag 3.4.0, for the
//! **default strict mode over any single shared sort order**. Compares two SAM files and reports, in
//! a [`SamComparisonMetric`] row, how many primary alignments match, differ, are unmapped on one or
//! both sides, are missing on one side, or disagree on duplicate marking; the tool prints
//! `SAM files match.` or `SAM files differ.` accordingly.
//!
//! Scope of this slice: inputs that share a sort order, dispatched on it to the coordinate,
//! queryname, or (any other value) unsorted comparison path. Each path matches reads by
//! `PrimaryAlignmentKey` rather than by position, so order within a coordinate (and, in the unsorted
//! path, order at all) does not matter; the per-type counts are commutative, so the port reproduces
//! htsjdk's totals without reproducing its exact `LinkedHashMap` iteration order. The headers are
//! compared field by field as `compareHeaders` does; the differences it prints are not reproduced,
//! only its verdict.
//!
//! Comparison options are carried in [`CompareOptions`], covering the mapping-quality features:
//! `LENIENT_LOW_MQ_ALIGNMENT` (count reads whose mapping quality is `<= LOW_MQ_THRESHOLD`, default 3,
//! in both files but mapped to different locations as matches), `LENIENT_UNKNOWN_MQ_ALIGNMENT` (the
//! same for the sentinel mapping quality 255), and `COMPARE_MQ` (emit a `## HISTOGRAM` of `mq1,mq2`
//! concordance counts for every matched pair), `LENIENT_HEADER` and `LENIENT_DUP` (the duplicate-set
//! swap logic). [`compare_sams`] runs the strict defaults; [`compare_sams_with_options`] takes an
//! explicit [`CompareOptions`]. Mismatched sort orders are deferred.
//!
//! Two Picard behaviours are reproduced deliberately. `alignmentsMatch` compares
//! `s1.getReadNegativeStrandFlag() == s1.getReadNegativeStrandFlag()` (the left record against
//! itself), so strand is effectively ignored: two reads at the same reference and start match even
//! on opposite strands. And when the left file has trailing records with no right counterpart,
//! htsjdk throws (it reads `it1.getCurrent()` after exhausting the iterator); the port instead
//! counts them as `MISSING_RIGHT`, a divergence on that malformed edge which the conformance corpus
//! does not cover.

use std::collections::BTreeMap;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use htsjdk_metrics::file::{Histogram, MetricBean, MetricsFile, Value};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// `SAMRecord.UNKNOWN_MAPPING_QUALITY`: the sentinel value for "no mapping quality available".
const UNKNOWN_MAPPING_QUALITY: u8 = 255;

/// `SAMComparisonArgumentCollection`: the comparison options this port covers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompareOptions {
    /// `LENIENT_LOW_MQ_ALIGNMENT`: count reads whose mapping quality is `<= low_mq_threshold` in
    /// both files but mapped to different locations as matches.
    pub lenient_low_mq_alignment: bool,
    /// `LENIENT_UNKNOWN_MQ_ALIGNMENT`: the same for reads with no mapping quality (value 255) in both.
    pub lenient_unknown_mq_alignment: bool,
    /// `LOW_MQ_THRESHOLD`: the ceiling for `lenient_low_mq_alignment`. Default 3.
    pub low_mq_threshold: i32,
    /// `COMPARE_MQ`: emit a `## HISTOGRAM` of `mq1,mq2` concordance counts for every matched pair.
    pub compare_mq: bool,
    /// `LENIENT_HEADER`: leave the program records and the sequences' species, assembly, MD5 and
    /// URI out of the header comparison.
    pub lenient_header: bool,
    /// `LENIENT_DUP`: a duplicate mark that differs only because the two files chose different
    /// representatives of one duplicate set is not counted.
    pub lenient_dup: bool,
}

impl Default for CompareOptions {
    fn default() -> Self {
        CompareOptions {
            lenient_low_mq_alignment: false,
            lenient_unknown_mq_alignment: false,
            low_mq_threshold: 3,
            compare_mq: false,
            lenient_header: false,
            lenient_dup: false,
        }
    }
}

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
#[derive(Debug, Clone, Default, PartialEq)]
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
    /// `SamComparison.mappingQualityHistogram`: accumulated `"mq1,mq2" -> count` bins for every
    /// matched pair, populated only when `COMPARE_MQ` is set. Not a metric column; carried here so
    /// [`write_report`] can emit it. Keyed in a `BTreeMap` so the bins come out in the same natural
    /// `String` order htsjdk's `TreeSet` comparator produces.
    pub mq_histogram: BTreeMap<String, i64>,
    /// Under `LENIENT_DUP`, the pairs whose duplicate marks differ, left record then right, kept
    /// for the duplicate-set pass instead of being counted. Not a metric column.
    pub differing_duplicates: Vec<(BamRecord, BamRecord)>,
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

/// `alignmentsMatch`: same reference and start (strand is not compared, reproducing htsjdk's
/// `s1 == s1` self-comparison), OR the lenient-MQ escape hatches: both mapping qualities
/// `<= low_mq_threshold` under `lenient_low_mq_alignment`, or both the unknown-MQ sentinel 255 under
/// `lenient_unknown_mq_alignment`.
fn alignments_match(
    h1: &SamHeader,
    s1: &BamRecord,
    h2: &SamHeader,
    s2: &BamRecord,
    opts: &CompareOptions,
) -> bool {
    (reference_name(h1, s1) == reference_name(h2, s2) && s1.alignment_start == s2.alignment_start)
        || (opts.lenient_low_mq_alignment
            && s1.mapping_quality as i32 <= opts.low_mq_threshold
            && s2.mapping_quality as i32 <= opts.low_mq_threshold)
        || (opts.lenient_unknown_mq_alignment
            && s1.mapping_quality == UNKNOWN_MAPPING_QUALITY
            && s2.mapping_quality == UNKNOWN_MAPPING_QUALITY)
}

/// `tallyAlignmentRecords`: `catalogDuplicateDifferences`, then `compareAlignmentRecords` +
/// `updateMetric`, then (under `COMPARE_MQ`) `compareAndUpdateMappingQualityConcordance`.
fn tally(
    h1: &SamHeader,
    s1: &BamRecord,
    h2: &SamHeader,
    s2: &BamRecord,
    opts: &CompareOptions,
    m: &mut SamComparisonMetric,
) {
    // `catalogDuplicateDifferences`: counted now, or kept for the lenient pass.
    if (s1.flags & DUPLICATE) != (s2.flags & DUPLICATE) {
        if opts.lenient_dup {
            m.differing_duplicates.push((s1.clone(), s2.clone()));
        } else {
            m.duplicate_markings_differ += 1;
        }
    }
    let s1_unmapped = s1.flags & READ_UNMAPPED != 0;
    let s2_unmapped = s2.flags & READ_UNMAPPED != 0;
    if s1_unmapped && s2_unmapped {
        m.unmapped_both += 1;
    } else if s1_unmapped {
        m.unmapped_left += 1;
    } else if s2_unmapped {
        m.unmapped_right += 1;
    } else if alignments_match(h1, s1, h2, s2, opts) {
        m.mappings_match += 1;
    } else {
        m.mappings_differ += 1;
    }
    if opts.compare_mq {
        // `String.format("%d,%d", s1.getMappingQuality(), s2.getMappingQuality())`.
        let key = format!("{},{}", s1.mapping_quality, s2.mapping_quality);
        *m.mq_histogram.entry(key).or_insert(0) += 1;
    }
}

/// `CompareSAMs`/`SamComparison` with the strict default options. `left_name`/`right_name` fill the
/// `LEFT_FILE`/`RIGHT_FILE` columns.
pub fn compare_sams(
    left_sam: &str,
    right_sam: &str,
    left_name: &str,
    right_name: &str,
) -> Result<SamComparisonMetric, CompareError> {
    compare_sams_with_options(
        left_sam,
        right_sam,
        left_name,
        right_name,
        &CompareOptions::default(),
    )
}

/// `CompareSAMs`/`SamComparison` with explicit comparison [`CompareOptions`].
pub fn compare_sams_with_options(
    left_sam: &str,
    right_sam: &str,
    left_name: &str,
    right_name: &str,
    opts: &CompareOptions,
) -> Result<SamComparisonMetric, CompareError> {
    let (h1, recs1) = read_sam_with(left_sam, ValidationStringency::Silent)?;
    let (h2, recs2) = read_sam_with(right_sam, ValidationStringency::Silent)?;
    compare_records(
        &h1,
        recs1.into_iter(),
        &h2,
        recs2.into_iter(),
        left_name,
        right_name,
        opts,
    )
}

/// [`compare_sams_with_options`] over records read as they are compared.
pub fn compare_records(
    h1: &SamHeader,
    recs1: impl Iterator<Item = BamRecord>,
    h2: &SamHeader,
    recs2: impl Iterator<Item = BamRecord>,
    left_name: &str,
    right_name: &str,
    opts: &CompareOptions,
) -> Result<SamComparisonMetric, CompareError> {
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

    let left = recs1.filter(|r| !is_secondary_or_supplementary(r));
    let right = recs2.filter(|r| !is_secondary_or_supplementary(r));

    let (headers_equal, dictionaries_differ) = compare_headers(h1, h2, opts.lenient_header);
    // A coordinate comparison needs one dictionary; with two it compares nothing and fails.
    let alignments_equal = match so1 {
        "queryname" => {
            compare_queryname(left, right, h1, h2, opts, &mut m);
            m.all_visited_equal()
        }
        "coordinate" if dictionaries_differ => false,
        "coordinate" => {
            compare_coordinate(left, right, h1, h2, opts, &mut m);
            m.all_visited_equal()
        }
        _ => {
            compare_unsorted(left, right, h1, h2, opts, &mut m);
            m.all_visited_equal()
        }
    };
    if opts.lenient_dup {
        let pairs = std::mem::take(&mut m.differing_duplicates);
        m.duplicate_markings_differ += lenient_duplicate_differences(h1, h2, &pairs);
    }
    m.are_equal = headers_equal && alignments_equal && m.duplicate_markings_differ == 0;
    Ok(m)
}

/// `compareHeaders`: whether the two headers agree, and whether their sequence dictionaries differ
/// (which stops a coordinate comparison). Each check reports and goes on, as the reference does.
fn compare_headers(h1: &SamHeader, h2: &SamHeader, lenient: bool) -> (bool, bool) {
    let hd = |h: &SamHeader, key: &str| h.attributes.get(key).map(str::to_string);
    let mut equal = hd(h1, "VN") == hd(h2, "VN");
    equal &= hd(h1, "CR") == hd(h2, "CR");
    equal &= hd(h1, "SO") == hd(h2, "SO");

    // `compareSequenceDictionaries`.
    let dictionaries_equal = if h1.sequences.len() != h2.sequences.len() {
        false
    } else {
        let mut all = true;
        for (s1, s2) in h1.sequences.iter().zip(&h2.sequences) {
            if s1.name != s2.name {
                all = false;
                continue;
            }
            all &= s1.length == s2.length;
            if !lenient {
                for key in ["SP", "AS", "M5", "UR"] {
                    all &= s1.attributes.get(key) == s2.attributes.get(key);
                }
            }
        }
        all
    };
    equal &= dictionaries_equal;

    // `compareReadGroups`.
    if h1.read_groups.len() != h2.read_groups.len() {
        equal = false;
    } else {
        for (g1, g2) in h1.read_groups.iter().zip(&h2.read_groups) {
            equal &= g1.id == g2.id;
            for key in ["SM", "LB", "DS", "PU", "PI", "CN", "DT", "PL"] {
                equal &= g1.attributes.get(key) == g2.attributes.get(key);
            }
        }
    }

    // `compareProgramRecords`, which LENIENT_HEADER skips.
    if !lenient {
        if h1.programs.len() != h2.programs.len() {
            equal = false;
        } else {
            for (p1, p2) in h1.programs.iter().zip(&h2.programs) {
                equal &= p1.id == p2.id;
                for key in ["VN", "CL"] {
                    equal &= p1.attributes.get(key) == p2.attributes.get(key);
                }
            }
        }
    }
    (equal, !dictionaries_equal)
}

/// A record as the duplicate comparator reads it.
fn duplicate_record(header: &SamHeader, rec: &BamRecord) -> crate::mark_duplicates::Record {
    use htsjdk_bam::tag::{Tag, TagValue};
    let group = match rec.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => header.read_groups.iter().position(|g| &g.id == id),
        _ => None,
    };
    let library = group
        .and_then(|index| header.read_groups[index].attributes.get("LB"))
        .unwrap_or("Unknown Library")
        .to_string();
    crate::mark_duplicates::Record {
        name: rec.read_name.clone(),
        flags: rec.flags,
        reference_index: rec.reference_index,
        alignment_start: rec.alignment_start,
        span: crate::mark_duplicates::Span::of(&rec.cigar),
        read_length: rec.read_bases.len() as u32,
        base_quality_score: crate::mark_duplicates::base_quality_score(&rec.base_qualities),
        mate_reference_index: rec.mate_reference_index,
        library: library.into(),
        read_group: group.map_or(-1, |index| index as i32),
        barcode: None,
        existing_dt: None,
        mate_span: match rec.tags.get(Tag::new(b"MC")) {
            Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text)
                .ok()
                .as_ref()
                .map(crate::mark_duplicates::Span::of),
            _ => None,
        },
        mate_alignment_start: rec.mate_alignment_start,
    }
}

/// `countLenientDuplicateMarkingDifferences`: the records whose marks differ are cut into duplicate
/// sets on each side. A left set's unmarked reads may swap with its marked ones; a right set's
/// unmarked read that finds such a partner, marked on the right, is matched with it, and neither
/// is counted. Every unmatched record is.
fn lenient_duplicate_differences(
    h1: &SamHeader,
    h2: &SamHeader,
    pairs: &[(BamRecord, BamRecord)],
) -> i64 {
    use crate::mark_duplicates::{Options, ScoringStrategy};
    use std::collections::{HashMap, HashSet};
    // `SAMRecordDuplicateComparator`'s own default strategy, over the mate's CIGAR.
    let options = Options {
        scoring: ScoringStrategy::TotalMappedReferenceLength,
        assume_mate_cigar: true,
        ..Options::default()
    };
    let left: Vec<_> = pairs.iter().map(|(l, _)| duplicate_record(h1, l)).collect();
    let right: Vec<_> = pairs.iter().map(|(_, r)| duplicate_record(h2, r)).collect();

    let mut swap_targets: HashMap<String, HashSet<String>> = HashMap::new();
    for set in crate::duplicate_set::duplicate_sets(&left, &options) {
        if set.len() > 1 {
            let marked: HashSet<String> = set
                .iter()
                .filter(|i| left[**i].flags & DUPLICATE != 0)
                .map(|i| left[*i].name.clone())
                .collect();
            for i in set.iter().filter(|i| left[**i].flags & DUPLICATE == 0) {
                swap_targets.insert(left[*i].name.clone(), marked.clone());
            }
        }
    }
    let mut matched: HashSet<String> = HashSet::new();
    let mut count = 0i64;
    for set in crate::duplicate_set::duplicate_sets(&right, &options) {
        if set.len() > 1 {
            let candidates: Vec<&str> = set
                .iter()
                .filter(|i| right[**i].flags & DUPLICATE != 0)
                .map(|i| right[*i].name.as_str())
                .filter(|name| swap_targets.contains_key(*name))
                .collect();
            for i in set.iter().filter(|i| right[**i].flags & DUPLICATE == 0) {
                let representative = right[*i].name.as_str();
                if matched.contains(representative) {
                    continue;
                }
                if let Some(partner) = candidates.iter().find(|candidate| {
                    !matched.contains(**candidate)
                        && swap_targets[**candidate].contains(representative)
                }) {
                    matched.insert((*partner).to_string());
                    matched.insert(representative.to_string());
                }
            }
        }
        count += set
            .iter()
            .filter(|i| !matched.contains(right[**i].name.as_str()))
            .count() as i64;
    }
    count
}

/// `compareQueryNameSortedAlignments`: a merge of the two name-ordered streams by `PrimaryAlignmentKey`.
fn compare_queryname(
    mut left: impl Iterator<Item = BamRecord>,
    right: impl Iterator<Item = BamRecord>,
    h1: &SamHeader,
    h2: &SamHeader,
    opts: &CompareOptions,
    m: &mut SamComparisonMetric,
) {
    let mut right = right.peekable();
    let mut next_left = left.next();
    while let Some(l) = &next_left {
        let Some(r) = right.peek() else {
            // htsjdk throws here; the port counts the trailing lefts (see the module note).
            m.missing_right += 1 + left.count() as i64;
            break;
        };
        match primary_alignment_key(l).cmp(&primary_alignment_key(r)) {
            std::cmp::Ordering::Less => {
                m.missing_right += 1;
                next_left = left.next();
            }
            std::cmp::Ordering::Greater => {
                m.missing_left += 1;
                right.next();
            }
            std::cmp::Ordering::Equal => {
                tally(h1, l, h2, r, opts, m);
                next_left = left.next();
                right.next();
            }
        }
    }
    m.missing_left += right.count() as i64;
}

/// `compareAlignmentCoordinates`: unmapped reads (no reference) sort last and equal to each other,
/// else order by reference index then alignment start.
fn compare_alignment_coordinates(left: (i32, i32), right: &BamRecord) -> std::cmp::Ordering {
    let lu = left.0 < 0;
    let ru = right.reference_index < 0;
    match (lu, ru) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) => left
            .0
            .cmp(&right.reference_index)
            .then(left.1.cmp(&right.alignment_start)),
    }
}

/// `compareCoordinateSortedAlignments`: within each coordinate, reads are matched by
/// `PrimaryAlignmentKey` rather than by position, so order within a coordinate does not matter. The
/// per-type counts are commutative, so the port matches htsjdk's totals without reproducing its
/// exact `LinkedHashMap` iteration order.
fn compare_coordinate(
    left: impl Iterator<Item = BamRecord>,
    right: impl Iterator<Item = BamRecord>,
    h1: &SamHeader,
    h2: &SamHeader,
    opts: &CompareOptions,
    m: &mut SamComparisonMetric,
) {
    use std::cmp::Ordering::{Equal, Greater};
    use std::collections::HashMap;

    let (mut left, mut right) = (left.peekable(), right.peekable());
    let mut left_unmatched: HashMap<(String, u8), BamRecord> = HashMap::new();
    let mut right_unmatched: HashMap<(String, u8), BamRecord> = HashMap::new();

    while let Some(anchor) = left.next() {
        if right.peek().is_none() {
            // Right exhausted: match remaining lefts against saved rights, else MISSING_RIGHT.
            for l in std::iter::once(anchor).chain(left.by_ref()) {
                match right_unmatched.remove(&primary_alignment_key(&l)) {
                    Some(r) => tally(h1, &l, h2, &r, opts, m),
                    None => m.missing_right += 1,
                }
            }
            break;
        }
        // Grab all lefts sharing this coordinate.
        let here = (anchor.reference_index, anchor.alignment_start);
        let mut left_here: HashMap<(String, u8), BamRecord> = HashMap::new();
        left_here.insert(primary_alignment_key(&anchor), anchor);
        while let Some(l) = left.next_if(|l| compare_alignment_coordinates(here, l) == Equal) {
            left_here.insert(primary_alignment_key(&l), l);
        }
        // Advance right past everything ordered before this coordinate, saving those rights.
        while let Some(r) = right.next_if(|r| compare_alignment_coordinates(here, r) == Greater) {
            right_unmatched.insert(primary_alignment_key(&r), r);
        }
        // Rights at this coordinate: match against left_here by key, else save.
        while let Some(r) = right.next_if(|r| compare_alignment_coordinates(here, r) == Equal) {
            let rk = primary_alignment_key(&r);
            match left_here.remove(&rk) {
                Some(l) => tally(h1, &l, h2, &r, opts, m),
                None => {
                    right_unmatched.insert(rk, r);
                }
            }
        }
        // Unmatched lefts at this coordinate carry over.
        left_unmatched.extend(left_here);
    }

    // Remaining rights: match against saved lefts, else MISSING_LEFT.
    for r in right {
        match left_unmatched.remove(&primary_alignment_key(&r)) {
            Some(l) => tally(h1, &l, h2, &r, opts, m),
            None => m.missing_left += 1,
        }
    }
    // Saved lefts: match against saved rights, else MISSING_RIGHT.
    for (k, l) in &left_unmatched {
        match right_unmatched.remove(k) {
            Some(r) => tally(h1, l, h2, &r, opts, m),
            None => m.missing_right += 1,
        }
    }
    m.missing_left += right_unmatched.len() as i64;
}

/// `compareUnsortedAlignments`: with no order assumed, index every left by `PrimaryAlignmentKey`,
/// match each right against it, and treat whatever is left over as missing. As above, the counts are
/// commutative, so a plain map matches htsjdk's totals.
fn compare_unsorted(
    left: impl Iterator<Item = BamRecord>,
    right: impl Iterator<Item = BamRecord>,
    h1: &SamHeader,
    h2: &SamHeader,
    opts: &CompareOptions,
    m: &mut SamComparisonMetric,
) {
    use std::collections::HashMap;

    let mut left_unmatched: HashMap<(String, u8), BamRecord> = HashMap::new();
    for l in left {
        left_unmatched.insert(primary_alignment_key(&l), l);
    }
    for r in right {
        match left_unmatched.remove(&primary_alignment_key(&r)) {
            Some(l) => tally(h1, &l, h2, &r, opts, m),
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
/// the caller supplies). One `SamComparisonMetric` row, followed by the mapping-quality concordance
/// `## HISTOGRAM` when `COMPARE_MQ` populated it. htsjdk adds the histogram whenever `COMPARE_MQ` is
/// set, but `MetricsFile.write` drops an empty histogram, so gating on non-empty is equivalent.
pub fn write_report(metric: &SamComparisonMetric) -> String {
    let mut mf = MetricsFile::new();
    mf.add_metric(metric);
    if !metric.mq_histogram.is_empty() {
        mf.histograms.push(Histogram {
            bin_label: "BIN".to_string(),
            value_label: "VALUE".to_string(),
            key_class: "java.lang.String".to_string(),
            bins: metric
                .mq_histogram
                .iter()
                .map(|(k, v)| (k.clone(), *v as f64))
                .collect(),
        });
    }
    mf.write()
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: &str = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n";

    /// A read whose `QUAL` is `*` has bases and no qualities: its length is `getReadLength`, the
    /// bases, which is what the mate-cigar window and `toString` read, and its quality score is
    /// the sum over none.
    #[test]
    fn a_read_without_qualities_keeps_its_length() {
        let rec = htsjdk_bam::text_parse::parse_line(
            "a\t0\tchr1\t10\t60\t6M\t*\t0\t0\tACGTAC\t*",
            |name| (name == "chr1").then_some(0),
        )
        .unwrap();
        assert!(rec.base_qualities.is_empty());
        let record = duplicate_record(&SamHeader::new(), &rec);
        assert_eq!(record.read_length, 6);
        assert_eq!(record.base_quality_score, 0);
        assert_eq!(
            crate::mate_cigar_duplicates::describe_with_contig(&record, Some("chr1")),
            "a 6b aligned to chr1:10-15."
        );
    }

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
    fn lenient_low_mq_counts_differing_positions_as_matches() {
        // Two reads at MQ=2 (<= LOW_MQ_THRESHOLD=3) at different positions: differ strictly, match
        // leniently.
        let l = format!("{H}a\t0\tchr1\t10\t2\t4M\t*\t0\t0\tACGT\tIIII\n");
        let r = format!("{H}a\t0\tchr1\t50\t2\t4M\t*\t0\t0\tACGT\tIIII\n");
        assert_eq!(compare_sams(&l, &r, "L", "R").unwrap().mappings_differ, 1);
        let opts = CompareOptions {
            lenient_low_mq_alignment: true,
            ..Default::default()
        };
        let m = compare_sams_with_options(&l, &r, "L", "R", &opts).unwrap();
        assert_eq!(m.mappings_match, 1);
        assert_eq!(m.mappings_differ, 0);
        assert!(m.are_equal);
    }

    #[test]
    fn lenient_unknown_mq_counts_differing_positions_as_matches() {
        let l = format!("{H}a\t0\tchr1\t10\t255\t4M\t*\t0\t0\tACGT\tIIII\n");
        let r = format!("{H}a\t0\tchr1\t50\t255\t4M\t*\t0\t0\tACGT\tIIII\n");
        assert_eq!(compare_sams(&l, &r, "L", "R").unwrap().mappings_differ, 1);
        let opts = CompareOptions {
            lenient_unknown_mq_alignment: true,
            ..Default::default()
        };
        let m = compare_sams_with_options(&l, &r, "L", "R", &opts).unwrap();
        assert_eq!(m.mappings_match, 1);
        assert!(m.are_equal);
    }

    #[test]
    fn compare_mq_records_a_concordance_histogram() {
        let l = format!("{H}a\t0\tchr1\t10\t40\t4M\t*\t0\t0\tACGT\tIIII\n");
        let r = format!("{H}a\t0\tchr1\t10\t30\t4M\t*\t0\t0\tACGT\tIIII\n");
        let opts = CompareOptions {
            compare_mq: true,
            ..Default::default()
        };
        let m = compare_sams_with_options(&l, &r, "L", "R", &opts).unwrap();
        assert_eq!(m.mq_histogram.get("40,30"), Some(&1));
        let report = write_report(&m);
        assert!(report.contains("## HISTOGRAM\tjava.lang.String\nBIN\tVALUE\n40,30\t1\n"));
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
