//! `UmiAwareMarkDuplicatesWithMateCigar`, end to end over a file's records.
//!
//! [`crate::umi_duplicates`] holds the pieces the conformance suite measured (the clustering, the
//! assigned UMI, the metrics arithmetic) over sets cut by position. The tool itself is
//! `SimpleMarkDuplicatesWithMateCigar` with htsjdk's `DuplicateSetIterator` wrapped in Picard's
//! `UmiAwareDuplicateSetIterator`, and what it writes depends on three things a position-cut model
//! does not have:
//!
//!  * **the duplicate sets are htsjdk's**: the whole file is re-sorted by
//!    `SAMRecordDuplicateComparator` (library, contig, unclipped 5' end, orientation, mate,
//!    mapped ends, score, name, end of pair) and cut wherever a record stops being
//!    `duplicateSetCompare`-equal to the set's representative. An unmapped mate placed beside its
//!    read lands in that read's set; an unmapped or secondary representative is a set of one;
//!  * **the UMI sub-sets come out of two Java `HashMap`s**: the UMIs are numbered in the iteration
//!    order of a `HashMap<String, Long>` and the sub-sets are emitted in the iteration order of a
//!    `HashMap<Integer, List>` keyed by union-find roots, and both orders decide which read a
//!    duplicate flag lands on and in what order the writer is fed;
//!  * **the flags only ever go up**: `DuplicateSet.sort` sets the flag on every record whose name
//!    differs from the representative's and clears it on the first record, once for the
//!    position's set and again for each UMI sub-set, and never clears anything else. A flag the
//!    input carried survives unless its record ends up first.
//!
//! The sub-sets are built with `new DuplicateSet()`, whose comparator is the default one
//! (`TOTAL_MAPPED_REFERENCE_LENGTH`), but `DuplicateScoringStrategy.computeDuplicateScore` caches
//! the first score it computes for a record in a transient attribute. Every record the sub-set
//! sort can reach the score step for was already scored by the outer sort (a comparison sort
//! compares every adjacent pair of its output, and the score step is only reached for records
//! that tie on everything before it), so the sub-sets are ordered by the USER's strategy and
//! this port scores once.
//!
//! The UMI metrics keep one quirk on purpose: the "first UMI seen" flag that fixes
//! `MEAN_UMI_LENGTH` belongs to the ITERATOR, not to the library, so a second library starts at a
//! mean length of zero and is refused as "UMIs of differing lengths were found.".
//!
//! Ported from `picard.sam.markduplicates.UmiAwareMarkDuplicatesWithMateCigar`,
//! `UmiAwareDuplicateSetIterator`, `UmiGraph`, `UmiUtil`, `UmiMetrics`,
//! `SimpleMarkDuplicatesWithMateCigar`, `picard.util.GraphUtils` (Picard 3.4.0) and
//! `htsjdk.samtools.DuplicateSetIterator`, `DuplicateSet`, `SAMRecordDuplicateComparator`,
//! `DuplicateScoringStrategy`, `SAMUtils`, `util.StringUtil` (htsjdk 4.2.0).

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::BTreeSet;

use htsjdk_bam::cigar::{Cigar, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::murmur3::Murmur3;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::java_hash_map::JavaHashMap;

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const SECONDARY: u16 = 0x100;
const VENDOR_FAILED: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// `SAMRecordDuplicateComparator.UNKNOWN_LIBRARY_STRING`.
pub const UNKNOWN_LIBRARY: &str = "Unknown Library";

/// `DuplicateScoringStrategy.ScoringStrategy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scoring {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
    Random,
}

/// The arguments the run reads.
#[derive(Debug, Clone)]
pub struct UmiArgs {
    pub scoring: Scoring,
    pub max_edit_distance_to_join: i32,
    pub umi_tag: String,
    pub molecular_identifier_tag: Option<String>,
    pub allow_missing_umis: bool,
    pub duplex_umi: bool,
    pub remove_duplicates: bool,
    /// `PROGRAM_RECORD_ID`: where set, every written record's `PG` is replaced by its chained id,
    /// and since this tool never sees a program id the chain is empty and the tag is removed.
    pub program_record_id: Option<String>,
}

/// An uncaught Java throwable: its class and its message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Thrown {
    pub class: &'static str,
    pub message: String,
}

impl Thrown {
    fn new(class: &'static str, message: impl Into<String>) -> Self {
        Thrown {
            class,
            message: message.into(),
        }
    }
    fn picard(message: impl Into<String>) -> Self {
        Self::new("picard.PicardException", message)
    }
}

/// One `UmiMetrics` row.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UmiMetricsRow {
    /// `null` where the read group has no `LB`.
    pub library: Option<String>,
    pub mean_umi_length: f64,
    pub observed_unique_umis: i64,
    pub inferred_unique_umis: i64,
    pub observed_base_errors: i64,
    pub duplicate_sets_ignoring_umi: i64,
    pub duplicate_sets_with_umi: i64,
    pub observed_umi_entropy: f64,
    pub inferred_umi_entropy: f64,
    pub umi_base_qualities: f64,
    pub pct_umi_with_n: f64,
}

/// `UmiMetrics` with its private accumulators.
#[derive(Debug, Clone, Default)]
struct Accumulator {
    row: UmiMetricsRow,
    observed: BTreeMap<String, f64>,
    inferred: BTreeMap<String, f64>,
    observed_umi_bases: i64,
    observed_with_n: i64,
    observed_without_n: i64,
}

impl Accumulator {
    /// `addUmiObservation`: the observed length counts the hyphen, unlike `MEAN_UMI_LENGTH`.
    fn add_observation(&mut self, observed: &str, inferred: &str) {
        *self.observed.entry(observed.to_string()).or_insert(0.0) += 1.0;
        *self.inferred.entry(inferred.to_string()).or_insert(0.0) += 1.0;
        self.observed_umi_bases += observed.encode_utf16().count() as i64;
        self.observed_without_n += 1;
    }

    /// `calculateDerivedFields`.
    fn finish(mut self) -> UmiMetricsRow {
        self.row.observed_unique_umis = self.observed.len() as i64;
        self.row.inferred_unique_umis = self.inferred.len() as i64;
        self.row.pct_umi_with_n = self.observed_with_n as f64
            / (self.observed_with_n as f64 + self.observed_without_n as f64);
        self.row.observed_umi_entropy = effective_number_of_bases(&self.observed);
        self.row.inferred_umi_entropy = effective_number_of_bases(&self.inferred);
        self.row.umi_base_qualities = htsjdk_bam::quality_util::phred_score_from_error_probability(
            self.row.observed_base_errors as f64 / self.observed_umi_bases as f64,
        ) as f64;
        self.row
    }
}

/// `UmiMetrics.effectiveNumberOfBases`: the entropy summed by `Collectors.summingDouble`, which is
/// Kahan-compensated (and, since JDK-8214761, subtracts the compensation at the end), over the
/// histogram's bins in key order, divided by `Math.log(4.0)`.
fn effective_number_of_bases(histogram: &BTreeMap<String, f64>) -> f64 {
    let total: f64 = histogram.values().fold(0.0, |sum, value| sum + value);
    let mut sum = 0.0f64;
    let mut compensation = 0.0f64;
    let mut simple = 0.0f64;
    for value in histogram.values() {
        let p = value / total;
        let term = -p * jmath::math::log(p);
        let tmp = term - compensation;
        let velvel = sum + tmp;
        compensation = (velvel - sum) - tmp;
        sum = velvel;
        simple += term;
    }
    let mut entropy = sum - compensation;
    if entropy.is_nan() && simple.is_infinite() {
        entropy = simple;
    }
    entropy / jmath::math::log(4.0)
}

/// What a successful run produced.
#[derive(Debug, Clone)]
pub struct UmiRun {
    /// The records in the order the writer receives them, before the writer's coordinate sort.
    pub written: Vec<BamRecord>,
    /// The `UmiMetrics` rows in the iteration order of the library map.
    pub metrics: Vec<UmiMetricsRow>,
}

/// Per-record values the comparator reads, computed once as htsjdk's transient attributes are.
struct Keys {
    library_id: i32,
    read_coordinate: i32,
    mate_coordinate: i32,
    score: i16,
    canonical_name: String,
}

fn flag(record: &BamRecord, bit: u16) -> bool {
    record.flags & bit != 0
}

fn paired_and_both_mapped(record: &BamRecord) -> bool {
    flag(record, PAIRED) && !flag(record, UNMAPPED) && !flag(record, MATE_UNMAPPED)
}

fn secondary_or_supplementary(record: &BamRecord) -> bool {
    flag(record, SECONDARY) || flag(record, SUPPLEMENTARY)
}

fn string_tag<'r>(record: &'r BamRecord, name: &str) -> Result<Option<&'r str>, Thrown> {
    let tag = tag_of(name)?;
    match record.tags.get(tag) {
        None => Ok(None),
        Some(TagValue::Str(value)) => Ok(Some(value)),
        Some(other) => Err(Thrown::new(
            "htsjdk.samtools.SAMException",
            format!(
                "Value for tag {name} is not a String: class {}",
                java_class(other)
            ),
        )),
    }
}

fn java_class(value: &TagValue) -> &'static str {
    match value {
        TagValue::Char(_) => "java.lang.Character",
        TagValue::Int(_) => "java.lang.Integer",
        TagValue::Float(_) => "java.lang.Float",
        TagValue::Str(_) => "java.lang.String",
        _ => "[B",
    }
}

fn tag_of(name: &str) -> Result<Tag, Thrown> {
    let bytes: [u8; 2] = name.as_bytes().try_into().map_err(|_| {
        Thrown::new(
            "java.lang.IllegalArgumentException",
            format!("String tag does not have length() == 2: {name}"),
        )
    })?;
    Ok(Tag::new(&bytes))
}

fn mate_cigar(record: &BamRecord) -> Option<Cigar> {
    match record.tags.get(Tag::new(b"MC")) {
        Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
        _ => None,
    }
}

/// `SAMRecord.getAlignmentEnd`, which is `0` for an unmapped read.
fn alignment_end(record: &BamRecord) -> i32 {
    if flag(record, UNMAPPED) {
        0
    } else {
        record.alignment_start + record.cigar.reference_length() as i32 - 1
    }
}

fn leading_clips(cigar: &Cigar) -> i32 {
    let mut clipped = 0;
    for element in &cigar.elements {
        match element.op {
            Op::S | Op::H => clipped += element.length as i32,
            _ => break,
        }
    }
    clipped
}

fn trailing_clips(cigar: &Cigar) -> i32 {
    let mut clipped = 0;
    for element in cigar.elements.iter().rev() {
        match element.op {
            Op::S | Op::H => clipped += element.length as i32,
            _ => break,
        }
    }
    clipped
}

fn unclipped_start(record: &BamRecord) -> i32 {
    record.alignment_start - leading_clips(&record.cigar)
}

fn unclipped_end(record: &BamRecord) -> i32 {
    alignment_end(record) + trailing_clips(&record.cigar)
}

/// `SAMRecord.toString`, which htsjdk appends to its mate-cigar refusal.
pub fn describe(record: &BamRecord, header: &SamHeader) -> String {
    let mut out = record.read_name.clone();
    if flag(record, PAIRED) {
        out.push_str(if flag(record, FIRST_OF_PAIR) {
            " 1/2"
        } else {
            " 2/2"
        });
    }
    out.push_str(&format!(" {}b", record.read_bases.len()));
    if flag(record, UNMAPPED) {
        out.push_str(" unmapped read.");
    } else {
        out.push_str(&format!(
            " aligned to {}:{}-{}.",
            contig(record, header).unwrap_or("null"),
            record.alignment_start,
            alignment_end(record)
        ));
    }
    out
}

/// `SAMRecord.getContig`: `null` for an unmapped read.
fn contig<'h>(record: &BamRecord, header: &'h SamHeader) -> Option<&'h str> {
    if flag(record, UNMAPPED) {
        return None;
    }
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.as_str())
}

/// `SAMUtils.getMateUnclippedStart` / `getMateUnclippedEnd`, which refuse a record without `MC`.
fn mate_five_prime(record: &BamRecord, header: &SamHeader) -> Result<i32, Thrown> {
    let Some(cigar) = mate_cigar(record) else {
        return Err(Thrown::new(
            "htsjdk.samtools.SAMException",
            format!(
                "Mate CIGAR (Tag MC) not found: {}",
                describe(record, header)
            ),
        ));
    };
    Ok(if flag(record, MATE_REVERSE) {
        let end = record.mate_alignment_start + cigar.reference_length() as i32 - 1;
        end + trailing_clips(&cigar)
    } else {
        record.mate_alignment_start - leading_clips(&cigar)
    })
}

/// The read group a record names, if the header has it.
fn read_group<'h>(
    record: &BamRecord,
    header: &'h SamHeader,
) -> Option<&'h htsjdk_bam::header::ReadGroup> {
    match record.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => header.read_groups.iter().find(|group| &group.id == id),
        _ => None,
    }
}

/// `SAMRecordDuplicateComparator.getLibraryName`.
fn comparator_library(record: &BamRecord, header: &SamHeader) -> String {
    read_group(record, header)
        .and_then(|group| group.attributes.get("LB"))
        .unwrap_or(UNKNOWN_LIBRARY)
        .to_string()
}

/// `DuplicateScoringStrategy.computeDuplicateScore(record, strategy, true)`.
fn duplicate_score(record: &BamRecord, scoring: Scoring) -> i16 {
    const CAP: i32 = i16::MAX as i32 / 2;
    let mut score: i16 = 0;
    match scoring {
        Scoring::SumOfBaseQualities => {
            let sum: i32 = record
                .base_qualities
                .iter()
                .filter(|quality| **quality >= 15)
                .map(|quality| *quality as i32)
                .sum();
            score = score.wrapping_add(sum.min(CAP) as i16);
        }
        Scoring::TotalMappedReferenceLength => {
            if !flag(record, UNMAPPED) {
                score = (record.cigar.reference_length() as i32).min(CAP) as i16;
            }
            if flag(record, PAIRED) && !flag(record, MATE_UNMAPPED) {
                // A mate with no `MC` is an NPE in the reference; every record that reaches this
                // with a mapped mate was already refused by the mate-coordinate check unless it is
                // itself unmapped, and the corpus gives those the tag.
                let mate = mate_cigar(record).map_or(0, |cigar| cigar.reference_length() as i32);
                score = score.wrapping_add(mate.min(CAP) as i16);
            }
        }
        Scoring::Random => {
            let hashed =
                Murmur3::new(1).hash_unencoded_chars(&record.read_name) & 0b11_1111_1111_1111;
            score = score.wrapping_add(hashed as i16);
            score = score.wrapping_sub(i16::MIN / 4);
        }
    }
    if flag(record, VENDOR_FAILED) {
        score = score.wrapping_add(i16::MIN / 2);
    }
    score
}

/// `SAMUtils.getCanonicalRecordName`: `<RG>:<name>` where the record names a read group.
fn canonical_name(record: &BamRecord) -> String {
    match record.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(group)) => format!("{group}:{}", record.read_name),
        _ => record.read_name.clone(),
    }
}

/// `String.compareTo`'s sign, over UTF-16 code units.
fn java_compare(a: &str, b: &str) -> Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

const FF: i32 = 0;
const FR: i32 = 1;
const F: i32 = 2;
const RF: i32 = 3;
const RR: i32 = 4;
const R: i32 = 5;

fn paired_orientation(record: &BamRecord) -> i32 {
    if paired_and_both_mapped(record) {
        match (flag(record, REVERSE), flag(record, MATE_REVERSE)) {
            (true, true) => RR,
            (true, false) => RF,
            (false, true) => FR,
            (false, false) => FF,
        }
    } else if flag(record, REVERSE) {
        R
    } else {
        F
    }
}

fn collapsed_orientation_compare(left: i32, right: i32) -> i32 {
    if left == F || left == R {
        if left == F && (right == F || right == FR || right == FF) {
            return 0;
        }
        if left == R && (right == R || right == RF || right == RR) {
            return 0;
        }
    } else if right == F || right == R {
        return -collapsed_orientation_compare(right, left);
    }
    left - right
}

fn has_unmapped_end(record: &BamRecord) -> bool {
    flag(record, UNMAPPED) || (flag(record, PAIRED) && flag(record, MATE_UNMAPPED))
}

fn has_mapped_end(record: &BamRecord) -> bool {
    !flag(record, UNMAPPED) || (flag(record, PAIRED) && !flag(record, MATE_UNMAPPED))
}

struct Comparator<'a> {
    records: &'a [BamRecord],
    keys: &'a [Keys],
}

impl Comparator<'_> {
    /// `fileOrderCompare(a, b, collapseOrientation, considerEnds)`, as a sign.
    fn file_order(&self, a: usize, b: usize, collapse: bool, consider_ends: bool) -> i32 {
        let (left, right) = (&self.records[a], &self.records[b]);
        let (lk, rk) = (&self.keys[a], &self.keys[b]);
        let mut cmp = lk.library_id - rk.library_id;
        if cmp == 0 {
            let (l, r) = (left.reference_index, right.reference_index);
            cmp = if l == -1 {
                if r == -1 {
                    0
                } else {
                    1
                }
            } else if r == -1 {
                -1
            } else {
                l - r
            };
        }
        if cmp == 0 {
            cmp = lk.read_coordinate - rk.read_coordinate;
        }
        if cmp == 0 {
            let (l, r) = (paired_orientation(left), paired_orientation(right));
            cmp = if collapse {
                collapsed_orientation_compare(l, r)
            } else {
                l - r
            };
        }
        if paired_and_both_mapped(left) && paired_and_both_mapped(right) {
            if cmp == 0 {
                cmp = left.mate_reference_index - right.mate_reference_index;
            }
            if cmp == 0 {
                cmp = lk.mate_coordinate - rk.mate_coordinate;
            }
        }
        if cmp == 0 {
            cmp = i32::from(!has_mapped_end(left)) - i32::from(!has_mapped_end(right));
        }
        if cmp == 0 && consider_ends {
            if flag(left, PAIRED) == flag(right, PAIRED) {
                cmp = i32::from(has_unmapped_end(left)) - i32::from(has_unmapped_end(right));
            } else {
                cmp = if flag(left, PAIRED) { -1 } else { 1 };
            }
        }
        cmp
    }

    fn duplicate_set_compare(&self, a: usize, b: usize) -> i32 {
        self.file_order(a, b, true, false)
    }

    /// `SAMRecordDuplicateComparator.compare`.
    fn compare(&self, a: usize, b: usize) -> Ordering {
        let file_order = self.file_order(a, b, false, true);
        if file_order != 0 {
            return file_order.cmp(&0);
        }
        let (left, right) = (&self.records[a], &self.records[b]);
        // `DuplicateScoringStrategy.compare`: paired first, then the score, then the canonical
        // name.
        if flag(left, PAIRED) != flag(right, PAIRED) {
            return if flag(left, PAIRED) {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        let scored = self.keys[b].score as i32 - self.keys[a].score as i32;
        if scored != 0 {
            return scored.cmp(&0);
        }
        let canonical = java_compare(&self.keys[a].canonical_name, &self.keys[b].canonical_name);
        if canonical != Ordering::Equal {
            return canonical;
        }
        let named = java_compare(&left.read_name, &right.read_name);
        if named != Ordering::Equal {
            return named;
        }
        if flag(left, PAIRED) && flag(right, PAIRED) {
            let l = i32::from(!flag(left, FIRST_OF_PAIR));
            let r = i32::from(!flag(right, FIRST_OF_PAIR));
            return l.cmp(&r);
        }
        Ordering::Equal
    }
}

/// `DuplicateSet`: the records it accepted, its representative, and whether `getRecords` still
/// has to sort.
struct DuplicateSet {
    records: Vec<usize>,
    representative: usize,
    needs_sorting: bool,
}

impl DuplicateSet {
    fn new() -> Self {
        DuplicateSet {
            records: Vec::new(),
            representative: 0,
            needs_sorting: false,
        }
    }

    /// `DuplicateSet.add`: `0` when the record joined, the comparison otherwise.
    fn add(&mut self, record: usize, comparator: &Comparator) -> i32 {
        if !self.records.is_empty() {
            let cmp = comparator.duplicate_set_compare(self.representative, record);
            if cmp != 0 {
                return cmp;
            }
            if comparator.compare(self.representative, record) == Ordering::Greater {
                self.representative = record;
            }
        } else {
            self.representative = record;
        }
        self.records.push(record);
        self.needs_sorting = true;
        0
    }

    /// `DuplicateSet.getRecords()`: sort, then set the duplicate flag on every mapped primary
    /// record whose name is not the representative's and clear it on the first.
    fn records(&mut self, comparator: &Comparator, flags: &mut [u16]) -> Vec<usize> {
        if self.needs_sorting && !self.records.is_empty() {
            if self.records.len() > 1 {
                self.records.sort_by(|a, b| comparator.compare(*a, *b));
            }
            let name = &comparator.records[self.representative].read_name;
            for index in &self.records {
                let record = &comparator.records[*index];
                if !flag(record, UNMAPPED)
                    && !secondary_or_supplementary(record)
                    && &record.read_name != name
                {
                    flags[*index] |= DUPLICATE;
                }
            }
            flags[self.records[0]] &= !DUPLICATE;
        }
        self.needs_sorting = false;
        self.records.clone()
    }
}

/// `UmiUtil.ReadStrand`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Strand {
    Top,
    Bottom,
    Unknown,
}

/// `UmiUtil.getStrand`.
fn strand(record: &BamRecord, header: &SamHeader) -> Result<Strand, Thrown> {
    if flag(record, UNMAPPED) {
        return Ok(Strand::Unknown);
    }
    if !flag(record, PAIRED) {
        return Err(Thrown::new(
            "java.lang.IllegalStateException",
            "Inappropriate call if not paired read",
        ));
    }
    if flag(record, MATE_UNMAPPED) {
        return Ok(Strand::Unknown);
    }
    let first = flag(record, FIRST_OF_PAIR);
    if record.reference_index != record.mate_reference_index {
        return Ok(
            if first == (record.reference_index < record.mate_reference_index) {
                Strand::Top
            } else {
                Strand::Bottom
            },
        );
    }
    let read = if flag(record, REVERSE) {
        unclipped_end(record)
    } else {
        unclipped_start(record)
    };
    let mate = mate_five_prime(record, header)?;
    Ok(if first == (read <= mate) {
        Strand::Top
    } else {
        Strand::Bottom
    })
}

/// Java's `String.split("-")`: no match leaves the string whole, and trailing empty strings go.
fn java_split_hyphen(text: &str) -> Vec<&str> {
    if !text.contains('-') {
        return vec![text];
    }
    let mut parts: Vec<&str> = text.split('-').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `UmiUtil.getTopStrandNormalizedUmi`.
fn normalized_umi(
    record: &BamRecord,
    umi_tag: &str,
    duplex: bool,
    header: &SamHeader,
) -> Result<Option<String>, Thrown> {
    let Some(umi) = string_tag(record, umi_tag)? else {
        return Ok(None);
    };
    if !umi.bytes().all(|base| b"ATCGNatcgn-".contains(&base)) {
        return Err(Thrown::picard(
            "UMI found with illegal characters.  UMIs must match the regular expression \
             ^[ATCGNatcgn-]*$.",
        ));
    }
    if !duplex {
        return Ok(Some(umi.to_string()));
    }
    let split = java_split_hyphen(umi);
    if split.len() != 2 {
        return Err(Thrown::picard(format!(
            "Duplex UMIs must be of the form X-Y where X and Y are equal length UMIs, for example \
             AT-GA.  Found UMI, {umi}"
        )));
    }
    Ok(Some(match strand(record, header)? {
        Strand::Bottom => format!("{}-{}", split[1], split[0]),
        _ => umi.to_string(),
    }))
}

/// `StringUtil.isWithinHammingDistance`, which refuses strings of different lengths.
fn within_hamming_distance(a: &str, b: &str, max: i32) -> Result<bool, Thrown> {
    let (a, b): (Vec<u16>, Vec<u16>) = (a.encode_utf16().collect(), b.encode_utf16().collect());
    if a.len() != b.len() {
        return Err(Thrown::new(
            "java.lang.IllegalArgumentException",
            "Attempted to determine if two strings of different length were within a specified \
             edit distance.",
        ));
    }
    let mut distance = 0;
    for (x, y) in a.iter().zip(&b) {
        if x != y {
            distance += 1;
            if distance > max {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// `StringUtil.hammingDistance`.
fn hamming_distance(a: &str, b: &str) -> Result<i64, Thrown> {
    let (a, b): (Vec<u16>, Vec<u16>) = (a.encode_utf16().collect(), b.encode_utf16().collect());
    if a.len() != b.len() {
        return Err(Thrown::new(
            "java.lang.IllegalArgumentException",
            format!(
                "Attempted to determine Hamming distance of strings with differing lengths. The \
                 first string has length {} and the second string has length {}.",
                a.len(),
                b.len()
            ),
        ));
    }
    Ok(a.iter().zip(&b).filter(|(x, y)| x != y).count() as i64)
}

/// `GraphUtils.Graph.findRepNode`, with its path compression.
fn find_rep(grouping: &mut [usize], mut node: usize) -> usize {
    let mut representative = node;
    while representative != grouping[representative] {
        representative = grouping[representative];
    }
    while node != representative {
        let next = grouping[node];
        grouping[node] = representative;
        node = next;
    }
    representative
}

/// `GraphUtils.Graph.cluster` over the UMI graph `UmiGraph` builds: node `i` per UMI, an edge for
/// every pair within the distance, neighbours in the order the double loop added them.
fn cluster(umis: &[String], max: i32) -> Result<Vec<usize>, Thrown> {
    let n = umis.len();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            if within_hamming_distance(&umis[i], &umis[j], max)? {
                if !neighbors[i].contains(&j) {
                    neighbors[i].push(j);
                }
                if !neighbors[j].contains(&i) {
                    neighbors[j].push(i);
                }
            }
        }
    }
    let mut grouping: Vec<usize> = (0..n).collect();
    for (i, list) in neighbors.iter().enumerate() {
        for j in list {
            // `joinNodes(cluster, j, i)`: the root of j is pointed at the root of i.
            let a = find_rep(&mut grouping, *j);
            let b = find_rep(&mut grouping, i);
            if a != b {
                grouping[a] = b;
            }
        }
    }
    Ok((0..n).map(|i| find_rep(&mut grouping, i)).collect())
}

/// A `java.util.HashMap<Integer, V>` kept for its iteration order. `Integer.hashCode` is the value
/// and `HashMap.hash` folds the high half down, which for these small non-negative keys is the
/// value itself.
struct IntHashMap<V> {
    table: Vec<Vec<(usize, V)>>,
    size: usize,
}

impl<V> IntHashMap<V> {
    fn new() -> Self {
        IntHashMap {
            table: Vec::new(),
            size: 0,
        }
    }
    fn spread(key: usize) -> usize {
        let h = key as u32;
        (h ^ (h >> 16)) as usize
    }
    fn get_mut(&mut self, key: usize) -> Option<&mut V> {
        if self.table.is_empty() {
            return None;
        }
        let index = Self::spread(key) & (self.table.len() - 1);
        self.table[index]
            .iter_mut()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v)
    }
    fn put(&mut self, key: usize, value: V) {
        if self.table.is_empty() {
            self.table = (0..16).map(|_| Vec::new()).collect();
        }
        let index = Self::spread(key) & (self.table.len() - 1);
        self.table[index].push((key, value));
        self.size += 1;
        if self.size > self.table.len() * 3 / 4 {
            let old = self.table.len();
            let mut grown: Vec<Vec<(usize, V)>> = (0..old * 2).map(|_| Vec::new()).collect();
            for (j, bucket) in std::mem::take(&mut self.table).into_iter().enumerate() {
                for (k, v) in bucket {
                    let high = Self::spread(k) & old != 0;
                    grown[if high { j + old } else { j }].push((k, v));
                }
            }
            self.table = grown;
        }
    }
    fn into_values(self) -> impl Iterator<Item = V> {
        self.table
            .into_iter()
            .flat_map(|bucket| bucket.into_iter().map(|(_, v)| v))
    }
}

/// `SimpleMarkDuplicatesWithMateCigar.doWork` driven by the UMI-aware iterator.
pub fn run(
    header: &SamHeader,
    mut records: Vec<BamRecord>,
    args: &UmiArgs,
) -> Result<UmiRun, Thrown> {
    // The comparator's library ids: "Unknown Library" and every header library, numbered from one
    // in sorted order.
    let mut names: BTreeSet<String> = BTreeSet::new();
    names.insert(UNKNOWN_LIBRARY.to_string());
    for group in &header.read_groups {
        if let Some(library) = group.attributes.get("LB") {
            names.insert(library.to_string());
        }
    }
    let ids: BTreeMap<String, i32> = names
        .into_iter()
        .enumerate()
        .map(|(position, name)| (name, position as i32 + 1))
        .collect();

    // The sort reaches the mate coordinate of `a[1]` first and of `a[0]` second; a missing `MC`
    // is refused as htsjdk computes it.
    if records.len() > 1 {
        let examination = [1usize, 0].into_iter().chain(2..records.len());
        for index in examination {
            let record = &records[index];
            if paired_and_both_mapped(record) && mate_cigar(record).is_none() {
                mate_five_prime(record, header)?;
            }
        }
    }
    let keys: Vec<Keys> = records
        .iter()
        .map(|record| Keys {
            library_id: ids[&comparator_library(record, header)],
            read_coordinate: if flag(record, REVERSE) {
                unclipped_end(record)
            } else {
                unclipped_start(record)
            },
            mate_coordinate: if paired_and_both_mapped(record) {
                mate_five_prime(record, header).unwrap_or(-1)
            } else {
                -1
            },
            score: duplicate_score(record, args.scoring),
            canonical_name: canonical_name(record),
        })
        .collect();

    let mut flags: Vec<u16> = records.iter().map(|record| record.flags).collect();
    // The tag values the run changes, applied to the records at the end.
    let mut umi_tag_value: Vec<Option<Option<String>>> = vec![None; records.len()];
    let mut molecular_identifier: Vec<Option<String>> = vec![None; records.len()];

    let snapshot = records.clone();
    let comparator = Comparator {
        records: &snapshot,
        keys: &keys,
    };

    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_by(|a, b| comparator.compare(*a, *b));

    // `DuplicateSetIterator`: cut a set wherever the next record is not comparable to the
    // representative, or the representative is unmapped or secondary.
    let mut outer_sets: Vec<DuplicateSet> = Vec::new();
    let mut current = DuplicateSet::new();
    for index in order {
        if current.records.is_empty() {
            current.add(index, &comparator);
            continue;
        }
        let representative = &snapshot[current.representative];
        if flag(representative, UNMAPPED) || secondary_or_supplementary(representative) {
            outer_sets.push(std::mem::replace(&mut current, DuplicateSet::new()));
            current.add(index, &comparator);
            continue;
        }
        let cmp = current.add(index, &comparator);
        if cmp > 0 {
            return Err(Thrown::new(
                "htsjdk.samtools.SAMException",
                "The input records were not sorted in duplicate order:",
            ));
        } else if cmp < 0 {
            outer_sets.push(std::mem::replace(&mut current, DuplicateSet::new()));
            current.add(index, &comparator);
        }
    }
    if !current.records.is_empty() {
        outer_sets.push(current);
    }

    let mut metrics: JavaHashMap<Accumulator> = JavaHashMap::new();
    let mut seen_first_read = false;
    let mut written: Vec<usize> = Vec::new();
    let umi_tag = tag_of(&args.umi_tag)?;

    // The current value of the UMI tag, as the run has left it.
    let current_umi = |index: usize, overrides: &[Option<Option<String>>]| -> Option<String> {
        match &overrides[index] {
            Some(value) => value.clone(),
            None => match snapshot[index].tags.get(umi_tag) {
                Some(TagValue::Str(value)) => Some(value.clone()),
                _ => None,
            },
        }
    };

    for mut outer in outer_sets {
        // `UmiGraph`'s constructor.
        let members = outer.records(&comparator, &mut flags);
        for index in &members {
            if string_tag(&snapshot[*index], &args.umi_tag)?.is_none() {
                if !args.allow_missing_umis {
                    return Err(Thrown::picard(format!(
                        "Read {} does not contain a UMI with the {} attribute.",
                        snapshot[*index].read_name, args.umi_tag
                    )));
                }
                umi_tag_value[*index] = Some(Some(String::new()));
            }
        }
        let normalize = |index: usize,
                         overrides: &[Option<Option<String>>]|
         -> Result<Option<String>, Thrown> {
            let mut record = snapshot[index].clone();
            match current_umi(index, overrides) {
                Some(value) => record.tags.insert(umi_tag, TagValue::Str(value)),
                None => record.tags.remove(umi_tag),
            }
            normalized_umi(&record, &args.umi_tag, args.duplex_umi, header)
        };
        let mut counts: JavaHashMap<i64> = JavaHashMap::new();
        let mut normalized: Vec<String> = Vec::with_capacity(members.len());
        for index in &members {
            let umi = normalize(*index, &umi_tag_value)?.unwrap_or_default();
            let count = counts.get(&umi).copied().unwrap_or(0);
            counts.put(&umi, count + 1);
            normalized.push(umi);
        }
        let umis: Vec<String> = counts.iter().map(|(k, _)| k.to_string()).collect();
        let count_of = |umi: &str| counts.get(umi).copied().unwrap_or(0);

        // `set.getRepresentative().getReadGroup().getLibrary()`.
        let representative = &snapshot[outer.representative];
        let Some(group) = read_group(representative, header) else {
            return Err(Thrown::new(
                "java.lang.NullPointerException",
                "Cannot invoke \"htsjdk.samtools.SAMReadGroupRecord.getLibrary()\" because the \
                 return value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null",
            ));
        };
        let library: Option<String> = group.attributes.get("LB").map(str::to_string);
        let library_key = library.clone().unwrap_or_default();
        if !metrics.contains_key(&library_key) {
            metrics.put(
                &library_key,
                Accumulator {
                    row: UmiMetricsRow {
                        library: library.clone(),
                        ..UmiMetricsRow::default()
                    },
                    ..Accumulator::default()
                },
            );
        }

        // `joinUmisIntoDuplicateSets`.
        let roots = cluster(&umis, args.max_edit_distance_to_join)?;
        let set_of = |umi: &str| -> usize {
            let position = umis.iter().position(|u| u == umi).unwrap_or(0);
            roots[position]
        };
        let mut lists: IntHashMap<Vec<usize>> = IntHashMap::new();
        for (position, index) in members.iter().enumerate() {
            let id = set_of(&normalized[position]);
            match lists.get_mut(id) {
                Some(list) => list.push(*index),
                None => lists.put(id, vec![*index]),
            }
        }
        let mut inferred: Vec<Option<String>> = vec![None; snapshot.len()];
        let mut sub_sets: Vec<DuplicateSet> = Vec::new();
        for list in lists.into_values() {
            let mut set = DuplicateSet::new();
            for index in &list {
                set.add(*index, &comparator);
            }
            let mut max_count = 0i64;
            let mut assigned: Option<String> = None;
            let mut fewest_n: Option<String> = None;
            let mut n_count = 0usize;
            for index in &list {
                let umi = normalize(*index, &umi_tag_value)?.unwrap_or_default();
                if umi.contains('N') {
                    let count = umi.matches('N').count();
                    if n_count == 0 || count < n_count {
                        n_count = count;
                        fewest_n = Some(umi);
                    }
                } else if count_of(&umi) > max_count {
                    max_count = count_of(&umi);
                    assigned = Some(umi);
                }
            }
            if assigned.is_none() {
                assigned = fewest_n;
            }
            for index in &list {
                let tag_now = current_umi(*index, &umi_tag_value);
                if args.allow_missing_umis && tag_now.as_deref() == Some("") {
                    umi_tag_value[*index] = Some(None);
                } else if let Some(tag) = &args.molecular_identifier_tag {
                    let record = &snapshot[*index];
                    let mut text = format!(
                        "{}:{}/{}",
                        contig(record, header).unwrap_or("null"),
                        if flag(record, REVERSE) {
                            record.alignment_start
                        } else {
                            record.mate_alignment_start
                        },
                        assigned.as_deref().unwrap_or("null")
                    );
                    if args.duplex_umi {
                        match strand(record, header)? {
                            Strand::Top => text.push_str("/A"),
                            Strand::Bottom => text.push_str("/B"),
                            Strand::Unknown => {}
                        }
                    }
                    let _ = tag_of(tag)?;
                    molecular_identifier[*index] = Some(text);
                }
                inferred[*index] = assigned.clone();
            }
            sub_sets.push(set);
        }

        // `UmiAwareDuplicateSetIterator.process`: the statistics, per sub-set record.
        let accumulator = metrics_entry(&mut metrics, &library_key);
        let mut emitted: Vec<Vec<usize>> = Vec::new();
        for set in &mut sub_sets {
            let members = set.records(&comparator, &mut flags);
            for index in &members {
                let Some(current) = normalize(*index, &umi_tag_value)? else {
                    continue;
                };
                if current.contains('N') {
                    accumulator.observed_with_n += 1;
                    continue;
                }
                let length = current.encode_utf16().count() - current.matches('-').count();
                if !seen_first_read {
                    accumulator.row.mean_umi_length = length as f64;
                    seen_first_read = true;
                } else if accumulator.row.mean_umi_length != length as f64 {
                    return Err(Thrown::picard("UMIs of differing lengths were found."));
                }
                let inferred_umi = inferred[*index].clone().unwrap_or_default();
                accumulator.row.observed_base_errors += hamming_distance(&current, &inferred_umi)?;
                accumulator.add_observation(&current, &inferred_umi);
            }
            emitted.push(members);
        }
        accumulator.row.duplicate_sets_with_umi += sub_sets.len() as i64;
        accumulator.row.duplicate_sets_ignoring_umi += 1;

        // `doWork`'s loop over the sub-sets: written unless removed.
        for members in emitted {
            for index in members {
                if !args.remove_duplicates || flags[index] & DUPLICATE == 0 {
                    written.push(index);
                }
            }
        }
    }

    // Apply what the run changed and hand the records over in the order they were written.
    for (index, record) in records.iter_mut().enumerate() {
        record.flags = flags[index];
        if let Some(value) = &umi_tag_value[index] {
            match value {
                Some(text) => record.tags.insert(umi_tag, TagValue::Str(text.clone())),
                None => record.tags.remove(umi_tag),
            }
        }
        if let (Some(tag), Some(text)) =
            (&args.molecular_identifier_tag, &molecular_identifier[index])
        {
            record
                .tags
                .insert(tag_of(tag)?, TagValue::Str(text.clone()));
        }
        if args.program_record_id.is_some() {
            record.tags.remove(Tag::new(b"PG"));
        }
    }
    let out: Vec<BamRecord> = written
        .iter()
        .map(|index| records[*index].clone())
        .collect();

    let rows = metrics
        .iter()
        .map(|(_, accumulator)| accumulator.clone().finish())
        .collect();
    Ok(UmiRun {
        written: out,
        metrics: rows,
    })
}

fn metrics_entry<'m>(metrics: &'m mut JavaHashMap<Accumulator>, key: &str) -> &'m mut Accumulator {
    metrics
        .get_mut(key)
        .expect("the library's metrics were created")
}
