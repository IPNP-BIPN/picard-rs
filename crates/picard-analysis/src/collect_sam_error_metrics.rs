//! `CollectSamErrorMetrics`: how often a base disagrees with the reference where it should not.
//!
//! The tool is not a mismatch counter. What makes it a quality estimate is everything it refuses
//! to count: a base at a site the sample is known to be polymorphic at, a base below a quality, a
//! read below a mapping quality, and the second observation of a pair that overlaps itself. And
//! the rate it reports is Bayesian rather than a ratio: the prior is a pseudo-count in phred
//! space, so a file with no errors at all reports a finite quality rather than an infinite one.
//!
//! One run writes one file per metric, named `<basename>.<suffix>`, and the suffix is the
//! calculator's own plus `_by_` plus the stratifier's. A stratifier splits the rows and nothing
//! else: the same bases, counted per bin.
//!
//! Ported from `picard.sam.SamErrorMetric.CollectSamErrorMetrics`,
//! `picard.sam.SamErrorMetric.ErrorMetric`, `picard.sam.SamErrorMetric.BaseErrorMetric`,
//! `picard.sam.SamErrorMetric.IndelErrorMetric`, `picard.sam.SamErrorMetric.OverlappingErrorMetric`,
//! `picard.sam.SamErrorMetric.SimpleErrorCalculator`, `picard.sam.SamErrorMetric.IndelErrorCalculator`,
//! `picard.sam.SamErrorMetric.OverlappingReadsErrorCalculator`,
//! `picard.sam.SamErrorMetric.BaseErrorAggregation` and
//! `picard.sam.SamErrorMetric.ReadBaseStratification` in Picard 3.4.0.

use std::collections::{BTreeMap, HashSet};

use crate::theoretical_sensitivity::JavaRandom;

/// The error probability a phred-scaled prior stands for: `PRIOR_Q` of 30 is one error in a
/// thousand.
pub fn prior_error(prior_q: i32) -> f64 {
    10f64.powf(-f64::from(prior_q) / 10.0)
}

/// `QualityUtil.getPhredScoreFromErrorProbability`, rounded the way Java rounds it.
pub fn phred_from_error_probability(probability: f64) -> i32 {
    // htsjdk's own function, not a local copy of its formula: `Math.log10` is correctly rounded and
    // `Math.round` is half up, and `f64::log10`/`f64::round` are neither.
    htsjdk_bam::quality_util::phred_score_from_error_probability(probability)
}

/// The quality of a count of errors out of a count of bases.
///
/// The prior is a pseudo-count: one prior's worth of error in the numerator and one whole base in
/// the denominator, so no errors at all still gives a finite number, and moving the prior moves
/// it.
pub fn q_score(errors: u64, total_bases: u64, prior_error: f64) -> i32 {
    phred_from_error_probability((errors as f64 + prior_error) / (total_bases as f64 + 1.0))
}

/// What a base is doing at the locus it was shown at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    Match,
    Insertion,
    Deletion,
}

/// One read, as much of it as the tool reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    /// One-based, on the single (possibly concatenated) contig the reads are placed on.
    pub start: i32,
    pub bases: Vec<u8>,
    /// Phred, not ASCII.
    pub qualities: Vec<u8>,
    pub flags: u16,
    /// One-based, on the same coordinate as `start`; zero when the mate has no position.
    pub mate_start: i32,
    pub cigar: Vec<(usize, char)>,
    pub mapping_quality: u8,
    /// The `RG` tag, which the read-group stratifier reads.
    pub read_group: String,
    /// The template length, which the insert-length stratifier reads.
    pub insert_size: i32,
}

impl Read {
    pub fn is_paired(&self) -> bool {
        self.flags & 0x1 != 0
    }
    pub fn is_first_of_pair(&self) -> bool {
        self.flags & 0x40 != 0
    }
    pub fn is_second_of_pair(&self) -> bool {
        self.flags & 0x80 != 0
    }
    pub fn is_unmapped(&self) -> bool {
        self.flags & 0x4 != 0
    }
    pub fn is_secondary(&self) -> bool {
        self.flags & 0x100 != 0
    }
    pub fn is_negative_strand(&self) -> bool {
        self.flags & 0x10 != 0
    }
}

/// One base of one read, at one locus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Observation {
    pub read: usize,
    /// The offset into the read's bases, which is what the cycle is counted from. A deletion is
    /// shown at the offset of the base BEFORE it, which is minus one for a read that starts with
    /// one.
    pub offset: i64,
    pub alignment: AlignmentType,
}

/// A locus and everything read over it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub position: i32,
    pub records: Vec<Observation>,
}

/// `SequenceUtil.isNoCall`.
pub fn is_no_call(base: u8) -> bool {
    matches!(base, b'N' | b'n' | b'.')
}

/// `SequenceUtil.basesEqual`, which is case-insensitive.
pub fn bases_equal(left: u8, right: u8) -> bool {
    left.eq_ignore_ascii_case(&right)
}

/// What a run was asked for.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub min_mapping_q: u8,
    pub min_base_q: u8,
    pub prior_q: i32,
    /// Zero is unlimited.
    pub max_loci: u64,
    /// One-based positions the sample is known to be polymorphic at.
    pub known_sites: Vec<i32>,
    /// The chance a locus is looked at at all; a locus that is not is neither counted nor
    /// skipped.
    pub probability: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            min_mapping_q: 20,
            min_base_q: 20,
            prior_q: 30,
            max_loci: 0,
            known_sites: Vec::new(),
            probability: 1.0,
        }
    }
}

/// What one locus holds before it is flattened into a [`Locus`]: the matched bases, then the
/// deletions, then the insertions, which is the order `addLocusBases` shows them in.
#[derive(Default)]
struct Pile {
    matches: Vec<Observation>,
    deletions: Vec<Observation>,
    insertions: Vec<Observation>,
}

fn consumes_read(operator: char) -> bool {
    matches!(operator, 'M' | 'I' | 'S' | '=' | 'X')
}

fn consumes_reference(operator: char) -> bool {
    matches!(operator, 'M' | 'D' | 'N' | '=' | 'X')
}

/// The pileup the tool walks: every locus a read covers, indels included, with the two thresholds
/// already applied.
///
/// Ported from htsjdk 4.2.0's `SamLocusIterator`, as `CollectSamErrorMetrics` configures it:
///
/// * secondary, supplementary and duplicate reads are filtered out, and so are unmapped reads and
///   reads below the mapping-quality cutoff;
/// * a matched base below the base-quality cutoff is not accumulated;
/// * an insertion is shown at the base BEFORE it, and is dropped when its first base is below the
///   cutoff; **the read offset only moves past an insertion that is kept**, so a dropped one
///   leaves every later offset of the read short by its length (an htsjdk quirk, ported);
/// * a deletion is shown at every reference position it spans, at the offset of the base before
///   it, and the base quality does not apply to it;
/// * a locus holding nothing at all is not emitted.
///
/// The thresholds drop observations rather than loci, one by read and one by base, which is why a
/// mismatch below `--MIN_BASE_Q` lowers the denominator instead of raising the error count.
pub fn pileup(reads: &[Read], options: &Options) -> Vec<Locus> {
    let mut loci: BTreeMap<i32, Pile> = BTreeMap::new();
    let cutoff = options.min_base_q;
    let passes = |read: &Read, offset: usize| {
        cutoff == 0
            || read.qualities.is_empty()
            || read.qualities.get(offset).is_none_or(|&q| q >= cutoff)
    };
    for (index, read) in reads.iter().enumerate() {
        if read.flags & (0x100 | 0x800 | 0x400) != 0
            || read.is_unmapped()
            || read.mapping_quality < options.min_mapping_q
        {
            continue;
        }
        // accumulateSamRecord: one observation per aligned base that meets the cutoff.
        let mut position = read.start;
        let mut offset = 0usize;
        for &(length, operator) in &read.cigar {
            if matches!(operator, 'M' | '=' | 'X') {
                for step in 0..length {
                    if passes(read, offset + step) {
                        loci.entry(position + step as i32)
                            .or_default()
                            .matches
                            .push(Observation {
                                read: index,
                                offset: (offset + step) as i64,
                                alignment: AlignmentType::Match,
                            });
                    }
                }
            }
            if consumes_read(operator) {
                offset += length;
            }
            if consumes_reference(operator) {
                position += length as i32;
            }
        }
        // accumulateIndels.
        let mut read_base: i64 = 0;
        let mut position = read.start;
        for &(length, operator) in &read.cigar {
            match operator {
                'I' => {
                    if passes(read, read_base.max(0) as usize) {
                        loci.entry(position - 1)
                            .or_default()
                            .insertions
                            .push(Observation {
                                read: index,
                                offset: read_base,
                                alignment: AlignmentType::Insertion,
                            });
                        read_base += length as i64;
                    }
                }
                'D' => {
                    for step in 0..length {
                        loci.entry(position + step as i32)
                            .or_default()
                            .deletions
                            .push(Observation {
                                read: index,
                                offset: read_base - 1,
                                alignment: AlignmentType::Deletion,
                            });
                    }
                    position += length as i32;
                }
                other => {
                    if consumes_read(other) {
                        read_base += length as i64;
                    }
                    if consumes_reference(other) {
                        position += length as i32;
                    }
                }
            }
        }
    }
    loci.into_iter()
        .map(|(position, pile)| {
            let mut records = pile.matches;
            records.extend(pile.deletions);
            records.extend(pile.insertions);
            Locus { position, records }
        })
        .collect()
}

/// The loci a run actually counts, in the order the tool's loop takes them.
///
/// Each locus first draws from `Random(42)` and is not looked at when the draw is above
/// `--PROBABILITY`; a known site is then skipped; the rest are counted, and `--MAX_LOCI` stops
/// the run once that many have been. So the cap counts what is left after the VCF and the draw
/// have taken their loci out.
pub fn processed_loci(loci: Vec<Locus>, options: &Options) -> Vec<Locus> {
    let mut random = JavaRandom::new(42);
    let mut kept = Vec::new();
    for locus in loci {
        if random.next_double() > options.probability {
            continue;
        }
        if options.known_sites.contains(&locus.position) {
            continue;
        }
        kept.push(locus);
        if options.max_loci != 0 && kept.len() as u64 >= options.max_loci {
            break;
        }
    }
    kept
}

/// The three calculators, by the suffix each one names its file with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Calculator {
    Error,
    OverlappingError,
    IndelError,
}

impl Calculator {
    pub fn suffix(&self) -> &'static str {
        match self {
            Calculator::Error => "error",
            Calculator::OverlappingError => "overlapping_error",
            Calculator::IndelError => "indel_error",
        }
    }

    /// `ErrorType.valueOf`.
    pub fn parse(name: &str) -> Option<Calculator> {
        match name {
            "ERROR" => Some(Calculator::Error),
            "OVERLAPPING_ERROR" => Some(Calculator::OverlappingError),
            "INDEL_ERROR" => Some(Calculator::IndelError),
            _ => None,
        }
    }
}

/// The stratifiers whose binning is ported.
///
/// Every stratifier's file suffix is in [`stratifier_suffix`]; these are the ones that also put a
/// base in a bin here, and [`Stratifier::parse`] names them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stratifier {
    All,
    BaseQuality,
    Cycle,
    GcContent,
    ReadDirection,
    ReadOrdinality,
    MappingQuality,
    ReadGroup,
    ReadBase,
    ReferenceBase,
    InsertLength,
    SoftClips,
}

impl Stratifier {
    /// `ReadBaseStratification.Stratifier.valueOf`, for the stratifiers this port bins.
    pub fn parse(name: &str) -> Option<Stratifier> {
        Some(match name {
            "ALL" => Stratifier::All,
            "BASE_QUALITY" => Stratifier::BaseQuality,
            "CYCLE" => Stratifier::Cycle,
            "GC_CONTENT" => Stratifier::GcContent,
            "READ_DIRECTION" => Stratifier::ReadDirection,
            "READ_ORDINALITY" => Stratifier::ReadOrdinality,
            "MAPPING_QUALITY" => Stratifier::MappingQuality,
            "READ_GROUP" => Stratifier::ReadGroup,
            "READ_BASE" => Stratifier::ReadBase,
            "REFERENCE_BASE" => Stratifier::ReferenceBase,
            "INSERT_LENGTH" => Stratifier::InsertLength,
            "SOFT_CLIPS" => Stratifier::SoftClips,
            _ => return None,
        })
    }
}

/// `ReadBaseStratification.Stratifier`, name to file suffix.
///
/// The suffix is not the name lower-cased: `GC_CONTENT` writes `gc`, and the two homopolymer
/// stratifiers name the reference base they are followed by.
pub fn stratifier_suffix(name: &str) -> Option<&'static str> {
    Some(match name {
        "ALL" => "all",
        "GC_CONTENT" => "gc",
        "READ_ORDINALITY" => "read_ordinality",
        "READ_BASE" => "read_base",
        "READ_DIRECTION" => "read_direction",
        "PAIR_ORIENTATION" => "pair_orientation",
        "PAIR_PROPERNESS" => "pair_proper",
        "REFERENCE_BASE" => "ref_base",
        "PRE_DINUC" => "pre_dinuc",
        "POST_DINUC" => "post_dinuc",
        "HOMOPOLYMER_LENGTH" => "homopolymer_length",
        "HOMOPOLYMER" => "homopolymer_and_following_ref_base",
        "BINNED_HOMOPOLYMER" => "binned_length_homopolymer_and_following_ref_base",
        "FLOWCELL_TILE" => "tile",
        "FLOWCELL_X" => "x",
        "FLOWCELL_Y" => "y",
        "READ_GROUP" => "read_group",
        "CYCLE" => "cycle",
        "BINNED_CYCLE" => "binned_cycle",
        "SOFT_CLIPS" => "softclipped_bases",
        "INSERT_LENGTH" => "insert_length",
        "BASE_QUALITY" => "base_quality",
        "MAPPING_QUALITY" => "mapping_quality",
        "MISMATCHES_IN_READ" => "mismatches_in_read",
        "ONE_BASE_PADDED_CONTEXT" => "one_base_padded_context",
        "TWO_BASE_PADDED_CONTEXT" => "two_base_padded_context",
        "CONSENSUS" => "consensus",
        "NS_IN_READ" => "ns_in_read",
        "INSERTIONS_IN_READ" => "cigar_elements_I_in_read",
        "DELETIONS_IN_READ" => "cigar_elements_D_in_read",
        "INDELS_IN_READ" => "indels_in_read",
        "INDEL_LENGTH" => "indel_length",
        _ => return None,
    })
}

/// The twenty-seven directives a run collects when it is not told otherwise.
pub const DEFAULT_ERROR_METRICS: [&str; 27] = [
    "ERROR",
    "ERROR:BASE_QUALITY",
    "ERROR:INSERT_LENGTH",
    "ERROR:GC_CONTENT",
    "ERROR:READ_DIRECTION",
    "ERROR:PAIR_ORIENTATION",
    "ERROR:HOMOPOLYMER",
    "ERROR:BINNED_HOMOPOLYMER",
    "ERROR:CYCLE",
    "ERROR:READ_ORDINALITY",
    "ERROR:READ_ORDINALITY:CYCLE",
    "ERROR:READ_ORDINALITY:HOMOPOLYMER",
    "ERROR:READ_ORDINALITY:GC_CONTENT",
    "ERROR:READ_ORDINALITY:PRE_DINUC",
    "ERROR:MAPPING_QUALITY",
    "ERROR:READ_GROUP",
    "ERROR:MISMATCHES_IN_READ",
    "ERROR:ONE_BASE_PADDED_CONTEXT",
    "OVERLAPPING_ERROR",
    "OVERLAPPING_ERROR:BASE_QUALITY",
    "OVERLAPPING_ERROR:INSERT_LENGTH",
    "OVERLAPPING_ERROR:READ_ORDINALITY",
    "OVERLAPPING_ERROR:READ_ORDINALITY:CYCLE",
    "OVERLAPPING_ERROR:READ_ORDINALITY:HOMOPOLYMER",
    "OVERLAPPING_ERROR:READ_ORDINALITY:GC_CONTENT",
    "OVERLAPPING_ERROR:READ_ORDINALITY:PRE_DINUC",
    "INDEL_ERROR",
];

/// The file suffix one directive writes.
///
/// Several stratifiers are folded into one from the left, and each fold joins with `_and_`, so
/// `ERROR:READ_ORDINALITY:CYCLE` writes `error_by_read_ordinality_and_cycle`.
pub fn aggregation_suffix(directive: &str) -> Option<String> {
    let mut terms = directive.split(':').map(str::trim);
    let calculator = Calculator::parse(terms.next()?)?;
    let stratifiers: Option<Vec<&str>> = terms.map(stratifier_suffix).collect();
    let stratifiers = stratifiers?;
    let joined = if stratifiers.is_empty() {
        "all".to_string()
    } else {
        stratifiers.join("_and_")
    };
    Some(format!("{}_by_{}", calculator.suffix(), joined))
}

/// What a list of directives is refused for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// Two directives that would write the same file.
    DuplicatedSuffix { suffix: String, class: String },
}

impl Refusal {
    pub fn message(&self) -> String {
        match self {
            Refusal::DuplicatedSuffix { suffix, class } => {
                format!("Duplicated suffix ({suffix}) found in aggregator {class}.")
            }
        }
    }
}

/// The suffixes a list of directives writes, or the duplicate that refuses it.
///
/// `--ERROR_METRICS` is a collection, and Picard's parser APPENDS to a collection rather than
/// replacing it. So naming a metric the default list already carries, without clearing the list
/// first, asks for the same file twice and the run is refused before a single locus is read.
pub fn suffixes(directives: &[String]) -> Result<Vec<String>, Refusal> {
    let mut seen: Vec<String> = Vec::new();
    for directive in directives {
        let Some(suffix) = aggregation_suffix(directive) else {
            continue;
        };
        if seen.contains(&suffix) {
            return Err(Refusal::DuplicatedSuffix {
                suffix,
                class: "class picard.sam.SamErrorMetric.BaseErrorAggregation".to_string(),
            });
        }
        seen.push(suffix);
    }
    Ok(seen)
}

/// How the two halves of a pair order two strata, which is how the tool's sorted set does.
#[derive(Debug, Clone, PartialEq)]
enum Order {
    /// An `Integer`, a `Byte` or a `Double`.
    Number(f64),
    /// An enum, by declaration order.
    Ordinal(u8),
    /// A `String` or a `Character`.
    Text(String),
}

/// One stratifier's answer for one base: what it sorts by, and what `toString` prints.
#[derive(Debug, Clone, PartialEq)]
pub struct Part {
    order: Order,
    text: String,
}

impl Part {
    fn number(value: i64) -> Part {
        Part {
            order: Order::Number(value as f64),
            text: value.to_string(),
        }
    }
    fn text(value: String) -> Part {
        Part {
            order: Order::Text(value.clone()),
            text: value,
        }
    }
}

fn compare_parts(left: &Part, right: &Part) -> std::cmp::Ordering {
    match (&left.order, &right.order) {
        (Order::Number(a), Order::Number(b)) => {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        }
        (Order::Ordinal(a), Order::Ordinal(b)) => a.cmp(b),
        (Order::Text(a), Order::Text(b)) => a.cmp(b),
        _ => std::cmp::Ordering::Equal,
    }
}

/// A stratum: one part per stratifier, which a `Pair` of pairs compares left to right and prints
/// joined by a comma.
#[derive(Debug, Clone, PartialEq)]
struct Key(Vec<Part>);

impl Key {
    fn text(&self) -> String {
        self.0
            .iter()
            .map(|part| part.text.as_str())
            .collect::<Vec<_>>()
            .join(",")
    }
}

impl Eq for Key {}

impl Ord for Key {
    fn cmp(&self, other: &Key) -> std::cmp::Ordering {
        for (left, right) in self.0.iter().zip(&other.0) {
            let order = compare_parts(left, right);
            if order != std::cmp::Ordering::Equal {
                return order;
            }
        }
        std::cmp::Ordering::Equal
    }
}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Key) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// `SequenceUtil.complement`, which leaves anything but a plain base as it is.
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        other => other,
    }
}

/// `ReadBaseStratification.stratifySequenceBase`: complemented on the reverse strand, upper case.
fn sequence_base(base: u8, reverse: bool) -> char {
    let base = if reverse { complement(base) } else { base };
    char::from(base.to_ascii_uppercase())
}

/// The bin one observation falls in, or nothing, which drops it.
pub fn stratify(
    stratifier: Stratifier,
    reads: &[Read],
    observation: &Observation,
    reference_base: u8,
) -> Option<Part> {
    let read = &reads[observation.read];
    let offset = observation.offset;
    let at = |index: i64| -> Option<u8> {
        usize::try_from(index)
            .ok()
            .and_then(|i| read.bases.get(i).copied())
    };
    Some(match stratifier {
        Stratifier::All => Part::text("all".to_string()),
        Stratifier::BaseQuality => Part::number(i64::from(
            *read.qualities.get(usize::try_from(offset).ok()?)?,
        )),
        Stratifier::Cycle => Part::number(cycle(read, offset) as i64),
        Stratifier::GcContent => {
            let value = gc_content(&read.bases);
            Part {
                order: Order::Number(value),
                text: format_double(value),
            }
        }
        Stratifier::ReadDirection => {
            if read.is_negative_strand() {
                Part {
                    order: Order::Ordinal(1),
                    text: "-".to_string(),
                }
            } else {
                Part {
                    order: Order::Ordinal(0),
                    text: "+".to_string(),
                }
            }
        }
        Stratifier::ReadOrdinality => {
            if !read.is_paired() {
                return None;
            }
            if read.is_first_of_pair() {
                Part {
                    order: Order::Ordinal(0),
                    text: "FIRST".to_string(),
                }
            } else {
                Part {
                    order: Order::Ordinal(1),
                    text: "SECOND".to_string(),
                }
            }
        }
        Stratifier::MappingQuality => Part::number(i64::from(read.mapping_quality)),
        Stratifier::ReadGroup => Part::text(read.read_group.clone()),
        Stratifier::ReadBase => {
            let base = at(offset)?;
            Part::text(sequence_base(base, read.is_negative_strand()).to_string())
        }
        Stratifier::ReferenceBase => {
            if is_no_call(reference_base) {
                return None;
            }
            Part::text(sequence_base(reference_base, read.is_negative_strand()).to_string())
        }
        Stratifier::InsertLength => Part::number(i64::from(
            (read.bases.len() as i32 * 10).min(read.insert_size.abs()),
        )),
        Stratifier::SoftClips => Part::number(
            read.cigar
                .iter()
                .filter(|(_, operator)| *operator == 'S')
                .map(|(length, _)| *length as i64)
                .sum(),
        ),
    })
}

/// The one-based cycle a base was read at, counted from whichever end the machine read from.
pub fn cycle(read: &Read, offset: i64) -> usize {
    let length = read.bases.len() as i64;
    let zero_based = if read.is_negative_strand() {
        length - offset - 1
    } else {
        offset
    };
    (zero_based + 1).max(0) as usize
}

/// The read's GC, rounded to whole percents and reported as a fraction.
///
/// `SequenceUtil.calculateGc` divides by every base of the read, an N included, and the stratifier
/// then takes `Math.round(100 * gc) / 100`.
pub fn gc_content(bases: &[u8]) -> f64 {
    let gc = bases
        .iter()
        .filter(|base| matches!(base, b'C' | b'G' | b'c' | b'g'))
        .count();
    let fraction = gc as f64 / bases.len() as f64;
    (100.0 * fraction).round() / 100.0
}

/// `Double.toString`, as far as the covariates go: a whole number keeps its `.0`, and the rest
/// print the shortest form that round-trips.
fn format_double(value: f64) -> String {
    let text = format!("{value}");
    if text.contains('.') {
        text
    } else {
        format!("{text}.0")
    }
}

/// `ReadBaseStratification.getIndelElement`: the cigar element an insertion or a deletion was
/// recorded for, found from the read offset it carries, or nothing when the offset points at no
/// element.
fn indel_element(read: &Read, observation: &Observation) -> Option<(char, usize)> {
    let offset = observation.offset;
    if read.cigar.is_empty() {
        return None;
    }
    if offset == -1 {
        return Some((read.cigar[0].1, read.cigar[0].0));
    }
    let mut read_position: i64 = 0;
    for &(length, operator) in &read.cigar {
        if read_position > offset + 1 {
            return None;
        }
        match observation.alignment {
            AlignmentType::Insertion => {
                if consumes_read(operator) && read_position == offset {
                    return Some((operator, length));
                }
            }
            AlignmentType::Deletion => {
                if read_position == offset + 1 {
                    return Some((operator, length));
                }
            }
            AlignmentType::Match => return None,
        }
        if consumes_read(operator) {
            read_position += length as i64;
        }
    }
    None
}

/// One row of an `error_by_*` table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseErrorMetric {
    pub covariate: String,
    pub total_bases: u64,
    pub error_bases: u64,
    pub q_score: i32,
}

/// One row of an `overlapping_error_by_*` table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OverlappingErrorMetric {
    pub covariate: String,
    pub total_bases: u64,
    pub bases_with_overlapping_reads: u64,
    /// The two reads disagree with the reference and agree with each other, which is the template
    /// differing from the reference rather than an error.
    pub disagrees_with_reference_only: u64,
    pub disagrees_with_reference_only_q: i32,
    /// The read disagrees with the reference and its mate agrees with it, which is an error in
    /// this read.
    pub disagrees_with_ref_and_mate: u64,
    pub disagrees_with_ref_and_mate_q: i32,
    pub three_ways_disagreement: u64,
    pub three_ways_disagreement_q: i32,
}

/// One row of an `indel_error_by_*` table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndelErrorMetric {
    pub covariate: String,
    pub total_bases: u64,
    pub insertions: u64,
    pub inserted_bases: u64,
    pub insertions_q: i32,
    pub deletions: u64,
    pub deleted_bases: u64,
    pub deletions_q: i32,
    /// The inserted and deleted bases together.
    pub error_bases: u64,
    /// Always zero: the indel metric derives its own fields and does not derive this one, so the
    /// column inherited from the base metric is written as it was initialised.
    pub q_score: i32,
}

/// A table, whichever metric it is a table of.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Table {
    Base(Vec<BaseErrorMetric>),
    Overlapping(Vec<OverlappingErrorMetric>),
    Indel(Vec<IndelErrorMetric>),
}

/// The counts one stratum accumulates, before they are turned into a row.
#[derive(Debug, Default, Clone)]
struct Counts {
    bases: u64,
    mismatches: u64,
    insertions: u64,
    inserted_bases: u64,
    deletions: u64,
    deleted_bases: u64,
    overlapping_bases: u64,
    both_disagree: u64,
    disagrees_with_ref_and_mate: u64,
    three_ways: u64,
}

/// `OverlappingReadsErrorCalculator.areReadsMates`, which is asked one way round only: the first
/// read's mate position must be the second's start, and nothing checks the reverse.
fn are_mates(left: &Read, right: &Read) -> bool {
    left.name == right.name
        && left.is_paired()
        && left.is_first_of_pair() != right.is_first_of_pair()
        && left.is_second_of_pair() != right.is_second_of_pair()
        && !left.is_unmapped()
        && !right.is_unmapped()
        && !left.is_secondary()
        && !right.is_secondary()
        && left.mate_start == right.start
}

/// Run one calculator, split by one stratifier, over a pileup.
///
/// The reference is the contig, one base per position from `reference_start`.
pub fn collect(
    reads: &[Read],
    reference: &[u8],
    reference_start: i32,
    loci: &[Locus],
    calculator: Calculator,
    stratifier: Stratifier,
    options: &Options,
) -> Table {
    collect_joint(
        reads,
        reference,
        reference_start,
        loci,
        calculator,
        &[stratifier],
        options,
    )
}

/// Run one calculator, split by several stratifiers at once.
///
/// `ERROR:READ_ORDINALITY:CYCLE` folds its stratifiers into pairs from the left, so a stratum is
/// the tuple of every stratifier's answer, ordered left to right and printed with commas, and a
/// base is dropped when any one of them answers nothing. No stratifier at all is the single
/// stratum `all`.
pub fn collect_joint(
    reads: &[Read],
    reference: &[u8],
    reference_start: i32,
    loci: &[Locus],
    calculator: Calculator,
    stratifiers: &[Stratifier],
    options: &Options,
) -> Table {
    let prior = prior_error(options.prior_q);
    let mut strata: BTreeMap<Key, Counts> = BTreeMap::new();
    let all = [Stratifier::All];
    let stratifiers = if stratifiers.is_empty() {
        &all[..]
    } else {
        stratifiers
    };
    // A deletion spans several loci, and the same record is shown at each of them; it is counted
    // once, at the first, and only for the loci that were counted at all.
    let mut seen_deletions: HashSet<(usize, i32)> = HashSet::new();

    for locus in loci {
        let reference_base = reference
            .get((locus.position - reference_start) as usize)
            .copied()
            .unwrap_or(b'N');
        for observation in &locus.records {
            if observation.alignment == AlignmentType::Deletion {
                let already = seen_deletions.contains(&(observation.read, locus.position - 1));
                seen_deletions.insert((observation.read, locus.position));
                if already {
                    continue;
                }
            }
            let parts: Option<Vec<Part>> = stratifiers
                .iter()
                .map(|&stratifier| stratify(stratifier, reads, observation, reference_base))
                .collect();
            let Some(parts) = parts else { continue };
            let counts = strata.entry(Key(parts)).or_default();
            let read = &reads[observation.read];
            let base = usize::try_from(observation.offset)
                .ok()
                .and_then(|i| read.bases.get(i))
                .copied()
                .unwrap_or(b'N');
            let element = match observation.alignment {
                AlignmentType::Match => None,
                _ => indel_element(read, observation),
            };

            // Every calculator counts its denominator the same way: matched bases that were
            // called, and the whole length of an insertion.
            match observation.alignment {
                AlignmentType::Match => {
                    if !is_no_call(base) {
                        counts.bases += 1;
                    }
                }
                AlignmentType::Insertion => {
                    if let Some((_, length)) = element {
                        counts.bases += length as u64;
                    }
                }
                AlignmentType::Deletion => {}
            }

            match calculator {
                Calculator::Error => {
                    if observation.alignment == AlignmentType::Match
                        && !is_no_call(base)
                        && !bases_equal(base, reference_base)
                    {
                        counts.mismatches += 1;
                    }
                }
                Calculator::IndelError => match observation.alignment {
                    AlignmentType::Insertion => {
                        counts.insertions += 1;
                        if let Some((_, length)) = element {
                            counts.inserted_bases += length as u64;
                        }
                    }
                    AlignmentType::Deletion => {
                        counts.deletions += 1;
                        if let Some((_, length)) = element {
                            counts.deleted_bases += length as u64;
                        }
                    }
                    AlignmentType::Match => {}
                },
                Calculator::OverlappingError => {
                    // The mate is looked for among the bases READ at the locus, not among its
                    // indels: the tool builds its name sets from `getRecordAndOffsets()`.
                    let mate = locus.records.iter().find(|other| {
                        other.alignment == AlignmentType::Match
                            && other.read != observation.read
                            && are_mates(read, &reads[other.read])
                    });
                    let Some(mate) = mate else { continue };
                    let mate_base = usize::try_from(mate.offset)
                        .ok()
                        .and_then(|i| reads[mate.read].bases.get(i))
                        .copied()
                        .unwrap_or(b'N');
                    if is_no_call(base) || is_no_call(mate_base) {
                        continue;
                    }
                    counts.overlapping_bases += 1;
                    if bases_equal(base, reference_base) {
                        continue;
                    }
                    if bases_equal(base, mate_base) {
                        counts.both_disagree += 1;
                    } else if bases_equal(mate_base, reference_base) {
                        counts.disagrees_with_ref_and_mate += 1;
                    } else {
                        counts.three_ways += 1;
                    }
                }
            }
        }
    }

    match calculator {
        Calculator::Error => Table::Base(
            strata
                .into_iter()
                // A stratum with no bases above the base-quality threshold is not a row of a
                // simple error table, which is why a file of poorly mapped reads writes a table
                // with no header at all.
                .filter(|(_, counts)| counts.bases != 0)
                .map(|(key, counts)| BaseErrorMetric {
                    covariate: key.text(),
                    total_bases: counts.bases,
                    error_bases: counts.mismatches,
                    q_score: q_score(counts.mismatches, counts.bases, prior),
                })
                .collect(),
        ),
        Calculator::OverlappingError => Table::Overlapping(
            strata
                .into_iter()
                .map(|(key, counts)| OverlappingErrorMetric {
                    covariate: key.text(),
                    total_bases: counts.bases,
                    bases_with_overlapping_reads: counts.overlapping_bases,
                    disagrees_with_reference_only: counts.both_disagree,
                    disagrees_with_reference_only_q: q_score(
                        counts.both_disagree,
                        counts.overlapping_bases,
                        prior,
                    ),
                    disagrees_with_ref_and_mate: counts.disagrees_with_ref_and_mate,
                    disagrees_with_ref_and_mate_q: q_score(
                        counts.disagrees_with_ref_and_mate,
                        counts.overlapping_bases,
                        prior,
                    ),
                    three_ways_disagreement: counts.three_ways,
                    three_ways_disagreement_q: q_score(
                        counts.three_ways,
                        counts.overlapping_bases,
                        prior,
                    ),
                })
                .collect(),
        ),
        Calculator::IndelError => Table::Indel(
            strata
                .into_iter()
                .map(|(key, counts)| IndelErrorMetric {
                    covariate: key.text(),
                    total_bases: counts.bases,
                    insertions: counts.insertions,
                    inserted_bases: counts.inserted_bases,
                    insertions_q: q_score(counts.insertions, counts.bases, prior),
                    deletions: counts.deletions,
                    deleted_bases: counts.deleted_bases,
                    deletions_q: q_score(counts.deletions, counts.bases, prior),
                    error_bases: counts.inserted_bases + counts.deleted_bases,
                    q_score: 0,
                })
                .collect(),
        ),
    }
}

/// The table as the metrics file writes it: a header and one line per row, tab separated.
///
/// A table with no rows writes nothing at all, header included.
pub fn render(table: &Table) -> String {
    let (header, rows): (&str, Vec<String>) = match table {
        Table::Base(rows) => (
            "ERROR_BASES\tQ_SCORE\tCOVARIATE\tTOTAL_BASES",
            rows.iter()
                .map(|row| {
                    format!(
                        "{}\t{}\t{}\t{}",
                        row.error_bases, row.q_score, row.covariate, row.total_bases
                    )
                })
                .collect(),
        ),
        Table::Overlapping(rows) => (
            "NUM_BASES_WITH_OVERLAPPING_READS\tNUM_DISAGREES_WITH_REFERENCE_ONLY\t\
             DISAGREES_WITH_REFERENCE_ONLY_Q\tNUM_DISAGREES_WITH_REF_AND_MATE\t\
             DISAGREES_WITH_REF_AND_MATE_ONLY_Q\tNUM_THREE_WAYS_DISAGREEMENT\t\
             THREE_WAYS_DISAGREEMENT_ONLY_Q\tCOVARIATE\tTOTAL_BASES",
            rows.iter()
                .map(|row| {
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        row.bases_with_overlapping_reads,
                        row.disagrees_with_reference_only,
                        row.disagrees_with_reference_only_q,
                        row.disagrees_with_ref_and_mate,
                        row.disagrees_with_ref_and_mate_q,
                        row.three_ways_disagreement,
                        row.three_ways_disagreement_q,
                        row.covariate,
                        row.total_bases
                    )
                })
                .collect(),
        ),
        Table::Indel(rows) => (
            "NUM_INSERTIONS\tNUM_INSERTED_BASES\tINSERTIONS_Q\tNUM_DELETIONS\tNUM_DELETED_BASES\t\
             DELETIONS_Q\tERROR_BASES\tQ_SCORE\tCOVARIATE\tTOTAL_BASES",
            rows.iter()
                .map(|row| {
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        row.insertions,
                        row.inserted_bases,
                        row.insertions_q,
                        row.deletions,
                        row.deleted_bases,
                        row.deletions_q,
                        row.error_bases,
                        row.q_score,
                        row.covariate,
                        row.total_bases
                    )
                })
                .collect(),
        ),
    };
    if rows.is_empty() {
        return String::new();
    }
    let mut text = String::from(header);
    for row in rows {
        text.push('\n');
        text.push_str(&row);
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(cigar: Vec<(usize, char)>, bases: &[u8], qualities: Vec<u8>) -> Read {
        Read {
            name: "r".to_string(),
            start: 11,
            bases: bases.to_vec(),
            qualities,
            flags: 0,
            mate_start: 0,
            cigar,
            mapping_quality: 60,
            read_group: "rg1".to_string(),
            insert_size: 0,
        }
    }

    fn at(loci: &[Locus], position: i32) -> &Locus {
        loci.iter()
            .find(|locus| locus.position == position)
            .expect("a locus")
    }

    /// `SequenceUtil.calculateGc` divides by every base, so an N counts against the fraction.
    #[test]
    fn gc_divides_by_every_base() {
        assert_eq!(gc_content(b"GCGCNNAAAA"), 0.4);
        assert_eq!(gc_content(b"GGGGGGGGGG"), 1.0);
        assert_eq!(gc_content(b"ACGTACGTAA"), 0.4);
    }

    /// An insertion is shown at the base before it, a deletion at every position it spans and at
    /// the offset of the base before it.
    #[test]
    fn indels_are_shown_where_htsjdk_shows_them() {
        let reads = [read(
            vec![(5, 'M'), (2, 'I'), (3, 'M'), (2, 'D'), (4, 'M')],
            b"ACGTACCACGTACG",
            vec![40; 14],
        )];
        let loci = pileup(&reads, &Options::default());
        let insertion = at(&loci, 15)
            .records
            .iter()
            .find(|o| o.alignment == AlignmentType::Insertion)
            .expect("the insertion rides the fifth base");
        assert_eq!(insertion.offset, 5);
        assert!(at(&loci, 16)
            .records
            .iter()
            .all(|o| o.alignment == AlignmentType::Match));
        for position in [19, 20] {
            let deletion = at(&loci, position)
                .records
                .iter()
                .find(|o| o.alignment == AlignmentType::Deletion)
                .expect("the deletion spans two positions");
            assert_eq!(deletion.offset, 9);
        }
    }

    /// An insertion whose first base is below the cutoff is dropped, and the read offset does not
    /// move past it, so the deletion after it is shown at an offset that is two short.
    #[test]
    fn a_dropped_insertion_leaves_the_offsets_short() {
        let mut qualities = vec![40; 14];
        qualities[5] = 2;
        let reads = [read(
            vec![(5, 'M'), (2, 'I'), (3, 'M'), (2, 'D'), (4, 'M')],
            b"ACGTACCACGTACG",
            qualities,
        )];
        let loci = pileup(&reads, &Options::default());
        assert!(at(&loci, 15)
            .records
            .iter()
            .all(|o| o.alignment != AlignmentType::Insertion));
        let deletion = at(&loci, 19)
            .records
            .iter()
            .find(|o| o.alignment == AlignmentType::Deletion)
            .expect("the deletion");
        assert_eq!(deletion.offset, 7);
    }

    /// A joint stratum prints its parts with commas and sorts them left to right, numbers by
    /// value, so cycle 10 comes after cycle 2 rather than before it.
    #[test]
    fn a_joint_stratum_sorts_part_by_part() {
        let mut first = read(vec![(12, 'M')], b"ACGTACGTACGT", vec![40; 12]);
        first.flags = 0x1 | 0x40;
        let mut second = first.clone();
        second.flags = 0x1 | 0x80;
        let reads = [first, second];
        let loci = pileup(&reads, &Options::default());
        let reference: Vec<u8> = (0..40).map(|i| b"ACGT"[i % 4]).collect();
        let table = collect_joint(
            &reads,
            &reference,
            1,
            &loci,
            Calculator::Error,
            &[Stratifier::ReadOrdinality, Stratifier::Cycle],
            &Options::default(),
        );
        let Table::Base(rows) = table else {
            panic!("a simple error table")
        };
        let covariates: Vec<&str> = rows.iter().map(|r| r.covariate.as_str()).collect();
        assert_eq!(covariates[0], "FIRST,1");
        assert_eq!(covariates[1], "FIRST,2");
        assert_eq!(covariates[9], "FIRST,10");
        assert_eq!(covariates[12], "SECOND,1");
        assert_eq!(rows.len(), 24);
    }

    /// The draw happens at every locus, so a probability of nothing keeps nothing, and the same
    /// seed keeps the same loci twice.
    #[test]
    fn the_draw_is_seeded() {
        let reads = [read(vec![(40, 'M')], &[b'A'; 40], vec![40; 40])];
        let loci = pileup(&reads, &Options::default());
        let none = Options {
            probability: 0.0,
            ..Options::default()
        };
        assert!(processed_loci(loci.clone(), &none).is_empty());
        let half = Options {
            probability: 0.5,
            ..Options::default()
        };
        let kept = processed_loci(loci.clone(), &half);
        assert_eq!(kept, processed_loci(loci.clone(), &half));
        assert!(!kept.is_empty() && kept.len() < loci.len());
    }
}
