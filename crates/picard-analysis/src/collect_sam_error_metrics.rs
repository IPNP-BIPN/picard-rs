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

use std::collections::{BTreeMap, HashMap, HashSet};

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
#[derive(Debug, Clone, PartialEq, Eq, Default)]
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
    /// The read's and its mate's contig, which `PAIR_PROPERNESS` compares.
    pub reference_index: i32,
    pub mate_reference_index: i32,
    /// `NM`, which `MISMATCHES_IN_READ` reads.
    pub nm: Attribute,
    /// Whether an `SA` tag is present, which `PAIR_PROPERNESS` asks.
    pub has_sa: bool,
    /// The consensus-calling tags `aD`, `bD` and `cD`, which `CONSENSUS` reads.
    pub first_strand_copies: Attribute,
    pub second_strand_copies: Attribute,
    pub both_strands_copies: Attribute,
}

/// A tag read through `SAMRecord.getIntegerAttribute`, which takes any `Number` and throws on
/// anything else.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Attribute {
    #[default]
    Absent,
    Int(i64),
    /// A float, kept as its bits: `Number.longValue` truncates it.
    Float(u32),
    /// Not a `Number`: the class `getClass` names.
    Other(&'static str),
}

impl Attribute {
    /// `getIntegerAttribute(tag)`: `None` when absent, the exception when not an integer.
    pub fn integer(&self, tag: &str) -> Result<Option<i32>, String> {
        let long = match *self {
            Attribute::Absent => return Ok(None),
            Attribute::Int(value) => value,
            Attribute::Float(bits) => f32::from_bits(bits) as i64,
            Attribute::Other(class) => {
                return Err(format!(
                    "java.lang.RuntimeException: Value for tag {tag} is not Number: class {class}"
                ))
            }
        };
        i32::try_from(long).map(Some).map_err(|_| {
            format!(
                "java.lang.RuntimeException: Value for tag {tag} is not in Integer range: {long}"
            )
        })
    }
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
    /// For each record, the record of its mate's base here, once [`Locus::find_mates`] has looked;
    /// empty before, and then every overlapping-error aggregator looks for itself.
    pub mates: Vec<Option<u32>>,
}

impl Locus {
    /// Looks up every record's mate once, for all the aggregators that ask.
    pub fn find_mates(&mut self, reads: &[Read]) {
        self.mates = mates(reads, &self.records);
    }
}

/// For each record, its mate's base at the same locus, the first in the locus's order.
///
/// The mate is looked for among the bases READ at the locus, not among its indels: the tool builds
/// its name sets from `getRecordAndOffsets()`. Only a read of the same name can be the mate, so the
/// bases are indexed by name once, and only when some read here can have a mate at all.
fn mates(reads: &[Read], records: &[Observation]) -> Vec<Option<u32>> {
    let pairs = |read: &Read| read.is_paired() && !read.is_unmapped() && !read.is_secondary();
    let mut mates = vec![None; records.len()];
    if !records.iter().any(|record| pairs(&reads[record.read])) {
        return mates;
    }
    let mut names: HashMap<&str, Vec<u32>> = HashMap::new();
    for (index, other) in records.iter().enumerate() {
        if other.alignment == AlignmentType::Match {
            names
                .entry(reads[other.read].name.as_str())
                .or_default()
                .push(index as u32);
        }
    }
    for (record, mate) in records.iter().zip(&mut mates) {
        let read = &reads[record.read];
        if !pairs(read) {
            continue;
        }
        *mate = names.get(read.name.as_str()).and_then(|candidates| {
            candidates.iter().copied().find(|&other| {
                let other = &records[other as usize];
                other.read != record.read && are_mates(read, &reads[other.read])
            })
        });
    }
    mates
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
    /// `LONG_HOMOPOLYMER`, the shortest homopolymer `BINNED_HOMOPOLYMER` calls long.
    pub long_homopolymer: i32,
    /// `LOCATION_BIN_SIZE`, what `FLOWCELL_X` and `FLOWCELL_Y` divide a coordinate by.
    pub location_bin_size: i32,
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
            long_homopolymer: 6,
            location_bin_size: 2500,
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
    let mut builder = PileupBuilder::new(options);
    for (index, read) in reads.iter().enumerate() {
        builder.add(index, read);
    }
    builder.finish()
}

/// [`pileup`] a read at a time, for a walk over a coordinate-sorted file that hands on each locus
/// as soon as no later read can reach it.
pub struct PileupBuilder {
    loci: BTreeMap<i32, Pile>,
    cutoff: u8,
    min_mapping_q: u8,
}

impl PileupBuilder {
    pub fn new(options: &Options) -> PileupBuilder {
        PileupBuilder {
            loci: BTreeMap::new(),
            cutoff: options.min_base_q,
            min_mapping_q: options.min_mapping_q,
        }
    }

    /// One read, the `index`-th, added to the loci it covers.
    pub fn add(&mut self, index: usize, read: &Read) {
        let cutoff = self.cutoff;
        let passes = |read: &Read, offset: usize| {
            cutoff == 0
                || read.qualities.is_empty()
                || read.qualities.get(offset).is_none_or(|&q| q >= cutoff)
        };
        if read.flags & (0x100 | 0x800 | 0x400) != 0
            || read.is_unmapped()
            || read.mapping_quality < self.min_mapping_q
        {
            return;
        }
        // accumulateSamRecord: one observation per aligned base that meets the cutoff.
        let mut position = read.start;
        let mut offset = 0usize;
        for &(length, operator) in &read.cigar {
            if matches!(operator, 'M' | '=' | 'X') {
                for step in 0..length {
                    if passes(read, offset + step) {
                        self.loci
                            .entry(position + step as i32)
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
                        self.loci
                            .entry(position - 1)
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
                        self.loci
                            .entry(position + step as i32)
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

    /// How many loci are held.
    pub fn held(&self) -> usize {
        self.loci.len()
    }

    /// The loci before `position`, which a read starting at `position + 1` or later cannot reach:
    /// it shows an insertion at the base before its first, and nothing earlier.
    pub fn take_before(&mut self, position: i32) -> Vec<Locus> {
        let rest = self.loci.split_off(&position);
        let done = std::mem::replace(&mut self.loci, rest);
        done.into_iter().map(flatten).collect()
    }

    /// Every locus still held.
    pub fn finish(self) -> Vec<Locus> {
        self.loci.into_iter().map(flatten).collect()
    }
}

/// A pile as the walk shows it: the matched bases, then the deletions, then the insertions.
fn flatten((position, pile): (i32, Pile)) -> Locus {
    let mut records = pile.matches;
    records.extend(pile.deletions);
    records.extend(pile.insertions);
    Locus {
        position,
        records,
        mates: Vec::new(),
    }
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

/// `ReadBaseStratification.Stratifier`, all thirty-two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stratifier {
    All,
    GcContent,
    ReadOrdinality,
    ReadBase,
    ReadDirection,
    PairOrientation,
    PairProperness,
    ReferenceBase,
    PreDinuc,
    PostDinuc,
    HomopolymerLength,
    Homopolymer,
    BinnedHomopolymer,
    FlowcellTile,
    FlowcellY,
    FlowcellX,
    ReadGroup,
    Cycle,
    BinnedCycle,
    SoftClips,
    InsertLength,
    BaseQuality,
    MappingQuality,
    MismatchesInRead,
    OneBasePaddedContext,
    TwoBasePaddedContext,
    Consensus,
    NsInRead,
    InsertionsInRead,
    DeletionsInRead,
    IndelsInRead,
    IndelLength,
}

impl Stratifier {
    /// `ReadBaseStratification.Stratifier.valueOf`.
    pub fn parse(name: &str) -> Option<Stratifier> {
        Some(match name {
            "ALL" => Stratifier::All,
            "GC_CONTENT" => Stratifier::GcContent,
            "READ_ORDINALITY" => Stratifier::ReadOrdinality,
            "READ_BASE" => Stratifier::ReadBase,
            "READ_DIRECTION" => Stratifier::ReadDirection,
            "PAIR_ORIENTATION" => Stratifier::PairOrientation,
            "PAIR_PROPERNESS" => Stratifier::PairProperness,
            "REFERENCE_BASE" => Stratifier::ReferenceBase,
            "PRE_DINUC" => Stratifier::PreDinuc,
            "POST_DINUC" => Stratifier::PostDinuc,
            "HOMOPOLYMER_LENGTH" => Stratifier::HomopolymerLength,
            "HOMOPOLYMER" => Stratifier::Homopolymer,
            "BINNED_HOMOPOLYMER" => Stratifier::BinnedHomopolymer,
            "FLOWCELL_TILE" => Stratifier::FlowcellTile,
            "FLOWCELL_Y" => Stratifier::FlowcellY,
            "FLOWCELL_X" => Stratifier::FlowcellX,
            "READ_GROUP" => Stratifier::ReadGroup,
            "CYCLE" => Stratifier::Cycle,
            "BINNED_CYCLE" => Stratifier::BinnedCycle,
            "SOFT_CLIPS" => Stratifier::SoftClips,
            "INSERT_LENGTH" => Stratifier::InsertLength,
            "BASE_QUALITY" => Stratifier::BaseQuality,
            "MAPPING_QUALITY" => Stratifier::MappingQuality,
            "MISMATCHES_IN_READ" => Stratifier::MismatchesInRead,
            "ONE_BASE_PADDED_CONTEXT" => Stratifier::OneBasePaddedContext,
            "TWO_BASE_PADDED_CONTEXT" => Stratifier::TwoBasePaddedContext,
            "CONSENSUS" => Stratifier::Consensus,
            "NS_IN_READ" => Stratifier::NsInRead,
            "INSERTIONS_IN_READ" => Stratifier::InsertionsInRead,
            "DELETIONS_IN_READ" => Stratifier::DeletionsInRead,
            "INDELS_IN_READ" => Stratifier::IndelsInRead,
            "INDEL_LENGTH" => Stratifier::IndelLength,
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

/// The twenty-six directives a run collects when it is not told otherwise.
pub const DEFAULT_ERROR_METRICS: [&str; 26] = [
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

/// How the two halves of a pair order two strata, which is how the tool's sorted set does, and
/// what `toString` prints for each.
#[derive(Debug, Clone)]
enum Order {
    /// An `Integer` or a `Byte`, which sorts as a double does.
    Integer(i64),
    /// A `Double`.
    Double(f64),
    /// An enum, by declaration order, printed by name.
    Ordinal(u8, &'static str),
    /// A `String` or a `Character`.
    Text(Text),
}

/// A `String` or a `Character`, held inline when it is short, which every one is but a long read
/// group: a base is counted without a single allocation.
#[derive(Debug, Clone)]
enum Text {
    Inline(u8, [u8; 22]),
    Heap(Box<str>),
}

impl Text {
    fn new(text: &str) -> Text {
        let bytes = text.as_bytes();
        if bytes.len() <= 22 {
            let mut inline = [0; 22];
            inline[..bytes.len()].copy_from_slice(bytes);
            Text::Inline(bytes.len() as u8, inline)
        } else {
            Text::Heap(text.into())
        }
    }

    fn as_bytes(&self) -> &[u8] {
        match self {
            Text::Inline(length, bytes) => &bytes[..*length as usize],
            Text::Heap(text) => text.as_bytes(),
        }
    }

    fn as_str(&self) -> &str {
        match self {
            // Only ever copied whole from a `&str`.
            Text::Inline(length, bytes) => {
                std::str::from_utf8(&bytes[..*length as usize]).unwrap_or_default()
            }
            Text::Heap(text) => text,
        }
    }
}

/// One stratifier's answer for one base: what it sorts by, and what `toString` prints.
#[derive(Debug, Clone)]
pub struct Part(Order);

impl Part {
    fn number(value: i64) -> Part {
        Part(Order::Integer(value))
    }
    fn text(value: &str) -> Part {
        Part(Order::Text(Text::new(value)))
    }
}

impl std::fmt::Display for Part {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            Order::Integer(value) => write!(f, "{value}"),
            Order::Double(value) => f.write_str(&format_double(*value)),
            Order::Ordinal(_, name) => f.write_str(name),
            Order::Text(text) => f.write_str(text.as_str()),
        }
    }
}

impl Ord for Part {
    fn cmp(&self, other: &Part) -> std::cmp::Ordering {
        let number = |order: &Order| match *order {
            Order::Integer(value) => Some(value as f64),
            Order::Double(value) => Some(value),
            _ => None,
        };
        match (&self.0, &other.0) {
            (Order::Ordinal(a, _), Order::Ordinal(b, _)) => a.cmp(b),
            (Order::Text(a), Order::Text(b)) => a.as_bytes().cmp(b.as_bytes()),
            (a, b) => match (number(a), number(b)) {
                (Some(a), Some(b)) => a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal),
                _ => std::cmp::Ordering::Equal,
            },
        }
    }
}

impl PartialOrd for Part {
    fn partial_cmp(&self, other: &Part) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Part {
    fn eq(&self, other: &Part) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for Part {}

/// Equal parts hash alike, as far as one stratifier's parts go: they are all of one kind, and
/// none is NaN.
impl std::hash::Hash for Part {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let number = |value: f64| if value == 0.0 { 0f64 } else { value }.to_bits();
        match &self.0 {
            Order::Integer(value) => number(*value as f64).hash(state),
            Order::Double(value) => number(*value).hash(state),
            Order::Ordinal(index, _) => index.hash(state),
            Order::Text(text) => text.as_bytes().hash(state),
        }
    }
}

impl Part {
    fn is_nan(&self) -> bool {
        matches!(self.0, Order::Double(value) if value.is_nan())
    }
}

/// A stratum, one part per stratifier, as a `Pair` of pairs compares them left to right (every
/// stratum of an aggregator has as many) and prints them joined by a comma.
fn stratum_text(parts: &[Part]) -> String {
    parts
        .iter()
        .map(Part::to_string)
        .collect::<Vec<_>>()
        .join(",")
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

/// An enum constant: it sorts by declaration order and prints its name (or, for `ReadDirection`,
/// its sign).
fn ordinal(index: u8, text: &'static str) -> Part {
    Part(Order::Ordinal(index, text))
}

fn character(c: char) -> Part {
    Part::text(c.encode_utf8(&mut [0; 4]))
}

/// A stratifier's parts, when it answered, added to the stratum being built.
fn answer<const N: usize>(parts: &mut Vec<Part>, answer: Option<[Part; N]>) -> bool {
    match answer {
        Some(answer) => {
            parts.extend(answer);
            true
        }
        None => false,
    }
}

/// `AbstractRecordAndOffset.getReadBase` / `getBaseQuality`: the array read at the offset, which
/// throws where it has no element (a deletion shown at offset -1, a read stored without qualities).
fn element(array: &[u8], offset: i64) -> Result<u8, String> {
    usize::try_from(offset)
        .ok()
        .and_then(|i| array.get(i).copied())
        .ok_or_else(|| {
            format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {offset} out of bounds for length {}",
                array.len()
            )
        })
}

/// The bin one observation falls in, added to `parts`: `Ok(false)` drops it, `Err` is the
/// exception the reference throws on it. A stratifier built from two others answers both parts,
/// and nothing when either answers nothing.
pub fn stratify(
    stratifier: Stratifier,
    reads: &[Read],
    observation: &Observation,
    reference_base: u8,
    options: &Options,
    parts: &mut Vec<Part>,
) -> Result<bool, String> {
    let read = &reads[observation.read];
    let offset = observation.offset;
    let reverse = read.is_negative_strand();
    // `stratifyReadBase`: the base `delta` cycles away in the read's own direction, or nothing
    // off either end.
    let read_char = |delta: i64| -> Option<char> {
        let requested = offset + if reverse { -delta } else { delta };
        let base = usize::try_from(requested)
            .ok()
            .and_then(|i| read.bases.get(i).copied())?;
        Some(sequence_base(base, reverse))
    };
    let read_base = |delta: i64| read_char(delta).map(character);
    // `stratifyReferenceBase`.
    let reference_char =
        || (!is_no_call(reference_base)).then(|| sequence_base(reference_base, reverse));
    let reference = || reference_char().map(character);
    // `stratifyHomopolymerLength`: the run of bases read before this one that equal the one read
    // just before it.
    let homopolymer_length = || -> Option<i64> {
        if is_no_call(reference_base) {
            return None;
        }
        let length = read.bases.len() as i64;
        if offset < 0 || offset >= length {
            return None;
        }
        let bases = &read.bases;
        let mut run = offset;
        if !reverse {
            loop {
                run -= 1;
                if run < 0 || bases[run as usize] != bases[(offset - 1) as usize] {
                    break;
                }
            }
            Some(offset - run - 1)
        } else {
            loop {
                run += 1;
                if run >= length || bases[run as usize] != bases[(offset + 1) as usize] {
                    break;
                }
            }
            Some(run - offset - 1)
        }
    };
    let pre_dinuc = || -> Option<[Part; 2]> { Some([read_base(-1)?, reference()?]) };
    let cigar_bases = |operator: char| -> i64 {
        read.cigar
            .iter()
            .filter(|(_, op)| *op == operator)
            .map(|(length, _)| *length as i64)
            .sum()
    };
    macro_rules! one {
        ($part:expr $(,)?) => {{
            parts.push($part);
            Ok(true)
        }};
    }

    match stratifier {
        Stratifier::All => one!(Part::text("all")),
        Stratifier::BaseQuality => one!(Part::number(i64::from(element(&read.qualities, offset)?))),
        Stratifier::Cycle => one!(Part::number(cycle(read, offset) as i64)),
        Stratifier::BinnedCycle => {
            let relative = cycle_signed(read, offset) as f64 / read.bases.len() as f64;
            let bin = (0..5u8).find(|&i| {
                let lower = [0.0, 0.2, 0.4, 0.6, 0.8][i as usize];
                let upper = [0.2, 0.4, 0.6, 0.8, 1.0][i as usize];
                relative >= lower && relative <= upper
            });
            match bin {
                Some(i) => one!(ordinal(
                    i,
                    ["QUINTILE_1", "QUINTILE_2", "QUINTILE_3", "QUINTILE_4", "QUINTILE_5"]
                        [i as usize]
                )),
                None => Err(format!(
                    "java.lang.IllegalArgumentException: Value for CycleBin must be between 0 and 1 \
                     (inclusive), found: {}",
                    java_g(relative)
                )),
            }
        }
        Stratifier::GcContent => {
            let value = gc_content(&read.bases);
            one!(Part(Order::Double(value)))
        }
        Stratifier::ReadDirection => one!(if reverse {
            ordinal(1, "-")
        } else {
            ordinal(0, "+")
        }),
        Stratifier::ReadOrdinality => {
            if !read.is_paired() {
                return Ok(false);
            }
            one!(if read.is_first_of_pair() {
                ordinal(0, "FIRST")
            } else {
                ordinal(1, "SECOND")
            })
        }
        Stratifier::PairOrientation => {
            // `PairOrientation.of`.
            if !read.is_paired() || read.is_unmapped() || read.flags & 0x8 != 0 {
                return Ok(false);
            }
            let mate_positive = read.flags & 0x20 == 0;
            let (first_positive, second_positive) = if read.is_first_of_pair() {
                (!reverse, mate_positive)
            } else {
                (mate_positive, !reverse)
            };
            one!(match (first_positive == second_positive, first_positive) {
                (true, true) => ordinal(3, "R1R2"),
                (true, false) => ordinal(2, "F1F2"),
                (false, true) => ordinal(0, "F1R2"),
                (false, false) => ordinal(1, "F2R1"),
            })
        }
        Stratifier::PairProperness => {
            // `ProperPaired.of`.
            let unmapped = read.is_unmapped();
            let mate_unmapped = read.flags & 0x8 != 0;
            let clipped = matches!(read.cigar.first(), Some((_, 'S' | 'H')))
                || matches!(read.cigar.last(), Some((_, 'S' | 'H')));
            one!(if read.is_paired()
                && !mate_unmapped
                && !unmapped
                && read.mate_reference_index != read.reference_index
            {
                ordinal(3, "DISCORDANT")
            } else if !unmapped && clipped && read.has_sa {
                ordinal(2, "CHIMERIC")
            } else if unmapped || read.is_paired() && mate_unmapped {
                ordinal(4, "UNKNOWN")
            } else if !read.is_paired() {
                // `getProperPairFlag` on a read that is not paired.
                return Err(
                    "java.lang.IllegalStateException: Inappropriate call if not paired read"
                        .to_string(),
                );
            } else if read.flags & 0x2 == 0 {
                ordinal(1, "IMPROPER")
            } else {
                ordinal(0, "PROPER")
            },)
        }
        Stratifier::MappingQuality => one!(Part::number(i64::from(read.mapping_quality))),
        Stratifier::ReadGroup => one!(Part::text(&read.read_group)),
        Stratifier::ReadBase => Ok(answer(parts, read_base(0).map(|part| [part]))),
        Stratifier::ReferenceBase => Ok(answer(parts, reference().map(|part| [part]))),
        Stratifier::PreDinuc => Ok(answer(parts, pre_dinuc())),
        Stratifier::PostDinuc => Ok(answer(parts, (|| Some([reference()?, read_base(1)?]))())),
        Stratifier::HomopolymerLength => Ok(answer(
            parts,
            homopolymer_length().map(|n| [Part::number(n)]),
        )),
        Stratifier::Homopolymer => Ok(answer(
            parts,
            (|| {
                let [before, reference] = pre_dinuc()?;
                Some([Part::number(homopolymer_length()?), before, reference])
            })(),
        )),
        Stratifier::BinnedHomopolymer => Ok(answer(
            parts,
            (|| {
                let length = homopolymer_length()?;
                let [before, reference] = pre_dinuc()?;
                let bin = if length < i64::from(options.long_homopolymer) {
                    ordinal(0, "SHORT_HOMOPOLYMER")
                } else {
                    ordinal(1, "LONG_HOMOPOLYMER")
                };
                Some([bin, before, reference])
            })(),
        )),
        Stratifier::FlowcellTile | Stratifier::FlowcellX | Stratifier::FlowcellY => {
            // A name the default parser cannot read leaves a fresh `PhysicalLocationInt` as it
            // was made: tile, x and y all -1.
            let location = crate::mark_duplicates::location(&read.name);
            let (tile, x, y) = if location.known {
                (i32::from(location.tile), location.x, location.y)
            } else {
                (-1, -1, -1)
            };
            let coordinate = match stratifier {
                Stratifier::FlowcellTile => return one!(Part::number(i64::from(tile))),
                Stratifier::FlowcellX => x,
                _ => y,
            };
            if options.location_bin_size == 0 {
                return Err("java.lang.ArithmeticException: / by zero".to_string());
            }
            one!(Part::number(i64::from(
                coordinate.wrapping_div(options.location_bin_size),
            )))
        }
        Stratifier::InsertLength => one!(Part::number(i64::from(
            (read.bases.len() as i32 * 10).min(read.insert_size.wrapping_abs()),
        ))),
        Stratifier::SoftClips => one!(Part::number(cigar_bases('S'))),
        Stratifier::MismatchesInRead => {
            let Some(mismatches) = read.nm.integer("NM")? else {
                return Ok(false);
            };
            // Not counting the base itself, which is compared byte for byte.
            let base = element(&read.bases, offset)?;
            one!(Part::number(i64::from(if base != reference_base {
                mismatches.wrapping_sub(1)
            } else {
                mismatches
            })))
        }
        Stratifier::OneBasePaddedContext | Stratifier::TwoBasePaddedContext => {
            let pad = if stratifier == Stratifier::OneBasePaddedContext {
                1
            } else {
                2
            };
            // `stratifySurroundingContext`: the reference base in the middle, read bases around
            // it, and nothing when a read base is missing. A missing reference base is the
            // string "null" in the middle, which is what `StringBuilder.append` writes for it.
            let mut text = String::new();
            for delta in -pad..=pad {
                if delta == 0 {
                    match reference_char() {
                        Some(base) => text.push(base),
                        None => text.push_str("null"),
                    }
                } else {
                    let Some(base) = read_char(delta) else {
                        return Ok(false);
                    };
                    text.push(base);
                }
            }
            one!(Part::text(&text))
        }
        Stratifier::Consensus => {
            // `ConsensusStratifier.stratify`.
            // Both strand tags or neither: one without the other is not even read.
            let (first, second) = if read.first_strand_copies != Attribute::Absent
                && read.second_strand_copies != Attribute::Absent
            {
                (
                    read.first_strand_copies.integer("aD")?.unwrap_or(0),
                    read.second_strand_copies.integer("bD")?.unwrap_or(0),
                )
            } else {
                (0, 0)
            };
            let both = read.both_strands_copies.integer("cD")?.unwrap_or(0);
            one!(if both == 1 {
                ordinal(0, "SIMPLEX_SINGLETON")
            } else if second == 0 && both > 1 {
                ordinal(1, "SIMPLEX_CONSENSUS")
            } else if first > 0 && second > 0 && (first == 1 || second == 1) {
                ordinal(2, "DUPLEX_SINGLETON")
            } else if first > 1 && second > 1 {
                ordinal(3, "DUPLEX_CONSENSUS")
            } else {
                ordinal(4, "UNKNOWN")
            })
        }
        Stratifier::NsInRead => one!(Part::number(
            read.bases.iter().filter(|&&base| is_no_call(base)).count() as i64,
        )),
        Stratifier::InsertionsInRead => one!(Part::number(cigar_bases('I'))),
        Stratifier::DeletionsInRead => one!(Part::number(cigar_bases('D'))),
        Stratifier::IndelsInRead => one!(Part::number(cigar_bases('I') + cigar_bases('D'))),
        Stratifier::IndelLength => {
            // `stratifyIndelLength`.
            let wanted = match observation.alignment {
                AlignmentType::Match => return one!(Part::number(0)),
                AlignmentType::Insertion => 'I',
                AlignmentType::Deletion => 'D',
            };
            match indel_element(read, observation) {
                None => Ok(false),
                Some((operator, _)) if operator != wanted => Err(
                    "java.lang.IllegalStateException: Wrong CIGAR operator for the given position."
                        .to_string(),
                ),
                Some((_, length)) => one!(Part::number(length as i64)),
            }
        }
    }
}

/// `stratifyCycle` without the floor [`cycle`] puts on it: `BINNED_CYCLE` divides the raw value.
fn cycle_signed(read: &Read, offset: i64) -> i64 {
    let length = read.bases.len() as i64;
    let zero_based = if read.is_negative_strand() {
        length - offset - 1
    } else {
        offset
    };
    zero_based + 1
}

/// `String.format("%g", value)`: six significant digits, scientific outside [1e-4, 1e6).
fn java_g(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    if value == 0.0 {
        return "0.00000".to_string();
    }
    let exponent = value.abs().log10().floor() as i32;
    if (-4..6).contains(&exponent) {
        format!("{:.*}", (5 - exponent).max(0) as usize, value)
    } else {
        let text = format!("{value:.5e}");
        let (mantissa, exp) = text.split_once('e').unwrap_or((&text, "0"));
        let exp: i32 = exp.parse().unwrap_or(0);
        format!(
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    }
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
    try_collect_joint(
        reads,
        reference,
        reference_start,
        loci,
        calculator,
        stratifiers,
        options,
    )
    .unwrap_or_else(|thrown| panic!("{}", thrown.message))
}

/// The exception a run throws, and where: the locus (in the order the run takes them) and the
/// observation within it. The tool hands every observation to every aggregator before it moves
/// on, so of several aggregators the one that throws first is the one at the earliest place.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Thrown {
    pub locus: usize,
    pub observation: usize,
    pub message: String,
}

/// [`collect_joint`], with the exception the reference throws on a base as `Err`.
pub fn try_collect_joint(
    reads: &[Read],
    reference: &[u8],
    reference_start: i32,
    loci: &[Locus],
    calculator: Calculator,
    stratifiers: &[Stratifier],
    options: &Options,
) -> Result<Table, Thrown> {
    let mut joint = Joint::new(calculator, stratifiers, options);
    joint.add(reads, reference, reference_start, loci, 0, 0, options)?;
    Ok(joint.finish())
}

/// One aggregator's state across the walk, for a run that hands it the loci a batch at a time:
/// [`try_collect_joint`] is one [`Joint::add`] of every locus.
pub struct Joint {
    calculator: Calculator,
    stratifiers: Vec<Stratifier>,
    prior: f64,
    // The strata as the tool's sorted map holds them, each the index of its counts.
    strata: BTreeMap<Vec<Part>, usize>,
    counts: Vec<Counts>,
    // The stratum each set of parts already found, so most bases skip the sorted map. Only while
    // no part is NaN: a NaN is equal to every number, and where the map then finds it depends on
    // the map's shape, so from then on every base asks the map itself.
    found: foldhash::HashMap<Vec<Part>, usize>,
    holds_nan: bool,
    // A deletion spans several loci, and the same record is shown at each of them; it is counted
    // once, at the first, and only for the loci that were counted at all.
    seen_deletions: HashSet<(usize, i32)>,
}

impl Joint {
    pub fn new(calculator: Calculator, stratifiers: &[Stratifier], options: &Options) -> Joint {
        Joint {
            calculator,
            stratifiers: if stratifiers.is_empty() {
                vec![Stratifier::All]
            } else {
                stratifiers.to_vec()
            },
            prior: prior_error(options.prior_q),
            strata: BTreeMap::new(),
            counts: Vec::new(),
            found: foldhash::HashMap::default(),
            holds_nan: false,
            seen_deletions: HashSet::new(),
        }
    }

    /// The next `loci` of the walk, the first of which is its `first_locus`-th: a throw names its
    /// place in the whole walk. `reads[0]` is the walk's `first_read`-th read, which is what lets
    /// a deletion be followed across batches whose reads are numbered from where each starts.
    #[allow(clippy::too_many_arguments)]
    pub fn add(
        &mut self,
        reads: &[Read],
        reference: &[u8],
        reference_start: i32,
        loci: &[Locus],
        first_locus: usize,
        first_read: usize,
        options: &Options,
    ) -> Result<(), Thrown> {
        let (calculator, stratifiers) = (self.calculator, &self.stratifiers[..]);
        let (strata, seen_deletions) = (&mut self.strata, &mut self.seen_deletions);
        let (all_counts, found, holds_nan) =
            (&mut self.counts, &mut self.found, &mut self.holds_nan);
        let mut parts = Vec::new();
        // Only a deletion seen at the locus before this batch's first can still be met again.
        if let Some(locus) = loci.first() {
            seen_deletions.retain(|&(_, position)| position >= locus.position - 1);
        }
        for (locus_index, locus) in loci.iter().enumerate() {
            let mut own_mates: Option<Vec<Option<u32>>> = None;
            let reference_base = reference
                .get((locus.position - reference_start) as usize)
                .copied()
                .unwrap_or(b'N');
            for (observation_index, observation) in locus.records.iter().enumerate() {
                let thrown = |message: String| Thrown {
                    locus: first_locus + locus_index,
                    observation: observation_index,
                    message,
                };
                if observation.alignment == AlignmentType::Deletion {
                    let read = first_read + observation.read;
                    let already = seen_deletions.contains(&(read, locus.position - 1));
                    seen_deletions.insert((read, locus.position));
                    if already {
                        continue;
                    }
                }
                // Every stratifier is asked, left to right, before any answer of nothing drops the
                // base: the pairs are built from the left and each asks both of its halves.
                // `addRecordAndOffset` reads the base quality of anything but a deletion before any
                // aggregator sees it, which throws for a read stored without qualities.
                if observation.alignment != AlignmentType::Deletion {
                    element(&reads[observation.read].qualities, observation.offset)
                        .map_err(thrown)?;
                }
                parts.clear();
                let mut answered = true;
                for &stratifier in stratifiers {
                    answered &= stratify(
                        stratifier,
                        reads,
                        observation,
                        reference_base,
                        options,
                        &mut parts,
                    )
                    .map_err(thrown)?;
                }
                if !answered {
                    continue;
                }
                let cached = !*holds_nan && !parts.iter().any(Part::is_nan);
                let known = if cached {
                    found.get(parts.as_slice()).copied()
                } else {
                    None
                };
                let index = known.unwrap_or_else(|| {
                    let index = match strata.get(parts.as_slice()) {
                        Some(&index) => index,
                        None => {
                            strata.insert(parts.clone(), all_counts.len());
                            all_counts.push(Counts::default());
                            *holds_nan |= !cached;
                            all_counts.len() - 1
                        }
                    };
                    if cached && !*holds_nan {
                        found.insert(parts.clone(), index);
                    }
                    index
                });
                let counts = &mut all_counts[index];
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
                        let found = if locus.mates.is_empty() {
                            &own_mates.get_or_insert_with(|| mates(reads, &locus.records))[..]
                        } else {
                            &locus.mates[..]
                        };
                        let Some(mate) = found[observation_index] else {
                            continue;
                        };
                        let mate = &locus.records[mate as usize];
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

        Ok(())
    }

    /// The table the walk has built.
    pub fn finish(self) -> Table {
        let (calculator, prior, counts) = (self.calculator, self.prior, self.counts);
        let strata = self
            .strata
            .into_iter()
            .map(|(key, index)| (key, counts[index].clone()));
        match calculator {
            Calculator::Error => Table::Base(
                strata
                    .into_iter()
                    // A stratum with no bases above the base-quality threshold is not a row of a
                    // simple error table, which is why a file of poorly mapped reads writes a table
                    // with no header at all.
                    .filter(|(_, counts)| counts.bases != 0)
                    .map(|(key, counts)| BaseErrorMetric {
                        covariate: stratum_text(&key),
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
                        covariate: stratum_text(&key),
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
                        covariate: stratum_text(&key),
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
            ..Read::default()
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
