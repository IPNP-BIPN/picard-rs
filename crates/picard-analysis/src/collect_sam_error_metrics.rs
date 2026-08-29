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

use std::collections::BTreeMap;

/// The error probability a phred-scaled prior stands for: `PRIOR_Q` of 30 is one error in a
/// thousand.
pub fn prior_error(prior_q: i32) -> f64 {
    10f64.powf(-f64::from(prior_q) / 10.0)
}

/// `QualityUtil.getPhredScoreFromErrorProbability`, rounded the way Java rounds it.
pub fn phred_from_error_probability(probability: f64) -> i32 {
    (-10.0 * probability.log10()).round() as i32
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
    /// One-based, on the single contig these fixtures use.
    pub start: i32,
    pub bases: Vec<u8>,
    /// Phred, not ASCII.
    pub qualities: Vec<u8>,
    pub flags: u16,
    pub mate_start: i32,
    pub cigar: Vec<(usize, char)>,
    pub mapping_quality: u8,
}

impl Read {
    pub fn is_paired(&self) -> bool {
        self.flags & 0x1 != 0
    }
    pub fn is_first_of_pair(&self) -> bool {
        self.flags & 0x40 != 0
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
    /// The offset into the read's bases, which is what the cycle is counted from.
    pub offset: usize,
    pub alignment: AlignmentType,
    /// The length of the cigar element an insertion or a deletion belongs to.
    pub indel_length: usize,
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
}

impl Default for Options {
    fn default() -> Self {
        Options {
            min_mapping_q: 20,
            min_base_q: 20,
            prior_q: 30,
            max_loci: 0,
            known_sites: Vec::new(),
        }
    }
}

/// The pileup the tool walks: every locus a read covers, indels included, with the two thresholds
/// already applied.
///
/// The thresholds drop observations rather than loci, one by read and one by base, which is why a
/// mismatch below `--MIN_BASE_Q` lowers the denominator instead of raising the error count.
pub fn pileup(reads: &[Read], options: &Options) -> Vec<Locus> {
    let mut loci: BTreeMap<i32, Vec<Observation>> = BTreeMap::new();
    for (index, read) in reads.iter().enumerate() {
        if read.mapping_quality < options.min_mapping_q || read.is_unmapped() {
            continue;
        }
        let mut position = read.start;
        let mut offset = 0usize;
        for &(length, operator) in &read.cigar {
            match operator {
                'M' | '=' | 'X' => {
                    for step in 0..length {
                        if read.qualities[offset + step] >= options.min_base_q {
                            loci.entry(position + step as i32)
                                .or_default()
                                .push(Observation {
                                    read: index,
                                    offset: offset + step,
                                    alignment: AlignmentType::Match,
                                    indel_length: 0,
                                });
                        }
                    }
                    position += length as i32;
                    offset += length;
                }
                'D' | 'N' => {
                    for step in 0..length {
                        loci.entry(position + step as i32)
                            .or_default()
                            .push(Observation {
                                read: index,
                                offset,
                                alignment: AlignmentType::Deletion,
                                indel_length: length,
                            });
                    }
                    position += length as i32;
                }
                'I' => {
                    loci.entry(position).or_default().push(Observation {
                        read: index,
                        offset,
                        alignment: AlignmentType::Insertion,
                        indel_length: length,
                    });
                    offset += length;
                }
                'S' => offset += length,
                _ => {}
            }
        }
    }
    loci.into_iter()
        .map(|(position, records)| Locus { position, records })
        .collect()
}

/// The loci a run actually counts: the known sites removed, and then the cap applied.
///
/// The order matters. A skipped locus is not a processed one, so `--MAX_LOCI` counts what is left
/// after the VCF has taken its sites out.
pub fn processed_loci(loci: Vec<Locus>, options: &Options) -> Vec<Locus> {
    let mut kept = Vec::new();
    for locus in loci {
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

/// The stratifiers whose binning is ported, which are the ones the goldens exercise.
///
/// Every stratifier's file suffix is in [`stratifier_suffix`]; these four are the ones that also
/// put a base in a bin here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stratifier {
    All,
    BaseQuality,
    Cycle,
    GcContent,
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
        "REFERENCE_BASE" => "reference_base",
        "PRE_DINUC" => "pre_dinuc",
        "POST_DINUC" => "post_dinuc",
        "HOMOPOLYMER_LENGTH" => "homopolymer_length",
        "HOMOPOLYMER" => "homopolymer_and_following_ref_base",
        "BINNED_HOMOPOLYMER" => "binned_length_homopolymer_and_following_ref_base",
        "FLOWCELL_TILE" => "tile",
        "FLOWCELL_X" => "flowcell_x",
        "FLOWCELL_Y" => "flowcell_y",
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
        "INSERTIONS_IN_READ" => "insertions_in_read",
        "DELETIONS_IN_READ" => "deletions_in_read",
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

/// The bin one observation falls in, or nothing, which drops it.
pub fn stratify(
    stratifier: Stratifier,
    reads: &[Read],
    observation: &Observation,
) -> Option<String> {
    let read = &reads[observation.read];
    Some(match stratifier {
        Stratifier::All => "all".to_string(),
        Stratifier::BaseQuality => read.qualities.get(observation.offset)?.to_string(),
        Stratifier::Cycle => cycle(read, observation.offset).to_string(),
        Stratifier::GcContent => format_double(gc_content(&read.bases)),
    })
}

/// The one-based cycle a base was read at, counted from whichever end the machine read from.
pub fn cycle(read: &Read, offset: usize) -> usize {
    if read.is_negative_strand() {
        read.bases.len() - offset - 1 + 1
    } else {
        offset + 1
    }
}

/// The read's GC, rounded to whole percents and reported as a fraction.
pub fn gc_content(bases: &[u8]) -> f64 {
    let counted = bases
        .iter()
        .filter(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
        .count();
    if counted == 0 {
        return 0.0;
    }
    let gc = bases
        .iter()
        .filter(|base| matches!(base.to_ascii_uppercase(), b'C' | b'G'))
        .count();
    (100.0 * gc as f64 / counted as f64).round() / 100.0
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

/// Whether two records are the two halves of one template.
fn are_mates(left: &Read, right: &Read) -> bool {
    left.name == right.name
        && left.is_paired()
        && left.is_first_of_pair() != right.is_first_of_pair()
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
    let prior = prior_error(options.prior_q);
    let mut strata: BTreeMap<Key, Counts> = BTreeMap::new();
    // A deletion spans several loci, and the same record is shown at each of them; it is counted
    // once, at the first.
    let mut seen_deletions: Vec<(usize, i32)> = Vec::new();

    for locus in loci {
        let reference_base = reference[(locus.position - reference_start) as usize];
        for observation in &locus.records {
            let Some(stratum) = stratify(stratifier, reads, observation) else {
                continue;
            };
            let counts = strata.entry(Key::new(&stratum)).or_default();
            let read = &reads[observation.read];
            let base = read.bases.get(observation.offset).copied().unwrap_or(b'N');

            // Every calculator counts its denominator the same way: matched bases that were
            // called, and the whole length of an insertion.
            match observation.alignment {
                AlignmentType::Match => {
                    if !is_no_call(base) {
                        counts.bases += 1;
                    }
                }
                AlignmentType::Insertion => counts.bases += observation.indel_length as u64,
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
                        counts.inserted_bases += observation.indel_length as u64;
                    }
                    AlignmentType::Deletion => {
                        let previous = (observation.read, locus.position - 1);
                        if !seen_deletions.contains(&previous) {
                            counts.deletions += 1;
                            counts.deleted_bases += observation.indel_length as u64;
                        }
                        seen_deletions.push((observation.read, locus.position));
                    }
                    AlignmentType::Match => {}
                },
                Calculator::OverlappingError => {
                    let mate = locus.records.iter().find(|other| {
                        other.read != observation.read
                            && are_mates(read, &reads[other.read])
                            && are_mates(&reads[other.read], read)
                    });
                    let Some(mate) = mate else { continue };
                    let mate_base = reads[mate.read]
                        .bases
                        .get(mate.offset)
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
                    covariate: key.text,
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
                    covariate: key.text,
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
                    covariate: key.text,
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

/// A stratum, ordered the way the tool's own sorted set orders it: numbers by value and everything
/// else by text.
#[derive(Debug, Clone, PartialEq)]
struct Key {
    text: String,
    number: Option<f64>,
}

impl Key {
    fn new(text: &str) -> Key {
        Key {
            text: text.to_string(),
            number: text.parse::<f64>().ok(),
        }
    }
}

impl Eq for Key {}

impl Ord for Key {
    fn cmp(&self, other: &Key) -> std::cmp::Ordering {
        match (self.number, other.number) {
            (Some(left), Some(right)) => left
                .partial_cmp(&right)
                .unwrap_or(std::cmp::Ordering::Equal),
            _ => self.text.cmp(&other.text),
        }
    }
}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Key) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
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
