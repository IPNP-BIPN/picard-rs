//! `UmiAwareMarkDuplicatesWithMateCigar`: duplicates split by the molecule they came from.
//!
//! A UMI is a barcode on the MOLECULE rather than on the library, so two reads at one position
//! carrying different UMIs came from different molecules and are not duplicates of each other.
//! This tool is [`crate::mate_cigar_duplicates`]'s simple sibling with that split applied inside
//! each duplicate set, and it writes a second metrics table about the UMIs it saw.
//!
//! # The split has a threshold
//!
//! Two UMIs a single base apart are one molecule and a sequencing error, not two molecules:
//! `MAX_EDIT_DISTANCE_TO_JOIN` is one by default, and the join is TRANSITIVE, so `ATCC`, `AACC`
//! and `AACG` end up in one set even though the first and the last are two apart. That is
//! union-find over a graph whose edges are the pairs within the distance, and it is why the
//! threshold cannot be applied pairwise.
//!
//! Each set is then ASSIGNED a UMI, which is the most frequent one in it, counted over the whole
//! position rather than over the set. A UMI containing `N` is never assigned while another choice
//! exists, and where every choice has one the fewest `N`s wins.
//!
//! # What is not ported
//!
//! The sets themselves come from htsjdk's `SAMRecordDuplicateComparator` in the reference, and
//! here they are cut on the same key `MarkDuplicates` uses for a fragment: library, barcode,
//! reference, 5' coordinate and orientation. Every fixture in the golden is a set of records at
//! one position, where the two agree; a file whose sets straddle positions would need the
//! comparator itself.
//!
//! `DUPLEX_UMI` is not ported either: a duplex UMI is two UMIs with a hyphen between them, and it
//! is normalised by the strand the read is on before anything else happens. No case measures it.
//!
//! Ported from `picard.sam.markduplicates.UmiAwareMarkDuplicatesWithMateCigar`,
//! `picard.sam.markduplicates.UmiGraph`, `picard.sam.markduplicates.UmiUtil`,
//! `picard.sam.markduplicates.UmiMetrics` and `htsjdk.samtools.util.QualityUtil` in Picard 3.4.0.

use std::collections::BTreeMap;

use crate::mark_duplicates::{mark, sets_by_position, Marking, Options, Record};
use crate::mate_cigar_duplicates::{needs_mate_cigar, Refusal, SortOrder};

/// The tool's own arguments.
#[derive(Debug, Clone)]
pub struct UmiOptions {
    pub base: Options,
    /// `MAX_EDIT_DISTANCE_TO_JOIN`, one by default.
    pub max_edit_distance_to_join: i32,
    /// `UMI_TAG_NAME`, `RX` by default.
    pub umi_tag: String,
    /// `MOLECULAR_IDENTIFIER_TAG`, which is written back where it is named.
    pub molecular_identifier_tag: Option<String>,
    /// `ALLOW_MISSING_UMIS`, false by default and documented as being for testing only.
    pub allow_missing_umis: bool,
}

impl Default for UmiOptions {
    fn default() -> Self {
        Self {
            base: Options::default(),
            max_edit_distance_to_join: 1,
            umi_tag: "RX".to_string(),
            molecular_identifier_tag: None,
            allow_missing_umis: false,
        }
    }
}

/// What the tool refuses a file for, beyond what its parent refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UmiRefusal {
    /// The parent's, for a file that is not coordinate sorted or has no mate cigar.
    Inherited(Refusal),
    /// `PicardException`, thrown by `UmiGraph` on the first record without a UMI.
    MissingUmi { read: String, tag: String },
    /// `PicardException`, thrown by `UmiUtil` on a UMI that is not made of bases.
    IllegalUmi,
}

impl UmiRefusal {
    pub fn message(&self) -> String {
        match self {
            UmiRefusal::Inherited(refusal) => refusal.message(),
            UmiRefusal::MissingUmi { read, tag } => {
                format!("Read {read} does not contain a UMI with the {tag} attribute.")
            }
            UmiRefusal::IllegalUmi => {
                "UMI found with illegal characters.  UMIs must match the regular expression \
                 ^[ATCGNatcgn-]*$."
                    .to_string()
            }
        }
    }

    pub fn exception(&self) -> &'static str {
        match self {
            UmiRefusal::Inherited(refusal) => refusal.exception(),
            _ => "picard.PicardException",
        }
    }
}

/// `StringUtil.isWithinHammingDistance`: the same length, and at most `distance` differences.
///
/// An `N` is a base like any other here, so `AANA` is one away from `AAAA` and joins it at the
/// default distance.
pub fn within_hamming_distance(left: &str, right: &str, distance: i32) -> bool {
    if left.len() != right.len() {
        return false;
    }
    let mut differences = 0;
    for (a, b) in left.bytes().zip(right.bytes()) {
        if a != b {
            differences += 1;
            if differences > distance {
                return false;
            }
        }
    }
    true
}

/// `UmiGraph.joinUmisIntoDuplicateSets`: union-find over the pairs within the distance.
///
/// The clusters are returned as one identifier per UMI, in the order the UMIs were given. The join
/// is transitive, which is the property a pairwise test would not have.
pub fn cluster(umis: &[String], max_edit_distance_to_join: i32) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..umis.len()).collect();
    fn find(parent: &mut [usize], index: usize) -> usize {
        let mut root = index;
        while parent[root] != root {
            root = parent[root];
        }
        let mut walk = index;
        while parent[walk] != root {
            let next = parent[walk];
            parent[walk] = root;
            walk = next;
        }
        root
    }
    for left in 0..umis.len() {
        for right in (left + 1)..umis.len() {
            if within_hamming_distance(&umis[left], &umis[right], max_edit_distance_to_join) {
                let a = find(&mut parent, left);
                let b = find(&mut parent, right);
                if a != b {
                    parent[a] = b;
                }
            }
        }
    }
    (0..umis.len())
        .map(|index| find(&mut parent, index))
        .collect()
}

/// The UMI a set is assigned: the most frequent one, counted over the whole position.
///
/// A UMI with an `N` is never assigned while another choice exists, and where every choice has one
/// the fewest `N`s wins. The count is `>` rather than `>=`, so the first of an equal-frequency set
/// is the one assigned.
pub fn assigned_umi(set: &[String], counts: &BTreeMap<String, u64>) -> Option<String> {
    let mut assigned: Option<String> = None;
    let mut max_count = 0u64;
    let mut fewest_n: Option<String> = None;
    let mut n_count = 0usize;
    for umi in set {
        let ns = umi.matches('N').count();
        if ns > 0 {
            if n_count == 0 || ns < n_count {
                n_count = ns;
                fewest_n = Some(umi.clone());
            }
        } else if counts.get(umi).copied().unwrap_or(0) > max_count {
            max_count = counts.get(umi).copied().unwrap_or(0);
            assigned = Some(umi.clone());
        }
    }
    assigned.or(fewest_n)
}

/// `UmiUtil.setMolecularIdentifier`: the contig, the OTHER end's position, and the assigned UMI.
///
/// The position is the record's own on a reverse read and its MATE's on a forward one, which is
/// the same position for both ends of a pair and is what makes the identifier shared.
pub fn molecular_identifier(
    contig: &str,
    alignment_start: i32,
    mate_alignment_start: i32,
    reverse_strand: bool,
    assigned: &str,
) -> String {
    let position = if reverse_strand {
        alignment_start
    } else {
        mate_alignment_start
    };
    format!("{contig}:{position}/{assigned}")
}

/// `UmiMetrics`, the row one library gets beside the duplication metrics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UmiMetrics {
    pub library: String,
    pub mean_umi_length: f64,
    pub observed_unique_umis: i64,
    pub inferred_unique_umis: i64,
    pub observed_base_errors: i64,
    pub duplicate_sets_ignoring_umi: i64,
    pub duplicate_sets_with_umi: i64,
    pub observed_umi_entropy: f64,
    pub inferred_umi_entropy: f64,
    pub umi_base_qualities: i32,
    pub percent_umi_with_n: f64,
}

/// `MathUtil.LOG_4_BASE_E`, which turns an entropy in nats into one in bases.
fn log_4() -> f64 {
    4f64.ln()
}

/// `UmiMetrics.effectiveNumberOfBases`: Shannon entropy over the counts, in base four.
pub fn effective_number_of_bases(counts: &BTreeMap<String, u64>) -> f64 {
    let total: u64 = counts.values().sum();
    if total == 0 {
        return 0.0;
    }
    let entropy: f64 = counts
        .values()
        .map(|count| {
            let p = *count as f64 / total as f64;
            -p * p.ln()
        })
        .sum();
    entropy / log_4()
}

/// `QualityUtil.getPhredScoreFromErrorProbability`, INCLUDING its overflow.
///
/// `Math.round` returns a `long`, and a probability of zero makes it `Long.MAX_VALUE`; the cast to
/// `int` keeps the low thirty-two bits of that, which is `-1`. Every fixture with no base errors
/// reports `-1` in the golden for exactly that reason, so the overflow is the measured behaviour
/// rather than an accident to be tidied away.
pub fn phred_from_error_probability(probability: f64) -> i32 {
    htsjdk_bam::quality_util::phred_score_from_error_probability(probability)
}

/// What one run decided: the marking, and the UMI metrics beside it.
#[derive(Debug, Clone, PartialEq)]
pub struct UmiMarking {
    pub marking: Marking,
    /// The identifier each record was assigned, where the tag was named.
    pub molecular_identifiers: Vec<Option<String>>,
    /// One row per library, in the order the libraries were first seen.
    pub umi_metrics: Vec<UmiMetrics>,
}

/// `UmiAwareMarkDuplicatesWithMateCigar.doWork`, over records already in memory.
pub fn mark_with_umis(
    records: &[Record],
    order: SortOrder,
    options: &UmiOptions,
) -> Result<UmiMarking, UmiRefusal> {
    if order != SortOrder::Coordinate {
        return Err(UmiRefusal::Inherited(Refusal::NotCoordinateSorted));
    }
    if let Some(record) = records
        .iter()
        .find(|record| needs_mate_cigar(record) && record.mate_cigar.is_none())
    {
        return Err(UmiRefusal::Inherited(Refusal::MateCigarNotFound {
            read: record.name.clone(),
        }));
    }
    for record in records {
        if record.barcode.is_none() && !options.allow_missing_umis {
            return Err(UmiRefusal::MissingUmi {
                read: record.name.clone(),
                tag: options.umi_tag.clone(),
            });
        }
        if let Some(umi) = &record.barcode {
            if !umi.bytes().all(|base| b"ATCGNatcgn-".contains(&base)) {
                return Err(UmiRefusal::IllegalUmi);
            }
        }
    }

    // The sets, and inside each of them the UMI clusters.
    let mut assigned: Vec<Option<String>> = vec![None; records.len()];
    let mut identifiers: Vec<Option<String>> = vec![None; records.len()];
    let mut observed: BTreeMap<String, u64> = BTreeMap::new();
    let mut inferred: BTreeMap<String, u64> = BTreeMap::new();
    let mut metrics = UmiMetrics {
        library: records
            .first()
            .map(|record| record.library.clone())
            .unwrap_or_default(),
        ..UmiMetrics::default()
    };
    let mut observed_bases = 0u64;
    let mut with_n = 0u64;
    let mut without_n = 0u64;
    let mut lengths: Vec<f64> = Vec::new();

    for set in sets_by_position(records, &options.base) {
        // The counts are over the position, which is what `umiCounts` holds.
        let mut counts: BTreeMap<String, u64> = BTreeMap::new();
        for index in &set {
            let umi = records[*index].barcode.clone().unwrap_or_default();
            *counts.entry(umi).or_insert(0) += 1;
        }
        let umis: Vec<String> = counts.keys().cloned().collect();
        let clusters = cluster(&umis, options.max_edit_distance_to_join);
        let mut subsets: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for (position, umi) in umis.iter().enumerate() {
            subsets
                .entry(clusters[position])
                .or_default()
                .push(umi.clone());
        }
        metrics.duplicate_sets_ignoring_umi += 1;
        metrics.duplicate_sets_with_umi += subsets.len() as i64;

        for members in subsets.values() {
            let choice = assigned_umi(members, &counts);
            for index in &set {
                let umi = records[*index].barcode.clone().unwrap_or_default();
                if !members.contains(&umi) {
                    continue;
                }
                assigned[*index] = choice.clone();
                if let (Some(tag), Some(choice)) = (&options.molecular_identifier_tag, &choice) {
                    let _ = tag;
                    identifiers[*index] = Some(molecular_identifier(
                        "chr1",
                        records[*index].alignment_start,
                        records[*index].mate_alignment_start,
                        records[*index].reverse_strand(),
                        choice,
                    ));
                }
                if umi.is_empty() {
                    continue;
                }
                if umi.contains('N') {
                    with_n += 1;
                    continue;
                }
                without_n += 1;
                lengths.push(umi.len() as f64);
                observed_bases += umi.len() as u64;
                if let Some(choice) = &choice {
                    metrics.observed_base_errors += umi
                        .bytes()
                        .zip(choice.bytes())
                        .filter(|(a, b)| a != b)
                        .count() as i64;
                    *inferred.entry(choice.clone()).or_insert(0) += 1;
                }
                *observed.entry(umi.clone()).or_insert(0) += 1;
            }
        }
    }

    metrics.mean_umi_length = lengths.first().copied().unwrap_or(0.0);
    metrics.observed_unique_umis = observed.len() as i64;
    metrics.inferred_unique_umis = inferred.len() as i64;
    metrics.observed_umi_entropy = effective_number_of_bases(&observed);
    metrics.inferred_umi_entropy = effective_number_of_bases(&inferred);
    metrics.umi_base_qualities =
        phred_from_error_probability(metrics.observed_base_errors as f64 / observed_bases as f64);
    // Not guarded: a run with no UMIs at all divides zero by zero, and the metrics file writes
    // the NaN that comes out as `?`. Guarding it to zero would report a percentage nobody
    // measured.
    metrics.percent_umi_with_n = with_n as f64 / (with_n + without_n) as f64;

    // The marking itself: the assigned UMI splits the position exactly the way a barcode does,
    // which is the mechanism `MarkDuplicates` already has.
    let mut barcoded: Vec<Record> = records.to_vec();
    for (record, umi) in barcoded.iter_mut().zip(&assigned) {
        record.barcode = umi.clone();
    }
    let mut base = options.base.clone();
    base.barcode_tag = Some(options.umi_tag.clone());
    Ok(UmiMarking {
        marking: mark(&barcoded, &base),
        molecular_identifiers: identifiers,
        umi_metrics: vec![metrics],
    })
}
