//! `UmiAwareMarkDuplicatesWithMateCigar`'s sets: htsjdk's duplicate sets, each cut again by UMI.
//!
//! [`crate::umi_duplicates`] cuts a position's reads by UMI over the sets `MarkDuplicates` cuts.
//! This is the tool as it is built: `SimpleMarkDuplicatesWithMateCigar` driving
//! `UmiAwareDuplicateSetIterator` over htsjdk's `DuplicateSetIterator`, and a good deal of what it
//! answers is in the order things are held in:
//!
//! * a set's records are first sorted and flagged by the tool's own comparator, with the scoring
//!   strategy it was given, and then regrouped by UMI into sets that are sorted and flagged AGAIN
//!   by a default `DuplicateSet`'s comparator, whose own scoring strategy is
//!   `TOTAL_MAPPED_REFERENCE_LENGTH` and which would score by it, except that
//!   `DuplicateScoringStrategy.computeDuplicateScore` stores a record's score on the record the
//!   first time ANY comparator asks, and every record of a set has been asked by the tool's
//!   comparator before the set is regrouped: the second pass sorts by the tool's scores, not its
//!   own. It only ever sets a flag (a set's first record is cleared), so a mate that is not first
//!   keeps what the first pass gave it;
//! * the distinct UMIs of a set are numbered in the order a `HashMap` iterates them, which is
//!   bucket order, and the sets they are joined into come out in the order a `HashMap` keyed by
//!   the union-find root iterates, which for small roots is the roots' own order;
//! * the metrics are kept per library in a `HashMap` filled by `computeIfAbsent`, written in the
//!   order it iterates, and their entropies are summed by `Collectors.summingDouble`, which is a
//!   compensated sum;
//! * the UMI metrics' MEAN_UMI_LENGTH is set from the first UMI without an N the iterator sees
//!   and compared, for every other UMI, against the metrics of ITS library: the second library's
//!   metrics still hold zero, so a file with two libraries throws "UMIs of differing lengths
//!   were found." at the second library's first UMI.
//!
//! Ported from `picard.sam.markduplicates.UmiAwareMarkDuplicatesWithMateCigar`,
//! `UmiAwareDuplicateSetIterator`, `UmiGraph`, `UmiUtil` and `UmiMetrics`, `picard.util.GraphUtils`
//! and htsjdk's `DuplicateSet` and `StringUtil.isWithinHammingDistance` in Picard 3.4.0 and
//! htsjdk 4.2.0.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};

use htsjdk_bam::cigar::Op;

use htsjdk_metrics::file::{MetricBean, Value};

use crate::duplicate_set::{compare, duplicate_sets, library_ids};
use crate::java_hash_map::JavaHashMap;
use crate::mark_duplicates::{Options, Record, ScoringStrategy};

/// `UmiUtil.ALLOWED_UMI`: `^[ATCGNatcgn-]*$`.
fn allowed_umi(umi: &str) -> bool {
    umi.bytes().all(|b| b"ATCGNatcgn-".contains(&b))
}

/// What `getStrand` can say of a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Top,
    Bottom,
    Unknown,
}

/// `SAMUtils.getMateUnclippedStart`, off the `MC` tag.
fn mate_unclipped_start(record: &Record) -> Option<i32> {
    let cigar = record.mate_cigar.as_ref()?;
    let mut start = record.mate_alignment_start;
    for element in &cigar.elements {
        match element.op {
            Op::S | Op::H => start -= element.length as i32,
            _ => break,
        }
    }
    Some(start)
}

/// `SAMUtils.getMateUnclippedEnd`, off the `MC` tag.
fn mate_unclipped_end(record: &Record) -> Option<i32> {
    let cigar = record.mate_cigar.as_ref()?;
    let mut end = record.mate_alignment_start + cigar.reference_length() as i32 - 1;
    for element in cigar.elements.iter().rev() {
        match element.op {
            Op::S | Op::H => end += element.length as i32,
            _ => break,
        }
    }
    Some(end)
}

/// `UmiUtil.getStrand`.
pub fn strand(record: &Record) -> Result<Strand, String> {
    if record.unmapped() {
        return Ok(Strand::Unknown);
    }
    // `getMateUnmappedFlag` asks a record that is not paired, which htsjdk refuses.
    if !record.paired() {
        return Err(
            "java.lang.IllegalStateException: Inappropriate call if not paired read".to_string(),
        );
    }
    if record.mate_unmapped() {
        return Ok(Strand::Unknown);
    }
    if record.reference_index != record.mate_reference_index {
        let first_before = record.reference_index < record.mate_reference_index;
        return Ok(if record.first_of_pair() == first_before {
            Strand::Top
        } else {
            Strand::Bottom
        });
    }
    let read_five_prime = if record.reverse_strand() {
        record.unclipped_end()
    } else {
        record.unclipped_start()
    };
    let missing = || {
        format!(
            "htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found: {}",
            record.name
        )
    };
    let mate_five_prime = if record.mate_reverse_strand() {
        mate_unclipped_end(record).ok_or_else(missing)?
    } else {
        mate_unclipped_start(record).ok_or_else(missing)?
    };
    Ok(
        if record.first_of_pair() == (read_five_prime <= mate_five_prime) {
            Strand::Top
        } else {
            Strand::Bottom
        },
    )
}

/// `UmiUtil.getTopStrandNormalizedUmi`, for a record that has a UMI.
pub fn normalized_umi(record: &Record, umi: &str, duplex: bool) -> Result<String, String> {
    if !allowed_umi(umi) {
        return Err(
            "picard.PicardException: UMI found with illegal characters.  UMIs must match the \
             regular expression ^[ATCGNatcgn-]*$."
                .to_string(),
        );
    }
    if !duplex {
        return Ok(umi.to_string());
    }
    // `String.split("-")` drops trailing empty strings.
    let mut parts: Vec<&str> = umi.split('-').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    if parts.len() != 2 {
        return Err(format!(
            "picard.PicardException: Duplex UMIs must be of the form X-Y where X and Y are equal \
             length UMIs, for example AT-GA.  Found UMI, {umi}"
        ));
    }
    Ok(match strand(record)? {
        Strand::Bottom => format!("{}-{}", parts[1], parts[0]),
        _ => umi.to_string(),
    })
}

/// `UmiUtil.getUmiLength`: the length with the dashes taken out.
fn umi_length(umi: &str) -> usize {
    umi.len() - umi.matches('-').count()
}

/// `StringUtil.isWithinHammingDistance`, which refuses strings of different lengths.
fn within_hamming_distance(a: &str, b: &str, max: i32) -> Result<bool, String> {
    if a.len() != b.len() {
        return Err(
            "java.lang.IllegalArgumentException: Attempted to determine if two strings of \
             different length were within a specified edit distance."
                .to_string(),
        );
    }
    let mut measured = 0;
    for (x, y) in a.bytes().zip(b.bytes()) {
        if x != y {
            measured += 1;
            if measured > max {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// `StringUtil.hammingDistance`.
fn hamming_distance(a: &str, b: &str) -> Result<i64, String> {
    if a.len() != b.len() {
        return Err(format!(
            "java.lang.IllegalArgumentException: Attempted to determine Hamming distance of \
             strings with differing lengths. The first string has length {} and the second string \
             has length {}.",
            a.len(),
            b.len()
        ));
    }
    Ok(a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count() as i64)
}

/// The order a `java.util.HashMap<Integer, _>` iterates keys that were `put` in this order:
/// buckets by `h ^ (h >>> 16)` masked to the table, a bucket in insertion order, the table growing
/// when the size passes three quarters of it.
fn integer_hash_order(first_seen: &[usize]) -> Vec<usize> {
    let mut table: Vec<Vec<usize>> = vec![Vec::new(); 16];
    let mut size = 0;
    for &key in first_seen {
        let h = key as u32;
        let spread = (h ^ (h >> 16)) as usize;
        let index = spread & (table.len() - 1);
        table[index].push(key);
        size += 1;
        if size > table.len() * 3 / 4 {
            let old = table.len();
            let mut grown: Vec<Vec<usize>> = vec![Vec::new(); old * 2];
            for (j, bucket) in table.into_iter().enumerate() {
                for k in bucket {
                    let h = k as u32;
                    let high = ((h ^ (h >> 16)) as usize) & old != 0;
                    grown[if high { j + old } else { j }].push(k);
                }
            }
            table = grown;
        }
    }
    table.into_iter().flatten().collect()
}

/// `Collectors.summingDouble`, JDK 17: Kahan's compensated sum, and the final subtraction.
fn summing_double(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut compensation, mut simple) = (0.0f64, 0.0f64, 0.0f64);
    for value in values {
        let tmp = value - compensation;
        let velvel = sum + tmp;
        compensation = (velvel - sum) - tmp;
        sum = velvel;
        simple += value;
    }
    let tmp = sum - compensation;
    if tmp.is_nan() && simple.is_infinite() {
        simple
    } else {
        tmp
    }
}

/// `UmiMetrics.effectiveNumberOfBases`: Shannon entropy over a `Histogram`'s bins, in base four.
fn effective_number_of_bases(histogram: &BTreeMap<String, f64>) -> f64 {
    let total: f64 = histogram.values().sum();
    let entropy = summing_double(histogram.values().map(|value| {
        let p = value / total;
        -p * jmath::math::log(p)
    }));
    entropy / jmath::math::log(4.0)
}

/// `UmiMetrics`.
#[derive(Debug, Clone, Default)]
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
    pub umi_base_qualities: f64,
    pub pct_umi_with_n: f64,
    observed_umis: BTreeMap<String, f64>,
    inferred_umis: BTreeMap<String, f64>,
    observed_umi_bases: i64,
    observed_umi_with_ns: i64,
    total_observed_umis_without_ns: i64,
}

impl UmiMetrics {
    fn new(library: &str) -> Self {
        UmiMetrics {
            library: library.to_string(),
            ..Default::default()
        }
    }

    /// `calculateDerivedFields`.
    fn calculate_derived_fields(&mut self) {
        self.observed_unique_umis = self.observed_umis.len() as i64;
        self.inferred_unique_umis = self.inferred_umis.len() as i64;
        self.pct_umi_with_n = self.observed_umi_with_ns as f64
            / (self.observed_umi_with_ns as f64 + self.total_observed_umis_without_ns as f64);
        self.observed_umi_entropy = effective_number_of_bases(&self.observed_umis);
        self.inferred_umi_entropy = effective_number_of_bases(&self.inferred_umis);
        self.umi_base_qualities = f64::from(
            htsjdk_bam::quality_util::phred_score_from_error_probability(
                self.observed_base_errors as f64 / self.observed_umi_bases as f64,
            ),
        );
    }

    fn add_umi_observation(&mut self, observed: &str, inferred: &str) {
        *self
            .observed_umis
            .entry(observed.to_string())
            .or_insert(0.0) += 1.0;
        *self
            .inferred_umis
            .entry(inferred.to_string())
            .or_insert(0.0) += 1.0;
        self.observed_umi_bases += observed.len() as i64;
        self.total_observed_umis_without_ns += 1;
    }
}

/// What a run decided.
#[derive(Debug, Clone)]
pub struct UmiAwareResult {
    /// The duplicate flag each record ends with.
    pub duplicate: Vec<bool>,
    /// The molecular identifier each record was given, where the tag was named.
    pub molecular_identifier: Vec<Option<String>>,
    /// Records whose (empty) UMI tag was taken off again, under `ALLOW_MISSING_UMIS`.
    pub umi_removed: Vec<bool>,
    /// The sets the iterator yielded, in order, each in the order its records are returned.
    pub sets: Vec<Vec<usize>>,
    /// One row per library, in the order the reference's map iterates them.
    pub metrics: Vec<UmiMetrics>,
}

/// The tool's own arguments.
#[derive(Debug, Clone)]
pub struct UmiAwareOptions {
    pub base: Options,
    pub max_edit_distance_to_join: i32,
    pub molecular_identifier_tag: bool,
    pub allow_missing_umis: bool,
    pub duplex_umi: bool,
}

/// `DuplicateSet.getRecords()` on a set that was filled by `add`: sorted by the comparator, the
/// flags reset, and the representative checked to be first.
fn sort_and_flag(
    set: &[usize],
    representative: usize,
    records: &[Record],
    library_of: &[i32],
    scoring: ScoringStrategy,
    duplicate: &mut [bool],
) -> Result<Vec<usize>, String> {
    let mut sorted = set.to_vec();
    if sorted.len() > 1 {
        sorted.sort_by(|a, b| {
            compare(
                &records[*a],
                &records[*b],
                library_of[*a],
                library_of[*b],
                scoring,
            )
        });
    }
    let name = &records[representative].name;
    for &index in &sorted {
        let record = &records[index];
        if !record.unmapped() && !record.secondary_or_supplementary() && record.name != *name {
            duplicate[index] = true;
        }
    }
    duplicate[sorted[0]] = false;
    if sorted[0] != representative {
        return Err(format!(
            "htsjdk.samtools.SAMException: BUG: the representative was not the first record after \
             sorting.\nFIRST: {}\nSECOND: {}",
            records[sorted[0]].name, records[representative].name
        ));
    }
    Ok(sorted)
}

/// The whole of it. `umis[i]` is record `i`'s UMI tag; `contigs` names the reference sequences.
pub fn run(
    records: &[Record],
    umis: &[Option<String>],
    contigs: &[String],
    options: &UmiAwareOptions,
) -> Result<UmiAwareResult, String> {
    let library_of = library_ids(records);
    let mut duplicate: Vec<bool> = records.iter().map(|r| r.flags & 0x400 != 0).collect();
    let mut molecular_identifier: Vec<Option<String>> = vec![None; records.len()];
    let mut umi_removed = vec![false; records.len()];
    let mut inferred: Vec<Option<String>> = vec![None; records.len()];
    let mut out_sets: Vec<Vec<usize>> = Vec::new();

    // The metrics map, filled by `computeIfAbsent`.
    let mut metrics_order: JavaHashMap<usize> = JavaHashMap::new();
    let mut metrics_index: HashMap<String, usize> = HashMap::new();
    let mut metrics: Vec<UmiMetrics> = Vec::new();
    let mut have_we_seen_first_read = false;

    for original in duplicate_sets(records, &options.base) {
        // `new DuplicateSet(comparator)`: the representative is the first of the sorted set,
        // which `duplicate_sets` already returns, and the flags are set here.
        let representative = original[0];
        let sorted = sort_and_flag(
            &original,
            representative,
            records,
            &library_of,
            options.base.scoring,
            &mut duplicate,
        )?;

        // UmiGraph's constructor: every record needs a UMI.
        let mut umi_of: Vec<String> = Vec::with_capacity(sorted.len());
        for &index in &sorted {
            match &umis[index] {
                Some(umi) => umi_of.push(umi.clone()),
                None if options.allow_missing_umis => umi_of.push(String::new()),
                None => {
                    return Err(format!(
                        "picard.PicardException: Read {} does not contain a UMI with the RX \
                         attribute.",
                        records[index].name
                    ))
                }
            }
        }
        let mut normalized: Vec<String> = Vec::with_capacity(sorted.len());
        for (position, &index) in sorted.iter().enumerate() {
            normalized.push(normalized_umi(
                &records[index],
                &umi_of[position],
                options.duplex_umi,
            )?);
        }

        // `groupingBy`: a HashMap filled by `computeIfAbsent`, whose key order numbers the UMIs.
        let mut order: JavaHashMap<usize> = JavaHashMap::new();
        let mut distinct: Vec<String> = Vec::new();
        let mut counts: HashMap<String, i64> = HashMap::new();
        for umi in &normalized {
            if !counts.contains_key(umi) {
                distinct.push(umi.clone());
                order.insert_front_if_absent(umi, distinct.len() - 1);
            }
            *counts.entry(umi.clone()).or_insert(0) += 1;
        }
        let umi_list: Vec<String> = order.iter().map(|(_, &at)| distinct[at].clone()).collect();
        let n = umi_list.len();

        // `joinUmisIntoDuplicateSets`: edges between UMIs within the distance, then union-find.
        let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for j in (i + 1)..n {
                if within_hamming_distance(
                    &umi_list[i],
                    &umi_list[j],
                    options.max_edit_distance_to_join,
                )? {
                    neighbors[i].push(j);
                    neighbors[j].push(i);
                }
            }
        }
        let mut cluster: Vec<usize> = (0..n).collect();
        fn find(grouping: &mut [usize], node: usize) -> usize {
            let mut representative = node;
            while representative != grouping[representative] {
                representative = grouping[representative];
            }
            let mut node = node;
            while node != representative {
                let next = grouping[node];
                grouping[node] = representative;
                node = next;
            }
            representative
        }
        for (i, adjacent) in neighbors.iter().enumerate() {
            for &j in adjacent {
                let rep_j = find(&mut cluster, j);
                let rep_i = find(&mut cluster, i);
                if rep_j != rep_i {
                    cluster[rep_j] = rep_i;
                }
            }
        }
        let mut set_id_of_umi: HashMap<String, usize> = HashMap::new();
        for (i, umi) in umi_list.iter().enumerate() {
            let root = find(&mut cluster, i);
            set_id_of_umi.insert(umi.clone(), root);
        }

        // Records grouped by set id, in the order the HashMap iterates the ids.
        let mut by_id: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut first_seen: Vec<usize> = Vec::new();
        for (position, _) in sorted.iter().enumerate() {
            let id = set_id_of_umi[&normalized[position]];
            if !by_id.contains_key(&id) {
                first_seen.push(id);
            }
            by_id.entry(id).or_default().push(position);
        }
        let iteration = integer_hash_order(&first_seen);

        let library = records[representative].library.clone();
        if !metrics_index.contains_key(&library) {
            metrics.push(UmiMetrics::new(&library));
            metrics_index.insert(library.clone(), metrics.len() - 1);
            metrics_order.insert_front_if_absent(&library, metrics.len() - 1);
        }
        let metric_at = metrics_index[&library];

        let mut subsets: Vec<Vec<usize>> = Vec::new();
        for id in &iteration {
            let positions = &by_id[id];
            // The assigned UMI: the most frequent without an N, else the one with the fewest Ns.
            let mut max_count = 0i64;
            let mut assigned: Option<String> = None;
            let mut fewest_n: Option<String> = None;
            let mut n_count = 0usize;
            for &position in positions {
                let umi = &normalized[position];
                if umi.contains('N') {
                    let count = umi.matches('N').count();
                    if n_count == 0 || count < n_count {
                        n_count = count;
                        fewest_n = Some(umi.clone());
                    }
                } else if counts[umi] > max_count {
                    max_count = counts[umi];
                    assigned = Some(umi.clone());
                }
            }
            let assigned = assigned.or(fewest_n).unwrap_or_default();

            let mut members: Vec<usize> = Vec::with_capacity(positions.len());
            for &position in positions {
                let index = sorted[position];
                if options.allow_missing_umis && umi_of[position].is_empty() {
                    umi_removed[index] = true;
                } else if options.molecular_identifier_tag {
                    // `getContig()` of an unmapped read is null, which concatenates as "null".
                    let contig = if records[index].reference_index < 0 {
                        "null"
                    } else {
                        contigs
                            .get(records[index].reference_index as usize)
                            .map(String::as_str)
                            .unwrap_or("null")
                    };
                    let record = &records[index];
                    let mut identifier = format!(
                        "{contig}:{}/{assigned}",
                        if record.reverse_strand() {
                            record.alignment_start
                        } else {
                            record.mate_alignment_start
                        }
                    );
                    if options.duplex_umi {
                        match strand(record)? {
                            Strand::Top => identifier.push_str("/A"),
                            Strand::Bottom => identifier.push_str("/B"),
                            Strand::Unknown => {}
                        }
                    }
                    molecular_identifier[index] = Some(identifier);
                }
                inferred[index] = Some(assigned.clone());
                members.push(index);
            }

            // `new DuplicateSet()`: records added one by one, the representative moving to any
            // record the DEFAULT comparator (total mapped reference length) puts first.
            let mut representative_of_subset = members[0];
            for &index in &members[1..] {
                if compare(
                    &records[representative_of_subset],
                    &records[index],
                    library_of[representative_of_subset],
                    library_of[index],
                    options.base.scoring,
                ) == Ordering::Greater
                {
                    representative_of_subset = index;
                }
            }
            let subset = sort_and_flag(
                &members,
                representative_of_subset,
                records,
                &library_of,
                options.base.scoring,
                &mut duplicate,
            )?;
            subsets.push(subset);
        }

        // The statistics, over each subset in turn.
        for subset in &subsets {
            for &index in subset {
                let position = sorted.iter().position(|s| *s == index).unwrap_or(0);
                if options.allow_missing_umis && umi_of[position].is_empty() {
                    continue; // the tag was taken off, so the UMI reads as null
                }
                let current = &normalized[position];
                let metric = &mut metrics[metric_at];
                if current.contains('N') {
                    metric.observed_umi_with_ns += 1;
                } else {
                    let length = umi_length(current) as f64;
                    if !have_we_seen_first_read {
                        metric.mean_umi_length = length;
                        have_we_seen_first_read = true;
                    } else if metric.mean_umi_length != length {
                        return Err(
                            "picard.PicardException: UMIs of differing lengths were found."
                                .to_string(),
                        );
                    }
                    let inferred_umi = inferred[index].clone().unwrap_or_default();
                    metric.observed_base_errors += hamming_distance(current, &inferred_umi)?;
                    metric.add_umi_observation(current, &inferred_umi);
                }
            }
        }
        let metric = &mut metrics[metric_at];
        metric.duplicate_sets_with_umi += subsets.len() as i64;
        metric.duplicate_sets_ignoring_umi += 1;
        out_sets.extend(subsets);
    }

    for metric in &mut metrics {
        metric.calculate_derived_fields();
    }
    let ordered: Vec<UmiMetrics> = metrics_order
        .iter()
        .map(|(_, &at)| metrics[at].clone())
        .collect();
    Ok(UmiAwareResult {
        duplicate,
        molecular_identifier,
        umi_removed,
        sets: out_sets,
        metrics: ordered,
    })
}

impl MetricBean for UmiMetrics {
    fn class_name(&self) -> &str {
        "picard.sam.markduplicates.UmiMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "LIBRARY",
            "MEAN_UMI_LENGTH",
            "OBSERVED_UNIQUE_UMIS",
            "INFERRED_UNIQUE_UMIS",
            "OBSERVED_BASE_ERRORS",
            "DUPLICATE_SETS_IGNORING_UMI",
            "DUPLICATE_SETS_WITH_UMI",
            "OBSERVED_UMI_ENTROPY",
            "INFERRED_UMI_ENTROPY",
            "UMI_BASE_QUALITIES",
            "PCT_UMI_WITH_N",
        ]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.library.clone()),
            Value::Double(self.mean_umi_length),
            Value::Long(self.observed_unique_umis),
            Value::Long(self.inferred_unique_umis),
            Value::Long(self.observed_base_errors),
            Value::Long(self.duplicate_sets_ignoring_umi),
            Value::Long(self.duplicate_sets_with_umi),
            Value::Double(self.observed_umi_entropy),
            Value::Double(self.inferred_umi_entropy),
            Value::Double(self.umi_base_qualities),
            Value::Double(self.pct_umi_with_n),
        ]
    }
}
