//! `MeanQualityByCycle` and `CollectBaseDistributionByCycle`.
//!
//! Ported from `picard.analysis.MeanQualityByCycle` with `picard.analysis.HistogramGenerator`,
//! and `picard.analysis.CollectBaseDistributionByCycle` with its private inner
//! `HistogramGenerator`, tag 3.4.0.
//!
//! The two are stratum-mates in the calibration gate, and they are ported **in one module on
//! purpose**: they share a cycle-indexing convention that is stated in neither of them, and
//! writing it once is exactly the amortisation the gate is trying to measure.
//!
//! ## The shared convention
//!
//! A cycle is **1-based** and, for a reverse-strand read, **counted from the other end**:
//!
//! ```text
//! cycle = if reverse { length - i } else { i + 1 }
//! ```
//!
//! Note the asymmetry: the forward form is `i + 1` and the reverse form is `length - i`, so a
//! reverse-strand read's first base is cycle `length` and its last is cycle 1. The arrays are
//! therefore sized `length + 1` and index 0 is never written, which is why every output loop
//! starts at 0 and relies on the count being zero there rather than skipping it explicitly.
//!
//! The second end's cycles are offset by `firstReadLength`, which is **the last cycle that had
//! any first-end data**, not the read length. On uniform input the two coincide; on a file of
//! mixed lengths they do not.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricBean, Value};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const READ_REVERSE: u16 = 0x10;
const SECOND_OF_PAIR: u16 = 0x80;
const VENDOR_FAILED: u16 = 0x200;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// Shared filter for both tools, from their `acceptRead`.
fn skip(rec: &BamRecord, pf_reads_only: bool, aligned_reads_only: bool) -> bool {
    (pf_reads_only && rec.flags & VENDOR_FAILED != 0)
        || (aligned_reads_only && rec.flags & READ_UNMAPPED != 0)
        || rec.flags & (SECONDARY | SUPPLEMENTARY) != 0
}

/// `cycle = rc ? length - i : i + 1`, the convention both tools index by.
fn cycle_of(reverse: bool, length: usize, i: usize) -> usize {
    if reverse {
        length - i
    } else {
        i + 1
    }
}

fn is_second_of_pair(rec: &BamRecord) -> bool {
    rec.flags & READ_PAIRED != 0 && rec.flags & SECOND_OF_PAIR != 0
}

/// `SAMUtils.fastqToPhred`, for the `OQ` tag.
fn fastq_to_phred(s: &str) -> Vec<i8> {
    s.chars()
        .map(|c| (c as u32 as u8).wrapping_sub(33) as i8)
        .collect()
}

// ---------------------------------------------------------------------------------------
// MeanQualityByCycle
// ---------------------------------------------------------------------------------------

/// `picard.analysis.HistogramGenerator`, reduced to what `MeanQualityByCycle` reads from it.
///
/// The class also accumulates error probabilities, which no path in this tool consults; they
/// exist for a flow-space method on the same class. Not porting them is a deliberate omission
/// rather than an oversight, and it is recorded here so a later reader does not have to
/// rediscover that the field is unused.
#[derive(Debug, Default)]
struct QualityByCycle {
    use_original_qualities: bool,
    first_totals: Vec<f64>,
    first_counts: Vec<i64>,
    second_totals: Vec<f64>,
    second_counts: Vec<i64>,
}

impl QualityByCycle {
    fn new(use_original_qualities: bool) -> Self {
        QualityByCycle {
            use_original_qualities,
            ..Default::default()
        }
    }

    fn ensure(&mut self, length: usize) {
        if length > self.first_totals.len() {
            self.first_totals.resize(length, 0.0);
            self.first_counts.resize(length, 0);
            self.second_totals.resize(length, 0.0);
            self.second_counts.resize(length, 0);
        }
    }

    fn add(&mut self, rec: &BamRecord) {
        // `getOriginalBaseQualities()` returns null when OQ is absent, and `addRecord` returns
        // immediately rather than falling back. That is the opposite of what
        // CollectQualityYieldMetrics does with the same tag.
        let quals: Vec<i8> = if self.use_original_qualities {
            match rec.tags.get(Tag::new(b"OQ")) {
                Some(TagValue::Str(s)) if !s.is_empty() => fastq_to_phred(s),
                _ => return,
            }
        } else {
            rec.base_qualities.iter().map(|&b| b as i8).collect()
        };

        let length = quals.len();
        let reverse = rec.flags & READ_REVERSE != 0;
        self.ensure(length + 1);
        let second = is_second_of_pair(rec);
        for (i, &q) in quals.iter().enumerate() {
            let cycle = cycle_of(reverse, length, i);
            if second {
                self.second_totals[cycle] += q as f64;
                self.second_counts[cycle] += 1;
            } else {
                self.first_totals[cycle] += q as f64;
                self.first_counts[cycle] += 1;
            }
        }
    }

    fn is_empty(&self) -> bool {
        self.first_totals.is_empty()
    }

    /// `getMeanQualityHistogram`.
    ///
    /// The first-end loop tests `firstReadTotalsByCycle[cycle] > 0`, on the **total** rather
    /// than the count, while the second-end loop tests the count. A cycle whose qualities are
    /// all zero therefore appears for a second end and not for a first.
    fn mean_quality_histogram(&self) -> OutHistogram {
        let label = if self.use_original_qualities {
            "MEAN_ORIGINAL_QUALITY"
        } else {
            "MEAN_QUALITY"
        };
        let mut bins: Vec<(String, f64)> = Vec::new();
        let mut first_read_length = 0usize;
        for cycle in 0..self.first_totals.len() {
            if self.first_totals[cycle] > 0.0 {
                bins.push((
                    cycle.to_string(),
                    self.first_totals[cycle] / self.first_counts[cycle] as f64,
                ));
                first_read_length = cycle;
            }
        }
        for i in 0..self.second_totals.len() {
            if self.second_counts[i] > 0 {
                let cycle = first_read_length + i;
                bins.push((
                    cycle.to_string(),
                    self.second_totals[i] / self.second_counts[i] as f64,
                ));
            }
        }
        // The Java accumulates into a Histogram keyed by cycle, so the bins come out in key
        // order and a repeated cycle would be summed, not duplicated.
        bins.sort_by_key(|(k, _)| k.parse::<i64>().unwrap_or(0));
        OutHistogram {
            bin_label: "CYCLE".to_string(),
            value_label: label.to_string(),
            key_class: "java.lang.Integer".to_string(),
            bins,
        }
    }
}

/// `MeanQualityByCycle`.
#[derive(Debug)]
pub struct MeanQualityByCycle {
    pub pf_reads_only: bool,
    pub aligned_reads_only: bool,
    q: QualityByCycle,
    oq: QualityByCycle,
}

impl Default for MeanQualityByCycle {
    fn default() -> Self {
        MeanQualityByCycle {
            pf_reads_only: false,
            aligned_reads_only: false,
            q: QualityByCycle::new(false),
            oq: QualityByCycle::new(true),
        }
    }
}

impl MeanQualityByCycle {
    pub fn accept(&mut self, rec: &BamRecord) {
        if skip(rec, self.pf_reads_only, self.aligned_reads_only) {
            return;
        }
        self.q.add(rec);
        self.oq.add(rec);
    }

    /// `finish`: the primary histogram always, the `OQ` one only if it saw anything.
    pub fn finish(&self) -> Vec<OutHistogram> {
        let mut out = vec![self.q.mean_quality_histogram()];
        if !self.oq.is_empty() {
            out.push(self.oq.mean_quality_histogram());
        }
        out
    }
}

// ---------------------------------------------------------------------------------------
// CollectBaseDistributionByCycle
// ---------------------------------------------------------------------------------------

/// `BaseDistributionByCycleMetrics`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BaseDistributionByCycleMetrics {
    pub read_end: i32,
    pub cycle: i32,
    pub pct_a: f64,
    pub pct_c: f64,
    pub pct_g: f64,
    pub pct_t: f64,
    pub pct_n: f64,
}

const BASE_DIST_COLUMNS: &[&str] = &[
    "READ_END", "CYCLE", "PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N",
];

impl MetricBean for BaseDistributionByCycleMetrics {
    fn class_name(&self) -> &str {
        // Extends MetricBase, not MultilevelMetrics, so there are no inherited columns.
        "picard.analysis.BaseDistributionByCycleMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        BASE_DIST_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Long(self.read_end as i64),
            Value::Long(self.cycle as i64),
            Value::Double(self.pct_a),
            Value::Double(self.pct_c),
            Value::Double(self.pct_g),
            Value::Double(self.pct_t),
            Value::Double(self.pct_n),
        ]
    }
}

/// `CollectBaseDistributionByCycle`.
#[derive(Debug, Default)]
pub struct CollectBaseDistributionByCycle {
    pub pf_reads_only: bool,
    pub aligned_reads_only: bool,
    max_length: usize,
    /// `[base][cycle]`, base indexed by [`base_to_int`].
    first_totals: [Vec<i64>; 5],
    first_counts: Vec<i64>,
    second_totals: [Vec<i64>; 5],
    second_counts: Vec<i64>,
    seen_second_end: bool,
}

/// The tool's own `baseToInt`: A, C, G, T, and **everything else** as 4.
///
/// Not just `N`: any byte that is not one of the four, including `=` and every IUPAC ambiguity
/// code, silently becomes `PCT_N`. So `PCT_N` means "not one of the four bases", not "N".
///
/// The case folding is **unreachable through BAM input**, and this is worth stating because a
/// sabotage test was written for it and found nothing. BAM stores bases as nibbles and decoding
/// always yields upper case (see htsjdk-rs decision 0008), so a lower-case base cannot survive
/// the file format. The branch matters only for SAM text input, which is not yet ported. The
/// ambiguity codes *do* survive, and sabotaging their mapping diverges all 7 corpus cases.
fn base_to_int(base: u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

impl CollectBaseDistributionByCycle {
    fn ensure(&mut self, length: usize) {
        if length > self.max_length {
            for i in 0..5 {
                self.first_totals[i].resize(length, 0);
                self.second_totals[i].resize(length, 0);
            }
            self.first_counts.resize(length, 0);
            self.second_counts.resize(length, 0);
            self.max_length = length;
        }
    }

    pub fn accept(&mut self, rec: &BamRecord) {
        if skip(rec, self.pf_reads_only, self.aligned_reads_only) {
            return;
        }
        let bases = &rec.read_bases;
        if bases.is_empty() {
            return;
        }
        let length = bases.len();
        let reverse = rec.flags & READ_REVERSE != 0;
        self.ensure(length + 1);
        let second = is_second_of_pair(rec);
        if second {
            self.seen_second_end = true;
        }
        for (i, &b) in bases.iter().enumerate() {
            let cycle = cycle_of(reverse, length, i);
            if second {
                self.second_totals[base_to_int(b)][cycle] += 1;
                self.second_counts[cycle] += 1;
            } else {
                self.first_totals[base_to_int(b)][cycle] += 1;
                self.first_counts[cycle] += 1;
            }
        }
    }

    /// `addToMetricsFile`.
    ///
    /// `firstReadLength` is the last cycle that had first-end data, and the second end's cycles
    /// are offset by it, exactly as in `MeanQualityByCycle`. The two tools reach that offset by
    /// separate code that happens to agree.
    pub fn finish(&self) -> Vec<BaseDistributionByCycleMetrics> {
        let mut out = Vec::new();
        let mut first_read_length = 0usize;
        for i in 0..self.max_length {
            if self.first_counts[i] != 0 {
                out.push(self.metric(1, i as i32, &self.first_totals, self.first_counts[i], i));
                first_read_length = i;
            }
        }
        if self.seen_second_end {
            for i in 0..self.max_length {
                if self.second_counts[i] != 0 {
                    out.push(self.metric(
                        2,
                        (i + first_read_length) as i32,
                        &self.second_totals,
                        self.second_counts[i],
                        i,
                    ));
                }
            }
        }
        out
    }

    fn metric(
        &self,
        read_end: i32,
        cycle: i32,
        totals: &[Vec<i64>; 5],
        count: i64,
        i: usize,
    ) -> BaseDistributionByCycleMetrics {
        let pct = |b: usize| 100.0 * totals[b][i] as f64 / count as f64;
        BaseDistributionByCycleMetrics {
            read_end,
            cycle,
            pct_a: pct(0),
            pct_c: pct(1),
            pct_g: pct(2),
            pct_t: pct(3),
            pct_n: pct(4),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.max_length == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

    fn read(bases: &[u8], quals: Vec<u8>, flags: u16) -> BamRecord {
        BamRecord {
            read_name: "r".into(),
            flags,
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                length: bases.len() as u32,
                op: Op::M,
            }]),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: bases.to_vec(),
            base_qualities: quals,
            tags: Default::default(),
        }
    }

    /// The convention is asymmetric: forward is `i + 1`, reverse is `length - i`. So a reverse
    /// read's first base is the LAST cycle, and its last base is cycle 1.
    #[test]
    fn the_cycle_convention_is_asymmetric() {
        assert_eq!(cycle_of(false, 5, 0), 1);
        assert_eq!(cycle_of(false, 5, 4), 5);
        assert_eq!(cycle_of(true, 5, 0), 5);
        assert_eq!(cycle_of(true, 5, 4), 1);
    }

    #[test]
    fn a_reverse_read_reverses_its_quality_profile() {
        let mut fwd = MeanQualityByCycle::default();
        fwd.accept(&read(b"ACGT", vec![10, 20, 30, 40], 0));
        let mut rev = MeanQualityByCycle::default();
        rev.accept(&read(b"ACGT", vec![10, 20, 30, 40], READ_REVERSE));

        let f: Vec<f64> = fwd.finish()[0].bins.iter().map(|(_, v)| *v).collect();
        let r: Vec<f64> = rev.finish()[0].bins.iter().map(|(_, v)| *v).collect();
        assert_eq!(f, vec![10.0, 20.0, 30.0, 40.0]);
        assert_eq!(r, vec![40.0, 30.0, 20.0, 10.0]);
    }

    /// `PCT_N` is "not one of the four bases", not "N". A `=` lands there too.
    #[test]
    fn pct_n_is_everything_that_is_not_acgt() {
        assert_eq!(base_to_int(b'N'), 4);
        assert_eq!(base_to_int(b'='), 4);
        for code in *b"MRSVWYHKDB" {
            assert_eq!(
                base_to_int(code),
                4,
                "IUPAC code {} is not A/C/G/T",
                code as char
            );
        }
    }

    /// The case folding exists and is unreachable through BAM, because BAM stores bases as
    /// nibbles that decode to upper case. Asserted here so the branch is covered by *something*
    /// while the conformance corpus honestly cannot reach it.
    #[test]
    fn case_folding_exists_but_no_bam_can_exercise_it() {
        assert_eq!(base_to_int(b'a'), base_to_int(b'A'));
        assert_eq!(base_to_int(b't'), base_to_int(b'T'));
        // What a BAM round trip actually yields:
        assert_eq!(
            htsjdk_bam::bases::compressed_bases_to_bytes(
                4,
                &htsjdk_bam::bases::bytes_to_compressed_bases(b"acgt").unwrap(),
                0
            ),
            b"ACGT",
            "lower case cannot survive the file format"
        );
    }

    #[test]
    fn base_percentages_are_out_of_a_hundred() {
        let mut c = CollectBaseDistributionByCycle::default();
        c.accept(&read(b"AACC", vec![30; 4], 0));
        c.accept(&read(b"AAGG", vec![30; 4], 0));
        let m = c.finish();
        assert_eq!(m.len(), 4);
        assert_eq!(m[0].pct_a, 100.0);
        assert_eq!(m[2].pct_c, 50.0);
        assert_eq!(m[2].pct_g, 50.0);
    }

    /// The `OQ` generator **returns without recording** when the tag is absent, rather than
    /// falling back to the primary qualities. That is the opposite of what
    /// `CollectQualityYieldMetrics` does with the same tag, and the two live in the same crate.
    #[test]
    fn the_oq_histogram_is_absent_rather_than_a_copy_of_the_primary_one() {
        let mut c = MeanQualityByCycle::default();
        c.accept(&read(b"ACGT", vec![30; 4], 0));
        assert_eq!(c.finish().len(), 1, "no OQ tag, so no second histogram");

        let mut r = read(b"ACGT", vec![30; 4], 0);
        r.tags.insert(Tag::new(b"OQ"), TagValue::Str("IIII".into()));
        let mut c2 = MeanQualityByCycle::default();
        c2.accept(&r);
        let out = c2.finish();
        assert_eq!(out.len(), 2);
        assert_eq!(out[1].value_label, "MEAN_ORIGINAL_QUALITY");
        assert_eq!(out[1].bins[0].1, 40.0);
    }

    #[test]
    fn secondary_and_supplementary_records_are_skipped_by_both_tools() {
        for flag in [SECONDARY, SUPPLEMENTARY] {
            let mut q = MeanQualityByCycle::default();
            q.accept(&read(b"ACGT", vec![30; 4], flag));
            assert!(q.finish()[0].bins.is_empty());

            let mut b = CollectBaseDistributionByCycle::default();
            b.accept(&read(b"ACGT", vec![30; 4], flag));
            assert!(b.finish().is_empty());
        }
    }

    /// The second end's cycles are offset by the last first-end cycle, not by the read length.
    #[test]
    fn the_second_end_is_offset_past_the_first() {
        let mut c = CollectBaseDistributionByCycle::default();
        c.accept(&read(b"AAAA", vec![30; 4], READ_PAIRED | 0x40));
        c.accept(&read(b"CCCC", vec![30; 4], READ_PAIRED | SECOND_OF_PAIR));
        let m = c.finish();
        let first: Vec<i32> = m
            .iter()
            .filter(|x| x.read_end == 1)
            .map(|x| x.cycle)
            .collect();
        let second: Vec<i32> = m
            .iter()
            .filter(|x| x.read_end == 2)
            .map(|x| x.cycle)
            .collect();
        assert_eq!(first, vec![1, 2, 3, 4]);
        assert_eq!(second, vec![5, 6, 7, 8]);
    }

    /// With no second end at all, the second-end block is skipped entirely rather than
    /// emitting zero rows.
    #[test]
    fn an_unpaired_file_emits_only_read_end_one() {
        let mut c = CollectBaseDistributionByCycle::default();
        c.accept(&read(b"ACGT", vec![30; 4], 0));
        assert!(c.finish().iter().all(|m| m.read_end == 1));
    }
}
