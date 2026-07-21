//! `CollectRnaSeqMetrics`, the per-unit collector.
//!
//! Ported from `picard.analysis.directed.RnaSeqMetricsCollector` (its inner
//! `PerUnitRnaSeqMetricsCollector`), `picard.analysis.RnaSeqMetrics`, and the setup in
//! `picard.analysis.CollectRnaSeqMetrics`, all at tag 3.4.0. Only the **ALL_READS** accumulation
//! level is ported: the corpora drive `METRIC_ACCUMULATION_LEVEL=[ALL_READS]`, so sample, library
//! and read group are null and the coverage prefix is `All_Reads.`. The multi-level fan-out
//! (`SAMRecordMultiLevelCollector`) is a separate symbol and is not claimed here.
//!
//! ## The one ordering subtlety, and what measurement said about it
//!
//! `computeCoverageMetrics` accumulates `normalized / transcriptCount` into the
//! `normalized_coverage` histogram **in floating point**, iterating the picked transcripts in Java
//! `HashMap` order (`RnaSeqMetricsCollector.java:373`). Floating-point addition is not associative,
//! so in principle those printed values depend on that iteration order, and reproducing it would
//! mean replaying `HashMap`'s bucket order (feasible, since `Gene.Transcript.hashCode` is
//! content-based, not identity-based).
//!
//! It was measured instead of assumed. Folding the differentiated-coverage corpus in three
//! different orders (the deterministic content order below, its reverse, and the `HashMap`
//! bucket order) yields the **same** `.rna_metrics` bytes: the last-ULP differences vanish under
//! `FormatUtil`'s formatting. So the fold order is unobservable at printed precision here, the
//! RnaSeq analogue of htsjdk-rs decision 0020, and this port does **not** claim to reproduce Java's
//! `HashMap` order. It folds in a single deterministic content order (so the output is reproducible
//! run to run despite the `OverlapDetector` iterating a Rust `HashMap`), and decision 0005 records
//! the residual risk: a corpus with many high-coverage transcripts could in principle make the fold
//! observable, at which point the exact `HashMap` order would have to be built and verified against
//! that corpus.
//!
//! The four `MEDIAN_*` metrics, by contrast, come from `Histogram.getMedian`, which sorts, so they
//! are order-independent regardless.

use std::collections::HashMap;

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::interval::Interval;
use htsjdk_bam::overlap::OverlapDetector;
use htsjdk_bam::record::BamRecord;
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricBean, Value};
use htsjdk_metrics::histogram::Histogram;

use crate::annotation::{Gene, LocusFunction};
use crate::insert_size::{pair_orientation, PairOrientation};

// SAM flag bits used here (htsjdk SAMFlag), named locally rather than pulling a dependency.
const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const FIRST_OF_PAIR: u16 = 0x40;
const NOT_PRIMARY_ALIGNMENT: u16 = 0x100;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

fn has(flags: u16, bit: u16) -> bool {
    flags & bit != 0
}

/// `RnaSeqMetricsCollector.StrandSpecificity`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandSpecificity {
    None,
    FirstReadTranscriptionStrand,
    SecondReadTranscriptionStrand,
}

/// The default of `RnaSeqMetricsCollector.defaultEndBiasBases`.
pub const DEFAULT_END_BIAS_BASES: i32 = 100;
/// `CollectRnaSeqMetrics.MINIMUM_LENGTH` default.
pub const DEFAULT_MINIMUM_LENGTH: i32 = 500;
/// `CollectRnaSeqMetrics.RRNA_FRAGMENT_PERCENTAGE` default.
pub const DEFAULT_RRNA_FRAGMENT_PERCENTAGE: f64 = 0.8;

/// `picard.analysis.RnaSeqMetrics`, fields in declaration order (which is column order).
///
/// `RIBOSOMAL_BASES` and `PCT_RIBOSOMAL_BASES` are the boxed `Long`/`Double` that are left null
/// when no ribosomal interval list is supplied; every other numeric field is a primitive that
/// starts at zero.
#[derive(Debug, Clone, Default)]
pub struct RnaSeqMetrics {
    pub pf_bases: i64,
    pub pf_aligned_bases: i64,
    pub ribosomal_bases: Option<i64>,
    pub coding_bases: i64,
    pub utr_bases: i64,
    pub intronic_bases: i64,
    pub intergenic_bases: i64,
    pub ignored_reads: i64,
    pub correct_strand_reads: i64,
    pub incorrect_strand_reads: i64,
    pub num_r1_transcript_strand_reads: i64,
    pub num_r2_transcript_strand_reads: i64,
    pub num_unexplained_reads: i64,
    pub pct_r1_transcript_strand_reads: f64,
    pub pct_r2_transcript_strand_reads: f64,
    pub pct_ribosomal_bases: Option<f64>,
    pub pct_coding_bases: f64,
    pub pct_utr_bases: f64,
    pub pct_intronic_bases: f64,
    pub pct_intergenic_bases: f64,
    pub pct_mrna_bases: f64,
    pub pct_usable_bases: f64,
    pub pct_correct_strand_reads: f64,
    pub median_cv_coverage: f64,
    pub median_5prime_bias: f64,
    pub median_3prime_bias: f64,
    pub median_5prime_to_3prime_bias: f64,
}

const COLUMNS: &[&str] = &[
    "PF_BASES",
    "PF_ALIGNED_BASES",
    "RIBOSOMAL_BASES",
    "CODING_BASES",
    "UTR_BASES",
    "INTRONIC_BASES",
    "INTERGENIC_BASES",
    "IGNORED_READS",
    "CORRECT_STRAND_READS",
    "INCORRECT_STRAND_READS",
    "NUM_R1_TRANSCRIPT_STRAND_READS",
    "NUM_R2_TRANSCRIPT_STRAND_READS",
    "NUM_UNEXPLAINED_READS",
    "PCT_R1_TRANSCRIPT_STRAND_READS",
    "PCT_R2_TRANSCRIPT_STRAND_READS",
    "PCT_RIBOSOMAL_BASES",
    "PCT_CODING_BASES",
    "PCT_UTR_BASES",
    "PCT_INTRONIC_BASES",
    "PCT_INTERGENIC_BASES",
    "PCT_MRNA_BASES",
    "PCT_USABLE_BASES",
    "PCT_CORRECT_STRAND_READS",
    "MEDIAN_CV_COVERAGE",
    "MEDIAN_5PRIME_BIAS",
    "MEDIAN_3PRIME_BIAS",
    "MEDIAN_5PRIME_TO_3PRIME_BIAS",
    // From MultilevelMetrics; HotSpot's getFields() puts the inherited public fields last.
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

impl MetricBean for RnaSeqMetrics {
    fn class_name(&self) -> &str {
        "picard.analysis.RnaSeqMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let long = |v: i64| Value::Long(v);
        let dbl = |v: f64| Value::Double(v);
        vec![
            long(self.pf_bases),
            long(self.pf_aligned_bases),
            self.ribosomal_bases.map(Value::Long).unwrap_or(Value::Null),
            long(self.coding_bases),
            long(self.utr_bases),
            long(self.intronic_bases),
            long(self.intergenic_bases),
            long(self.ignored_reads),
            long(self.correct_strand_reads),
            long(self.incorrect_strand_reads),
            long(self.num_r1_transcript_strand_reads),
            long(self.num_r2_transcript_strand_reads),
            long(self.num_unexplained_reads),
            dbl(self.pct_r1_transcript_strand_reads),
            dbl(self.pct_r2_transcript_strand_reads),
            self.pct_ribosomal_bases
                .map(Value::Double)
                .unwrap_or(Value::Null),
            dbl(self.pct_coding_bases),
            dbl(self.pct_utr_bases),
            dbl(self.pct_intronic_bases),
            dbl(self.pct_intergenic_bases),
            dbl(self.pct_mrna_bases),
            dbl(self.pct_usable_bases),
            dbl(self.pct_correct_strand_reads),
            dbl(self.median_cv_coverage),
            dbl(self.median_5prime_bias),
            dbl(self.median_3prime_bias),
            dbl(self.median_5prime_to_3prime_bias),
            // SAMPLE, LIBRARY, READ_GROUP are null at the ALL_READS level.
            Value::Null,
            Value::Null,
            Value::Null,
        ]
    }
}

/// The per-unit collector for the ALL_READS level.
pub struct RnaSeqMetricsCollector<'a> {
    seq_names: &'a [String],
    gene_overlap: &'a OverlapDetector<Gene>,
    ribosomal_overlap: &'a OverlapDetector<Interval>,
    ignored_sequence_indices: &'a [i32],
    minimum_length: i32,
    strand_specificity: StrandSpecificity,
    rrna_fragment_percentage: f64,
    end_bias_bases: i32,

    metrics: RnaSeqMetrics,
    /// Coverage arrays keyed by the identity of the `Transcript` they belong to, expressed as
    /// `(gene pointer, transcript index)` since the genes live in `gene_overlap` for the whole run.
    coverage_by_transcript: HashMap<(usize, usize), Vec<i32>>,
}

impl<'a> RnaSeqMetricsCollector<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        seq_names: &'a [String],
        gene_overlap: &'a OverlapDetector<Gene>,
        ribosomal_overlap: &'a OverlapDetector<Interval>,
        ribosomal_initial_value: Option<i64>,
        ignored_sequence_indices: &'a [i32],
        minimum_length: i32,
        strand_specificity: StrandSpecificity,
        rrna_fragment_percentage: f64,
        end_bias_bases: i32,
    ) -> Self {
        RnaSeqMetricsCollector {
            seq_names,
            gene_overlap,
            ribosomal_overlap,
            ignored_sequence_indices,
            minimum_length,
            strand_specificity,
            rrna_fragment_percentage,
            end_bias_bases,
            metrics: RnaSeqMetrics {
                ribosomal_bases: ribosomal_initial_value,
                ..Default::default()
            },
            coverage_by_transcript: HashMap::new(),
        }
    }

    fn reference_name(&self, index: i32) -> &str {
        &self.seq_names[index as usize]
    }

    /// `getNumAlignedBases`: the summed length of the alignment blocks.
    fn num_aligned_bases(rec: &BamRecord) -> i64 {
        alignment_blocks(&rec.cigar, rec.alignment_start)
            .iter()
            .map(|b| b.length as i64)
            .sum()
    }

    /// `PerUnitRnaSeqMetricsCollector.acceptRecord`.
    pub fn accept(&mut self, rec: &BamRecord) {
        let flags = rec.flags;

        if has(flags, READ_FAILS_VENDOR_QUALITY) {
            return;
        }

        // PF bases counts unmapped reads, but not non-primary alignments.
        if !has(flags, NOT_PRIMARY_ALIGNMENT) {
            self.metrics.pf_bases += rec.read_length() as i64;
        }

        // A primary, mapped read on an ignored sequence is counted and dropped.
        if !has(flags, READ_UNMAPPED)
            && !has(flags, NOT_PRIMARY_ALIGNMENT)
            && self.ignored_sequence_indices.contains(&rec.reference_index)
        {
            self.metrics.ignored_reads += 1;
            return;
        }

        // From here, secondary and unmapped reads are done with.
        if has(flags, NOT_PRIMARY_ALIGNMENT) || has(flags, READ_UNMAPPED) {
            return;
        }

        let contig = self.reference_name(rec.reference_index).to_string();
        let read_start = rec.alignment_start;
        let read_end = rec.alignment_end();

        // The fragment interval: the whole template for a properly-paired read, else the read.
        let fragment: Option<(i32, i32)> = if !has(flags, READ_PAIRED) {
            Some((read_start, read_end))
        } else if has(flags, MATE_UNMAPPED) || rec.reference_index != rec.mate_reference_index {
            None
        } else {
            let fragment_start = read_start.min(rec.mate_alignment_start);
            // CoordMath.getEnd(start, length) = start + length - 1.
            let fragment_end = fragment_start + rec.inferred_insert_size.abs() - 1;
            Some((fragment_start, fragment_end))
        };

        if let Some((frag_start, frag_end)) = fragment {
            let frag_len = frag_end - frag_start + 1; // Interval.length()
            let mut intersection_length = 0;
            for interval in self
                .ribosomal_overlap
                .get_overlaps(&contig, frag_start, frag_end)
            {
                // getIntersectionLength on same-contig overlapping intervals.
                let this_len = interval.end.min(frag_end) - interval.start.max(frag_start) + 1;
                intersection_length = intersection_length.max(this_len);
            }
            if intersection_length as f64 / frag_len as f64 >= self.rrna_fragment_percentage {
                let aligned = Self::num_aligned_bases(rec);
                // RIBOSOMAL_BASES is non-null whenever a ribosomal list is present, which is the
                // only way this branch is reached.
                self.metrics.ribosomal_bases =
                    Some(self.metrics.ribosomal_bases.unwrap_or(0) + aligned);
                self.metrics.pf_aligned_bases += aligned;
                return;
            }
        }

        let overlapping_genes: Vec<&Gene> = self
            .gene_overlap
            .get_overlaps(&contig, read_start, read_end);
        let blocks = alignment_blocks(&rec.cigar, rec.alignment_start);
        let mut overlaps_exon = false;

        for block in &blocks {
            let len = block.length as usize;
            let mut locus_functions = vec![LocusFunction::Intergenic; len];

            for gene in &overlapping_genes {
                let gene_ptr = *gene as *const Gene as usize;
                for (tx_idx, transcript) in gene.transcripts.iter().enumerate() {
                    transcript.assign_locus_function_for_range(
                        block.reference_start,
                        &mut locus_functions,
                    );
                    // collectCoverageStatistics is always true for CollectRnaSeqMetrics.
                    let coverage = self
                        .coverage_by_transcript
                        .entry((gene_ptr, tx_idx))
                        .or_insert_with(|| vec![0; transcript.length() as usize]);
                    transcript.add_coverage_counts(
                        block.reference_start,
                        block.reference_start + block.length - 1,
                        coverage,
                    );
                }
            }

            for lf in &locus_functions {
                self.metrics.pf_aligned_bases += 1;
                match lf {
                    LocusFunction::Intergenic => self.metrics.intergenic_bases += 1,
                    LocusFunction::Intronic => self.metrics.intronic_bases += 1,
                    LocusFunction::Utr => {
                        self.metrics.utr_bases += 1;
                        overlaps_exon = true;
                    }
                    LocusFunction::Coding => {
                        self.metrics.coding_bases += 1;
                        overlaps_exon = true;
                    }
                    LocusFunction::Ribosomal => {
                        // Unreachable via transcripts (they assign at most CODING); kept for parity.
                        if let Some(rb) = self.metrics.ribosomal_bases.as_mut() {
                            *rb += 1;
                        }
                    }
                }
            }
        }

        // Strandedness is charged per read, and only when the read hits exactly one gene's exon.
        if !has(flags, SUPPLEMENTARY_ALIGNMENT) && overlaps_exon && overlapping_genes.len() == 1 {
            let gene = overlapping_genes[0];
            let negative_transcription_strand = gene.negative_strand;
            let read_one_or_unpaired = !has(flags, READ_PAIRED) || has(flags, FIRST_OF_PAIR);
            let negative_read_strand = has(flags, READ_NEGATIVE_STRAND);

            if self.strand_specificity != StrandSpecificity::None {
                let strands_agree = negative_read_strand == negative_transcription_strand;
                let first_expected_to_agree =
                    self.strand_specificity == StrandSpecificity::FirstReadTranscriptionStrand;
                let this_read_expected_to_agree = read_one_or_unpaired == first_expected_to_agree;
                if strands_agree == this_read_expected_to_agree {
                    self.metrics.correct_strand_reads += 1;
                } else {
                    self.metrics.incorrect_strand_reads += 1;
                }
            }

            if read_one_or_unpaired {
                let proper_orientation: bool;
                let left_most: i32;
                let right_most: i32;
                if has(flags, READ_PAIRED) {
                    if has(flags, MATE_UNMAPPED) {
                        proper_orientation = false;
                        left_most = 0;
                        right_most = 0;
                    } else {
                        // No MC tag is set by the corpora, so the mate length falls back to the
                        // read length, matching SAMUtils.getMateCigar == null.
                        let mate_reference_length =
                            mate_cigar_reference_length(rec).unwrap_or(rec.read_length() as i32);
                        let mate_alignment_end =
                            rec.mate_alignment_start + mate_reference_length - 1;
                        proper_orientation = pair_orientation(rec) == PairOrientation::Fr;
                        left_most = read_start.min(rec.mate_alignment_start);
                        right_most = read_end.max(mate_alignment_end);
                    }
                } else {
                    proper_orientation = true;
                    left_most = read_start;
                    right_most = read_end;
                }

                // CoordMath.encloses(gene.start, gene.end, leftMost, rightMost).
                let enclosed = left_most >= gene.start && right_most <= gene.end;
                if proper_orientation && enclosed {
                    if negative_read_strand == negative_transcription_strand {
                        self.metrics.num_r1_transcript_strand_reads += 1;
                    } else {
                        self.metrics.num_r2_transcript_strand_reads += 1;
                    }
                } else {
                    self.metrics.num_unexplained_reads += 1;
                }
            }
        }
    }

    /// `finish` then `addMetricsToFile`: closes the ratios, computes the coverage metrics, and
    /// returns the metric row plus the (possibly empty) normalized-coverage histogram.
    pub fn finish(mut self) -> (RnaSeqMetrics, OutHistogram) {
        // finish()
        if self.metrics.pf_aligned_bases > 0 {
            let aligned = self.metrics.pf_aligned_bases as f64;
            if let Some(rb) = self.metrics.ribosomal_bases {
                self.metrics.pct_ribosomal_bases = Some(rb as f64 / aligned);
            }
            self.metrics.pct_coding_bases = self.metrics.coding_bases as f64 / aligned;
            self.metrics.pct_utr_bases = self.metrics.utr_bases as f64 / aligned;
            self.metrics.pct_intronic_bases = self.metrics.intronic_bases as f64 / aligned;
            self.metrics.pct_intergenic_bases = self.metrics.intergenic_bases as f64 / aligned;
            self.metrics.pct_mrna_bases =
                self.metrics.pct_coding_bases + self.metrics.pct_utr_bases;
            self.metrics.pct_usable_bases = (self.metrics.coding_bases + self.metrics.utr_bases)
                as f64
                / self.metrics.pf_bases as f64;
        }
        if self.metrics.correct_strand_reads > 0 || self.metrics.incorrect_strand_reads > 0 {
            self.metrics.pct_correct_strand_reads = self.metrics.correct_strand_reads as f64
                / (self.metrics.correct_strand_reads + self.metrics.incorrect_strand_reads) as f64;
        }
        let reads_examined = self.metrics.num_r1_transcript_strand_reads
            + self.metrics.num_r2_transcript_strand_reads;
        if reads_examined > 0 {
            self.metrics.pct_r1_transcript_strand_reads =
                self.metrics.num_r1_transcript_strand_reads as f64 / reads_examined as f64;
            self.metrics.pct_r2_transcript_strand_reads =
                self.metrics.num_r2_transcript_strand_reads as f64 / reads_examined as f64;
        }

        let histogram = self.compute_coverage_metrics();
        (self.metrics, histogram)
    }

    /// `computeCoverageMetrics`.
    fn compute_coverage_metrics(&mut self) -> OutHistogram {
        let picked = self.pick_transcripts();
        let transcript_count = picked.len() as f64;

        let mut cvs = Histogram::new("", "");
        let mut five_prime_skews = Histogram::new("", "");
        let mut three_prime_skews = Histogram::new("", "");
        let mut five_to_three_skews = Histogram::new("", "");
        let mut normalized = Histogram::new("normalized_position", "All_Reads.normalized_coverage");

        for tx in &picked {
            // promote(int[]) then reverse when the gene is on the negative strand.
            let mut coverage: Vec<f64> = tx.coverage.iter().map(|&c| c as f64).collect();
            if tx.negative_strand {
                coverage.reverse();
            }
            let mean = mathutil_mean(&coverage, 0, coverage.len());

            let stdev = mathutil_stddev(&coverage, mean);
            cvs.increment(stdev / mean);

            let five_prime = mathutil_mean(&coverage, 0, self.end_bias_bases as usize);
            let three_prime = mathutil_mean(
                &coverage,
                coverage.len() - self.end_bias_bases as usize,
                coverage.len(),
            );
            five_prime_skews.increment(five_prime / mean);
            three_prime_skews.increment(three_prime / mean);
            five_to_three_skews.increment(mathutil_divide(five_prime, three_prime));

            let last_index = (coverage.len() - 1) as f64;
            for percent in 0..=100i32 {
                let p = percent as f64 / 100.0;
                let start = ((last_index * (p - 0.005)).max(0.0)) as i32;
                let end = ((last_index * (p + 0.005)).min(last_index)) as i32;
                let length = end - start + 1;
                let mut sum = 0.0;
                for i in start..=end {
                    sum += coverage[i as usize];
                }
                let norm = (sum / length as f64) / mean;
                normalized.increment_by(percent as f64, norm / transcript_count);
            }
        }

        self.metrics.median_cv_coverage = cvs.median();
        self.metrics.median_5prime_bias = five_prime_skews.median();
        self.metrics.median_3prime_bias = three_prime_skews.median();
        self.metrics.median_5prime_to_3prime_bias = five_to_three_skews.median();

        // Convert the htsjdk Histogram<Integer> into the printable form: keys sorted numerically
        // and rendered as plain integers.
        let mut bins: Vec<(i64, f64)> = normalized.bins().map(|(k, v)| (k as i64, v)).collect();
        bins.sort_by_key(|(k, _)| *k);
        OutHistogram {
            bin_label: "normalized_position".to_string(),
            value_label: "All_Reads.normalized_coverage".to_string(),
            key_class: "java.lang.Integer".to_string(),
            bins: bins.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        }
    }

    /// `pickTranscripts`, returning the picked transcripts already in Java `HashMap` iteration
    /// order (which is the order `computeCoverageMetrics` folds them in).
    fn pick_transcripts(&self) -> Vec<PickedTranscript> {
        // bestPerGene: the highest-mean-coverage qualifying transcript of each gene.
        struct Best {
            coverage: Vec<i32>,
            negative_strand: bool,
            mean: f64,
            /// A deterministic content sort key: the four coordinates then the name. Used only to
            /// fix a reproducible fold order, not to reproduce Java's HashMap order (measured
            /// unobservable, see the module docs and decision 0003).
            sort_key: (i32, i32, i32, i32, String),
        }
        let mut best_per_gene: Vec<Best> = Vec::new();

        let bias = self.minimum_length.max(self.end_bias_bases);
        for gene in self.gene_overlap.get_all() {
            let gene_ptr = gene as *const Gene as usize;
            let mut best: Option<Best> = None;
            for (tx_idx, tx) in gene.transcripts.iter().enumerate() {
                let cov = match self.coverage_by_transcript.get(&(gene_ptr, tx_idx)) {
                    Some(c) => c,
                    None => continue, // uncovered transcript: absent from the coverage map
                };
                if tx.length() < bias {
                    continue;
                }
                let mean = mathutil_mean_ints(cov);
                if mean < 1.0 {
                    continue;
                }
                let better = match &best {
                    None => true,
                    Some(b) => mean > b.mean,
                };
                if better {
                    best = Some(Best {
                        coverage: cov.clone(),
                        negative_strand: gene.negative_strand,
                        mean,
                        sort_key: (
                            tx.transcription_start,
                            tx.transcription_end,
                            tx.coding_start,
                            tx.coding_end,
                            tx.name.clone(),
                        ),
                    });
                }
            }
            if let Some(b) = best {
                best_per_gene.push(b);
            }
        }

        // The 1000th-best coverage is the floor; everything at or above it is kept.
        let mut coverages: Vec<f64> = best_per_gene.iter().map(|b| b.mean).collect();
        coverages.sort_by(|a, b| a.total_cmp(b));
        let min = if coverages.is_empty() {
            0.0
        } else {
            coverages[coverages.len().saturating_sub(1001)]
        };

        let mut kept: Vec<Best> = best_per_gene
            .into_iter()
            .filter(|b| b.mean >= min)
            .collect();
        kept.sort_by(|a, b| a.sort_key.cmp(&b.sort_key));
        kept.into_iter()
            .map(|b| PickedTranscript {
                coverage: b.coverage,
                negative_strand: b.negative_strand,
            })
            .collect()
    }
}

struct PickedTranscript {
    coverage: Vec<i32>,
    negative_strand: bool,
}

// --- MathUtil, ported exactly (naive sequential sums) ---

/// `MathUtil.mean(in, start, stop)`.
fn mathutil_mean(input: &[f64], start: usize, stop: usize) -> f64 {
    let mut total = 0.0;
    for &v in &input[start..stop] {
        total += v;
    }
    total / (stop - start) as f64
}

fn mathutil_mean_ints(input: &[i32]) -> f64 {
    // MathUtil.mean(MathUtil.promote(cov), 0, cov.length): promote then sum.
    let mut total = 0.0;
    for &v in input {
        total += v as f64;
    }
    total / input.len() as f64
}

/// `MathUtil.stddev(in, mean)` = `stddev(in, 0, in.length, mean)`.
fn mathutil_stddev(input: &[f64], mean: f64) -> f64 {
    let mut total = 0.0;
    for &v in input {
        total += v * v;
    }
    ((total / input.len() as f64) - (mean * mean)).sqrt()
}

/// `MathUtil.divide`: guarded against a (near-)zero denominator.
fn mathutil_divide(numerator: f64, denominator: f64) -> f64 {
    if (0.0 - denominator).abs() > 0.000001 {
        numerator / denominator
    } else {
        0.0
    }
}

/// `SAMUtils.getMateCigar`: reads the "MC" tag's reference length, or `None` when absent.
fn mate_cigar_reference_length(rec: &BamRecord) -> Option<i32> {
    use htsjdk_bam::tag::{Tag, TagValue};
    let tag = Tag::new(b"MC");
    match rec.tags.get(tag) {
        Some(TagValue::Str(s)) => htsjdk_bam::text_parse::parse_cigar(s)
            .ok()
            .map(|c| c.reference_length() as i32),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn divide_guards_a_zero_denominator() {
        assert_eq!(mathutil_divide(1.0, 0.0), 0.0);
        assert_eq!(mathutil_divide(6.0, 3.0), 2.0);
    }

    #[test]
    fn mean_and_stddev_are_the_naive_sequential_forms() {
        let xs = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(mathutil_mean(&xs, 0, 4), 2.5);
        // sqrt(mean(x^2) - mean^2) = sqrt(7.5 - 6.25) = sqrt(1.25).
        assert!((mathutil_stddev(&xs, 2.5) - 1.25_f64.sqrt()).abs() < 1e-12);
    }
}
