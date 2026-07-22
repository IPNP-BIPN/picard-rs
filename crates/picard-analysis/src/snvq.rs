//! `CollectQualityYieldMetricsSNVQ`.
//!
//! Ported from `picard.analysis.CollectQualityYieldMetricsSNVQ` and its inner
//! `QualityYieldMetrics` (a `MergeableMetricBase`), tag 3.4.0. A fan-out of the QualityYield
//! archetype: it counts base-quality yield at the 20/30/40 thresholds exactly as
//! `CollectQualityYieldMetrics` does, and adds the same counts for the **SNV qualities**, the
//! per-alternate-base qualities stored in the tags `qa`, `qc`, `qg`, `qt` (q + a lowercased base
//! of `ACGT`), FASTQ encoded.
//!
//! It is entirely integer counting plus ratio derivation, with no transcendental math, so it is
//! byte-identical without any of the `gatk-jmath` machinery. The default output (this port) is the
//! metrics row only; `INCLUDE_BQ_HISTOGRAM` adds `SeriesStats`-backed histograms and is a separate
//! surface not ported here.
//!
//! The SNVQ rule for one base position: for each of the four bases `ACGT` that is **not** the read
//! base at that position, the corresponding tag's quality is one SNVQ observation. A non-`ACGT`
//! base (an `N`) is unequal to all four, so it contributes four SNVQ observations.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, Value};

const SECONDARY_ALIGNMENT: u16 = 0x100;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

/// `SNVQ_BASE_ORDER = "ACGT"`, and the tag for each is `q` + its lowercased letter.
const BASE_ORDER: [u8; 4] = *b"ACGT";
const SNVQ_TAGS: [&[u8; 2]; 4] = [b"qa", b"qc", b"qg", b"qt"];

/// `CollectQualityYieldMetricsSNVQ$QualityYieldMetrics`, fields in declaration order.
#[derive(Debug, Clone, Default)]
pub struct SnvqMetrics {
    pub total_reads: i64,
    pub pf_reads: i64,
    pub read_length: i64,
    pub total_bases: i64,
    pub pf_bases: i64,
    pub q20_bases: i64,
    pub pf_q20_bases: i64,
    pub q30_bases: i64,
    pub pf_q30_bases: i64,
    pub q40_bases: i64,
    pub pf_q40_bases: i64,
    pub pct_q20_bases: f64,
    pub pct_q30_bases: f64,
    pub pct_q40_bases: f64,
    pub pct_pf_q20_bases: f64,
    pub pct_pf_q30_bases: f64,
    pub pct_pf_q40_bases: f64,
    pub total_snvq: i64,
    pub pf_snvq: i64,
    pub q20_snvq: i64,
    pub pf_q20_snvq: i64,
    pub q30_snvq: i64,
    pub pf_q30_snvq: i64,
    pub q40_snvq: i64,
    pub pf_q40_snvq: i64,
    pub pct_q20_snvq: f64,
    pub pct_q30_snvq: f64,
    pub pct_q40_snvq: f64,
    pub pct_pf_q20_snvq: f64,
    pub pct_pf_q30_snvq: f64,
    pub pct_pf_q40_snvq: f64,
}

const COLUMNS: &[&str] = &[
    "TOTAL_READS",
    "PF_READS",
    "READ_LENGTH",
    "TOTAL_BASES",
    "PF_BASES",
    "Q20_BASES",
    "PF_Q20_BASES",
    "Q30_BASES",
    "PF_Q30_BASES",
    "Q40_BASES",
    "PF_Q40_BASES",
    "PCT_Q20_BASES",
    "PCT_Q30_BASES",
    "PCT_Q40_BASES",
    "PCT_PF_Q20_BASES",
    "PCT_PF_Q30_BASES",
    "PCT_PF_Q40_BASES",
    "TOTAL_SNVQ",
    "PF_SNVQ",
    "Q20_SNVQ",
    "PF_Q20_SNVQ",
    "Q30_SNVQ",
    "PF_Q30_SNVQ",
    "Q40_SNVQ",
    "PF_Q40_SNVQ",
    "PCT_Q20_SNVQ",
    "PCT_Q30_SNVQ",
    "PCT_Q40_SNVQ",
    "PCT_PF_Q20_SNVQ",
    "PCT_PF_Q30_SNVQ",
    "PCT_PF_Q40_SNVQ",
];

impl MetricBean for SnvqMetrics {
    fn class_name(&self) -> &str {
        "picard.analysis.CollectQualityYieldMetricsSNVQ$QualityYieldMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Long(self.total_reads),
            Value::Long(self.pf_reads),
            Value::Long(self.read_length),
            Value::Long(self.total_bases),
            Value::Long(self.pf_bases),
            Value::Long(self.q20_bases),
            Value::Long(self.pf_q20_bases),
            Value::Long(self.q30_bases),
            Value::Long(self.pf_q30_bases),
            Value::Long(self.q40_bases),
            Value::Long(self.pf_q40_bases),
            Value::Double(self.pct_q20_bases),
            Value::Double(self.pct_q30_bases),
            Value::Double(self.pct_q40_bases),
            Value::Double(self.pct_pf_q20_bases),
            Value::Double(self.pct_pf_q30_bases),
            Value::Double(self.pct_pf_q40_bases),
            Value::Long(self.total_snvq),
            Value::Long(self.pf_snvq),
            Value::Long(self.q20_snvq),
            Value::Long(self.pf_q20_snvq),
            Value::Long(self.q30_snvq),
            Value::Long(self.pf_q30_snvq),
            Value::Long(self.q40_snvq),
            Value::Long(self.pf_q40_snvq),
            Value::Double(self.pct_q20_snvq),
            Value::Double(self.pct_q30_snvq),
            Value::Double(self.pct_q40_snvq),
            Value::Double(self.pct_pf_q20_snvq),
            Value::Double(self.pct_pf_q30_snvq),
            Value::Double(self.pct_pf_q40_snvq),
        ]
    }
}

/// The collector. `INCLUDE_SECONDARY_ALIGNMENTS`/`INCLUDE_SUPPLEMENTAL_ALIGNMENTS` default false,
/// as in the tool; the default (`ALTERNATE_QUALITY_ATTRIBUTE = null`) quality source is the read's
/// own base qualities.
#[derive(Default)]
pub struct SnvqCollector {
    include_secondary: bool,
    include_supplemental: bool,
    metrics: SnvqMetrics,
}

impl SnvqCollector {
    pub fn new(include_secondary: bool, include_supplemental: bool) -> Self {
        SnvqCollector {
            include_secondary,
            include_supplemental,
            metrics: SnvqMetrics::default(),
        }
    }

    fn tag_bytes<'a>(rec: &'a BamRecord, name: &[u8; 2]) -> &'a [u8] {
        match rec.tags.get(Tag::new(name)) {
            Some(TagValue::Str(s)) => s.as_bytes(),
            _ => panic!("CollectQualityYieldMetricsSNVQ requires the {:?} tag", name),
        }
    }

    pub fn accept(&mut self, rec: &BamRecord) {
        if !self.include_secondary && rec.flags & SECONDARY_ALIGNMENT != 0 {
            return;
        }
        if !self.include_supplemental && rec.flags & SUPPLEMENTARY_ALIGNMENT != 0 {
            return;
        }

        let m = &mut self.metrics;
        let length = rec.read_length() as i64;
        m.total_reads += 1;
        m.total_bases += length;

        let is_pf = rec.flags & READ_FAILS_VENDOR_QUALITY == 0;
        if is_pf {
            m.pf_reads += 1;
            m.pf_bases += length;
        }

        let quals = &rec.base_qualities;
        let bases = &rec.read_bases;
        assert_eq!(
            quals.len(),
            bases.len(),
            "quality string length does not match bases string"
        );
        let snvq: [&[u8]; 4] = [
            Self::tag_bytes(rec, SNVQ_TAGS[0]),
            Self::tag_bytes(rec, SNVQ_TAGS[1]),
            Self::tag_bytes(rec, SNVQ_TAGS[2]),
            Self::tag_bytes(rec, SNVQ_TAGS[3]),
        ];

        for (read_position, &qb) in quals.iter().enumerate() {
            let qual = qb as i32;
            // Base-quality yield: Q40 implies Q30 implies Q20.
            if qual >= 40 {
                m.q20_bases += 1;
                m.q30_bases += 1;
                m.q40_bases += 1;
            } else if qual >= 30 {
                m.q20_bases += 1;
                m.q30_bases += 1;
            } else if qual >= 20 {
                m.q20_bases += 1;
            }
            if is_pf {
                if qual >= 40 {
                    m.pf_q20_bases += 1;
                    m.pf_q30_bases += 1;
                    m.pf_q40_bases += 1;
                } else if qual >= 30 {
                    m.pf_q20_bases += 1;
                    m.pf_q30_bases += 1;
                } else if qual >= 20 {
                    m.pf_q20_bases += 1;
                }
            }

            let base = bases[read_position];
            for i in 0..BASE_ORDER.len() {
                if base != BASE_ORDER[i] {
                    // fastqToPhred: the FASTQ character minus 33.
                    let q = snvq[i][read_position] as i32 - 33;
                    m.total_snvq += 1;
                    if is_pf {
                        m.pf_snvq += 1;
                    }
                    if q >= 40 {
                        m.q20_snvq += 1;
                        m.q30_snvq += 1;
                        m.q40_snvq += 1;
                        if is_pf {
                            m.pf_q20_snvq += 1;
                            m.pf_q30_snvq += 1;
                            m.pf_q40_snvq += 1;
                        }
                    } else if q >= 30 {
                        m.q20_snvq += 1;
                        m.q30_snvq += 1;
                        if is_pf {
                            m.pf_q20_snvq += 1;
                            m.pf_q30_snvq += 1;
                        }
                    } else if q >= 20 {
                        m.q20_snvq += 1;
                        if is_pf {
                            m.pf_q20_snvq += 1;
                        }
                    }
                }
            }
        }
    }

    /// `calculateDerivedFields`: the ratios and the floored mean read length.
    pub fn finish(mut self) -> SnvqMetrics {
        let m = &mut self.metrics;
        m.read_length = if m.total_reads == 0 {
            0
        } else {
            m.total_bases / m.total_reads // (int) truncation
        };

        if m.total_bases != 0 {
            let t = m.total_bases as f64;
            m.pct_q20_bases = m.q20_bases as f64 / t;
            m.pct_q30_bases = m.q30_bases as f64 / t;
            m.pct_q40_bases = m.q40_bases as f64 / t;
        }
        if m.pf_bases != 0 {
            let t = m.pf_bases as f64;
            m.pct_pf_q20_bases = m.pf_q20_bases as f64 / t;
            m.pct_pf_q30_bases = m.pf_q30_bases as f64 / t;
            m.pct_pf_q40_bases = m.pf_q40_bases as f64 / t;
        }
        if m.total_snvq != 0 {
            let t = m.total_snvq as f64;
            m.pct_q20_snvq = m.q20_snvq as f64 / t;
            m.pct_q30_snvq = m.q30_snvq as f64 / t;
            m.pct_q40_snvq = m.q40_snvq as f64 / t;
        }
        if m.pf_snvq != 0 {
            let t = m.pf_snvq as f64;
            m.pct_pf_q20_snvq = m.pf_q20_snvq as f64 / t;
            m.pct_pf_q30_snvq = m.pf_q30_snvq as f64 / t;
            m.pct_pf_q40_snvq = m.pf_q40_snvq as f64 / t;
        }
        self.metrics
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_length_is_the_floored_mean() {
        // 16 bases over 3 reads floors to 5.
        let m = SnvqMetrics {
            total_reads: 3,
            total_bases: 16,
            ..Default::default()
        };
        let c = SnvqCollector {
            metrics: m,
            ..Default::default()
        };
        assert_eq!(c.finish().read_length, 5);
    }

    #[test]
    fn a_zero_denominator_leaves_the_ratio_zero() {
        let c = SnvqCollector::default();
        let m = c.finish();
        assert_eq!(m.pct_q20_bases, 0.0);
        assert_eq!(m.read_length, 0);
    }
}
