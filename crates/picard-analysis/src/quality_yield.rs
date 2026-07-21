//! `CollectQualityYieldMetrics`.
//!
//! Ported from `picard.analysis.CollectQualityYieldMetrics`, tag 3.4.0: the
//! `QualityYieldMetricsCollector` inner class and the `QualityYieldMetrics` bean.
//!
//! This is the first member of the metrics-collector archetype, which covers 57 of the 311
//! tools. It is ported at full price so the second and third can be measured against it.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, Value};

/// `SAMFlag.NOT_PRIMARY_ALIGNMENT`.
pub const SECONDARY_ALIGNMENT: u16 = 0x100;
/// `SAMFlag.READ_FAILS_VENDOR_QUALITY_CHECK`.
pub const VENDOR_QUALITY_CHECK_FAILED: u16 = 0x200;
/// `SAMFlag.SUPPLEMENTARY_ALIGNMENT`.
pub const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

/// `CollectQualityYieldMetrics$QualityYieldMetrics`.
///
/// Field order is the declaration order in the Java source, because that is what HotSpot's
/// `Class.getFields()` returns and therefore the column order of the output file.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct QualityYieldMetrics {
    pub total_reads: i64,
    pub pf_reads: i64,
    pub read_length: i32,
    pub total_bases: i64,
    pub pf_bases: i64,
    pub q20_bases: i64,
    pub pf_q20_bases: i64,
    pub q30_bases: i64,
    pub pf_q30_bases: i64,
    pub q20_equivalent_yield: i64,
    pub pf_q20_equivalent_yield: i64,
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
    "Q20_EQUIVALENT_YIELD",
    "PF_Q20_EQUIVALENT_YIELD",
];

impl MetricBean for QualityYieldMetrics {
    fn class_name(&self) -> &str {
        // The `$` is Java's inner-class separator and appears verbatim in the file.
        "picard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }

    fn values(&self) -> Vec<Value> {
        vec![
            Value::Long(self.total_reads),
            Value::Long(self.pf_reads),
            Value::Long(self.read_length as i64),
            Value::Long(self.total_bases),
            Value::Long(self.pf_bases),
            Value::Long(self.q20_bases),
            Value::Long(self.pf_q20_bases),
            Value::Long(self.q30_bases),
            Value::Long(self.pf_q30_bases),
            Value::Long(self.q20_equivalent_yield),
            Value::Long(self.pf_q20_equivalent_yield),
        ]
    }
}

/// The tool's arguments, with Picard's defaults.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `USE_ORIGINAL_QUALITIES`, default **true**.
    pub use_original_qualities: bool,
    pub include_secondary_alignments: bool,
    pub include_supplemental_alignments: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            use_original_qualities: true,
            include_secondary_alignments: false,
            include_supplemental_alignments: false,
        }
    }
}

/// `QualityYieldMetricsCollector`.
#[derive(Debug, Clone)]
pub struct QualityYieldMetricsCollector {
    options: Options,
    metrics: QualityYieldMetrics,
}

/// `SAMUtils.fastqToPhred`: one Phred score per character, offset 33.
///
/// The cast to `byte` is Java's, and it is why a character above 160 wraps into a negative
/// score rather than being rejected.
fn fastq_to_phred(s: &str) -> Vec<i8> {
    s.chars().map(|c| (c as u32 as u8).wrapping_sub(33) as i8).collect()
}

impl QualityYieldMetricsCollector {
    pub fn new(options: Options) -> Self {
        QualityYieldMetricsCollector {
            options,
            metrics: QualityYieldMetrics::default(),
        }
    }

    /// `QualityYieldMetricsCollector.acceptRecord`.
    pub fn accept(&mut self, rec: &BamRecord) {
        if !self.options.include_secondary_alignments && rec.flags & SECONDARY_ALIGNMENT != 0 {
            return;
        }
        if !self.options.include_supplemental_alignments
            && rec.flags & SUPPLEMENTARY_ALIGNMENT != 0
        {
            return;
        }

        let length = rec.read_length() as i64;
        self.metrics.total_reads += 1;
        self.metrics.total_bases += length;

        let is_pf_read = rec.flags & VENDOR_QUALITY_CHECK_FAILED == 0;
        if is_pf_read {
            self.metrics.pf_reads += 1;
            self.metrics.pf_bases += length;
        }

        // `getOriginalBaseQualities()` returns null when OQ is absent or empty, and the
        // collector then falls back to the primary qualities.
        let original: Option<Vec<i8>> = if self.options.use_original_qualities {
            match rec.tags.get(Tag::new(b"OQ")) {
                Some(TagValue::Str(s)) if !s.is_empty() => Some(fastq_to_phred(s)),
                _ => None,
            }
        } else {
            None
        };
        // Java iterates `byte[]` widened to `int`, so a stored quality above 127 is negative
        // here rather than large. Reproduced by keeping the scores signed.
        let quals: Vec<i8> = match original {
            Some(oq) => oq,
            None => rec.base_qualities.iter().map(|&b| b as i8).collect(),
        };

        for qual in quals {
            let qual = qual as i64;
            self.metrics.q20_equivalent_yield += qual;
            if qual >= 30 {
                self.metrics.q20_bases += 1;
                self.metrics.q30_bases += 1;
            } else if qual >= 20 {
                self.metrics.q20_bases += 1;
            }
            if is_pf_read {
                self.metrics.pf_q20_equivalent_yield += qual;
                if qual >= 30 {
                    self.metrics.pf_q20_bases += 1;
                    self.metrics.pf_q30_bases += 1;
                } else if qual >= 20 {
                    self.metrics.pf_q20_bases += 1;
                }
            }
        }
    }

    /// `QualityYieldMetricsCollector.finish` plus `calculateDerivedFields`.
    ///
    /// Both divisions are integer divisions, truncating. `Q20_EQUIVALENT_YIELD` is a sum of
    /// quality scores divided by 20, and `READ_LENGTH` is the mean read length floored, so a
    /// file of 99-base and 101-base reads reports 100.
    pub fn finish(&mut self) {
        self.metrics.q20_equivalent_yield /= 20;
        self.metrics.pf_q20_equivalent_yield /= 20;
        self.metrics.read_length = if self.metrics.total_reads == 0 {
            0
        } else {
            (self.metrics.total_bases / self.metrics.total_reads) as i32
        };
    }

    pub fn metrics(&self) -> &QualityYieldMetrics {
        &self.metrics
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

    fn read(quals: Vec<u8>, flags: u16) -> BamRecord {
        let n = quals.len();
        BamRecord {
            read_name: "r".into(),
            flags,
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                length: n as u32,
                op: Op::M,
            }]),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: vec![b'A'; n],
            base_qualities: quals,
            tags: Default::default(),
        }
    }

    fn collect(records: &[BamRecord], options: Options) -> QualityYieldMetrics {
        let mut c = QualityYieldMetricsCollector::new(options);
        for r in records {
            c.accept(r);
        }
        c.finish();
        c.metrics().clone()
    }

    /// A quality of exactly 20 counts toward Q20 and not Q30; exactly 30 counts toward both.
    #[test]
    fn the_thresholds_are_inclusive_at_twenty_and_thirty() {
        let m = collect(&[read(vec![19, 20, 29, 30], 0)], Options::default());
        assert_eq!(m.q20_bases, 3, "20, 29 and 30 are all >= 20");
        assert_eq!(m.q30_bases, 1, "only 30 is >= 30");
    }

    /// The equivalent yield is the quality sum divided by 20, truncating.
    #[test]
    fn the_equivalent_yield_is_an_integer_division() {
        // 39 = one base of quality 39; 39/20 = 1.
        let m = collect(&[read(vec![39], 0)], Options::default());
        assert_eq!(m.q20_equivalent_yield, 1);
        let m = collect(&[read(vec![40], 0)], Options::default());
        assert_eq!(m.q20_equivalent_yield, 2);
    }

    /// `READ_LENGTH` is the mean floored, not rounded.
    #[test]
    fn the_read_length_is_a_floored_mean() {
        let m = collect(
            &[read(vec![30; 99], 0), read(vec![30; 101], 0)],
            Options::default(),
        );
        assert_eq!(m.read_length, 100);
        let m = collect(
            &[read(vec![30; 10], 0), read(vec![30; 11], 0)],
            Options::default(),
        );
        assert_eq!(m.read_length, 10, "21 / 2 truncates to 10");
    }

    #[test]
    fn an_empty_input_reports_zero_read_length_rather_than_dividing() {
        let m = collect(&[], Options::default());
        assert_eq!(m.read_length, 0);
        assert_eq!(m.total_reads, 0);
    }

    /// A vendor-failed read counts toward TOTAL but not PF, and its qualities are excluded
    /// from every PF counter.
    #[test]
    fn vendor_failed_reads_count_toward_total_only() {
        let m = collect(
            &[
                read(vec![30; 10], 0),
                read(vec![30; 10], VENDOR_QUALITY_CHECK_FAILED),
            ],
            Options::default(),
        );
        assert_eq!((m.total_reads, m.pf_reads), (2, 1));
        assert_eq!((m.total_bases, m.pf_bases), (20, 10));
        assert_eq!((m.q30_bases, m.pf_q30_bases), (20, 10));
    }

    #[test]
    fn secondary_and_supplementary_are_excluded_by_default() {
        let records = [
            read(vec![30; 10], 0),
            read(vec![30; 10], SECONDARY_ALIGNMENT),
            read(vec![30; 10], SUPPLEMENTARY_ALIGNMENT),
        ];
        assert_eq!(collect(&records, Options::default()).total_reads, 1);
        assert_eq!(
            collect(
                &records,
                Options {
                    include_secondary_alignments: true,
                    ..Options::default()
                }
            )
            .total_reads,
            2
        );
        assert_eq!(
            collect(
                &records,
                Options {
                    include_supplemental_alignments: true,
                    ..Options::default()
                }
            )
            .total_reads,
            2
        );
    }

    /// `USE_ORIGINAL_QUALITIES` defaults to **true**, so a record carrying `OQ` is measured on
    /// its original scores rather than its current ones. Getting the default backwards would
    /// silently measure the wrong column of data.
    #[test]
    fn the_oq_tag_is_preferred_by_default() {
        let mut r = read(vec![10; 4], 0);
        // "IIII" is Phred 40 at offset 33.
        r.tags.insert(Tag::new(b"OQ"), TagValue::Str("IIII".into()));

        let with_oq = collect(std::slice::from_ref(&r), Options::default());
        assert_eq!(with_oq.q30_bases, 4, "OQ says 40, so all four are >= 30");

        let without = collect(
            std::slice::from_ref(&r),
            Options {
                use_original_qualities: false,
                ..Options::default()
            },
        );
        assert_eq!(without.q30_bases, 0, "the primary qualities are 10");
    }

    /// An empty `OQ` is treated as absent, matching `getOriginalBaseQualities`.
    #[test]
    fn an_empty_oq_falls_back_to_the_primary_qualities() {
        let mut r = read(vec![30; 4], 0);
        r.tags.insert(Tag::new(b"OQ"), TagValue::Str(String::new()));
        assert_eq!(collect(&[r], Options::default()).q30_bases, 4);
    }

    #[test]
    fn fastq_to_phred_uses_offset_thirty_three() {
        assert_eq!(fastq_to_phred("!"), vec![0]);
        assert_eq!(fastq_to_phred("I"), vec![40]);
        assert_eq!(fastq_to_phred("#$%"), vec![2, 3, 4]);
    }

    /// The column order is the Java declaration order, which is what the output file uses.
    #[test]
    fn the_columns_are_in_declaration_order() {
        let m = QualityYieldMetrics::default();
        assert_eq!(m.columns()[0], "TOTAL_READS");
        assert_eq!(m.columns()[2], "READ_LENGTH");
        assert_eq!(m.columns().len(), m.values().len());
    }
}
