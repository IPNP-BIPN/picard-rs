//! `QualityScoreDistribution`.
//!
//! Ported from `picard.analysis.QualityScoreDistribution`, tag 3.4.0.
//!
//! The smallest tool in the metrics stratum, and the only one so far whose output has **no
//! metric rows at all**: the entire body is one or two histogram tables. That makes it the case
//! that exercises two things the other suites never reach.
//!
//! **The histogram key class is `java.lang.Byte`, not `Integer`.** The counts are incremented
//! with `qHisto.increment((byte) i, ...)`, so the `## HISTOGRAM` line names `java.lang.Byte` and
//! the keys sort as bytes.
//!
//! **The two histograms have different key sets.** The `OQ` histogram is emitted only when it is
//! non-empty, and when it is, its qualities need not overlap the primary ones at all. htsjdk
//! builds the combined key set as a `TreeSet` over the *first non-empty* histogram's comparator,
//! so the union is re-sorted rather than concatenated. A file whose Q histogram holds only
//! quality 30 and whose OQ histogram holds 2, 3, 40 and 41 comes out as `2 3 30 40 41`, with the
//! Q key in the middle. That is the case htsjdk-rs's key-union fix was written for, and here it
//! arrives through a real tool rather than a synthetic test.
//!
//! Two filters worth stating because they are asymmetric:
//!
//!   * secondary **and** supplementary records are dropped, by `isSecondaryOrSupplementary`;
//!   * vendor-failed and unmapped reads are **kept** by default, and dropped only when
//!     `PF_READS_ONLY` or `ALIGNED_READS_ONLY` is set. So the default distribution includes
//!     reads that most other collectors in this crate exclude.
//!
//! And the no-call rule is on the **base**, not the quality: `INCLUDE_NO_CALLS` defaults to
//! false, so the quality at an `N` position is discarded however good it is.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::is_no_call;
use htsjdk_bam::tag::{Tag, TagValue};

const READ_UNMAPPED: u16 = 0x4;
const SECONDARY: u16 = 0x100;
const VENDOR_FAILED: u16 = 0x200;
const SUPPLEMENTARY: u16 = 0x800;

/// The tool's arguments, with Picard's defaults.
#[derive(Debug, Clone, Copy, Default)]
pub struct Options {
    pub aligned_reads_only: bool,
    pub pf_reads_only: bool,
    pub include_no_calls: bool,
}

/// `QualityScoreDistribution`, accumulating the two 128-entry count arrays.
///
/// The arrays are 128 long because htsjdk indexes them with a Java `byte`, which is signed: a
/// quality above 127 would index negatively and throw rather than wrap. BAM qualities are
/// 0..93 and `OQ` characters are 33..126, so the range is never approached, but the size is
/// htsjdk's and is kept.
pub struct QualityScoreDistribution {
    options: Options,
    q_counts: [i64; 128],
    oq_counts: [i64; 128],
}

impl QualityScoreDistribution {
    pub fn new(options: Options) -> Self {
        QualityScoreDistribution {
            options,
            q_counts: [0; 128],
            oq_counts: [0; 128],
        }
    }

    /// `acceptRead`.
    pub fn accept(&mut self, rec: &BamRecord) {
        if self.options.pf_reads_only && rec.flags & VENDOR_FAILED != 0 {
            return;
        }
        if self.options.aligned_reads_only && rec.flags & READ_UNMAPPED != 0 {
            return;
        }
        if rec.flags & SECONDARY != 0 || rec.flags & SUPPLEMENTARY != 0 {
            return;
        }

        let bases = &rec.read_bases;
        let quals = &rec.base_qualities;
        let oq = original_base_qualities(rec);

        // The loop bound is the *qualities'* length while the no-call test indexes the bases.
        // htsjdk would throw on a record whose SEQ is `*` and whose QUAL is not; reproduced by
        // the bounds check here, which is the same outcome without the stack trace.
        for (i, &qual) in quals.iter().enumerate() {
            let is_call = match bases.get(i) {
                Some(&b) => !is_no_call(b),
                None => break,
            };
            if self.options.include_no_calls || is_call {
                self.q_counts[qual as usize] += 1;
                if let Some(oq) = &oq {
                    if let Some(&v) = oq.get(i) {
                        self.oq_counts[v as usize] += 1;
                    }
                }
            }
        }
    }

    /// `finish`: the `QUALITY` histogram of `COUNT_OF_Q`, and `COUNT_OF_OQ` when non-empty.
    ///
    /// Returns `(bins_of_q, bins_of_oq)` as `(quality, count)` pairs in ascending quality, which
    /// is what `Histogram<Byte>`'s natural ordering gives. An empty OQ list means the second
    /// histogram is not added to the file at all, which is a different output shape and not an
    /// empty column.
    pub fn finish(&self) -> (Bins, Bins) {
        let collect = |counts: &[i64; 128]| -> Bins {
            counts
                .iter()
                .enumerate()
                .filter(|(_, &c)| c > 0)
                .map(|(i, &c)| (i as u8, c as f64))
                .collect()
        };
        (collect(&self.q_counts), collect(&self.oq_counts))
    }
}

/// One histogram's bins: `(quality, count)` in ascending quality.
pub type Bins = Vec<(u8, f64)>;

/// `SAMRecord.getOriginalBaseQualities`: the `OQ` tag, decoded from ASCII by subtracting 33.
///
/// Returns `None` when the tag is absent, which is what makes the second histogram appear or not.
pub fn original_base_qualities(rec: &BamRecord) -> Option<Vec<u8>> {
    match rec.tags.get(Tag::new(b"OQ")) {
        Some(TagValue::Str(s)) => Some(s.bytes().map(|b| b.wrapping_sub(33)).collect()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
    use htsjdk_bam::tag::Tags;

    fn rec(bases: &[u8], quals: &[u8], flags: u16, oq: Option<&str>) -> BamRecord {
        let mut tags = Tags::new();
        if let Some(oq) = oq {
            tags.insert(Tag::new(b"OQ"), TagValue::Str(oq.to_string()));
        }
        BamRecord {
            read_name: "r".to_string(),
            flags,
            reference_index: 0,
            alignment_start: 1,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                length: bases.len() as u32,
                op: Op::M,
            }]),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: bases.to_vec(),
            base_qualities: quals.to_vec(),
            tags,
        }
    }

    #[test]
    fn qualities_are_counted_per_value() {
        let mut d = QualityScoreDistribution::new(Options::default());
        d.accept(&rec(b"ACGT", &[30, 30, 20, 10], 0, None));
        let (q, oq) = d.finish();
        assert_eq!(q, [(10, 1.0), (20, 1.0), (30, 2.0)]);
        assert!(oq.is_empty(), "no OQ tag, no second histogram");
    }

    /// The no-call rule looks at the base, so a good quality at an `N` is discarded.
    #[test]
    fn no_call_bases_are_excluded_by_default() {
        let mut d = QualityScoreDistribution::new(Options::default());
        d.accept(&rec(b"ANNT", &[30, 5, 6, 30], 0, None));
        assert_eq!(
            d.finish().0,
            [(30, 2.0)],
            "the 5 and 6 are gone with the Ns"
        );

        let mut d = QualityScoreDistribution::new(Options {
            include_no_calls: true,
            ..Options::default()
        });
        d.accept(&rec(b"ANNT", &[30, 5, 6, 30], 0, None));
        assert_eq!(d.finish().0, [(5, 1.0), (6, 1.0), (30, 2.0)]);
    }

    #[test]
    fn secondary_and_supplementary_are_dropped() {
        let mut d = QualityScoreDistribution::new(Options::default());
        d.accept(&rec(b"ACGT", &[30; 4], SECONDARY, None));
        d.accept(&rec(b"ACGT", &[31; 4], SUPPLEMENTARY, None));
        assert!(d.finish().0.is_empty());
    }

    /// Vendor-failed and unmapped reads are kept by default, which most collectors here do not
    /// do, and dropped only on request.
    #[test]
    fn failed_and_unmapped_reads_are_kept_unless_asked_otherwise() {
        let mut d = QualityScoreDistribution::new(Options::default());
        d.accept(&rec(b"ACGT", &[11; 4], VENDOR_FAILED, None));
        d.accept(&rec(b"ACGT", &[12; 4], READ_UNMAPPED, None));
        assert_eq!(d.finish().0, [(11, 4.0), (12, 4.0)]);

        let mut d = QualityScoreDistribution::new(Options {
            pf_reads_only: true,
            ..Options::default()
        });
        d.accept(&rec(b"ACGT", &[11; 4], VENDOR_FAILED, None));
        d.accept(&rec(b"ACGT", &[12; 4], READ_UNMAPPED, None));
        assert_eq!(d.finish().0, [(12, 4.0)]);

        let mut d = QualityScoreDistribution::new(Options {
            aligned_reads_only: true,
            ..Options::default()
        });
        d.accept(&rec(b"ACGT", &[11; 4], VENDOR_FAILED, None));
        d.accept(&rec(b"ACGT", &[12; 4], READ_UNMAPPED, None));
        assert_eq!(d.finish().0, [(11, 4.0)]);
    }

    /// The OQ tag is ASCII-33, and its keys need not overlap the primary qualities at all.
    #[test]
    fn the_oq_histogram_has_its_own_keys() {
        let mut d = QualityScoreDistribution::new(Options::default());
        // Qualities 2, 3, 40, 41 encoded as ASCII + 33.
        let oq: String = [2u8, 3, 40, 41].iter().map(|&v| (v + 33) as char).collect();
        d.accept(&rec(b"ACGT", &[30; 4], 0, Some(&oq)));
        let (q, oq_bins) = d.finish();
        assert_eq!(q, [(30, 4.0)]);
        assert_eq!(oq_bins, [(2, 1.0), (3, 1.0), (40, 1.0), (41, 1.0)]);
        // The union the writer must produce puts 30 between 3 and 40, which a concatenation of
        // the two key lists would not.
        let mut union: Vec<u8> = q.iter().chain(&oq_bins).map(|(k, _)| *k).collect();
        union.sort();
        assert_eq!(union, [2, 3, 30, 40, 41]);
    }

    #[test]
    fn quality_zero_is_a_bin_like_any_other() {
        let mut d = QualityScoreDistribution::new(Options::default());
        d.accept(&rec(b"ACGT", &[0, 0, 30, 30], 0, None));
        assert_eq!(d.finish().0, [(0, 2.0), (30, 2.0)]);
    }
}
