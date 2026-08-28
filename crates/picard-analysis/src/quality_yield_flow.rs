//! `CollectQualityYieldMetricsFlow`: the flow-based cousin of `CollectQualityYieldMetrics`, which
//! counts flows rather than bases.
//!
//! Reading the file and turning a read into a flow-based one are not ported: the key and the
//! per-flow error probabilities come out of the tp and t0 matrices. What is ported is which reads
//! reach the tally, how a flow's quality is derived from its error probability, and the metrics
//! the tally produces.
//!
//! Ported from `picard.analysis.CollectQualityYieldMetricsFlow` in Picard 3.4.0.

/// `CollectQualityYieldMetricsFlow.MIN_QUAL`.
pub const MIN_QUAL: i64 = 0;
/// `CollectQualityYieldMetricsFlow.MAX_QUAL`.
pub const MAX_QUAL: i64 = 100;
/// `CollectQualityYieldMetricsFlow.CYCLE_LENGTH`, the flows one histogram cycle holds.
pub const CYCLE_LENGTH: usize = 4;
/// `acceptRecord`, on a read whose read group is not a flow platform.
pub const NOT_A_FLOW_PLATFORM_MESSAGE: &str = "Reads should originate from a flow based platform";

/// One record, reduced to what `acceptRecord` reads before the flow conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Record<'a> {
    pub read_length: usize,
    pub secondary: bool,
    pub supplementary: bool,
    pub fails_vendor_quality: bool,
    /// The flow qualities, which the reference derives from the read's error probabilities.
    pub flow_qualities: &'a [u8],
}

/// What `acceptRecord` does with one record before it looks at any flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// Returned before `TOTAL_READS` is touched, so the read is counted NOWHERE.
    Skipped,
    /// Counted in `TOTAL_READS` and nowhere else, its flows included.
    Total,
    /// Counted everywhere.
    Counted,
}

/// The three early returns of `acceptRecord`, in their own order.
///
/// A read of no bases, a secondary read that is not included and a supplementary read that is not
/// included all return BEFORE `TOTAL_READS` is incremented, so they are counted nowhere at all. A
/// read that fails vendor quality returns AFTER it, so it is counted in `TOTAL_READS` and left out
/// of `PF_READS`. The two include arguments are independent: naming one leaves the other's
/// records out.
pub fn outcome(record: &Record, include_secondary: bool, include_supplementary: bool) -> Outcome {
    if record.read_length == 0
        || (!include_secondary && record.secondary)
        || (!include_supplementary && record.supplementary)
    {
        return Outcome::Skipped;
    }
    if record.fails_vendor_quality {
        return Outcome::Total;
    }
    Outcome::Counted
}

/// `getFlowQualities`' per-flow conversion: the phred of the error probability, rounded and
/// clamped.
///
/// A probability of exactly zero is the special case, and it answers `MAX_QUAL` rather than the
/// infinity the logarithm would give.
pub fn flow_quality(error_probability: f64) -> u8 {
    if error_probability == 0.0 {
        return MAX_QUAL as u8;
    }
    let q = (-10.0 * error_probability.log10()).round() as i64;
    q.clamp(MIN_QUAL, MAX_QUAL) as u8
}

/// The metrics one run writes, before the derived fields are filled in.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Metrics {
    pub total_reads: i64,
    pub pf_reads: i64,
    pub pf_flows: i64,
    pub pf_q20_flows: i64,
    pub pf_q30_flows: i64,
    /// The SUM of the flow qualities while the tally runs, divided by twenty at the end.
    pub pf_q20_equivalent_yield: i64,
}

impl Metrics {
    /// `MEAN_PF_READ_NUMBER_OF_FLOWS`, an INTEGER division that truncates.
    pub fn mean_pf_read_number_of_flows(&self) -> i32 {
        if self.pf_reads == 0 {
            0
        } else {
            (self.pf_flows / self.pf_reads) as i32
        }
    }

    /// `PCT_PF_Q20_FLOWS`, which is a fraction and not a percentage.
    pub fn pct_pf_q20_flows(&self) -> f64 {
        if self.pf_flows == 0 {
            0.0
        } else {
            self.pf_q20_flows as f64 / self.pf_flows as f64
        }
    }

    /// `PCT_PF_Q30_FLOWS`, likewise.
    pub fn pct_pf_q30_flows(&self) -> f64 {
        if self.pf_flows == 0 {
            0.0
        } else {
            self.pf_q30_flows as f64 / self.pf_flows as f64
        }
    }
}

/// The whole tally: `acceptRecord` over every record, then the division `finish` does.
///
/// `PF_Q20_FLOWS` counts every flow at 20 or over, the 30s included, because the branch that
/// increments `PF_Q30_FLOWS` increments it too. It is therefore never smaller than
/// `PF_Q30_FLOWS`. `PF_Q20_EQUIVALENT_YIELD` is not a count of anything: it is the sum of the
/// qualities divided by twenty, so it moves with the qualities and not only with the flows.
pub fn collect(
    records: &[Record],
    include_secondary: bool,
    include_supplementary: bool,
) -> Metrics {
    let mut metrics = Metrics::default();
    for record in records {
        match outcome(record, include_secondary, include_supplementary) {
            Outcome::Skipped => continue,
            Outcome::Total => {
                metrics.total_reads += 1;
            }
            Outcome::Counted => {
                metrics.total_reads += 1;
                metrics.pf_reads += 1;
                metrics.pf_flows += record.flow_qualities.len() as i64;
                for quality in record.flow_qualities {
                    let quality = i64::from(*quality);
                    metrics.pf_q20_equivalent_yield += quality;
                    if quality >= 30 {
                        metrics.pf_q20_flows += 1;
                        metrics.pf_q30_flows += 1;
                    } else if quality >= 20 {
                        metrics.pf_q20_flows += 1;
                    }
                }
            }
        }
    }
    metrics.pf_q20_equivalent_yield /= 20;
    metrics
}

/// `Math.ceil((float) quals.length / CYCLE_LENGTH)`, the histogram's cycle count for one read.
pub fn cycle_count(flows: usize) -> usize {
    flows.div_ceil(CYCLE_LENGTH)
}
