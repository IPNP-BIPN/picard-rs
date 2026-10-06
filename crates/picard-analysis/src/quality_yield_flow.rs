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
    // `Math.round(-10 * Math.log10(p))`: the Java rounding and the Java logarithm.
    let q = jmath::math::round(-10.0 * jmath::math::log10(error_probability));
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

// ---------------------------------------------------------------------------------------------
// The tool itself: reads in, the metrics and the histograms out.
// ---------------------------------------------------------------------------------------------

use std::collections::{BTreeMap, HashMap};

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricBean, Value};

use crate::flow_based::{FlowArguments, FlowBasedRead, FlowReadGroupInfo};
use crate::series_stats::SeriesStats;

const COLUMNS: &[&str] = &[
    "TOTAL_READS",
    "PF_READS",
    "MEAN_PF_READ_NUMBER_OF_FLOWS",
    "PF_FLOWS",
    "PF_Q20_FLOWS",
    "PCT_PF_Q20_FLOWS",
    "PF_Q30_FLOWS",
    "PCT_PF_Q30_FLOWS",
    "PF_Q20_EQUIVALENT_YIELD",
];

impl MetricBean for Metrics {
    fn class_name(&self) -> &str {
        "picard.analysis.CollectQualityYieldMetricsFlow$QualityYieldMetricsFlow"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Long(self.total_reads),
            Value::Long(self.pf_reads),
            Value::Long(i64::from(self.mean_pf_read_number_of_flows())),
            Value::Long(self.pf_flows),
            Value::Long(self.pf_q20_flows),
            Value::Double(self.pct_pf_q20_flows()),
            Value::Long(self.pf_q30_flows),
            Value::Double(self.pct_pf_q30_flows()),
            Value::Long(self.pf_q20_equivalent_yield),
        ]
    }
}

/// `QualityYieldMetricsCollectorFlow`, with the histogram state the tool keeps beside it.
pub struct FlowCollector {
    include_secondary: bool,
    include_supplementary: bool,
    include_histogram: bool,
    arguments: FlowArguments,
    metrics: Metrics,
    quality_histogram: BTreeMap<i64, i64>,
    flow_quality_stats: Vec<SeriesStats>,
    /// `FlowBasedKeyCodec.readGroupInfo`, filled on first use of each read group.
    infos: HashMap<String, FlowReadGroupInfo>,
}

impl FlowCollector {
    pub fn new(
        include_secondary: bool,
        include_supplementary: bool,
        include_histogram: bool,
        arguments: FlowArguments,
    ) -> Self {
        FlowCollector {
            include_secondary,
            include_supplementary,
            include_histogram,
            arguments,
            metrics: Metrics::default(),
            quality_histogram: BTreeMap::new(),
            flow_quality_stats: Vec::new(),
            infos: HashMap::new(),
        }
    }

    /// `acceptRecord`. An `Err` is the Java exception, class and message, as the run prints it.
    pub fn accept(&mut self, header: &SamHeader, record: &BamRecord) -> Result<(), String> {
        let view = Record {
            read_length: record.read_length(),
            secondary: record.flags & 0x100 != 0,
            supplementary: record.flags & 0x800 != 0,
            fails_vendor_quality: record.flags & 0x200 != 0,
            flow_qualities: &[],
        };
        match outcome(&view, self.include_secondary, self.include_supplementary) {
            Outcome::Skipped => return Ok(()),
            Outcome::Total => {
                self.metrics.total_reads += 1;
                return Ok(());
            }
            Outcome::Counted => self.metrics.total_reads += 1,
        }

        let group_id = match record.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(id)) => id.clone(),
            _ => String::new(),
        };
        if !self.infos.contains_key(&group_id) {
            let group = header
                .read_groups
                .iter()
                .find(|g| g.id == group_id)
                .ok_or_else(|| {
                    "java.lang.NullPointerException: Cannot invoke \
                     \"htsjdk.samtools.SAMReadGroupRecord.getReadGroupId()\" because the return \
                     value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null"
                        .to_string()
                })?;
            self.infos
                .insert(group_id.clone(), FlowReadGroupInfo::new(group)?);
        }
        let info = &self.infos[&group_id];
        if !info.is_flow_platform {
            return Err(format!(
                "picard.PicardException: {NOT_A_FLOW_PLATFORM_MESSAGE}"
            ));
        }
        let flow_order = info
            .flow_order
            .clone()
            .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
        let read = FlowBasedRead::new(record, &flow_order, info.max_class, &self.arguments)?;

        self.metrics.pf_reads += 1;
        self.metrics.pf_flows += read.key().len() as i64;
        let quals = flow_qualities(&read);

        let cycles = if self.include_histogram {
            cycle_count(quals.len())
        } else {
            0
        };
        let mut cycle_qual_count = vec![0i32; cycles];
        let mut cycle_qual_sum = vec![0i32; cycles];
        for (flow, &qual) in quals.iter().enumerate() {
            let qual = i64::from(qual);
            self.metrics.pf_q20_equivalent_yield += qual;
            if qual >= 30 {
                self.metrics.pf_q20_flows += 1;
                self.metrics.pf_q30_flows += 1;
            } else if qual >= 20 {
                self.metrics.pf_q20_flows += 1;
            }
            if self.include_histogram {
                *self.quality_histogram.entry(qual).or_insert(0) += 1;
                let cycle = flow / CYCLE_LENGTH;
                cycle_qual_count[cycle] += 1;
                cycle_qual_sum[cycle] += qual as i32;
            }
        }
        if self.include_histogram {
            while self.flow_quality_stats.len() < cycles {
                self.flow_quality_stats.push(SeriesStats::new());
            }
            let negative = record.flags & 0x10 != 0;
            for cycle in 0..cycles {
                let id = if !negative { cycle } else { cycles - 1 - cycle };
                self.flow_quality_stats[id]
                    .add(f64::from(cycle_qual_sum[cycle]) / f64::from(cycle_qual_count[cycle]));
            }
        }
        Ok(())
    }

    /// `finish` and `addMetricsToFile`: the equivalent yield is divided by twenty here.
    pub fn finish(mut self) -> (Metrics, Vec<Histogram>) {
        self.metrics.pf_q20_equivalent_yield /= 20;
        let mut histograms = Vec::new();
        if self.include_histogram {
            histograms.push(integer_histogram(
                "QUAL_COUNT",
                self.quality_histogram
                    .iter()
                    .map(|(k, v)| (*k, *v as f64))
                    .collect(),
            ));
            let stats = &self.flow_quality_stats;
            let series = |label: &str, value: &dyn Fn(&SeriesStats) -> f64| {
                integer_histogram(
                    label,
                    stats
                        .iter()
                        .enumerate()
                        .map(|(i, s)| (i as i64, value(s)))
                        .collect(),
                )
            };
            histograms.push(series("MEAN_CYCLE_QUAL", &|s| s.mean()));
            histograms.push(series("MEDIAN_CYCLE_QUAL", &|s| s.median()));
            histograms.push(series("Q25_CYCLE_QUAL", &|s| s.percentile(25.0)));
            histograms.push(series("Q75_CYCLE_QUAL", &|s| s.percentile(75.0)));
        }
        (self.metrics, histograms)
    }
}

/// A `Histogram<Integer>` named `KEY` / `label`, bins already in key order.
pub fn integer_histogram(label: &str, bins: Vec<(i64, f64)>) -> Histogram {
    Histogram {
        bin_label: "KEY".to_string(),
        value_label: label.to_string(),
        key_class: "java.lang.Integer".to_string(),
        bins: bins.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
    }
}

/// `getFlowQualities` over `computeErrorProb`.
fn flow_qualities(read: &FlowBasedRead) -> Vec<u8> {
    let key = read.key();
    let max_hmer = read.max_hmer();
    let mut column = vec![0.0f64; max_hmer as usize + 1];
    let mut result = Vec::with_capacity(key.len());
    for (i, &call) in key.iter().enumerate() {
        let mut sum = 0.0;
        for (j, cell) in column.iter_mut().enumerate() {
            *cell = read.prob(i, j as i32);
            sum += *cell;
        }
        for cell in column.iter_mut() {
            *cell /= sum;
        }
        let error = 1.0 - column[call.min(max_hmer) as usize];
        result.push(flow_quality(error));
    }
    result
}
