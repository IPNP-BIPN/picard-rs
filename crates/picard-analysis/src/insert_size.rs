//! `CollectInsertSizeMetrics`.
//!
//! Ported from `picard.analysis.CollectInsertSizeMetrics`, `InsertSizeMetrics` and
//! `picard.analysis.directed.InsertSizeMetricsCollector`, tag 3.4.0.
//!
//! First member of the `histogram` + `multi_level` + `single_pass` stratum, chosen by the
//! measurement in `docs/decisions/0001`. It is ported at full price so `CollectRnaSeqMetrics`
//! and `CollectAlignmentSummaryMetrics` can be measured against it.

use htsjdk_bam::record::BamRecord;
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricBean, Value};
use htsjdk_metrics::histogram::Histogram;

/// `SamPairUtil.PairOrientation`, in declaration order.
///
/// The order matters twice over. `InsertSizeMetricsCollector` stores the histograms in an
/// `EnumMap` and **puts them in FR, TANDEM, RF order** while the map iterates in *declaration*
/// order, FR, RF, TANDEM. So the rows of the output file are not in the order the code writes
/// them, and a port that preserved insertion order would emit them transposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PairOrientation {
    Fr,
    Rf,
    Tandem,
}

impl PairOrientation {
    pub const DECLARATION_ORDER: [PairOrientation; 3] = [
        PairOrientation::Fr,
        PairOrientation::Rf,
        PairOrientation::Tandem,
    ];

    /// `FormatUtil.format(Enum)` is `value.name()`.
    pub fn name(self) -> &'static str {
        match self {
            PairOrientation::Fr => "FR",
            PairOrientation::Rf => "RF",
            PairOrientation::Tandem => "TANDEM",
        }
    }

    /// The suffix `PerUnitInsertSizeMetricsCollector` gives the histogram's value label.
    fn label_suffix(self) -> &'static str {
        match self {
            PairOrientation::Fr => "fr_count",
            PairOrientation::Rf => "rf_count",
            PairOrientation::Tandem => "tandem_count",
        }
    }
}

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// `SamPairUtil.getPairOrientation`.
///
/// Two reads on the same strand are `TANDEM`. Otherwise the five-prime position of the forward
/// read is compared against that of the reverse read, and the *reverse* read's position is
/// derived from the alignment end or from the insert size depending on which mate is in hand.
pub fn pair_orientation(rec: &BamRecord) -> PairOrientation {
    let read_reverse = rec.flags & READ_REVERSE != 0;
    let mate_reverse = rec.flags & MATE_REVERSE != 0;
    if read_reverse == mate_reverse {
        return PairOrientation::Tandem;
    }
    let positive_five_prime = if read_reverse {
        rec.mate_alignment_start as i64
    } else {
        rec.alignment_start as i64
    };
    let negative_five_prime = if read_reverse {
        rec.alignment_end() as i64
    } else {
        rec.alignment_start as i64 + rec.inferred_insert_size as i64
    };
    if positive_five_prime < negative_five_prime {
        PairOrientation::Fr
    } else {
        PairOrientation::Rf
    }
}

/// `InsertSizeMetrics`, with the column order HotSpot's `Class.getFields()` produces.
///
/// `InsertSizeMetrics extends MultilevelMetrics`, and the golden shows the **declared** fields
/// come first, with the inherited `SAMPLE`, `LIBRARY` and `READ_GROUP` last. `getFields()` is
/// documented as returning fields in no particular order, so like the ordering in decision
/// 0009 this is a property of the reference implementation rather than of the language.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InsertSizeMetrics {
    pub median_insert_size: f64,
    pub mode_insert_size: f64,
    pub median_absolute_deviation: f64,
    pub min_insert_size: i32,
    pub max_insert_size: i32,
    pub mean_insert_size: f64,
    pub standard_deviation: f64,
    pub read_pairs: i64,
    pub pair_orientation: Option<PairOrientation>,
    pub width_of_10_percent: i32,
    pub width_of_20_percent: i32,
    pub width_of_30_percent: i32,
    pub width_of_40_percent: i32,
    pub width_of_50_percent: i32,
    pub width_of_60_percent: i32,
    pub width_of_70_percent: i32,
    pub width_of_80_percent: i32,
    pub width_of_90_percent: i32,
    pub width_of_95_percent: i32,
    pub width_of_99_percent: i32,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

const COLUMNS: &[&str] = &[
    "MEDIAN_INSERT_SIZE",
    "MODE_INSERT_SIZE",
    "MEDIAN_ABSOLUTE_DEVIATION",
    "MIN_INSERT_SIZE",
    "MAX_INSERT_SIZE",
    "MEAN_INSERT_SIZE",
    "STANDARD_DEVIATION",
    "READ_PAIRS",
    "PAIR_ORIENTATION",
    "WIDTH_OF_10_PERCENT",
    "WIDTH_OF_20_PERCENT",
    "WIDTH_OF_30_PERCENT",
    "WIDTH_OF_40_PERCENT",
    "WIDTH_OF_50_PERCENT",
    "WIDTH_OF_60_PERCENT",
    "WIDTH_OF_70_PERCENT",
    "WIDTH_OF_80_PERCENT",
    "WIDTH_OF_90_PERCENT",
    "WIDTH_OF_95_PERCENT",
    "WIDTH_OF_99_PERCENT",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

impl MetricBean for InsertSizeMetrics {
    fn class_name(&self) -> &str {
        "picard.analysis.InsertSizeMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }

    fn values(&self) -> Vec<Value> {
        let text = |o: &Option<String>| match o {
            Some(s) => Value::Str(s.clone()),
            None => Value::Null,
        };
        vec![
            Value::Double(self.median_insert_size),
            Value::Double(self.mode_insert_size),
            Value::Double(self.median_absolute_deviation),
            Value::Long(self.min_insert_size as i64),
            Value::Long(self.max_insert_size as i64),
            Value::Double(self.mean_insert_size),
            Value::Double(self.standard_deviation),
            Value::Long(self.read_pairs),
            match self.pair_orientation {
                Some(p) => Value::Str(p.name().to_string()),
                None => Value::Null,
            },
            Value::Long(self.width_of_10_percent as i64),
            Value::Long(self.width_of_20_percent as i64),
            Value::Long(self.width_of_30_percent as i64),
            Value::Long(self.width_of_40_percent as i64),
            Value::Long(self.width_of_50_percent as i64),
            Value::Long(self.width_of_60_percent as i64),
            Value::Long(self.width_of_70_percent as i64),
            Value::Long(self.width_of_80_percent as i64),
            Value::Long(self.width_of_90_percent as i64),
            Value::Long(self.width_of_95_percent as i64),
            Value::Long(self.width_of_99_percent as i64),
            text(&self.sample),
            text(&self.library),
            text(&self.read_group),
        ]
    }
}

/// The tool's arguments, with Picard's defaults.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `MINIMUM_PCT`, default 0.05. Declared `float` in Java, so the comparison it drives runs
    /// at `f32` precision widened to `f64`.
    pub minimum_pct: f32,
    /// `DEVIATIONS`, default 10.
    pub deviations: f64,
    pub histogram_width: Option<i32>,
    pub min_histogram_width: Option<i32>,
    pub include_duplicates: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            minimum_pct: 0.05,
            deviations: 10.0,
            histogram_width: None,
            min_histogram_width: None,
            include_duplicates: false,
        }
    }
}

/// `InsertSizeMetricsCollector` restricted to the ALL_READS level, which is the default.
pub struct InsertSizeMetricsCollector {
    options: Options,
    histograms: [Histogram; 3],
}

impl InsertSizeMetricsCollector {
    pub fn new(options: Options) -> Self {
        // The prefix is "All_Reads." when no sample, library or read group is set, which is
        // the only level the default METRIC_ACCUMULATION_LEVEL produces.
        let mk = |o: PairOrientation| {
            Histogram::new("insert_size", &format!("All_Reads.{}", o.label_suffix()))
        };
        InsertSizeMetricsCollector {
            options,
            histograms: [
                mk(PairOrientation::Fr),
                mk(PairOrientation::Rf),
                mk(PairOrientation::Tandem),
            ],
        }
    }

    /// `InsertSizeMetricsCollector.acceptRecord`.
    ///
    /// The filter that surprises: `getFirstOfPairFlag()` returning true is a **rejection**, so
    /// only the second read of each pair is counted. Counting both would double every number
    /// while leaving the distribution's shape, the median and every width intact, so it would
    /// look right.
    pub fn accept(&mut self, rec: &BamRecord) {
        if rec.flags & READ_PAIRED == 0
            || rec.flags & READ_UNMAPPED != 0
            || rec.flags & MATE_UNMAPPED != 0
            || rec.flags & FIRST_OF_PAIR != 0
            || rec.flags & (SECONDARY | SUPPLEMENTARY) != 0
            || (rec.flags & DUPLICATE != 0 && !self.options.include_duplicates)
            || rec.inferred_insert_size == 0
        {
            return;
        }
        let insert_size = rec.inferred_insert_size.unsigned_abs() as f64;
        let orientation = pair_orientation(rec);
        self.histograms[orientation as usize].increment(insert_size);
    }

    /// `PerUnitInsertSizeMetricsCollector.addMetricsToFile`.
    ///
    /// Returns one metric and one histogram per orientation that clears `MINIMUM_PCT`, in
    /// `EnumMap` order.
    pub fn finish(mut self) -> Vec<(InsertSizeMetrics, OutHistogram)> {
        let total_inserts: f64 = self.histograms.iter().map(|h| h.count()).sum();
        if total_inserts == 0.0 {
            return Vec::new();
        }

        let mut out = Vec::new();
        for orientation in PairOrientation::DECLARATION_ORDER {
            let histogram = &mut self.histograms[orientation as usize];
            let total = histogram.count();
            // The threshold is a `float` in Java, widened for the comparison.
            if total < total_inserts * self.options.minimum_pct as f64 {
                continue;
            }

            let mut m = InsertSizeMetrics {
                pair_orientation: Some(orientation),
                ..Default::default()
            };

            if !histogram.is_empty() {
                m.read_pairs = total as i64;
                m.max_insert_size = histogram.max().unwrap_or(0.0) as i32;
                m.min_insert_size = histogram.min().unwrap_or(0.0) as i32;
                m.median_insert_size = histogram.median();
                m.mode_insert_size = histogram.mode().unwrap_or(0.0);
                m.median_absolute_deviation = histogram.median_absolute_deviation();

                widths(histogram, total, &mut m);
            }

            // `trimByWidth` mutates the histogram in place, and the *trimmed* one is what is
            // written to the file. The mean and standard deviation are taken after the trim,
            // so they describe the core of the distribution rather than all of it.
            histogram.trim_by_width(width_to_trim_to(&m, &self.options));
            if !histogram.is_empty() {
                m.mean_insert_size = histogram.mean();
                m.standard_deviation = histogram.standard_deviation();
            }

            let bins = histogram
                .bins()
                .map(|(k, v)| (format!("{}", k as i64), v))
                .collect();
            out.push((
                m,
                OutHistogram {
                    bin_label: histogram.bin_label.clone(),
                    value_label: histogram.value_label.clone(),
                    key_class: "java.lang.Integer".to_string(),
                    bins,
                },
            ));
        }
        out
    }
}

/// The widening scan around the median that fills the `WIDTH_OF_n_PERCENT` fields.
///
/// It walks outward one bin at a time from the median, accumulating coverage, and records the
/// first distance at which each threshold is reached. `distance` is `(int)(high - low) + 1`,
/// so it counts bins inclusive of both ends, and the loop starts with `low == high` at the
/// median, giving a minimum width of 1 rather than 0.
fn widths(histogram: &Histogram, total: f64, m: &mut InsertSizeMetrics) {
    const THRESHOLDS: [f64; 11] = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99];
    let mut widths = [0i32; 11];
    let median = histogram.median();
    let (min, max) = (
        histogram.min().unwrap_or(0.0),
        histogram.max().unwrap_or(0.0),
    );
    let mut covered = 0.0;
    let mut low = median;
    let mut high = median;

    while low >= min - 1.0 || high <= max + 1.0 {
        // `histogram.get((int) low)` truncates toward zero, so a median of 300.5 probes 300.
        if let Some(v) = histogram.get(low as i64 as f64) {
            covered += v;
        }
        if low != high {
            if let Some(v) = histogram.get(high as i64 as f64) {
                covered += v;
            }
        }
        let percent_covered = covered / total;
        let distance = (high - low) as i32 + 1;

        // The eleven near-identical Java lines, kept as a table so the thresholds stay
        // readable side by side. Each field is set once, at the first distance that reaches
        // its threshold, which is what `== 0` guards.
        for (threshold, slot) in THRESHOLDS.iter().zip(widths.iter_mut()) {
            if percent_covered >= *threshold && *slot == 0 {
                *slot = distance;
            }
        }

        low -= 1.0;
        high += 1.0;
    }

    let [w10, w20, w30, w40, w50, w60, w70, w80, w90, w95, w99] = widths;
    m.width_of_10_percent = w10;
    m.width_of_20_percent = w20;
    m.width_of_30_percent = w30;
    m.width_of_40_percent = w40;
    m.width_of_50_percent = w50;
    m.width_of_60_percent = w60;
    m.width_of_70_percent = w70;
    m.width_of_80_percent = w80;
    m.width_of_90_percent = w90;
    m.width_of_95_percent = w95;
    m.width_of_99_percent = w99;
}

/// `PerUnitInsertSizeMetricsCollector.getWidthToTrimTo`.
fn width_to_trim_to(m: &InsertSizeMetrics, options: &Options) -> i32 {
    match options.histogram_width {
        Some(w) => w,
        None => {
            let calculated =
                (m.median_insert_size + options.deviations * m.median_absolute_deviation) as i32;
            match options.min_histogram_width {
                Some(min) => min.max(calculated),
                None => calculated,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

    fn pair(start: i32, insert: i32, orientation: PairOrientation, extra: u16) -> Vec<BamRecord> {
        let len = 50u32;
        (0..2)
            .map(|mate| {
                let first = mate == 0;
                let (first_rev, second_rev) = match orientation {
                    PairOrientation::Fr => (false, true),
                    PairOrientation::Rf => (true, false),
                    PairOrientation::Tandem => (false, false),
                };
                let this_rev = if first { first_rev } else { second_rev };
                let mate_rev = if first { second_rev } else { first_rev };
                let mut flags = READ_PAIRED | 0x2 | extra;
                flags |= if first { FIRST_OF_PAIR } else { 0x80 };
                if this_rev {
                    flags |= READ_REVERSE;
                }
                if mate_rev {
                    flags |= MATE_REVERSE;
                }
                let right = start + insert - len as i32;
                BamRecord {
                    read_name: "p".into(),
                    flags,
                    reference_index: 0,
                    alignment_start: if first { start } else { right },
                    mapping_quality: 60,
                    cigar: Cigar::new(vec![CigarElement {
                        length: len,
                        op: Op::M,
                    }]),
                    mate_reference_index: 0,
                    mate_alignment_start: if first { right } else { start },
                    inferred_insert_size: if first { insert } else { -insert },
                    read_bases: vec![b'A'; len as usize],
                    base_qualities: vec![30; len as usize],
                    tags: Default::default(),
                }
            })
            .collect()
    }

    fn collect(records: &[BamRecord], options: Options) -> Vec<(InsertSizeMetrics, OutHistogram)> {
        let mut c = InsertSizeMetricsCollector::new(options);
        for r in records {
            c.accept(r);
        }
        c.finish()
    }

    /// Only the second read of a pair is counted. Counting both would double every number and
    /// leave the shape, the median and every width intact, so it would look right.
    #[test]
    fn only_the_second_read_of_a_pair_is_counted() {
        let records = pair(1000, 300, PairOrientation::Fr, 0);
        assert_eq!(records.len(), 2);
        let out = collect(&records, Options::default());
        assert_eq!(out[0].0.read_pairs, 1, "one pair, not two reads");
    }

    #[test]
    fn the_insert_size_is_the_absolute_value() {
        let out = collect(&pair(1000, 300, PairOrientation::Fr, 0), Options::default());
        assert_eq!(out[0].0.median_insert_size, 300.0);
        assert_eq!(out[0].0.min_insert_size, 300);
    }

    #[test]
    fn duplicates_are_excluded_by_default_and_included_on_request() {
        let mut records = pair(1000, 300, PairOrientation::Fr, 0);
        records.extend(pair(2000, 900, PairOrientation::Fr, DUPLICATE));
        assert_eq!(collect(&records, Options::default())[0].0.read_pairs, 1);
        assert_eq!(
            collect(
                &records,
                Options {
                    include_duplicates: true,
                    ..Options::default()
                }
            )[0]
            .0
            .read_pairs,
            2
        );
    }

    #[test]
    fn secondary_and_supplementary_records_are_skipped() {
        for flag in [SECONDARY, SUPPLEMENTARY] {
            let records = pair(1000, 300, PairOrientation::Fr, flag);
            assert!(collect(&records, Options::default()).is_empty());
        }
    }

    /// The orientations come out in enum declaration order, FR then RF then TANDEM, which is
    /// not the order the Java puts them into the map.
    #[test]
    fn orientations_are_emitted_in_declaration_order() {
        let mut records = Vec::new();
        for i in 0..30 {
            records.extend(pair(1000 + i * 500, 300, PairOrientation::Tandem, 0));
            records.extend(pair(2000 + i * 500, 300, PairOrientation::Rf, 0));
            records.extend(pair(3000 + i * 500, 300, PairOrientation::Fr, 0));
        }
        let out = collect(&records, Options::default());
        let seen: Vec<&str> = out
            .iter()
            .map(|(m, _)| m.pair_orientation.unwrap().name())
            .collect();
        assert_eq!(seen, vec!["FR", "RF", "TANDEM"]);
    }

    /// An orientation holding less than MINIMUM_PCT of the data is dropped entirely.
    #[test]
    fn a_rare_orientation_is_dropped_below_the_minimum_percentage() {
        let mut records = Vec::new();
        for i in 0..100 {
            records.extend(pair(1000 + i * 500, 300, PairOrientation::Fr, 0));
        }
        records.extend(pair(900_000, 500, PairOrientation::Rf, 0));

        let out = collect(&records, Options::default());
        assert_eq!(out.len(), 1, "1 of 101 is below the 5% default");

        // With MINIMUM_PCT=0 the test `total >= totalInserts * 0` passes for every
        // orientation, including TANDEM with a count of zero, so all three produce a metric
        // row. The `isEmpty()` guard inside only skips filling the statistics; it does not
        // skip `addMetric`. So an empty orientation yields a row of zeros.
        let all = collect(
            &records,
            Options {
                minimum_pct: 0.0,
                ..Options::default()
            },
        );
        assert_eq!(
            all.len(),
            3,
            "MINIMUM_PCT=0 keeps even the empty orientation"
        );
        assert_eq!(all[2].0.read_pairs, 0, "TANDEM saw nothing");
        assert!(all[2].1.bins.is_empty(), "and its histogram is empty");
    }

    /// The width scan starts at the median with low == high, so the narrowest width is 1.
    #[test]
    fn the_narrowest_width_is_one_bin_not_zero() {
        let mut records = Vec::new();
        for i in 0..50 {
            records.extend(pair(1000 + i * 500, 300, PairOrientation::Fr, 0));
        }
        let out = collect(&records, Options::default());
        assert_eq!(out[0].0.width_of_10_percent, 1);
        assert_eq!(out[0].0.width_of_99_percent, 1, "all mass is in one bin");
    }

    /// The mean and standard deviation are taken *after* the trim, so a long tail does not
    /// drag them. The trimmed histogram is also what is written to the file.
    #[test]
    fn the_mean_is_computed_after_the_trim() {
        let mut records = Vec::new();
        for i in 0..200 {
            records.extend(pair(1000 + i * 500, 300, PairOrientation::Fr, 0));
        }
        // One extreme outlier, far beyond median + 10 * MAD.
        records.extend(pair(900_000, 50_000, PairOrientation::Fr, 0));

        let out = collect(&records, Options::default());
        let (m, h) = &out[0];
        assert_eq!(m.max_insert_size, 50_000, "MAX is from the untrimmed data");
        assert_eq!(
            m.mean_insert_size, 300.0,
            "the mean is from the trimmed data, so the outlier is gone"
        );
        assert!(
            !h.bins.iter().any(|(k, _)| k == "50000"),
            "the written histogram is the trimmed one"
        );
    }

    #[test]
    fn an_explicit_histogram_width_overrides_the_computed_one() {
        let mut records = Vec::new();
        for i in 0..200 {
            records.extend(pair(1000 + i * 500, 300, PairOrientation::Fr, 0));
        }
        records.extend(pair(900_000, 5_000, PairOrientation::Fr, 0));
        let out = collect(
            &records,
            Options {
                histogram_width: Some(10_000),
                ..Options::default()
            },
        );
        assert!(
            out[0].1.bins.iter().any(|(k, _)| k == "5000"),
            "a width of 10000 keeps the 5000 bin"
        );
    }

    #[test]
    fn an_empty_input_produces_nothing() {
        assert!(collect(&[], Options::default()).is_empty());
    }

    #[test]
    fn the_column_order_puts_inherited_fields_last() {
        let m = InsertSizeMetrics::default();
        assert_eq!(m.columns()[0], "MEDIAN_INSERT_SIZE");
        assert_eq!(
            &m.columns()[20..],
            &["SAMPLE", "LIBRARY", "READ_GROUP"],
            "MultilevelMetrics' fields are inherited and come last"
        );
        assert_eq!(m.columns().len(), m.values().len());
    }

    #[test]
    fn a_null_multilevel_field_is_the_empty_string() {
        let m = InsertSizeMetrics::default();
        let v = m.values();
        assert_eq!(v[20].format(), "");
        assert_eq!(v[21].format(), "");
        assert_eq!(v[22].format(), "");
    }
}
