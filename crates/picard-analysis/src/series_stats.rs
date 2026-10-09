//! `picard.util.SeriesStats`: the running summary the flow and SNVQ yield collectors keep per
//! cycle or per read position.
//!
//! The mean is the sum over the count; the percentile is read off a sorted map of the distinct
//! values with their multiplicities, indexed by `(int)(count * percentile / 100)`. The values are
//! doubles, so the keys are ordered by `total_cmp`, which agrees with `Double.compareTo` on every
//! value these collectors produce (finite, never NaN).
//!
//! Ported from `picard.util.SeriesStats` in Picard 3.4.0.

/// One series.
#[derive(Debug, Clone, Default)]
pub struct SeriesStats {
    last: f64,
    count: i64,
    sum: f64,
    /// The distinct values in ascending order, each with how many times it was added.
    bins: Vec<(f64, i64)>,
}

impl SeriesStats {
    pub fn new() -> Self {
        SeriesStats {
            last: f64::NAN,
            ..Default::default()
        }
    }

    /// `add`.
    pub fn add(&mut self, v: f64) {
        self.last = v;
        self.sum += v;
        self.count += 1;
        match self.bins.binary_search_by(|(k, _)| k.total_cmp(&v)) {
            Ok(at) => self.bins[at].1 += 1,
            Err(at) => self.bins.insert(at, (v, 1)),
        }
    }

    /// `getMean`: NaN for an empty series.
    pub fn mean(&self) -> f64 {
        if self.count != 0 {
            self.sum / self.count as f64
        } else {
            f64::NAN
        }
    }

    /// `getMedian`.
    pub fn median(&self) -> f64 {
        self.percentile(50.0)
    }

    /// `getPercentile`: the last value for a single observation, otherwise the first bin whose
    /// cumulative range holds `(int)(count * percentile / 100)`, else the highest.
    pub fn percentile(&self, percentile: f64) -> f64 {
        if self.count == 0 {
            return f64::NAN;
        }
        if self.count == 1 {
            return self.last;
        }
        let percentile_index = (self.count as f64 * percentile / 100.0) as i64;
        let mut index = 0i64;
        for (value, size) in &self.bins {
            if percentile_index >= index && percentile_index < index + size {
                return *value;
            }
            index += size;
        }
        self.bins.last().map(|(v, _)| *v).unwrap_or(f64::NAN)
    }
}
