//! `picard.util.SeriesStats`: the running summary the flow and SNVQ yield collectors keep per
//! cycle or per read position.
//!
//! The mean is the sum over the count; the percentile is read off a sorted map of the distinct
//! values with their multiplicities, indexed by `(int)(count * precentile / 100)` (the reference's
//! own spelling of the parameter). The values are doubles, so the keys are ordered by `total_cmp`,
//! which agrees with `Double.compareTo` on every value these collectors produce (finite, never
//! NaN).
//!
//! Ported from `picard.util.SeriesStats` in Picard 3.5.0; the file is unchanged since 3.4.0.
//!
//! Of the reference's members, `add`, `getMean`, `getMedian` and `getPercentile` are here;
//! `getMin`, `getMax`, `getStd` and `getUniq` are not.
//!
//! # Reading this file without reading Rust
//!
//! - `f64` is Java's `double` and `i64` Java's `long`.
//! - `x as i64` on an `f64` truncates toward zero and saturates at the type's range, NaN giving
//!   zero, which is what Java's `(int)` / `(long)` narrowing of a `double` does.
//! - `binary_search_by` returns `Ok(index)` when it finds the value and `Err(index)` with the
//!   place it would be inserted when it does not; `match` takes the two cases apart.
//! - `Option` is "a value or nothing"; `map` and `unwrap_or` work on it as on Java's `Optional`.

/// One series.
///
/// Java's `TreeMap<Double, AtomicLong>` becomes a sorted vector of pairs: Rust's `f64` has no
/// total order of its own (NaN), so it cannot key a `BTreeMap`, and a vector kept sorted by
/// `total_cmp` gives the same ascending walk the reference's `entrySet()` does.
#[derive(Debug, Clone, Default)]
pub struct SeriesStats {
    /// `last`, the most recent value, which `getPercentile` returns for a one-value series.
    last: f64,
    /// `count`; an `int` in the reference, wider here, which is the same below 2^31 values.
    count: i64,
    sum: f64,
    /// The distinct values in ascending order, each with how many times it was added.
    bins: Vec<(f64, i64)>,
}

impl SeriesStats {
    /// The reference's field initialisers: `last = Double.NaN`, the rest zero and empty.
    /// `..Default::default()` fills the fields not named from the derived all-zero value.
    pub fn new() -> Self {
        SeriesStats {
            last: f64::NAN,
            ..Default::default()
        }
    }

    /// `add`.
    ///
    /// Picard 3.5.0, `src/main/java/picard/util/SeriesStats.java:44-48`:
    ///
    /// ```java
    /// if ( bins.containsKey(v) ) {
    ///     bins.get(v).incrementAndGet();
    /// } else {
    ///     bins.put(v, new AtomicLong(1));
    /// }
    /// ```
    ///
    /// The `TreeMap` finds keys by `Double.compareTo`, so `0.0` and `-0.0` are separate bins;
    /// `total_cmp` separates them too.
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
    /// cumulative range holds `(int)(count * precentile / 100)`, else the highest.
    ///
    /// Picard 3.5.0, `src/main/java/picard/util/SeriesStats.java:89-100`:
    ///
    /// ```java
    /// int percentileIndex = (int)(count * precentile / 100);
    /// int index = 0;
    /// for (Map.Entry<Double, AtomicLong> entry : bins.entrySet() ) {
    ///     long binSize = entry.getValue().get();
    ///     if ( percentileIndex >= index && (percentileIndex < (index + binSize)) ) {
    ///         return entry.getKey();
    ///     }
    ///     index += binSize;
    /// }
    ///
    /// // if here, we need the highest entry
    /// return bins.lastKey();
    /// ```
    ///
    /// The index is the floor of a product, not an interpolation, so the median of an even-sized
    /// series is its upper middle value.
    pub fn percentile(&self, percentile: f64) -> f64 {
        if self.count == 0 {
            return f64::NAN;
        }
        if self.count == 1 {
            return self.last;
        }
        let percentile_index = (self.count as f64 * percentile / 100.0) as i64;
        let mut index = 0i64;
        // `&self.bins` iterates the pairs by reference; `*value` reads the `f64` out of one.
        for (value, size) in &self.bins {
            if percentile_index >= index && percentile_index < index + size {
                return *value;
            }
            index += size;
        }
        // `bins.lastKey()`; reached only when the index falls outside the series (a percentile
        // of 100 or more, or a negative index), and the series has at least two values here, so
        // the `NaN` fallback is never taken.
        self.bins.last().map(|(v, _)| *v).unwrap_or(f64::NAN)
    }
}
