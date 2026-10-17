//! Two position indexes several tools need where a scan per record was quadratic.
//!
//! [`Spans`] answers whether a stretch overlaps any of a set of intervals, which is what an
//! interval-list query or `ByIntervalListVariantContextIterator` decides record by record.
//! [`PositionSet`] is `DbSnpBitSetUtil`: one bit per position of each contig, so a whole-genome
//! dbSNP costs a bitset per contig rather than an entry per base.

use std::hash::Hash;

/// Intervals grouped by contig, each group sorted by start and carrying the furthest end so far.
/// A stretch overlaps one of them when the furthest end among those starting at or before its
/// end reaches its start.
pub struct Spans<K> {
    by_contig: foldhash::HashMap<K, Vec<(i64, i64)>>,
}

impl<K: Hash + Eq> Spans<K> {
    /// The index of `(contig, start, end)` intervals, both ends inclusive.
    pub fn new(intervals: impl IntoIterator<Item = (K, i64, i64)>) -> Spans<K> {
        let mut by_contig: foldhash::HashMap<K, Vec<(i64, i64)>> = Default::default();
        for (contig, start, end) in intervals {
            by_contig.entry(contig).or_default().push((start, end));
        }
        for spans in by_contig.values_mut() {
            spans.sort_unstable();
            let mut furthest = i64::MIN;
            for span in spans.iter_mut() {
                furthest = furthest.max(span.1);
                span.1 = furthest;
            }
        }
        Spans { by_contig }
    }

    /// Whether `start..=end` on `contig` overlaps any interval.
    pub fn overlaps<Q>(&self, contig: &Q, start: i64, end: i64) -> bool
    where
        K: std::borrow::Borrow<Q>,
        Q: Hash + Eq + ?Sized,
    {
        self.by_contig.get(contig).is_some_and(|spans| {
            let before = spans.partition_point(|&(s, _)| s <= end);
            before > 0 && spans[before - 1].1 >= start
        })
    }
}

/// `DbSnpBitSetUtil`'s `sequenceToBitSet`: a bitset per contig, grown to the largest position set.
#[derive(Default)]
pub struct PositionSet {
    by_contig: foldhash::HashMap<String, Vec<u64>>,
}

impl PositionSet {
    /// `bits.set(i, true)` for every `i` from `start` to `end`, inclusive.
    pub fn set_range(&mut self, contig: &str, start: i64, end: i64) {
        if end < start || end < 0 {
            return;
        }
        let bits = match self.by_contig.get_mut(contig) {
            Some(bits) => bits,
            None => self.by_contig.entry(contig.to_string()).or_default(),
        };
        let words = (end as usize) / 64 + 1;
        if bits.len() < words {
            bits.resize(words, 0);
        }
        for position in start.max(0) as usize..=end as usize {
            bits[position / 64] |= 1 << (position % 64);
        }
    }

    /// `isDbSnpSite(contig, position)`.
    pub fn contains(&self, contig: &str, position: i64) -> bool {
        position >= 0
            && self.by_contig.get(contig).is_some_and(|bits| {
                bits.get(position as usize / 64)
                    .is_some_and(|word| word & (1 << (position as usize % 64)) != 0)
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_stretch_overlaps_an_interval_that_started_long_before() {
        let spans = Spans::new([("c", 10, 1000), ("c", 20, 30), ("d", 5, 5)]);
        assert!(spans.overlaps("c", 500, 600));
        assert!(spans.overlaps("c", 1000, 1200));
        assert!(!spans.overlaps("c", 1001, 1200));
        assert!(!spans.overlaps("c", 1, 9));
        assert!(spans.overlaps("d", 5, 5));
        assert!(!spans.overlaps("e", 5, 5));
    }

    #[test]
    fn a_range_sets_every_position_in_it_and_no_other() {
        let mut set = PositionSet::default();
        set.set_range("chr1", 62, 130);
        assert!(!set.contains("chr1", 61));
        assert!((62..=130).all(|p| set.contains("chr1", p)));
        assert!(!set.contains("chr1", 131));
        assert!(!set.contains("chr2", 62));
        assert!(!set.contains("chr1", 1_000_000));
    }
}
