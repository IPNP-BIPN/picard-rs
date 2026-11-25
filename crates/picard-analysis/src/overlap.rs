//! Two position indexes several tools need where a scan per record was quadratic.
//!
//! [`Spans`] answers whether a stretch overlaps any of a set of intervals, which is what an
//! interval-list query or `ByIntervalListVariantContextIterator` decides record by record.
//! [`PositionSet`] is `DbSnpBitSetUtil`: one bit per position of each contig, so a whole-genome
//! dbSNP costs a bitset per contig rather than an entry per base.
//!
//! ## What `Spans` stands in for (Picard 3.5.0, htsjdk 5.0.0)
//!
//! - `CollectVariantCallingMetrics`, the INPUT filter: an `OverlapDetector` over the uniqued
//!   target intervals, asked once per record (`picard/vcf/processor/VariantIteratorProducer.java`
//!   lines 126-128 and 266):
//!
//!   ```java
//!   final List<Interval> uniques = intervals.uniqued(false).getIntervals();
//!   this.intervalsOfInterestDetector = new OverlapDetector<>(0, 0);
//!   intervalsOfInterestDetector.addAll(uniques, uniques);
//!   ```
//!
//!   and each record is kept when `intervalsOfInterestDetector.getOverlaps(new Interval(...))`
//!   over `vc.getContig(), vc.getStart(), vc.getEnd()` is not empty.
//! - `CollectVariantCallingMetrics`'s dbSNP load and `CreateExtendedIlluminaManifest`'s dbSNP
//!   pass: `ByIntervalListVariantContextIterator` (`picard/vcf/`, lines 33, 51 and 53), which
//!   runs `this.reader.query(interval.getContig(), interval.getStart(), interval.getEnd())` once
//!   per interval of `intervals.uniqued()`, and drops a record already returned for the previous
//!   interval:
//!
//!   ```java
//!   null == previousInterval || !overlapsInterval(ctx, previousInterval)
//!   ```
//!
//! - `CollectIndependentReplicateMetrics`: an indexed BAM query over the het sites
//!   (`CollectIndependentReplicateMetrics.java:245`), `in.query(QueryInterval.optimizeIntervals(
//!   ...), false)`, where `false` asks for overlapping reads rather than contained ones.
//!
//! ## Why the answers are the same
//!
//! Every one of those decides closed overlap, both ends inclusive: an interval `[s, e]` and a
//! stretch `[start, end]` overlap when `s <= end && start <= e`. htsjdk's `IntervalTree`, under
//! `OverlapDetector.getOverlaps`, makes that test in `minOverlapper` (`IntervalTree.java:258-259`)
//! and again in `getNextOverlapper` (`IntervalTree.java:635`):
//!
//! ```java
//! } while (node != null && !(node.mStart <= end && start <= node.mEnd));
//! ```
//!
//! the BAM query filter (`BAMQueryMultipleIntervalsIteratorFilter.java:57-58`, either test true
//! means no overlap):
//!
//! ```java
//! else if (intervalEnd < record.getAlignmentStart()) return IntervalComparison.BEFORE;
//! else if (alignmentEnd < interval.start) return IntervalComparison.AFTER;
//! ```
//!
//! and an indexed VCF query: a Tribble index stops or skips on `f.getStart() > end` and
//! `f.getEnd() < start` (`TribbleIndexedFeatureReader.java:534` and `:542`), and a tabix index,
//! queried with `start - 1, end` (`TabixFeatureReader.java:140`), keeps a line when
//! `intv.end > beg && intv.beg < end` (`TabixReader.java:450`), the same test in 0-based
//! half-open coordinates.
//!
//! The Java merges first (`uniqued`, `optimizeIntervals`) and `Spans` does not. That changes no
//! answer: intervals are merged only when they overlap or abut, so a merged interval covers
//! exactly its members' positions, and a stretch overlaps it exactly when it overlaps a member.
//! The Java's per-interval queries plus the previous-interval check return each record once;
//! `Spans` is asked once per record, so it keeps the same records, once each.
//!
//! ## Where they are not the same
//!
//! - **A zero-length interval** (`end == start - 1`). `Spans` keeps it and applies the same two
//!   comparisons, so `(10, 9)` overlaps `5..=20` but not `10..=10`. `OverlapDetector.addLhs`
//!   drops it (`OverlapDetector.java:86`, `if (start <= end) { // Don't put in sequences that
//!   have no overlappable bases`), and a tabix query of it finds nothing (`TabixReader.java:194`,
//!   `if (beg >= end) return 0;`), so there it overlaps nothing. Of the callers, only
//!   `CollectVariantCallingMetrics` reads its intervals from a user's interval file
//!   (TARGET_INTERVALS); the manifest loci end at `start + length`, never before their start.
//! - **A stretch with `end < start`.** `OverlapDetector.getOverlaps` answers empty for it
//!   (`OverlapDetector.java:191`, `if (start > end) {`); `Spans`, like the BAM filter above,
//!   just applies the two comparisons.
//! - **`QueryInterval.end <= 0`** means "to the end of the reference" in the BAM filter
//!   (`BAMQueryMultipleIntervalsIteratorFilter.java:46`). Here an end of 0 is an ordinary end.
//!   The replicate tool builds its intervals from `vc.getStart(), vc.getEnd()`
//!   (`CollectIndependentReplicateMetrics.java:606`), so this differs only for a record that
//!   ends at position 0 or below.
//!
//! ## Reading this file without reading Rust
//!
//! - `Spans<K>` with `K: Hash + Eq` is a generic class whose key needs `hashCode` and
//!   `equals`; the callers key by contig name or by reference index.
//! - `impl IntoIterator<Item = T>` as a parameter accepts anything iterable, like `Iterable<T>`.
//! - `map.entry(k).or_default()` is `map.computeIfAbsent(k, k -> new ArrayList<>())`.
//! - `x.is_some_and(|v| ...)` is "present and the predicate holds", Java's
//!   `opt.filter(v -> ...).isPresent()`. `|v| ...` is a lambda.
//! - `as usize` converts a signed `i64` to the unsigned index type; it is only applied here to
//!   values already checked to be non-negative, so nothing wraps.
//! - `foldhash::HashMap` is a `HashMap` with a different hash function. These maps are only
//!   looked up, never iterated, so their order cannot reach any output.

use std::hash::Hash;

/// Intervals grouped by contig, each group sorted by start and carrying the furthest end so far.
/// A stretch overlaps one of them when the furthest end among those starting at or before its
/// end reaches its start.
///
/// That is the closed-overlap test of the module note, asked of all the intervals at once: the
/// intervals with `s <= end` are a prefix of the sorted list, and one of them has `e >= start`
/// exactly when the largest `e` in that prefix does. Keeping the running maximum is what lets
/// the index skip the merging the Java does: an interval nested inside a longer earlier one is
/// still answered for by the longer one's end.
pub struct Spans<K> {
    by_contig: foldhash::HashMap<K, Vec<(i64, i64)>>,
}

impl<K: Hash + Eq> Spans<K> {
    /// The index of `(contig, start, end)` intervals, both ends inclusive.
    ///
    /// No interval is dropped, merged or checked, including a zero-length one; see the module
    /// note for where that differs from `OverlapDetector`.
    pub fn new(intervals: impl IntoIterator<Item = (K, i64, i64)>) -> Spans<K> {
        let mut by_contig: foldhash::HashMap<K, Vec<(i64, i64)>> = Default::default();
        for (contig, start, end) in intervals {
            by_contig.entry(contig).or_default().push((start, end));
        }
        for spans in by_contig.values_mut() {
            // Tuples sort by start, then end. "Unstable" means equal elements may swap, which
            // cannot matter: two equal tuples are the same interval.
            spans.sort_unstable();
            // From here on `.1` (the tuple's second field) is no longer the interval's own end
            // but the largest end of it and every interval sorted before it.
            let mut furthest = i64::MIN;
            for span in spans.iter_mut() {
                furthest = furthest.max(span.1);
                span.1 = furthest;
            }
        }
        Spans { by_contig }
    }

    /// Whether `start..=end` on `contig` overlaps any interval.
    ///
    /// A contig with no intervals answers false, as `OverlapDetector.getOverlaps` answers an
    /// empty set when it has no tree for the contig (`OverlapDetector.java:185-186`).
    ///
    /// The `where` clause lets a `Spans<String>` be asked with a `&str` (`K: Borrow<Q>`), so
    /// a lookup does not build a `String`; `?Sized` admits `str`, which has no fixed size.
    pub fn overlaps<Q>(&self, contig: &Q, start: i64, end: i64) -> bool
    where
        K: std::borrow::Borrow<Q>,
        Q: Hash + Eq + ?Sized,
    {
        self.by_contig.get(contig).is_some_and(|spans| {
            // `partition_point` is a binary search: the number of leading intervals with
            // `s <= end`, valid because the list is sorted by `s`. `|&(s, _)|` unpacks the
            // tuple and ignores its second field.
            let before = spans.partition_point(|&(s, _)| s <= end);
            // The last of those carries the furthest end among them.
            before > 0 && spans[before - 1].1 >= start
        })
    }
}

/// `DbSnpBitSetUtil`'s `sequenceToBitSet`: a bitset per contig, grown to the largest position set.
///
/// ```java
/// private final Map<String, BitSet> sequenceToBitSet = new HashMap<>();
/// ```
///
/// Position `p` is bit `p % 64` of word `p / 64`. The Java allocates a contig's bitset on the
/// first record seen there, whatever its type, sized from the sequence dictionary; this
/// allocates on the first range set and grows to fit. Only membership is ever read, so the
/// size is not observable.
#[derive(Default)]
pub struct PositionSet {
    by_contig: foldhash::HashMap<String, Vec<u64>>,
}

impl PositionSet {
    /// `bits.set(i, true)` for every `i` from `start` to `end`, inclusive.
    ///
    /// `DbSnpBitSetUtil.java:253`, run by `loadVcf` for each SNP or indel record:
    ///
    /// ```java
    /// for (int i = kv.getStart(); i <= kv.getEnd(); i++)  bits.set(i, true);
    /// ```
    ///
    /// `end < start` sets nothing, as the loop runs zero times. A negative position is skipped
    /// here, where `BitSet.set` would throw `IndexOutOfBoundsException`.
    pub fn set_range(&mut self, contig: &str, start: i64, end: i64) {
        if end < start || end < 0 {
            return;
        }
        // Look up by `&str` first: `entry` needs an owned `String` key, and building one for
        // every record of a contig already present would be wasted.
        let bits = match self.by_contig.get_mut(contig) {
            Some(bits) => bits,
            None => self.by_contig.entry(contig.to_string()).or_default(),
        };
        let words = (end as usize) / 64 + 1;
        if bits.len() < words {
            bits.resize(words, 0);
        }
        // `a..=b` is the inclusive range; `1 << n` is a `u64` here because it is or-ed into one.
        for position in start.max(0) as usize..=end as usize {
            bits[position / 64] |= 1 << (position % 64);
        }
    }

    /// `isDbSnpSite(contig, position)`.
    ///
    /// `DbSnpBitSetUtil.java:267-269`:
    ///
    /// ```java
    /// return sequenceToBitSet.get(sequenceName) != null &&
    ///         pos <= sequenceToBitSet.get(sequenceName).length() &&
    ///         sequenceToBitSet.get(sequenceName).get(pos);
    /// ```
    ///
    /// A position past the last word answers false, as a bit past `length()` is unset in the
    /// Java. A negative position answers false, where `BitSet.get` would throw.
    pub fn contains(&self, contig: &str, position: i64) -> bool {
        position >= 0
            && self.by_contig.get(contig).is_some_and(|bits| {
                // `Vec::get` is a bounds-checked read: `None` past the end instead of a throw.
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
