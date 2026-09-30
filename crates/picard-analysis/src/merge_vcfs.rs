//! Port of the parts of `picard.vcf.MergeVcfs` (Picard 3.4.0) that are not the shared VCF round
//! trip: htsjdk's `MergingIterator`, which decides the order of records that tie.
//!
//! # A tie is broken by a binary heap, not by the input order
//!
//! `MergingIterator` keeps one `ComparableIterator` per input in a `java.util.PriorityQueue`,
//! ordered by `VariantContextComparator` over each iterator's next record -- contig index, then
//! start, and nothing else. Two records at one position compare equal, and a `PriorityQueue` is a
//! binary heap whose `siftUp` stops at an equal parent and whose `siftDown` prefers the left child
//! of two equal ones, so which of them comes out first depends on where the heap happened to put
//! each iterator. Measured on the oracle with two inputs meeting at three positions: the first
//! input's record came first at two of them and the second's at the third. The heap below is
//! `PriorityQueue`'s own `offer` and `poll`, step for step, so the order is the same.
//!
//! # An unsorted input is noticed only when it goes backwards
//!
//! `next()` compares each record with the one returned before it -- from whichever input -- and
//! throws `IllegalStateException` when it is smaller. Everything before that record has already
//! been written.

use std::cmp::Ordering;

/// `java.util.PriorityQueue` over indices, with a comparison that can throw.
///
/// Only `offer` and `poll` are ported, which is all `MergingIterator` calls.
pub struct JavaPriorityQueue {
    queue: Vec<usize>,
}

impl JavaPriorityQueue {
    pub fn new() -> Self {
        Self { queue: Vec::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    /// `offer`: append, then `siftUpComparable`.
    pub fn offer<E>(
        &mut self,
        item: usize,
        compare: &mut impl FnMut(usize, usize) -> Result<Ordering, E>,
    ) -> Result<(), E> {
        let mut k = self.queue.len();
        self.queue.push(item);
        while k > 0 {
            let parent = (k - 1) >> 1;
            let e = self.queue[parent];
            if compare(item, e)? != Ordering::Less {
                break;
            }
            self.queue[k] = e;
            k = parent;
        }
        self.queue[k] = item;
        Ok(())
    }

    /// `poll`: take the head, move the last element to it, then `siftDownComparable`.
    pub fn poll<E>(
        &mut self,
        compare: &mut impl FnMut(usize, usize) -> Result<Ordering, E>,
    ) -> Result<Option<usize>, E> {
        if self.queue.is_empty() {
            return Ok(None);
        }
        let result = self.queue[0];
        let x = self.queue.pop().expect("not empty");
        let n = self.queue.len();
        if n > 0 {
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = (k << 1) + 1;
                let mut c = self.queue[child];
                let right = child + 1;
                if right < n && compare(c, self.queue[right])? == Ordering::Greater {
                    child = right;
                    c = self.queue[child];
                }
                if compare(x, c)? != Ordering::Greater {
                    break;
                }
                self.queue[k] = c;
                k = child;
            }
            self.queue[k] = x;
        }
        Ok(Some(result))
    }
}

impl Default for JavaPriorityQueue {
    fn default() -> Self {
        Self::new()
    }
}

/// Why a merge stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeFailure<E> {
    /// The comparator threw.
    Compare(E),
    /// `next()`'s order check: the exception line.
    NotSorted(String),
}

/// `MergingIterator` drained: the inputs' items in merged order, as `(input, position)` pairs.
///
/// `compare` is the comparator over two items; `comparator_class` is the class name the order
/// check's message carries. On a failure, what was returned before it is returned beside it,
/// because the tool has written those records by then.
pub fn merge_sorted<T, E>(
    inputs: &[Vec<T>],
    comparator_class: &str,
    compare: impl Fn(&T, &T) -> Result<Ordering, E>,
) -> (Vec<(usize, usize)>, Option<MergeFailure<E>>) {
    let mut positions = vec![0usize; inputs.len()];
    let mut queue = JavaPriorityQueue::new();
    let mut out: Vec<(usize, usize)> = Vec::new();

    macro_rules! peek_compare {
        ($positions:expr) => {
            |a: usize, b: usize| compare(&inputs[a][$positions[a]], &inputs[b][$positions[b]])
        };
    }

    for (input, items) in inputs.iter().enumerate() {
        // `addIfNotEmpty`: an exhausted iterator is closed rather than queued.
        if !items.is_empty() {
            if let Err(e) = queue.offer(input, &mut peek_compare!(positions)) {
                return (out, Some(MergeFailure::Compare(e)));
            }
        }
    }

    let mut last: Option<(usize, usize)> = None;
    loop {
        let polled = match queue.poll(&mut peek_compare!(positions)) {
            Ok(Some(input)) => input,
            Ok(None) => break,
            Err(e) => return (out, Some(MergeFailure::Compare(e))),
        };
        let next = (polled, positions[polled]);
        positions[polled] += 1;
        if let Some((li, lp)) = last {
            match compare(&inputs[li][lp], &inputs[next.0][next.1]) {
                Ok(Ordering::Greater) => {
                    return (
                        out,
                        Some(MergeFailure::NotSorted(format!(
                            "java.lang.IllegalStateException: The elements of the input \
                             Iterators are not sorted according to the comparator \
                             {comparator_class}"
                        ))),
                    )
                }
                Ok(_) => {}
                Err(e) => return (out, Some(MergeFailure::Compare(e))),
            }
        }
        if positions[polled] < inputs[polled].len() {
            if let Err(e) = queue.offer(polled, &mut peek_compare!(positions)) {
                return (out, Some(MergeFailure::Compare(e)));
            }
        }
        out.push(next);
        last = Some(next);
    }
    (out, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ties_follow_the_heap_rather_than_the_input_order() {
        // The two corpus files MergeVcfs is measured on, as (position, name): chr2 is 10000 on.
        let a = vec![
            (100, "a100"),
            (200, "a200"),
            (550, "a550"),
            (10100, "a2:100"),
        ];
        let b = vec![
            (100, "b100"),
            (150, "b150"),
            (550, "b550"),
            (700, "b700"),
            (10050, "b2:50"),
            (10100, "b2:100"),
        ];
        let inputs = [a, b];
        let (order, failure) = merge_sorted(&inputs, "C", |x: &(i32, &str), y: &(i32, &str)| {
            Ok::<_, ()>(x.0.cmp(&y.0))
        });
        assert!(failure.is_none());
        let names: Vec<&str> = order.iter().map(|&(i, p)| inputs[i][p].1).collect();
        // What the oracle wrote: a first at chr1:100 and chr2:100, b first at chr1:550.
        assert_eq!(
            names,
            ["a100", "b100", "b150", "a200", "b550", "a550", "b700", "b2:50", "a2:100", "b2:100"]
        );
    }

    #[test]
    fn an_input_that_goes_backwards_stops_the_merge() {
        let (order, failure) = merge_sorted(&[vec![30, 10]], "C", |a: &i32, b: &i32| {
            Ok::<_, ()>(a.cmp(b))
        });
        assert_eq!(order.len(), 1);
        assert!(matches!(failure, Some(MergeFailure::NotSorted(_))));
    }
}
