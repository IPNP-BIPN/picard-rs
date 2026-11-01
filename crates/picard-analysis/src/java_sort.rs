//! `Arrays.parallelSort(T[], from, to, Comparator)`, for a comparator that is not transitive.
//!
//! Over a total preorder every stable sort gives the one answer, so the callers keep Rust's sort
//! there. `SAMRecordQueryNameComparator` is one only while every paired record says it is exactly
//! one of first and second of pair ([`orders_loosely`]). Otherwise the answer is whatever the
//! algorithm's comparisons make of it: Rust's depends on the size of the element and may panic, and
//! Java's may throw. At [`MIN_ARRAY_SORT_GRAN`] records or fewer Java runs `java.util.TimSort`
//! whatever the machine, and [`sort`] is that, comparison for comparison, down to the exception.
//!
//! The elements sorted are indices, so what a record weighs decides nothing.

use std::cmp::Ordering;

/// `Arrays.MIN_ARRAY_SORT_GRAN`.
const MIN_ARRAY_SORT_GRAN: usize = 1 << 13;
/// `TimSort.MIN_MERGE`.
const MIN_MERGE: usize = 32;
/// `TimSort.MIN_GALLOP`.
const MIN_GALLOP: usize = 7;

/// `IllegalArgumentException("Comparison method violates its general contract!")`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContractViolated;

impl std::fmt::Display for ContractViolated {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(
            "java.lang.IllegalArgumentException: Comparison method violates its general contract!",
        )
    }
}

impl std::error::Error for ContractViolated {}

impl From<ContractViolated> for std::io::Error {
    fn from(violated: ContractViolated) -> Self {
        std::io::Error::other(violated)
    }
}

/// Whether a record with these flags can make `SAMRecordQueryNameComparator` intransitive: paired,
/// and flagged both first and second of pair or neither.
pub fn orders_loosely(flags: u16) -> bool {
    const PAIRED: u16 = 0x1;
    const FIRST: u16 = 0x40;
    const SECOND: u16 = 0x80;
    flags & PAIRED != 0 && (flags & FIRST != 0) == (flags & SECOND != 0)
}

/// Sorts `items` as `Arrays.parallelSort` does.
pub fn sort<T>(
    items: &mut Vec<T>,
    mut compare: impl FnMut(&T, &T) -> Ordering,
) -> Result<(), ContractViolated> {
    let order = {
        let items = &*items;
        sorted_indices(items.len(), |a, b| {
            compare(&items[a as usize], &items[b as usize])
        })?
    };
    let mut slots: Vec<Option<T>> = items.drain(..).map(Some).collect();
    items.extend(
        order
            .into_iter()
            .map(|i| slots[i as usize].take().expect("an index appears once")),
    );
    Ok(())
}

fn sorted_indices(
    n: usize,
    mut compare: impl FnMut(u32, u32) -> Ordering,
) -> Result<Vec<u32>, ContractViolated> {
    let mut order: Vec<u32> = (0..n as u32).collect();
    if n <= MIN_ARRAY_SORT_GRAN {
        tim_sort(&mut order, &mut compare)?;
    } else {
        // shortcut: above the grain Java's answer depends on the common pool's parallelism, so
        // there is none to match. This one is stable, does not depend on the machine, and neither
        // panics nor refuses, where Java at some parallelisms throws.
        merge_sort(&mut order, &mut compare);
    }
    Ok(order)
}

fn merge_sort(a: &mut [u32], compare: &mut impl FnMut(u32, u32) -> Ordering) {
    if a.len() < 2 {
        return;
    }
    let (left, right) = a.split_at_mut(a.len() / 2);
    merge_sort(left, compare);
    merge_sort(right, compare);
    let mut merged = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0, 0);
    while i < left.len() && j < right.len() {
        if compare(right[j], left[i]) == Ordering::Less {
            merged.push(right[j]);
            j += 1;
        } else {
            merged.push(left[i]);
            i += 1;
        }
    }
    merged.extend_from_slice(&left[i..]);
    merged.extend_from_slice(&right[j..]);
    a.copy_from_slice(&merged);
}

/// `TimSort.sort`, from the JDK 17 the oracle runs.
fn tim_sort(
    a: &mut [u32],
    compare: &mut impl FnMut(u32, u32) -> Ordering,
) -> Result<(), ContractViolated> {
    let n = a.len();
    if n < 2 {
        return Ok(());
    }
    if n < MIN_MERGE {
        let run = count_run_and_make_ascending(a, 0, n, compare);
        binary_sort(a, 0, n, run, compare);
        return Ok(());
    }
    let mut state = TimSort {
        a,
        compare,
        min_gallop: MIN_GALLOP,
        runs: Vec::new(),
    };
    let min_run = min_run_length(n);
    let (mut lo, mut remaining) = (0, n);
    while remaining != 0 {
        let mut run = count_run_and_make_ascending(state.a, lo, n, state.compare);
        if run < min_run {
            let force = remaining.min(min_run);
            binary_sort(state.a, lo, lo + force, lo + run, state.compare);
            run = force;
        }
        state.runs.push((lo, run));
        state.merge_collapse()?;
        lo += run;
        remaining -= run;
    }
    state.merge_force_collapse()
}

fn min_run_length(mut n: usize) -> usize {
    let mut r = 0;
    while n >= MIN_MERGE {
        r |= n & 1;
        n >>= 1;
    }
    n + r
}

fn binary_sort(
    a: &mut [u32],
    lo: usize,
    hi: usize,
    start: usize,
    compare: &mut impl FnMut(u32, u32) -> Ordering,
) {
    for start in start.max(lo + 1)..hi {
        let pivot = a[start];
        let (mut left, mut right) = (lo, start);
        while left < right {
            let mid = (left + right) >> 1;
            if compare(pivot, a[mid]) == Ordering::Less {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        a.copy_within(left..start, left + 1);
        a[left] = pivot;
    }
}

fn count_run_and_make_ascending(
    a: &mut [u32],
    lo: usize,
    hi: usize,
    compare: &mut impl FnMut(u32, u32) -> Ordering,
) -> usize {
    let mut run_hi = lo + 1;
    if run_hi == hi {
        return 1;
    }
    run_hi += 1;
    if compare(a[lo + 1], a[lo]) == Ordering::Less {
        while run_hi < hi && compare(a[run_hi], a[run_hi - 1]) == Ordering::Less {
            run_hi += 1;
        }
        a[lo..run_hi].reverse();
    } else {
        while run_hi < hi && compare(a[run_hi], a[run_hi - 1]) != Ordering::Less {
            run_hi += 1;
        }
    }
    run_hi - lo
}

fn gallop_left(
    key: u32,
    a: &[u32],
    base: usize,
    len: usize,
    hint: usize,
    compare: &mut impl FnMut(u32, u32) -> Ordering,
) -> usize {
    let (mut last_ofs, mut ofs) = (0isize, 1isize);
    let hint = hint as isize;
    let at = |i: isize| a[(base as isize + i) as usize];
    if compare(key, at(hint)) == Ordering::Greater {
        let max_ofs = len as isize - hint;
        while ofs < max_ofs && compare(key, at(hint + ofs)) == Ordering::Greater {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        last_ofs += hint;
        ofs += hint;
    } else {
        let max_ofs = hint + 1;
        while ofs < max_ofs && compare(key, at(hint - ofs)) != Ordering::Greater {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        (last_ofs, ofs) = (hint - ofs, hint - last_ofs);
    }
    last_ofs += 1;
    while last_ofs < ofs {
        let m = last_ofs + ((ofs - last_ofs) >> 1);
        if compare(key, at(m)) == Ordering::Greater {
            last_ofs = m + 1;
        } else {
            ofs = m;
        }
    }
    ofs as usize
}

fn gallop_right(
    key: u32,
    a: &[u32],
    base: usize,
    len: usize,
    hint: usize,
    compare: &mut impl FnMut(u32, u32) -> Ordering,
) -> usize {
    let (mut last_ofs, mut ofs) = (0isize, 1isize);
    let hint = hint as isize;
    let at = |i: isize| a[(base as isize + i) as usize];
    if compare(key, at(hint)) == Ordering::Less {
        let max_ofs = hint + 1;
        while ofs < max_ofs && compare(key, at(hint - ofs)) == Ordering::Less {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        (last_ofs, ofs) = (hint - ofs, hint - last_ofs);
    } else {
        let max_ofs = len as isize - hint;
        while ofs < max_ofs && compare(key, at(hint + ofs)) != Ordering::Less {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        last_ofs += hint;
        ofs += hint;
    }
    last_ofs += 1;
    while last_ofs < ofs {
        let m = last_ofs + ((ofs - last_ofs) >> 1);
        if compare(key, at(m)) == Ordering::Less {
            ofs = m;
        } else {
            last_ofs = m + 1;
        }
    }
    ofs as usize
}

struct TimSort<'a, F> {
    a: &'a mut [u32],
    compare: &'a mut F,
    min_gallop: usize,
    /// `runBase` and `runLen`.
    runs: Vec<(usize, usize)>,
}

impl<F: FnMut(u32, u32) -> Ordering> TimSort<'_, F> {
    fn merge_collapse(&mut self) -> Result<(), ContractViolated> {
        while self.runs.len() > 1 {
            let len = |i: usize| self.runs[i].1;
            let mut n = self.runs.len() - 2;
            if (n > 0 && len(n - 1) <= len(n) + len(n + 1))
                || (n > 1 && len(n - 2) <= len(n) + len(n - 1))
            {
                if len(n - 1) < len(n + 1) {
                    n -= 1;
                }
            } else if len(n) > len(n + 1) {
                break;
            }
            self.merge_at(n)?;
        }
        Ok(())
    }

    fn merge_force_collapse(&mut self) -> Result<(), ContractViolated> {
        while self.runs.len() > 1 {
            let mut n = self.runs.len() - 2;
            if n > 0 && self.runs[n - 1].1 < self.runs[n + 1].1 {
                n -= 1;
            }
            self.merge_at(n)?;
        }
        Ok(())
    }

    fn merge_at(&mut self, i: usize) -> Result<(), ContractViolated> {
        let (mut base1, mut len1) = self.runs[i];
        let (base2, mut len2) = self.runs[i + 1];
        self.runs[i].1 = len1 + len2;
        self.runs.remove(i + 1);
        let k = gallop_right(self.a[base2], self.a, base1, len1, 0, self.compare);
        base1 += k;
        len1 -= k;
        if len1 == 0 {
            return Ok(());
        }
        len2 = gallop_left(
            self.a[base1 + len1 - 1],
            self.a,
            base2,
            len2,
            len2 - 1,
            self.compare,
        );
        if len2 == 0 {
            return Ok(());
        }
        if len1 <= len2 {
            self.merge_lo(base1, len1, base2, len2)
        } else {
            self.merge_hi(base1, len1, base2, len2)
        }
    }

    fn merge_lo(
        &mut self,
        base1: usize,
        mut len1: usize,
        base2: usize,
        mut len2: usize,
    ) -> Result<(), ContractViolated> {
        let a = &mut *self.a;
        let compare = &mut *self.compare;
        let tmp = a[base1..base1 + len1].to_vec();
        let (mut cursor1, mut cursor2, mut dest) = (0, base2, base1);
        a[dest] = a[cursor2];
        dest += 1;
        cursor2 += 1;
        len2 -= 1;
        if len2 == 0 {
            a[dest..dest + len1].copy_from_slice(&tmp[cursor1..cursor1 + len1]);
            return Ok(());
        }
        if len1 == 1 {
            a.copy_within(cursor2..cursor2 + len2, dest);
            a[dest + len2] = tmp[cursor1];
            return Ok(());
        }
        let mut min_gallop = self.min_gallop as isize;
        'outer: loop {
            let (mut count1, mut count2) = (0usize, 0usize);
            loop {
                if compare(a[cursor2], tmp[cursor1]) == Ordering::Less {
                    a[dest] = a[cursor2];
                    dest += 1;
                    cursor2 += 1;
                    count2 += 1;
                    count1 = 0;
                    len2 -= 1;
                    if len2 == 0 {
                        break 'outer;
                    }
                } else {
                    a[dest] = tmp[cursor1];
                    dest += 1;
                    cursor1 += 1;
                    count1 += 1;
                    count2 = 0;
                    len1 -= 1;
                    if len1 == 1 {
                        break 'outer;
                    }
                }
                if ((count1 | count2) as isize) >= min_gallop {
                    break;
                }
            }
            loop {
                count1 = gallop_right(a[cursor2], &tmp, cursor1, len1, 0, compare);
                if count1 != 0 {
                    a[dest..dest + count1].copy_from_slice(&tmp[cursor1..cursor1 + count1]);
                    dest += count1;
                    cursor1 += count1;
                    len1 -= count1;
                    if len1 <= 1 {
                        break 'outer;
                    }
                }
                a[dest] = a[cursor2];
                dest += 1;
                cursor2 += 1;
                len2 -= 1;
                if len2 == 0 {
                    break 'outer;
                }
                count2 = gallop_left(tmp[cursor1], a, cursor2, len2, 0, compare);
                if count2 != 0 {
                    a.copy_within(cursor2..cursor2 + count2, dest);
                    dest += count2;
                    cursor2 += count2;
                    len2 -= count2;
                    if len2 == 0 {
                        break 'outer;
                    }
                }
                a[dest] = tmp[cursor1];
                dest += 1;
                cursor1 += 1;
                len1 -= 1;
                if len1 == 1 {
                    break 'outer;
                }
                min_gallop -= 1;
                if count1 < MIN_GALLOP && count2 < MIN_GALLOP {
                    break;
                }
            }
            min_gallop = min_gallop.max(0) + 2;
        }
        self.min_gallop = min_gallop.max(1) as usize;
        match len1 {
            1 => {
                a.copy_within(cursor2..cursor2 + len2, dest);
                a[dest + len2] = tmp[cursor1];
                Ok(())
            }
            0 => Err(ContractViolated),
            _ => {
                a[dest..dest + len1].copy_from_slice(&tmp[cursor1..cursor1 + len1]);
                Ok(())
            }
        }
    }

    fn merge_hi(
        &mut self,
        base1: usize,
        mut len1: usize,
        base2: usize,
        mut len2: usize,
    ) -> Result<(), ContractViolated> {
        let a = &mut *self.a;
        let compare = &mut *self.compare;
        let tmp = a[base2..base2 + len2].to_vec();
        // Each of these can step to one before the start of what it walks.
        let (mut cursor1, mut cursor2, mut dest) = (
            (base1 + len1) as isize - 1,
            len2 as isize - 1,
            (base2 + len2) as isize - 1,
        );
        let copy_run_of_a = |a: &mut [u32], cursor1: isize, dest: isize, count: usize| {
            let from = (cursor1 + 1) as usize;
            a.copy_within(from..from + count, (dest + 1) as usize);
        };
        a[dest as usize] = a[cursor1 as usize];
        dest -= 1;
        cursor1 -= 1;
        len1 -= 1;
        if len1 == 0 {
            let to = (dest + 1) as usize - len2;
            a[to..to + len2].copy_from_slice(&tmp[..len2]);
            return Ok(());
        }
        if len2 == 1 {
            dest -= len1 as isize;
            cursor1 -= len1 as isize;
            copy_run_of_a(a, cursor1, dest, len1);
            a[dest as usize] = tmp[cursor2 as usize];
            return Ok(());
        }
        let mut min_gallop = self.min_gallop as isize;
        'outer: loop {
            let (mut count1, mut count2) = (0usize, 0usize);
            loop {
                if compare(tmp[cursor2 as usize], a[cursor1 as usize]) == Ordering::Less {
                    a[dest as usize] = a[cursor1 as usize];
                    dest -= 1;
                    cursor1 -= 1;
                    count1 += 1;
                    count2 = 0;
                    len1 -= 1;
                    if len1 == 0 {
                        break 'outer;
                    }
                } else {
                    a[dest as usize] = tmp[cursor2 as usize];
                    dest -= 1;
                    cursor2 -= 1;
                    count2 += 1;
                    count1 = 0;
                    len2 -= 1;
                    if len2 == 1 {
                        break 'outer;
                    }
                }
                if ((count1 | count2) as isize) >= min_gallop {
                    break;
                }
            }
            loop {
                count1 =
                    len1 - gallop_right(tmp[cursor2 as usize], a, base1, len1, len1 - 1, compare);
                if count1 != 0 {
                    dest -= count1 as isize;
                    cursor1 -= count1 as isize;
                    len1 -= count1;
                    copy_run_of_a(a, cursor1, dest, count1);
                    if len1 == 0 {
                        break 'outer;
                    }
                }
                a[dest as usize] = tmp[cursor2 as usize];
                dest -= 1;
                cursor2 -= 1;
                len2 -= 1;
                if len2 == 1 {
                    break 'outer;
                }
                count2 = len2 - gallop_left(a[cursor1 as usize], &tmp, 0, len2, len2 - 1, compare);
                if count2 != 0 {
                    dest -= count2 as isize;
                    cursor2 -= count2 as isize;
                    len2 -= count2;
                    let from = (cursor2 + 1) as usize;
                    let to = (dest + 1) as usize;
                    a[to..to + count2].copy_from_slice(&tmp[from..from + count2]);
                    if len2 <= 1 {
                        break 'outer;
                    }
                }
                a[dest as usize] = a[cursor1 as usize];
                dest -= 1;
                cursor1 -= 1;
                len1 -= 1;
                if len1 == 0 {
                    break 'outer;
                }
                min_gallop -= 1;
                if count1 < MIN_GALLOP && count2 < MIN_GALLOP {
                    break;
                }
            }
            min_gallop = min_gallop.max(0) + 2;
        }
        self.min_gallop = min_gallop.max(1) as usize;
        match len2 {
            1 => {
                dest -= len1 as isize;
                cursor1 -= len1 as isize;
                copy_run_of_a(a, cursor1, dest, len1);
                a[dest as usize] = tmp[cursor2 as usize];
                Ok(())
            }
            0 => Err(ContractViolated),
            _ => {
                let to = (dest + 1) as usize - len2;
                a[to..to + len2].copy_from_slice(&tmp[..len2]);
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// What `SAMRecordQueryNameComparator` reads of two records with one name and no `HI`.
    fn by_flags(a: u16, b: u16) -> Ordering {
        let (a_paired, b_paired) = (a & 1 != 0, b & 1 != 0);
        if a_paired || b_paired {
            if !a_paired {
                return Ordering::Greater;
            }
            if !b_paired {
                return Ordering::Less;
            }
            if a & 0x40 != 0 && b & 0x80 != 0 {
                return Ordering::Less;
            }
            if a & 0x80 != 0 && b & 0x40 != 0 {
                return Ordering::Greater;
            }
        }
        [0x10, 0x100, 0x800]
            .iter()
            .map(|bit| (a & bit != 0).cmp(&(b & bit != 0)))
            .find(|cmp| *cmp != Ordering::Equal)
            .unwrap_or(Ordering::Equal)
    }

    #[test]
    fn only_the_flags_called_loose_break_the_comparators_order() {
        let bits = [0x1, 0x10, 0x40, 0x80, 0x100, 0x800];
        let every: Vec<u16> = (0..1u16 << bits.len())
            .map(|set| {
                (0..bits.len())
                    .filter(|i| set >> i & 1 != 0)
                    .map(|i| bits[i])
                    .sum()
            })
            .collect();
        let transitive = |flags: &[u16]| {
            flags.iter().all(|&a| {
                flags.iter().all(|&b| {
                    flags.iter().all(|&c| {
                        by_flags(a, b) == Ordering::Greater
                            || by_flags(b, c) == Ordering::Greater
                            || by_flags(a, c) != Ordering::Greater
                    })
                })
            })
        };
        let tight: Vec<u16> = every
            .iter()
            .copied()
            .filter(|f| !orders_loosely(*f))
            .collect();
        assert!(transitive(&tight));
        assert!(!transitive(&every));
    }

    #[test]
    fn over_a_total_order_it_is_a_stable_sort() {
        let mut state = 11u32;
        for n in [0, 1, 2, 31, 32, 33, 1000, 8192, 8193, 30000] {
            let mut items: Vec<(u32, usize)> = (0..n)
                .map(|i| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    // Long runs in places, for the gallops, and many equal keys, for stability.
                    (
                        if i % 1000 < 400 {
                            i as u32 / 8
                        } else {
                            state >> 24
                        },
                        i,
                    )
                })
                .collect();
            let mut expected = items.clone();
            expected.sort_by_key(|item| item.0);
            sort(&mut items, |a, b| a.0.cmp(&b.0)).unwrap();
            assert_eq!(items, expected, "{n}");
        }
    }

    #[test]
    fn the_exception_is_the_one_java_throws() {
        assert_eq!(
            ContractViolated.to_string(),
            "java.lang.IllegalArgumentException: Comparison method violates its general contract!"
        );
    }
}
