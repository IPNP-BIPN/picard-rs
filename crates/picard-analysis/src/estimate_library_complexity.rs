//! `EstimateLibraryComplexity`: which pairs a group calls duplicates, and what the histogram gives.
//!
//! Reading the file and sorting the pairs are not ported. What is ported is the quality check that
//! admits a pair, the rule that calls two pairs duplicates, the floor that keeps a bin out of the
//! metrics, and the way the estimate is asked for.
//!
//! Ported from `picard.sam.markduplicates.EstimateLibraryComplexity`,
//! `picard.sam.markduplicates.ElcIdenticalBasesDuplicatesFinder` and
//! `picard.sam.DuplicationMetrics` in Picard 3.4.0.

use crate::jumping_library::estimate_library_size;

/// `MIN_IDENTICAL_BASES`: the prefix of BOTH ends that a group agrees on.
pub const DEFAULT_MIN_IDENTICAL_BASES: usize = 5;
/// `MAX_DIFF_RATE`: how far the rest may differ.
pub const DEFAULT_MAX_DIFF_RATE: f64 = 0.03;
/// `MIN_MEAN_QUALITY`: the floor a pair'send must clear.
pub const DEFAULT_MIN_MEAN_QUALITY: i32 = 20;
/// `MIN_GROUP_COUNT`: how many groups a bin needs before the METRICS count it.
pub const DEFAULT_MIN_GROUP_COUNT: i64 = 2;

/// `passesQualityCheck`.
///
/// The mean is an INTEGER division over the read's length, so qualities averaging nineteen and a
/// half are dropped at twenty. A read shorter than the seed fails outright, and an `N` anywhere in
/// the seed fails whatever the qualities say.
pub fn passes_quality_check(
    bases: &[u8],
    qualities: &[u8],
    seed_length: usize,
    minimum_quality: i32,
    max_read_length: usize,
) -> bool {
    if bases.len() < seed_length {
        return false;
    }
    if bases[..seed_length]
        .iter()
        .any(|base| base.eq_ignore_ascii_case(&b'N'))
    {
        return false;
    }
    let read_length = if max_read_length == 0 {
        bases.len()
    } else {
        bases.len().min(max_read_length)
    };
    let total: i32 = qualities[..read_length].iter().map(|q| i32::from(*q)).sum();
    total / read_length as i32 >= minimum_quality
}

/// Whether two pairs land in the same group: the first `seed_length` bases of BOTH ends agree.
pub fn same_group(left: (&[u8], &[u8]), right: (&[u8], &[u8]), seed_length: usize) -> bool {
    left.0.len() >= seed_length
        && right.0.len() >= seed_length
        && left.1.len() >= seed_length
        && right.1.len() >= seed_length
        && left.0[..seed_length] == right.0[..seed_length]
        && left.1[..seed_length] == right.1[..seed_length]
}

/// `ElcIdenticalBasesDuplicatesFinder.matches`.
///
/// The comparison starts AT the seed, because the grouping has already settled it, so a difference
/// inside the prefix is not an error: it is what put the two pairs in different groups. The
/// allowance is a rate over the two ends' compared lengths together, floored, and the ends are
/// compared over the SHORTER of the two.
pub fn matches(
    left: (&[u8], &[u8]),
    right: (&[u8], &[u8]),
    seed_length: usize,
    max_diff_rate: f64,
    max_read_length: usize,
) -> bool {
    let truncate = |length: usize| {
        if max_read_length == 0 {
            length
        } else {
            length.min(max_read_length)
        }
    };
    let read_one = truncate(left.0.len().min(right.0.len()));
    let read_two = truncate(left.1.len().min(right.1.len()));
    let max_errors = ((read_one + read_two) as f64 * max_diff_rate).floor() as i64;
    let mut errors = 0;
    for (a, b) in [
        (&left.0[..read_one], &right.0[..read_one]),
        (&left.1[..read_two], &right.1[..read_two]),
    ] {
        for index in seed_length..a.len() {
            if a[index] != b[index] {
                errors += 1;
                if errors > max_errors {
                    return false;
                }
            }
        }
    }
    true
}

/// The metrics one library's histogram gives.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Metrics {
    pub read_pairs_examined: i64,
    pub read_pair_duplicates: i64,
    pub read_pair_optical_duplicates: i64,
    pub percent_duplication: f64,
    pub estimated_library_size: Option<i64>,
}

/// The metrics a library's duplicate-set histogram produces.
///
/// A bin holding fewer than `min_group_count` groups is dropped HERE and nowhere else: the
/// histogram file still carries it, which is why a single duplicate pair reports nothing examined
/// beside a histogram that says there were two.
///
/// `bins` maps a duplicate-set size to how many sets of that size a library has, and `optical` to
/// how many of those were optical duplicates.
pub fn metrics(bins: &[(i64, i64, i64)], min_group_count: i64) -> Metrics {
    let mut out = Metrics::default();
    for (size, groups, optical) in bins {
        if *groups >= min_group_count {
            out.read_pairs_examined += size * groups;
            out.read_pair_duplicates += (size - 1) * groups;
            out.read_pair_optical_duplicates += optical;
        }
    }
    out.percent_duplication = if out.read_pairs_examined == 0 {
        0.0
    } else {
        (out.read_pair_duplicates * 2) as f64 / (out.read_pairs_examined * 2) as f64
    };
    // `calculateDerivedFields` takes the OPTICAL duplicates off the pair count before it estimates,
    // so a library whose every duplicate is optical has nothing left to estimate from.
    out.estimated_library_size = estimate_library_size(
        out.read_pairs_examined - out.read_pair_optical_duplicates,
        out.read_pairs_examined - out.read_pair_duplicates,
    );
    out
}

use crate::mark_duplicates::{optical_duplicates, Location, ReadEnds};
use std::collections::BTreeMap;

/// `MAX_GROUP_RATIO`'s default, the multiple of the expected group size past which a group is
/// dropped with a warning rather than searched.
pub const DEFAULT_MAX_GROUP_RATIO: i64 = 500;

/// One template as the tool holds it: both ends in READ order, and where the cluster was.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairedRead {
    pub read1: Vec<u8>,
    pub read2: Vec<u8>,
    pub library: String,
    pub read_group: i32,
    pub location: Location,
    /// `PairedReadSequenceWithBarcodes`: the three tags' `String.hashCode()`, or zero for a tag
    /// that was not asked for or not present. They are compared for equality and nothing else.
    pub barcodes: (i32, i32, i32),
}

/// `PairedReadComparator`: the first `seed` bases of read one, then of read two.
///
/// The comparison is a BYTE SUBTRACTION in the reference, over `byte`, which is signed there. Every
/// base is ASCII, so the sign never shows; the ordering is the unsigned one either way.
pub fn seed_order(left: &PairedRead, right: &PairedRead, seed: usize) -> std::cmp::Ordering {
    for index in 0..seed {
        match left.read1[index].cmp(&right.read1[index]) {
            std::cmp::Ordering::Equal => {}
            other => return other,
        }
    }
    for index in 0..seed {
        match left.read2[index].cmp(&right.read2[index]) {
            std::cmp::Ordering::Equal => {}
            other => return other,
        }
    }
    std::cmp::Ordering::Equal
}

/// `getNextGroup`: the run of pairs whose seeds equal the run's FIRST pair, which is what the sort
/// has already put together.
pub fn groups(pairs: &[PairedRead], seed: usize) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let mut start = 0;
    while start < pairs.len() {
        let mut end = start + 1;
        while end < pairs.len()
            && same_group(
                (&pairs[start].read1, &pairs[start].read2),
                (&pairs[end].read1, &pairs[end].read2),
                seed,
            )
        {
            end += 1;
        }
        out.push(start..end);
        start = end;
    }
    out
}

/// How many pairs [`PairSorter`] holds before it spills, when `MAX_RECORDS_IN_RAM` is not given.
/// shortcut: a fixed count (about 1 GB of 2x150 pairs), where Picard sizes it from the heap; raise
/// it when runs of this size are too many for a whole genome.
pub const DEFAULT_MAX_PAIRS_IN_RAM: usize = 2_000_000;

/// How many runs are merged at once; past it the runs so far are merged into one first, which keeps
/// the order since they all precede the runs still to come (as `sort_sam`'s collection does).
const MAX_OPEN_RUNS: usize = 256;

/// The `SortingCollection<PairedReadSequence>` with `PairedReadComparator`: holds `max_in_ram`
/// pairs, spills each full batch to `tmp_dir` as a stably sorted run, and merges the runs with
/// equal pairs taken from the earliest run first. That is a stable sort of the whole input, so the
/// order does not depend on the batch size.
pub struct PairSorter {
    seed: usize,
    max_in_ram: usize,
    tmp_dir: std::path::PathBuf,
    ram: Vec<PairedRead>,
    runs: Vec<std::fs::File>,
}

/// A pair's sort key: the first `seed` bases of each end, which compare as [`seed_order`] does.
fn seed_key(pair: &PairedRead, seed: usize) -> Vec<u8> {
    [&pair.read1[..seed], &pair.read2[..seed]].concat()
}

impl PairSorter {
    pub fn new(seed: usize, max_in_ram: usize, tmp_dir: std::path::PathBuf) -> Self {
        PairSorter {
            seed,
            max_in_ram: max_in_ram.max(1),
            tmp_dir,
            ram: Vec::new(),
            runs: Vec::new(),
        }
    }

    /// `add`: spill first when the buffer is already full.
    pub fn add(&mut self, pair: PairedRead) -> std::io::Result<()> {
        if self.ram.len() == self.max_in_ram {
            self.spill()?;
        }
        self.ram.push(pair);
        Ok(())
    }

    fn spill(&mut self) -> std::io::Result<()> {
        let seed = self.seed;
        self.ram.sort_by(|l, r| seed_order(l, r, seed));
        let run = write_run(&self.tmp_dir, self.ram.drain(..).map(Ok))?;
        self.runs.push(run);
        if self.runs.len() == MAX_OPEN_RUNS {
            let merged = PairMerge::new(seed, std::mem::take(&mut self.runs))?;
            let run = write_run(&self.tmp_dir, merged)?;
            self.runs.push(run);
        }
        Ok(())
    }

    /// `doneAdding` and `iterator`: the pairs in order, from RAM when nothing spilled.
    pub fn into_sorted(
        mut self,
    ) -> std::io::Result<Box<dyn Iterator<Item = std::io::Result<PairedRead>>>> {
        let seed = self.seed;
        if self.runs.is_empty() {
            self.ram.sort_by(|l, r| seed_order(l, r, seed));
            return Ok(Box::new(self.ram.into_iter().map(Ok)));
        }
        if !self.ram.is_empty() {
            self.spill()?;
        }
        Ok(Box::new(PairMerge::new(seed, self.runs)?))
    }
}

/// Writes pairs to a fresh file in `dir`, unlinked at once so nothing is left behind, and rewound.
fn write_run(
    dir: &std::path::Path,
    pairs: impl Iterator<Item = std::io::Result<PairedRead>>,
) -> std::io::Result<std::fs::File> {
    use std::io::{Seek, Write};
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = dir.join(format!(
        "elc-sortingcollection.{}.{n}.tmp",
        std::process::id()
    ));
    let file = std::fs::File::options()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)?;
    // shortcut: unlinking an open file keeps it readable on Unix only; Windows is not a target.
    std::fs::remove_file(&path)?;
    let mut out = std::io::BufWriter::new(file);
    for pair in pairs {
        let p = pair?;
        for bytes in [&p.read1, &p.read2, p.library.as_bytes()] {
            out.write_all(&(bytes.len() as u32).to_le_bytes())?;
            out.write_all(bytes)?;
        }
        out.write_all(&p.read_group.to_le_bytes())?;
        out.write_all(&p.location.tile.to_le_bytes())?;
        out.write_all(&p.location.x.to_le_bytes())?;
        out.write_all(&p.location.y.to_le_bytes())?;
        out.write_all(&[u8::from(p.location.known)])?;
        for b in [p.barcodes.0, p.barcodes.1, p.barcodes.2] {
            out.write_all(&b.to_le_bytes())?;
        }
    }
    let mut file = out.into_inner().map_err(|e| e.into_error())?;
    file.seek(std::io::SeekFrom::Start(0))?;
    Ok(file)
}

/// One spilled run, read back pair by pair; `None` at its end.
fn read_pair(reader: &mut impl std::io::BufRead) -> std::io::Result<Option<PairedRead>> {
    if reader.fill_buf()?.is_empty() {
        return Ok(None);
    }
    fn take<const N: usize>(r: &mut impl std::io::Read) -> std::io::Result<[u8; N]> {
        let mut b = [0u8; N];
        r.read_exact(&mut b)?;
        Ok(b)
    }
    let mut bytes = || -> std::io::Result<Vec<u8>> {
        let mut v = vec![0u8; u32::from_le_bytes(take(reader)?) as usize];
        reader.read_exact(&mut v)?;
        Ok(v)
    };
    let read1 = bytes()?;
    let read2 = bytes()?;
    let library = String::from_utf8(bytes()?).map_err(std::io::Error::other)?;
    Ok(Some(PairedRead {
        read1,
        read2,
        library,
        read_group: i32::from_le_bytes(take(reader)?),
        location: Location {
            tile: i16::from_le_bytes(take(reader)?),
            x: i32::from_le_bytes(take(reader)?),
            y: i32::from_le_bytes(take(reader)?),
            known: take::<1>(reader)?[0] != 0,
        },
        barcodes: (
            i32::from_le_bytes(take(reader)?),
            i32::from_le_bytes(take(reader)?),
            i32::from_le_bytes(take(reader)?),
        ),
    }))
}

/// The runs merged; of equal pairs the one from the earliest run comes first
/// (`MergingIterator`'s tie-break).
struct PairMerge {
    seed: usize,
    runs: Vec<std::io::BufReader<std::fs::File>>,
    heap: std::collections::BinaryHeap<std::cmp::Reverse<(Vec<u8>, usize, Slot)>>,
}

/// A pair in the heap, ignored by its ordering: the key and the run index decide.
struct Slot(PairedRead);
impl PartialEq for Slot {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}
impl Eq for Slot {}
impl PartialOrd for Slot {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Slot {
    fn cmp(&self, _: &Self) -> std::cmp::Ordering {
        std::cmp::Ordering::Equal
    }
}

impl PairMerge {
    fn new(seed: usize, files: Vec<std::fs::File>) -> std::io::Result<Self> {
        let mut merge = PairMerge {
            seed,
            runs: files.into_iter().map(std::io::BufReader::new).collect(),
            heap: Default::default(),
        };
        for run in 0..merge.runs.len() {
            merge.refill(run)?;
        }
        Ok(merge)
    }

    fn refill(&mut self, run: usize) -> std::io::Result<()> {
        if let Some(pair) = read_pair(&mut self.runs[run])? {
            let key = seed_key(&pair, self.seed);
            self.heap.push(std::cmp::Reverse((key, run, Slot(pair))));
        }
        Ok(())
    }
}

impl Iterator for PairMerge {
    type Item = std::io::Result<PairedRead>;

    fn next(&mut self) -> Option<Self::Item> {
        let std::cmp::Reverse((_, run, Slot(pair))) = self.heap.pop()?;
        Some(self.refill(run).map(|()| pair))
    }
}

/// What one library's search produced: how many groups of each size, and how many of those groups'
/// members were optical duplicates.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LibraryHistograms {
    pub duplication: BTreeMap<i64, i64>,
    pub optical: BTreeMap<i64, i64>,
}

/// `OPTICAL_DUPLICATE_PIXEL_DISTANCE` and `MAX_OPTICAL_DUPLICATE_SET_SIZE`: how close two reads
/// must be to be one cluster, and how large a set the finder looks at.
#[derive(Debug, Clone, Copy)]
pub struct OpticalLimits {
    pub distance: i32,
    pub max_set_size: i64,
}

/// `ElcIdenticalBasesDuplicatesFinder.searchDuplicates` and `fillHistogram`.
///
/// Each pair not already claimed takes every later pair that matches it, and the SIZE of what it
/// took -- itself included -- is the bin incremented. A pair that took nobody increments bin one,
/// and that bin is what `MIN_GROUP_COUNT` then decides about.
///
/// The optical count is per MEMBER and not per group: the flags come back one per member of the
/// set, and every flag set increments the same bin again. The keeper handed to the finder is the
/// pair that claimed the others, which sits LAST in the list the finder is given.
pub fn search_duplicates(
    group: &[&PairedRead],
    seed: usize,
    max_diff_rate: f64,
    max_read_length: usize,
    optical: OpticalLimits,
    use_barcodes: bool,
    histograms: &mut LibraryHistograms,
) {
    let mut claimed = vec![false; group.len()];
    for left in 0..group.len() {
        if claimed[left] {
            continue;
        }
        let mut dupes: Vec<usize> = Vec::new();
        for right in (left + 1)..group.len() {
            if claimed[right] {
                continue;
            }
            let same_barcodes = !use_barcodes || group[left].barcodes == group[right].barcodes;
            if same_barcodes
                && matches(
                    (&group[left].read1, &group[left].read2),
                    (&group[right].read1, &group[right].read2),
                    seed,
                    max_diff_rate,
                    max_read_length,
                )
            {
                dupes.push(right);
                claimed[right] = true;
            }
        }
        if dupes.is_empty() {
            *histograms.duplication.entry(1).or_insert(0) += 1;
            continue;
        }
        // `dupes.add(prs)`: the claiming pair is appended, so it is the LAST element and the
        // keeper at the same time.
        dupes.push(left);
        let size = dupes.len() as i64;
        *histograms.duplication.entry(size).or_insert(0) += 1;
        let ends: Vec<ReadEnds> = dupes
            .iter()
            .map(|index| ReadEnds {
                location: group[*index].location,
                read_group: group[*index].read_group,
                ..ReadEnds {
                    library: "".into(),
                    read1_reference_index: -1,
                    read1_coordinate: 0,
                    read2_reference_index: -1,
                    read2_coordinate: 0,
                    orientation: 0,
                    orientation_for_optical_duplicates: 0,
                    read1_index_in_file: 0,
                    read2_index_in_file: 0,
                    score: 0,
                    location: Location::default(),
                    read_group: -1,
                    barcode: None,
                    is_optical_duplicate: false,
                }
            })
            .collect();
        // `findOpticalDuplicates` returns nothing for a set longer than its maximum.
        let flags = if ends.len() as i64 > optical.max_set_size {
            vec![false; ends.len()]
        } else {
            optical_duplicates(&ends, Some(ends.len() - 1), optical.distance)
        };
        for flag in flags {
            if flag {
                *histograms.optical.entry(size).or_insert(0) += 1;
            }
        }
    }
}

#[cfg(test)]
mod sorter_tests {
    use super::*;

    #[test]
    fn spilled_runs_come_back_in_the_in_memory_order() {
        // Few distinct seeds, so most pairs tie and only stability orders them.
        let mut state = 7u64;
        let pairs: Vec<PairedRead> = (0..500)
            .map(|i| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let base = |s: u64| b"ACGT"[(s % 4) as usize];
                PairedRead {
                    read1: vec![base(state >> 33), base(state >> 40), b'A', b'C'],
                    read2: vec![base(state >> 50), b'T', b'G', b'G'],
                    library: format!("lib{}", i % 3),
                    read_group: i % 5,
                    location: Location {
                        tile: (i % 7) as i16,
                        x: i,
                        y: -i,
                        known: i % 2 == 0,
                    },
                    barcodes: (i, i + 1, -i),
                }
            })
            .collect();
        let mut expected = pairs.clone();
        expected.sort_by(|l, r| seed_order(l, r, 2));
        for max_in_ram in [1, 7, 499, 500, 10_000] {
            let mut sorter = PairSorter::new(2, max_in_ram, std::env::temp_dir());
            for pair in pairs.clone() {
                sorter.add(pair).unwrap();
            }
            let got: Vec<PairedRead> = sorter.into_sorted().unwrap().map(Result::unwrap).collect();
            assert_eq!(got, expected, "max_in_ram {max_in_ram}");
        }
    }
}
