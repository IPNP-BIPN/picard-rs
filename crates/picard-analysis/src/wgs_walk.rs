//! The walk behind `CollectWgsMetrics`, `CollectRawWgsMetrics` and
//! `CollectWgsMetricsWithNonZeroCoverage`: a locus iterator, the record filters it pulls through,
//! and the two collectors that turn its pileups into the histograms `WgsMetrics` is derived from.
//!
//! [`crate::collect_wgs_metrics`] holds the accounting rules in isolation. This is the machinery
//! that applies them to a file, ported line for line because the machinery leaks into the answer
//! in three places a cleaner design would lose:
//!
//! * **The filters count what the iterator has PULLED, not what the walk used.** The counting
//!   filters sit under a `FilteringSamIterator` that prefetches its next passing record, inside a
//!   `PeekableIterator` that prefetches again. With `STOP_AFTER` the walk ends mid-file, and the
//!   `PCT_EXC_*` columns then count every record up to two passing records past the last one the
//!   pileup saw. So the locus iterator is ported as the pull machine it is.
//! * **The fast iterator files a block's edges relative to the accumulator's head**, which is
//!   one locus earlier than the read whenever the previous read started one base before it. The
//!   fast collector reads positions from the edge itself and is mostly indifferent, but which
//!   locus processes an edge decides what `STOP_AFTER` cuts and what the overlap test sees.
//! * **The accumulation cap counts the head locus**, whichever locus that happens to be.
//!
//! Ported from htsjdk 4.2.0 `AbstractLocusIterator`, `SamLocusIterator`, `EdgeReadIterator`,
//! `EdgingRecordAndOffset`, `IntervalListReferenceSequenceMask`,
//! `WholeGenomeReferenceSequenceMask`, `FilteringSamIterator`, `PeekableIterator`,
//! `SamRecordIntervalIteratorFactory`, and Picard 3.4.0 `CollectWgsMetrics.WgsMetricsCollector`,
//! `FastWgsMetricsCollector`, `CounterManager`, `WgsMetricsProcessorImpl` and the `Counting*Filter`s.

use std::collections::{HashMap, HashSet, VecDeque};

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::record::BamRecord;

use crate::adapter::AdapterUtility;

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const SECONDARY: u16 = 0x100;
const VENDOR_FAILED: u16 = 0x200;
const DUPLICATE: u16 = 0x400;

/// A Java exception, as `Exception in thread "main"` prints it: `class: message`.
pub type Thrown = String;

/// The arguments the walk reads.
#[derive(Debug, Clone)]
pub struct WalkOptions {
    pub minimum_mapping_quality: i32,
    pub minimum_base_quality: i32,
    pub coverage_cap: i32,
    pub locus_accumulation_cap: i32,
    pub stop_after: i64,
    pub count_unpaired: bool,
    /// Picks the iterator: `EdgeReadIterator` when set, `SamLocusIterator` otherwise.
    pub use_fast_algorithm: bool,
    /// Picks the collector. `CollectWgsMetrics` pairs it with the iterator;
    /// `CollectWgsMetricsWithNonZeroCoverage` always uses the slow one, whatever the iterator.
    pub fast_collector: bool,
    pub read_length: i32,
}

/// How the records reach the locus iterator.
pub enum Source<'a> {
    /// `samReader.iterator()`: every record, unmapped ones included.
    WholeFile,
    /// An interval query. `indexed` is `samReader.hasIndex()`, which decides between the index's
    /// overlap query and a filtering scan, and the two disagree about unmapped reads that sit at
    /// their mate's position: the query returns them, the scan does not.
    Intervals {
        intervals: &'a [(i32, i32, i32)],
        indexed: bool,
    },
}

/// What the walk produces: the collector's arrays and the filters' base counts.
#[derive(Debug, Clone, Default)]
pub struct WalkResult {
    pub high_quality_depth: Vec<i64>,
    pub unfiltered_depth: Vec<i64>,
    pub unfiltered_baseq: Vec<i64>,
    pub excluded_by_baseq: i64,
    pub excluded_by_overlap: i64,
    pub excluded_by_capping: i64,
    pub excluded_by_adapter: i64,
    pub excluded_by_mapq: i64,
    pub excluded_by_dupe: i64,
    pub excluded_by_pairing: i64,
}

/// `SequenceUtil.isNoCall`.
pub fn is_no_call(base: u8) -> bool {
    base == b'N' || base == b'n' || base == b'.'
}

fn aligned_bases(record: &BamRecord) -> i64 {
    alignment_blocks(&record.cigar, record.alignment_start)
        .iter()
        .map(|b| b.length as i64)
        .sum()
}

// ------------------------------------------------------------------------------------------------
// The record stream: base iterator, counting filters, and the two lookaheads.
// ------------------------------------------------------------------------------------------------

struct Filters {
    adapter: AdapterUtility,
    minimum_mapping_quality: i32,
    count_unpaired: bool,
    adapter_bases: i64,
    mapq_bases: i64,
    dupe_bases: i64,
    pair_bases: i64,
}

impl Filters {
    /// `AggregateFilter.filterOut` over `[Secondary, Adapter, MapQ, Duplicate, (Paired)]`: the
    /// first filter that takes a record ends the chain, so a record is counted once at most.
    fn filter_out(&mut self, record: &BamRecord) -> bool {
        if record.flags & SECONDARY != 0 {
            return true;
        }
        if self.adapter.is_adapter(record) {
            self.adapter_bases += aligned_bases(record);
            return true;
        }
        if (record.mapping_quality as i32) < self.minimum_mapping_quality {
            self.mapq_bases += aligned_bases(record);
            return true;
        }
        if record.flags & DUPLICATE != 0 {
            self.dupe_bases += aligned_bases(record);
            return true;
        }
        if !self.count_unpaired && (record.flags & PAIRED == 0 || record.flags & MATE_UNMAPPED != 0)
        {
            self.pair_bases += aligned_bases(record);
            return true;
        }
        false
    }
}

/// `PeekableIterator(FilteringSamIterator(base, filters))`.
struct SamStream<'a> {
    records: &'a [BamRecord],
    base: Vec<usize>,
    position: usize,
    filtering_next: Option<usize>,
    outer: Option<usize>,
    filters: Filters,
}

impl<'a> SamStream<'a> {
    fn new(records: &'a [BamRecord], base: Vec<usize>, filters: Filters) -> Self {
        let mut stream = SamStream {
            records,
            base,
            position: 0,
            filtering_next: None,
            outer: None,
            filters,
        };
        // FilteringSamIterator's constructor prefetches; PeekableIterator's advances once.
        stream.filtering_next = stream.next_passing();
        stream.outer = stream.filtering_advance();
        stream
    }

    fn next_passing(&mut self) -> Option<usize> {
        while self.position < self.base.len() {
            let index = self.base[self.position];
            self.position += 1;
            if !self.filters.filter_out(&self.records[index]) {
                return Some(index);
            }
        }
        None
    }

    fn filtering_advance(&mut self) -> Option<usize> {
        let result = self.filtering_next?;
        self.filtering_next = self.next_passing();
        Some(result)
    }

    fn peek(&self) -> Option<usize> {
        self.outer
    }

    fn next(&mut self) {
        self.outer = self.filtering_advance();
    }
}

/// The records an interval query hands the locus iterator, in file order.
fn interval_base(
    records: &[BamRecord],
    intervals: &[(i32, i32, i32)],
    indexed: bool,
) -> Vec<usize> {
    let mut out = Vec::new();
    if indexed {
        // `BAMQueryMultipleIntervalsIteratorFilter.compareIntervalToRecord`, OVERLAPPING: an
        // unmapped read placed at its mate's position ends where it starts.
        for (index, record) in records.iter().enumerate() {
            if record.reference_index < 0 {
                continue;
            }
            let end = if record.flags & UNMAPPED != 0 && record.alignment_start != 0 {
                record.alignment_start
            } else {
                record.alignment_end()
            };
            let overlaps = intervals.iter().any(|&(seq, start, stop)| {
                seq == record.reference_index && stop >= record.alignment_start && end >= start
            });
            if overlaps {
                out.push(index);
            }
        }
        return out;
    }
    // `StopAfterFilteringIterator` over an `IntervalFilter`.
    let (stop_seq, stop_pos) = match intervals.last() {
        Some(&(seq, _, end)) => (seq, end),
        None => (-1, -1),
    };
    let mut current = 0usize;
    for (index, record) in records.iter().enumerate() {
        if record.reference_index == -1 || record.reference_index > stop_seq {
            break;
        }
        if record.reference_index == stop_seq && record.alignment_start > stop_pos {
            break;
        }
        while current < intervals.len()
            && (intervals[current].0 < record.reference_index
                || (intervals[current].0 == record.reference_index
                    && intervals[current].2 < record.alignment_start))
        {
            current += 1;
        }
        let keep = current < intervals.len()
            && intervals[current].0 == record.reference_index
            && intervals[current].1 <= record.alignment_end();
        if keep {
            out.push(index);
        }
    }
    out
}

// ------------------------------------------------------------------------------------------------
// The reference-sequence mask.
// ------------------------------------------------------------------------------------------------

enum Mask {
    Whole {
        lengths: Vec<i32>,
    },
    Intervals {
        per_sequence: Vec<Vec<(i32, i32)>>,
        current: i32,
        last_sequence: i32,
        last_position: i32,
    },
}

impl Mask {
    fn whole(lengths: Vec<i32>) -> Mask {
        Mask::Whole { lengths }
    }

    fn intervals(sequences: usize, intervals: &[(i32, i32, i32)]) -> Mask {
        let mut per_sequence = vec![Vec::new(); sequences];
        for &(seq, start, end) in intervals {
            if seq >= 0 && (seq as usize) < sequences {
                per_sequence[seq as usize].push((start, end));
            }
        }
        let (last_sequence, last_position) = match intervals.last() {
            Some(&(seq, _, end)) => (seq, end),
            None => (-1, 0),
        };
        Mask::Intervals {
            per_sequence,
            current: -1,
            last_sequence,
            last_position,
        }
    }

    fn load(&mut self, sequence: i32) -> Result<(), Thrown> {
        if let Mask::Intervals { current, .. } = self {
            if sequence < *current {
                return Err(format!(
                    "java.lang.IllegalArgumentException: Cannot look at an earlier sequence.  Current: {current}; requested: {sequence}"
                ));
            }
            if sequence > *current {
                *current = sequence;
            }
        }
        Ok(())
    }

    fn get(&mut self, sequence: i32, position: i32) -> Result<bool, Thrown> {
        match self {
            Mask::Whole { lengths } => {
                if sequence < 0 {
                    return Err(format!(
                        "java.lang.IllegalArgumentException: Negative sequence index {sequence}"
                    ));
                }
                match lengths.get(sequence as usize) {
                    None => Ok(false),
                    Some(&length) => Ok(position <= length),
                }
            }
            Mask::Intervals { .. } => {
                self.load(sequence)?;
                if let Mask::Intervals { per_sequence, .. } = self {
                    let ranges = per_sequence
                        .get(sequence as usize)
                        .map(Vec::as_slice)
                        .unwrap_or(&[]);
                    Ok(ranges.iter().any(|&(s, e)| s <= position && position <= e))
                } else {
                    unreachable!()
                }
            }
        }
    }

    fn next_position(&mut self, sequence: i32, position: i32) -> Result<i32, Thrown> {
        match self {
            Mask::Whole { .. } => {
                if self.get(sequence, position + 1)? {
                    Ok(position + 1)
                } else {
                    Ok(-1)
                }
            }
            Mask::Intervals { .. } => {
                self.load(sequence)?;
                if let Mask::Intervals { per_sequence, .. } = self {
                    let ranges = per_sequence
                        .get(sequence as usize)
                        .map(Vec::as_slice)
                        .unwrap_or(&[]);
                    let from = position + 1;
                    let mut best = -1;
                    for &(s, e) in ranges {
                        if e >= from {
                            let candidate = s.max(from);
                            if best == -1 || candidate < best {
                                best = candidate;
                            }
                        }
                    }
                    Ok(best)
                } else {
                    unreachable!()
                }
            }
        }
    }

    fn max_sequence(&self) -> i32 {
        match self {
            Mask::Whole { lengths } => lengths.len() as i32 - 1,
            Mask::Intervals { last_sequence, .. } => *last_sequence,
        }
    }

    fn max_position(&self) -> i32 {
        match self {
            Mask::Whole { lengths } => *lengths.last().unwrap_or(&0),
            Mask::Intervals { last_position, .. } => *last_position,
        }
    }
}

// ------------------------------------------------------------------------------------------------
// The locus iterator.
// ------------------------------------------------------------------------------------------------

/// One `RecordAndOffset` (slow) or `EdgingRecordAndOffset` (fast).
#[derive(Debug, Clone, Copy)]
pub struct Entry {
    pub record: usize,
    pub offset: i32,
    /// The edge fields, which the slow iterator leaves at zero.
    pub length: i32,
    pub ref_pos: i32,
    pub begin: bool,
    /// For a BEGIN, its own id; for an END, the id of the BEGIN it closes.
    pub id: usize,
}

#[derive(Debug)]
pub struct Locus {
    pub sequence: i32,
    pub position: i32,
    pub entries: Vec<Entry>,
    id: usize,
}

/// `LocusComparator.compare`, which returns a difference rather than a sign.
fn compare(a: (i32, i32), b: (i32, i32)) -> i32 {
    let reference = a.0.wrapping_sub(b.0);
    if reference == 0 {
        a.1.wrapping_sub(b.1)
    } else {
        reference
    }
}

struct EdgeIntervals<'a> {
    intervals: &'a [(i32, i32, i32)],
    /// Index of `currentInterval`, and of the `PeekableIterator`'s next element.
    current: Option<usize>,
    peek: usize,
}

struct LocusIterator<'a> {
    fast: bool,
    records: &'a [BamRecord],
    sam: SamStream<'a>,
    mask: Mask,
    complete: VecDeque<Locus>,
    accumulator: VecDeque<Locus>,
    include_non_pf_reads: bool,
    mapping_quality_cutoff: i32,
    max_reads_per_locus: i32,
    enforced_limit: bool,
    last_sequence: i32,
    last_position: i32,
    finished_aligned_reads: bool,
    next_id: usize,
    edge_intervals: Option<EdgeIntervals<'a>>,
}

impl<'a> LocusIterator<'a> {
    fn new_locus(&mut self, sequence: i32, position: i32) -> Locus {
        self.next_id += 1;
        Locus {
            sequence,
            position,
            entries: Vec::new(),
            id: self.next_id,
        }
    }

    fn sam_has_more(&self) -> bool {
        !self.finished_aligned_reads && self.sam.peek().is_some()
    }

    fn has_remaining_mask_bases(&mut self) -> Result<bool, Thrown> {
        let max = self.mask.max_sequence();
        if self.last_sequence < max {
            return Ok(true);
        }
        if self.last_sequence == max {
            let next = self
                .mask
                .next_position(self.last_sequence, self.last_position)?;
            return Ok(self.last_position < next);
        }
        Ok(false)
    }

    fn has_next(&mut self) -> Result<bool, Thrown> {
        while self.complete.is_empty()
            && (!self.accumulator.is_empty()
                || self.sam_has_more()
                || self.has_remaining_mask_bases()?)
        {
            if let Some(locus) = self.next()? {
                self.complete.push_front(locus);
            }
        }
        Ok(!self.complete.is_empty())
    }

    fn next(&mut self) -> Result<Option<Locus>, Thrown> {
        while self.complete.is_empty() && self.sam_has_more() {
            let index = self.sam.peek().expect("sam_has_more");
            let record = &self.records[index];
            if record.reference_index == -1 {
                self.finished_aligned_reads = true;
                continue;
            }
            if record.flags & UNMAPPED != 0
                || (record.mapping_quality as i32) < self.mapping_quality_cutoff
                || (!self.include_non_pf_reads && record.flags & VENDOR_FAILED != 0)
            {
                self.sam.next();
                continue;
            }
            let start = (record.reference_index, record.alignment_start);
            while !self.accumulator.is_empty()
                && (compare(
                    (self.accumulator[0].sequence, self.accumulator[0].position),
                    start,
                ) < -1
                    || self.accumulator[0].sequence != start.0)
            {
                let first = self.accumulator[0].id;
                self.populate_complete_queue(start)?;
                if let Some(locus) = self.complete.pop_front() {
                    return Ok(Some(locus));
                }
                if !self.accumulator.is_empty() && self.accumulator[0].id == first {
                    return Err("htsjdk.samtools.SAMException: Stuck in infinite loop".to_string());
                }
            }
            if let Some(head) = self.accumulator.front() {
                if head.sequence != record.reference_index
                    || record.alignment_start - head.position > 1
                {
                    return Err("java.lang.IllegalStateException: Accumulator should be empty or aligned with current or previous SAMRecord".to_string());
                }
            }
            if !self.surpassed_accumulation_threshold() {
                if self.fast {
                    self.accumulate_edges(index)?;
                } else {
                    self.accumulate_bases(index);
                }
            }
            self.sam.next();
        }

        let end = (i32::MAX, i32::MAX);
        if self.complete.is_empty() && !self.sam_has_more() {
            while !self.accumulator.is_empty() {
                self.populate_complete_queue(end)?;
                if let Some(locus) = self.complete.pop_front() {
                    return Ok(Some(locus));
                }
            }
        }
        if let Some(locus) = self.complete.pop_front() {
            return Ok(Some(locus));
        }
        let after_last = (self.mask.max_sequence(), self.mask.max_position() + 1);
        self.create_next_uncovered_locus(after_last)
    }

    fn surpassed_accumulation_threshold(&mut self) -> bool {
        let surpasses = match self.accumulator.front() {
            Some(head) => head.entries.len() as i64 >= self.max_reads_per_locus as i64,
            None => false,
        };
        if surpasses && !self.enforced_limit {
            self.enforced_limit = true;
        }
        surpasses
    }

    fn create_next_uncovered_locus(&mut self, stop: (i32, i32)) -> Result<Option<Locus>, Thrown> {
        while self.last_sequence <= stop.0 && self.last_sequence <= self.mask.max_sequence() {
            if self.last_sequence == stop.0 && self.last_position.wrapping_add(1) >= stop.1 {
                return Ok(None);
            }
            let next_bit = self
                .mask
                .next_position(self.last_sequence, self.last_position)?;
            if next_bit == -1 {
                if self.last_sequence == stop.0 {
                    self.last_position = stop.1;
                    return Ok(None);
                }
                self.last_sequence += 1;
                self.last_position = 0;
            } else if self.last_sequence < stop.0 || next_bit < stop.1 {
                self.last_position = next_bit;
                let locus = self.new_locus(self.last_sequence, self.last_position);
                return Ok(Some(locus));
            } else {
                return Ok(None);
            }
        }
        Ok(None)
    }

    fn populate_complete_queue(&mut self, stop: (i32, i32)) -> Result<(), Thrown> {
        // removeSkippedRegion
        let mut skipped = 0;
        while skipped < self.accumulator.len()
            && self.accumulator[skipped].entries.is_empty()
            && compare(
                (
                    self.accumulator[skipped].sequence,
                    self.accumulator[skipped].position,
                ),
                stop,
            ) < 0
        {
            skipped += 1;
        }
        self.accumulator.drain(..skipped);

        let Some(head) = self.accumulator.front() else {
            return Ok(());
        };
        let head_locus = (head.sequence, head.position);
        if compare(stop, head_locus) <= 0 {
            return Ok(());
        }
        if let Some(zero) = self.create_next_uncovered_locus(head_locus)? {
            self.complete.push_back(zero);
            return Ok(());
        }
        let locus = self.accumulator.pop_front().expect("non-empty");
        let (sequence, position) = (locus.sequence, locus.position);
        if self.mask.get(sequence, position)? {
            self.complete.push_back(locus);
        }
        self.last_sequence = sequence;
        self.last_position = position;
        Ok(())
    }

    /// `SamLocusIterator.accumulateSamRecord`, with the quality cutoff at zero (no check).
    fn accumulate_bases(&mut self, index: usize) {
        let record = &self.records[index];
        let sequence = record.reference_index;
        let alignment_start = record.alignment_start;
        let alignment_end = record.alignment_end();
        let alignment_length = alignment_end - alignment_start;
        let acc_index_where_read_starts = match self.accumulator.front() {
            None => 0,
            Some(head) => alignment_start - head.position,
        };
        let new_loci =
            acc_index_where_read_starts + alignment_length - self.accumulator.len() as i32;
        let mut i = 0;
        while i <= new_loci {
            let locus = self.new_locus(sequence, alignment_end - new_loci + i);
            self.accumulator.push_back(locus);
            i += 1;
        }
        let head = self.accumulator[0].position;
        for block in alignment_blocks(&record.cigar, record.alignment_start) {
            let block_start = block.reference_start - head;
            for i in 0..block.length {
                let read_offset = block.read_start + i - 1;
                self.accumulator[(block_start + i) as usize]
                    .entries
                    .push(Entry {
                        record: index,
                        offset: read_offset,
                        length: 1,
                        ref_pos: 0,
                        begin: true,
                        id: 0,
                    });
            }
        }
    }

    /// `EdgeReadIterator.advanceCurrentIntervalAndCheckIfIntervalContainsRead`.
    fn interval_contains_read(&mut self, record: &BamRecord) -> bool {
        let Some(edge) = self.edge_intervals.as_mut() else {
            return false;
        };
        let Some(mut current) = edge.current else {
            return false;
        };
        let read = (
            record.reference_index,
            record.alignment_start,
            record.alignment_end(),
        );
        while edge.peek < edge.intervals.len() {
            let peek = edge.intervals[edge.peek];
            let mut order = read.0 - peek.0;
            if order == 0 {
                order = read.1 - peek.1;
            }
            if order == 0 {
                order = read.2 - peek.2;
            }
            if order > 0 {
                current = edge.peek;
                edge.peek += 1;
            } else {
                break;
            }
        }
        edge.current = Some(current);
        let interval = edge.intervals[current];
        interval.0 == read.0 && interval.1 <= read.1 && interval.2 >= read.2
    }

    /// `EdgeReadIterator.accumulateSamRecord`.
    fn accumulate_edges(&mut self, index: usize) -> Result<(), Thrown> {
        let records = self.records;
        let record = &records[index];
        let need_intervals = self.edge_intervals.is_some() && !self.interval_contains_read(record);
        let sequence = record.reference_index;
        for block in alignment_blocks(&record.cigar, record.alignment_start) {
            let offset_in_read = block.read_start - 1;
            let block_ref_start = block.reference_start;
            if self.accumulator.is_empty() {
                let locus = self.new_locus(sequence, record.alignment_start);
                self.accumulator.push_back(locus);
            }
            let next_position = self.accumulator[0].position + self.accumulator.len() as i32;
            if next_position != self.accumulator.back().expect("non-empty").position + 1 {
                return Err("java.lang.IllegalStateException: The accumulator has gotten into a funk. Cannot continue".to_string());
            }
            let mut locus_position = next_position;
            while locus_position <= block_ref_start + block.length {
                let locus = self.new_locus(sequence, locus_position);
                self.accumulator.push_back(locus);
                locus_position += 1;
            }
            let off_start = block_ref_start - record.alignment_start;
            let off_end = off_start + block.length;
            let mut pieces: Vec<(i32, i32, i32, i32, i32)> = Vec::new();
            if need_intervals {
                let query_end = block_ref_start + block.length;
                let intervals = self
                    .edge_intervals
                    .as_ref()
                    .expect("need_intervals")
                    .intervals;
                for &(seq, start, end) in intervals {
                    if seq != sequence || end < block_ref_start || start > query_end {
                        continue;
                    }
                    let start_in_block = if block_ref_start < start {
                        start - block_ref_start
                    } else {
                        0
                    };
                    let actual_start = off_start + start_in_block;
                    let ref_end_block = block_ref_start + block.length;
                    let actual_end = off_end
                        - if ref_end_block > end {
                            ref_end_block - end - 1
                        } else {
                            0
                        };
                    pieces.push((
                        actual_start,
                        actual_end,
                        offset_in_read + start_in_block,
                        actual_end - actual_start,
                        block_ref_start + start_in_block,
                    ));
                }
            } else {
                pieces.push((
                    off_start,
                    off_end,
                    offset_in_read,
                    off_end - off_start,
                    block_ref_start,
                ));
            }
            for (at_start, at_end, offset, length, ref_pos) in pieces {
                if length as usize > record.read_length() {
                    return Err("java.lang.IllegalArgumentException: Block length cannot be larger than whole read length".to_string());
                }
                self.next_id += 1;
                let id = self.next_id;
                let begin = Entry {
                    record: index,
                    offset,
                    length,
                    ref_pos,
                    begin: true,
                    id,
                };
                let end = Entry {
                    begin: false,
                    ..begin
                };
                self.accumulator[at_start as usize].entries.push(begin);
                self.accumulator[at_end as usize].entries.push(end);
            }
        }
        Ok(())
    }
}

// ------------------------------------------------------------------------------------------------
// The collectors.
// ------------------------------------------------------------------------------------------------

/// `CounterManager`, with its two counters.
struct CounterManager {
    array_length: i32,
    read_length: i32,
    offset: i32,
    pileup: Vec<i32>,
    unfiltered: Vec<i32>,
}

const COUNTER_BOUNDS: &str = "is out of counter bounds. Possible cause of exception can be wrong READ_LENGTH parameter (much smaller than actual read length)";

impl CounterManager {
    fn new(array_length: i32, read_length: i32) -> Self {
        let length = array_length.max(0) as usize;
        CounterManager {
            array_length,
            read_length,
            offset: 0,
            pileup: vec![0; length],
            unfiltered: vec![0; length],
        }
    }

    fn check_out_of_bounds(&mut self, position: i32) -> Result<(), Thrown> {
        if position - self.offset + self.read_length >= self.array_length {
            if position - self.offset < self.array_length {
                if position < self.offset {
                    return Err("java.lang.IllegalArgumentException: Position in the reference sequence is lesser than offset.".to_string());
                }
                let skip = (position - self.offset) as usize;
                let length = self.array_length as usize;
                for array in [&mut self.pileup, &mut self.unfiltered] {
                    array.copy_within(skip..length, 0);
                    for value in &mut array[length - skip..length] {
                        *value = 0;
                    }
                }
                self.offset = position;
            } else {
                self.clear();
                self.offset = position;
            }
        }
        Ok(())
    }

    fn clear(&mut self) {
        self.pileup.iter_mut().for_each(|v| *v = 0);
        self.unfiltered.iter_mut().for_each(|v| *v = 0);
        self.offset = 0;
    }

    fn slot(&self, index: i32) -> Result<usize, Thrown> {
        let at = index - self.offset;
        if at < 0 || at >= self.array_length {
            return Err(format!(
                "java.lang.ArrayIndexOutOfBoundsException: The requested index {index} {COUNTER_BOUNDS}"
            ));
        }
        Ok(at as usize)
    }
}

struct Collector<'a> {
    records: &'a [BamRecord],
    references: &'a [Vec<u8>],
    options: WalkOptions,
    result: WalkResult,
    counter: i64,
    // fast-only state
    previous_sequence: i32,
    counters: CounterManager,
    read_names: Option<HashMap<&'a str, Vec<Entry>>>,
}

impl<'a> Collector<'a> {
    fn is_time_to_stop(&self, processed: i64) -> bool {
        self.options.stop_after > 0 && processed > self.options.stop_after - 1
    }

    fn reference_base_n(&self, sequence: i32, position: i32) -> Result<bool, Thrown> {
        let bases = &self.references[sequence as usize];
        match bases.get((position - 1) as usize) {
            Some(&b) => Ok(is_no_call(b)),
            None => Err(format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {} out of bounds for length {}",
                position - 1,
                bases.len()
            )),
        }
    }

    /// `CollectWgsMetrics.WgsMetricsCollector.addInfo`.
    fn add_info_slow(&mut self, locus: &Locus, reference_n: bool) -> Result<(), Thrown> {
        if reference_n {
            return Ok(());
        }
        // The slow collector handed the fast iterator's edges: the enhanced `for` casts each
        // element to `SamLocusIterator.RecordAndOffset`, and the first one fails.
        if self.options.use_fast_algorithm {
            if let Some(first) = locus.entries.first() {
                let class = if first.begin {
                    "htsjdk.samtools.util.EdgingRecordAndOffset$StartEdgingRecordAndOffset"
                } else {
                    "htsjdk.samtools.util.EdgingRecordAndOffset$EndEdgingRecordAndOffset"
                };
                return Err(format!(
                    "java.lang.ClassCastException: class {class} cannot be cast to class htsjdk.samtools.util.SamLocusIterator$RecordAndOffset ({class} and htsjdk.samtools.util.SamLocusIterator$RecordAndOffset are in unnamed module of loader 'app')"
                ));
            }
        }
        let cap = self.options.coverage_cap;
        let mut names: HashSet<&str> = HashSet::with_capacity(locus.entries.len());
        let mut pileup = 0i32;
        let mut unfiltered = 0i32;
        for entry in &locus.entries {
            let record = &self.records[entry.record];
            let quality = record.base_qualities[entry.offset as usize] as i8 as i32;
            if quality <= 2 {
                self.result.excluded_by_baseq += 1;
                continue;
            }
            if unfiltered < cap {
                self.result.unfiltered_baseq[quality as usize] += 1;
                unfiltered += 1;
            }
            if quality < self.options.minimum_base_quality
                || is_no_call(record.read_bases[entry.offset as usize])
            {
                self.result.excluded_by_baseq += 1;
                continue;
            }
            if !names.insert(record.read_name.as_str()) {
                self.result.excluded_by_overlap += 1;
                continue;
            }
            pileup += 1;
        }
        let high_quality = pileup.min(cap);
        if high_quality < pileup {
            self.result.excluded_by_capping += (pileup - cap) as i64;
        }
        self.result.high_quality_depth[high_quality as usize] += 1;
        self.result.unfiltered_depth[unfiltered as usize] += 1;
        Ok(())
    }

    /// `FastWgsMetricsCollector.addInfo`.
    fn add_info_fast(&mut self, locus: &Locus, reference_n: bool) -> Result<(), Thrown> {
        if self.read_names.is_none() {
            self.read_names = Some(HashMap::new());
        }
        if self.previous_sequence != locus.sequence {
            self.read_names.as_mut().expect("set").clear();
            self.counters.clear();
            self.previous_sequence = locus.sequence;
        }
        self.counters.check_out_of_bounds(locus.position)?;
        let records = self.records;
        for entry in &locus.entries {
            let name = records[entry.record].read_name.as_str();
            if entry.begin {
                self.process_record(locus.sequence, entry, name)?;
            } else {
                let names = self.read_names.as_mut().expect("set");
                if let Some(set) = names.get_mut(name) {
                    if set.len() == 1 {
                        names.remove(name);
                    } else {
                        set.retain(|e| e.id != entry.id);
                    }
                }
            }
        }
        if !reference_n {
            let cap = self.options.coverage_cap;
            let slot = self.counters.slot(locus.position)?;
            let pileup = self.counters.pileup[slot];
            let high_quality = pileup.min(cap);
            if high_quality < pileup {
                self.result.excluded_by_capping += (pileup - cap) as i64;
            }
            self.result.high_quality_depth[high_quality as usize] += 1;
            let slot = self.counters.slot(locus.position)?;
            let unfiltered = self.counters.unfiltered[slot];
            self.result.unfiltered_depth[unfiltered as usize] += 1;
        }
        Ok(())
    }

    fn process_record(
        &mut self,
        sequence: i32,
        entry: &Entry,
        name: &'a str,
    ) -> Result<(), Thrown> {
        let mut processed = self.counter;
        let set = self
            .read_names
            .as_mut()
            .expect("set")
            .remove(name)
            .unwrap_or_default();
        // The set is in the map while the record is processed; it is held here and put back.
        let record = &self.records[entry.record];
        let cap = self.options.coverage_cap;
        let mut outcome = Ok(());
        for i in 0..entry.length {
            let index = i + entry.ref_pos;
            if self.reference_base_n(sequence, index)? {
                continue;
            }
            let quality = record.base_qualities[(i + entry.offset) as usize] as i8 as i32;
            if quality <= 2 {
                self.result.excluded_by_baseq += 1;
            } else {
                let slot = match self.counters.slot(index) {
                    Ok(slot) => slot,
                    Err(e) => {
                        outcome = Err(e);
                        break;
                    }
                };
                if self.counters.unfiltered[slot] < cap {
                    self.result.unfiltered_baseq[quality as usize] += 1;
                    self.counters.unfiltered[slot] += 1;
                }
                if quality < self.options.minimum_base_quality
                    || is_no_call(record.read_bases[(i + entry.offset) as usize])
                {
                    self.result.excluded_by_baseq += 1;
                } else {
                    let mut low = 0usize;
                    for other in &set {
                        if index - other.ref_pos >= other.length {
                            low += 1;
                            continue;
                        }
                        let quals = &self.records[other.record].base_qualities;
                        let relative = index - other.ref_pos + other.offset;
                        if relative < 0 || relative as usize >= quals.len() {
                            outcome = Err("java.lang.IllegalArgumentException: The requested position is not covered by this StartEdgingRecordAndOffset object. ".to_string());
                            break;
                        }
                        if (quals[relative as usize] as i8 as i32)
                            < self.options.minimum_base_quality
                        {
                            low += 1;
                        }
                    }
                    if outcome.is_err() {
                        break;
                    }
                    if set.len() - low > 0 {
                        self.result.excluded_by_overlap += 1;
                    } else {
                        self.counters.pileup[slot] += 1;
                    }
                }
            }
            processed += 1;
            if self.is_time_to_stop(processed) {
                break;
            }
        }
        outcome?;
        let mut set = set;
        set.push(*entry);
        self.read_names.as_mut().expect("set").insert(name, set);
        Ok(())
    }
}

/// `WgsMetricsProcessorImpl.processFile` over a fresh locus iterator and collector.
///
/// `references` holds each dictionary sequence's bases, in dictionary order; `lengths` the
/// dictionary's lengths, which the whole-genome mask reads.
pub fn walk(
    records: &[BamRecord],
    references: &[Vec<u8>],
    lengths: &[i32],
    source: Source<'_>,
    options: &WalkOptions,
) -> Result<WalkResult, Thrown> {
    let base: Vec<usize> = match &source {
        Source::WholeFile => (0..records.len()).collect(),
        Source::Intervals { intervals, indexed } => interval_base(records, intervals, *indexed),
    };
    let filters = Filters {
        adapter: AdapterUtility::with_defaults(),
        minimum_mapping_quality: options.minimum_mapping_quality,
        count_unpaired: options.count_unpaired,
        adapter_bases: 0,
        mapq_bases: 0,
        dupe_bases: 0,
        pair_bases: 0,
    };
    let (mask, edge_intervals) = match &source {
        Source::WholeFile => (Mask::whole(lengths.to_vec()), None),
        Source::Intervals { intervals, .. } => (
            Mask::intervals(lengths.len(), intervals),
            Some(EdgeIntervals {
                intervals,
                current: if intervals.is_empty() { None } else { Some(0) },
                peek: 1,
            }),
        ),
    };
    let mut iterator = LocusIterator {
        fast: options.use_fast_algorithm,
        records,
        sam: SamStream::new(records, base, filters),
        mask,
        complete: VecDeque::new(),
        accumulator: VecDeque::new(),
        include_non_pf_reads: false,
        mapping_quality_cutoff: 0,
        max_reads_per_locus: if options.use_fast_algorithm {
            i32::MAX
        } else {
            options.locus_accumulation_cap
        },
        enforced_limit: false,
        last_sequence: 0,
        last_position: 0,
        finished_aligned_reads: false,
        next_id: 0,
        edge_intervals: if options.use_fast_algorithm {
            edge_intervals
        } else {
            None
        },
    };
    let cap = options.coverage_cap.max(0) as usize;
    let mut collector = Collector {
        records,
        references,
        options: options.clone(),
        result: WalkResult {
            high_quality_depth: vec![0; cap + 1],
            unfiltered_depth: vec![0; cap + 1],
            unfiltered_baseq: vec![0; 127],
            ..WalkResult::default()
        },
        counter: 0,
        previous_sequence: -1,
        counters: CounterManager::new(options.read_length.wrapping_mul(2000), options.read_length),
        read_names: None,
    };

    let mut counter: i64 = 0;
    while iterator.has_next()? {
        let locus = iterator.next()?.expect("has_next");
        let reference_n = collector.reference_base_n(locus.sequence, locus.position)?;
        if options.fast_collector {
            collector.add_info_fast(&locus, reference_n)?;
        } else {
            collector.add_info_slow(&locus, reference_n)?;
        }
        if reference_n {
            continue;
        }
        counter += 1;
        if collector.is_time_to_stop(counter) {
            break;
        }
        collector.counter = counter;
    }

    let mut result = collector.result;
    let filters = &iterator.sam.filters;
    result.excluded_by_adapter = filters.adapter_bases;
    result.excluded_by_mapq = filters.mapq_bases;
    result.excluded_by_dupe = filters.dupe_bases;
    result.excluded_by_pairing = filters.pair_bases;
    Ok(result)
}

/// `IntervalList.uniqued()` over `(sequence, start, end)`: sorted, then overlapping AND abutting
/// intervals merged.
pub fn uniqued(mut intervals: Vec<(i32, i32, i32)>) -> Vec<(i32, i32, i32)> {
    intervals.sort();
    let mut out: Vec<(i32, i32, i32)> = Vec::new();
    for interval in intervals {
        if let Some(last) = out.last_mut() {
            if last.0 == interval.0 && interval.1 <= last.2 + 1 {
                last.2 = last.2.max(interval.2);
                continue;
            }
        }
        out.push(interval);
    }
    out
}

// ------------------------------------------------------------------------------------------------
// WgsMetrics.
// ------------------------------------------------------------------------------------------------

/// The `WgsMetrics` columns, in declaration order.
pub const WGS_COLUMNS: [&str; 32] = [
    "GENOME_TERRITORY",
    "MEAN_COVERAGE",
    "SD_COVERAGE",
    "MEDIAN_COVERAGE",
    "MAD_COVERAGE",
    "PCT_EXC_ADAPTER",
    "PCT_EXC_MAPQ",
    "PCT_EXC_DUPE",
    "PCT_EXC_UNPAIRED",
    "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP",
    "PCT_EXC_CAPPED",
    "PCT_EXC_TOTAL",
    "PCT_1X",
    "PCT_5X",
    "PCT_10X",
    "PCT_15X",
    "PCT_20X",
    "PCT_25X",
    "PCT_30X",
    "PCT_40X",
    "PCT_50X",
    "PCT_60X",
    "PCT_70X",
    "PCT_80X",
    "PCT_90X",
    "PCT_100X",
    "FOLD_80_BASE_PENALTY",
    "FOLD_90_BASE_PENALTY",
    "FOLD_95_BASE_PENALTY",
    "HET_SNP_SENSITIVITY",
    "HET_SNP_Q",
];

/// One `WgsMetrics` row: `GENOME_TERRITORY` and then every double, in column order.
#[derive(Debug, Clone)]
pub struct WgsMetricsRow {
    pub genome_territory: i64,
    /// The 31 double columns from `MEAN_COVERAGE` to `HET_SNP_Q`.
    pub doubles: Vec<f64>,
}

fn histogram_of(array: &[i64]) -> htsjdk_metrics::histogram::Histogram {
    let mut histogram = htsjdk_metrics::histogram::Histogram::new("coverage", "count");
    for (i, &value) in array.iter().enumerate() {
        histogram.increment_by(i as f64, value as f64);
    }
    histogram
}

const LOG_ODDS_THRESHOLD: f64 = 3.0;

/// `CollectWgsMetrics.generateWgsMetrics` (the counting overload) and the `WgsMetrics`
/// constructor's `calculateDerivedFields`.
pub fn wgs_metrics(
    high_quality_depth: &[i64],
    unfiltered_depth: &[i64],
    unfiltered_baseq: &[i64],
    walk: &WalkResult,
    coverage_cap: i32,
    sample_size: i32,
) -> Result<WgsMetricsRow, Thrown> {
    let high_quality = histogram_of(high_quality_depth);
    let total = high_quality.sum();
    let total_with_excludes = total
        + walk.excluded_by_dupe as f64
        + walk.excluded_by_adapter as f64
        + walk.excluded_by_mapq as f64
        + walk.excluded_by_pairing as f64
        + walk.excluded_by_baseq as f64
        + walk.excluded_by_overlap as f64
        + walk.excluded_by_capping as f64;
    let pct = |count: i64| {
        if total_with_excludes == 0.0 {
            0.0
        } else {
            count as f64 / total_with_excludes
        }
    };
    let pct_total = if total_with_excludes == 0.0 {
        0.0
    } else {
        (total_with_excludes - total) / total_with_excludes
    };

    if sample_size <= 0 {
        return Err("picard.PicardException: Sample size is required when a baseQ histogram is given when deriving metrics.".to_string());
    }
    let mut depth_array = vec![0i64; coverage_cap as usize + 1];
    for (id, value) in high_quality.bins() {
        let depth = (id as i32).min(coverage_cap);
        depth_array[depth as usize] += value as i64;
    }
    let territory = high_quality.sum_of_values() as i64;
    let mean = high_quality.mean();
    let sd = high_quality.standard_deviation();
    let median = high_quality.median();
    let mad = high_quality.median_absolute_deviation();
    let at_least = |depth: usize| {
        let sum: i64 = depth_array.iter().skip(depth).sum();
        sum as f64 / territory as f64
    };
    let (fold80, fold90, fold95) = if high_quality.count() > 0.0 {
        let p = |q: f64| high_quality.percentile(q).unwrap_or(f64::NAN);
        (mean / p(0.2), mean / p(0.1), mean / p(0.05))
    } else {
        (0.0, 0.0, 0.0)
    };
    let depth_distribution = crate::theoretical_sensitivity::normalize(
        &unfiltered_depth
            .iter()
            .map(|&v| v as f64)
            .collect::<Vec<_>>(),
    );
    let baseq_distribution = crate::theoretical_sensitivity::normalize(
        &unfiltered_baseq
            .iter()
            .map(|&v| v as f64)
            .collect::<Vec<_>>(),
    );
    let sensitivity = crate::theoretical_sensitivity::het_snp_sensitivity(
        &depth_distribution,
        &baseq_distribution,
        sample_size,
        LOG_ODDS_THRESHOLD,
    )?;
    let q = crate::theoretical_sensitivity::phred_from_error_probability(1.0 - sensitivity);

    let mut doubles = vec![
        mean,
        sd,
        median,
        mad,
        pct(walk.excluded_by_adapter),
        pct(walk.excluded_by_mapq),
        pct(walk.excluded_by_dupe),
        pct(walk.excluded_by_pairing),
        pct(walk.excluded_by_baseq),
        pct(walk.excluded_by_overlap),
        pct(walk.excluded_by_capping),
        pct_total,
    ];
    for depth in [1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100] {
        doubles.push(at_least(depth));
    }
    doubles.extend([fold80, fold90, fold95, sensitivity, q as f64]);
    Ok(WgsMetricsRow {
        genome_territory: territory,
        doubles,
    })
}
