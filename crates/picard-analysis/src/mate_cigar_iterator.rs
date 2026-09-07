//! `MarkDuplicatesWithMateCigarIterator`: the sliding window `MarkDuplicatesWithMateCigar` marks
//! with, and the `MarkQueue` inside it.
//!
//! [`crate::mark_duplicates`] sorts every read end in the file and cuts sets on equal keys, which
//! needs two passes and the whole file in hand. This one makes ONE pass, and it can because the
//! `MC` tag tells it where the mate lies without waiting for the mate: a read end is complete the
//! moment it is read.
//!
//! The consequence is that duplicate marking happens in a WINDOW rather than over a sorted whole,
//! and the window is what makes the two tools disagree. `MarkQueue` holds at most one end per
//! comparable position; a second end at that position is compared against it there and then, and
//! the loser is a duplicate immediately. `MarkDuplicates` instead gathers every end at a position
//! and picks the best of the whole set at the end. Where three or more ends share a position and
//! the arrival order is not the score order, the two can keep a different read -- which is what
//! the covering array found on `read0208`.
//!
//! Three things decide when the queue is drained, and all three are the reference's:
//!
//!   * the minimum distance is `max(2 * <first read's length>, 100)`, set from the FIRST mapped
//!     read of the file and never revised;
//!   * an end leaves the queue when the current read's coordinate is more than that distance past
//!     the queue's head, or when the reference index changes, or at the end of the file;
//!   * leaving the queue is what makes an end a NON-duplicate, and taking a paired end out drags
//!     out every fragment and unpaired read at the same position as duplicates of it.
//!
//! Ported from `picard.sam.markduplicates.MarkDuplicatesWithMateCigarIterator`,
//! `picard.sam.markduplicates.util.MarkQueue` and
//! `picard.sam.markduplicates.util.ReadEndsForMateCigar` at tag 3.4.0.

use std::cmp::Ordering;

use crate::mark_duplicates::{
    duplicate_score_with, location, orientation_byte, Location, Options, Record, ScoringStrategy,
    F, FF, FR, R,
};

/// One read end, as `ReadEndsForMateCigar` builds it.
///
/// `read2_coordinate` comes from the mate's `MC` tag rather than from the mate itself, which is
/// the whole reason one pass is enough.
#[derive(Debug, Clone)]
struct End {
    library_id: i32,
    read1_reference_index: i32,
    read1_coordinate: i32,
    read2_reference_index: i32,
    read2_coordinate: i32,
    orientation: u8,
    orientation_for_optical: u8,
    has_unmapped: i32,
    /// The record's index in the file, which is how a decision reaches the output.
    index: usize,
    paired: bool,
    first_of_pair: bool,
    /// `SAMUtils.getCanonicalRecordName`, which is the read name: the tie-break that makes the
    /// choice deterministic when two ends score alike.
    name: String,
    score: i16,
    location: Location,
    read_group: i32,
    /// The index of this end's location set, for optical duplicate tracking.
    location_set: Option<usize>,
}

impl End {
    fn is_paired(&self) -> bool {
        self.read2_reference_index != -1
    }
}

/// `MarkQueue.MarkQueueComparator`: the comparator BOTH sets are keyed on, and the reason the
/// queue holds one end per position. It stops at the coordinates: two ends that differ only in
/// score are the same key, which is what makes "is there already a comparable record here?" mean
/// "is this position taken?".
fn compare_queue(left: &End, right: &End) -> Ordering {
    left.library_id
        .cmp(&right.library_id)
        .then(left.read1_reference_index.cmp(&right.read1_reference_index))
        .then(left.read1_coordinate.cmp(&right.read1_coordinate))
        // Reversed, to get pairs first, on the order `ReadEnds` defines.
        .then(right.orientation.cmp(&left.orientation))
        .then(left.read2_reference_index.cmp(&right.read2_reference_index))
        .then(left.read2_coordinate.cmp(&right.read2_coordinate))
}

/// `MarkQueue.ReadEndsMCComparator`: NOT a key, only a verdict. It is what `MarkQueue.add` calls
/// to decide which of two ends at one position is the better, and its tail -- paired first, then
/// the higher score, then the read name -- is the whole of that decision.
fn compare_ends(left: &End, right: &End, scoring: ScoringStrategy, records: &[Record]) -> Ordering {
    let coarse = left
        .library_id
        .cmp(&right.library_id)
        .then(left.read1_reference_index.cmp(&right.read1_reference_index))
        .then(left.read1_coordinate.cmp(&right.read1_coordinate))
        .then(right.orientation.cmp(&left.orientation));
    if coarse != Ordering::Equal {
        return coarse;
    }
    // Unpaired goes first, on `ReadEnds`'s own order.
    if left.is_paired() != right.is_paired() {
        return if left.is_paired() {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    let rest = left
        .has_unmapped
        .cmp(&right.has_unmapped)
        .then(left.read2_reference_index.cmp(&right.read2_reference_index))
        .then(left.read2_coordinate.cmp(&right.read2_coordinate));
    if rest != Ordering::Equal {
        return rest;
    }
    // `DuplicateScoringStrategy.compare(lhs, rhs, strategy, assumeMateCigar = true)`: paired
    // before unpaired, then the HIGHER score first, then the read name.
    let (left_record, right_record) = (&records[left.index], &records[right.index]);
    if left_record.paired() != right_record.paired() {
        return if left_record.paired() {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    let by_score = duplicate_score_with(right_record, scoring, true).cmp(&duplicate_score_with(
        left_record,
        scoring,
        true,
    ));
    if by_score != Ordering::Equal {
        return by_score;
    }
    left.name.cmp(&right.name)
}

/// A `TreeSet` under a comparator: a sorted vector whose keys are unique under that comparator, so
/// that "contains" and "the one comparable entry" mean what they mean in the Java.
struct SortedEnds {
    entries: Vec<End>,
}

impl SortedEnds {
    fn new() -> Self {
        SortedEnds {
            entries: Vec::new(),
        }
    }

    fn find(&self, wanted: &End, compare: &dyn Fn(&End, &End) -> Ordering) -> Result<usize, usize> {
        self.entries
            .binary_search_by(|entry| compare(entry, wanted))
    }

    fn get(&self, wanted: &End, compare: &dyn Fn(&End, &End) -> Ordering) -> Option<&End> {
        self.find(wanted, compare).ok().map(|at| &self.entries[at])
    }

    fn insert(&mut self, end: End, compare: &dyn Fn(&End, &End) -> Ordering) {
        match self.find(&end, compare) {
            // `TreeSet.add` on an entry the comparator calls equal keeps the one already there.
            Ok(_) => {}
            Err(at) => self.entries.insert(at, end),
        }
    }

    fn remove(&mut self, wanted: &End, compare: &dyn Fn(&End, &End) -> Ordering) -> Option<End> {
        match self.find(wanted, compare) {
            Ok(at) => Some(self.entries.remove(at)),
            Err(_) => None,
        }
    }

    fn poll_first(&mut self) -> Option<End> {
        if self.entries.is_empty() {
            None
        } else {
            Some(self.entries.remove(0))
        }
    }

    fn first(&self) -> Option<&End> {
        self.entries.first()
    }

    fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// What one pass decided: the duplicate flag and the optical flag of every record.
pub struct Decisions {
    pub duplicate: Vec<bool>,
    pub optical: Vec<bool>,
}

/// The whole of `MarkDuplicatesWithMateCigar`'s marking, over a coordinate-sorted file.
///
/// `skipped` are the records the caller has already decided to leave out (a pair with no `MC`
/// under `SKIP_PAIRS_WITH_NO_MATE_CIGAR`): they reach the output unmarked and take no part.
pub fn mark_with_queue(records: &[Record], options: &Options, skipped: &[usize]) -> Decisions {
    let mut duplicate = vec![false; records.len()];
    let mut optical = vec![false; records.len()];

    // `LibraryIdGenerator`: an id per library, in the order the libraries are first seen.
    let mut libraries: Vec<String> = Vec::new();
    let library_id = |libraries: &mut Vec<String>, name: &str| -> i32 {
        match libraries.iter().position(|known| known == name) {
            Some(at) => at as i32,
            None => {
                libraries.push(name.to_string());
                (libraries.len() - 1) as i32
            }
        }
    };

    // Both sets are keyed on the coarse comparator, so each position holds one end; the fine one
    // decides which end that is.
    let better = |left: &End, right: &End| compare_ends(left, right, options.scoring, records);

    let mut non_duplicate = SortedEnds::new();
    let mut other_end = SortedEnds::new();
    // Every location set ever made, so that an end can point at one by index.
    let mut location_sets: Vec<Vec<End>> = Vec::new();
    let mut minimum_distance: i32 = -1;
    let mut reference_index: i32 = -1;

    for (index, record) in records.iter().enumerate() {
        if skipped.contains(&index) {
            continue;
        }
        if record.unmapped() || record.secondary_or_supplementary() {
            // Not considered for duplicate marking at all; the queue never sees them.
            continue;
        }

        // "If not already set, this sets the minimum distance to twice the read length, or 100,
        // whichever is larger" -- from the FIRST mapped record, and never revised.
        if minimum_distance == -1 {
            minimum_distance = std::cmp::max(2 * record.qualities.len() as i32, 100);
        }

        let current = build_end(
            records,
            index,
            library_id(&mut libraries, &record.library),
            options,
        );

        // Drain everything the window has left behind, which is what makes those ends
        // non-duplicates.
        drain(
            &mut non_duplicate,
            &mut other_end,
            &mut location_sets,
            &mut duplicate,
            &mut optical,
            options,
            records,
            false,
            Some(&current),
            minimum_distance,
            reference_index,
        );
        reference_index = current.read1_reference_index;

        add(
            &mut non_duplicate,
            &mut other_end,
            &mut location_sets,
            &mut duplicate,
            &better,
            current,
        );
    }

    // End of file: everything still in the queue leaves it.
    drain(
        &mut non_duplicate,
        &mut other_end,
        &mut location_sets,
        &mut duplicate,
        &mut optical,
        options,
        records,
        true,
        None,
        minimum_distance,
        reference_index,
    );

    Decisions { duplicate, optical }
}

/// `new ReadEndsForMateCigar(...)`.
fn build_end(records: &[Record], index: usize, library_id: i32, options: &Options) -> End {
    let record = &records[index];
    let mut end = End {
        library_id,
        read1_reference_index: record.reference_index,
        read1_coordinate: record.five_prime_coordinate(),
        read2_reference_index: -1,
        read2_coordinate: -1,
        orientation: if record.reverse_strand() { R } else { F },
        orientation_for_optical: 0,
        has_unmapped: 0,
        index,
        paired: record.paired(),
        first_of_pair: record.first_of_pair(),
        name: record.name.clone(),
        score: duplicate_score_with(record, options.scoring, true),
        location: Location::default(),
        read_group: -1,
        location_set: None,
    };

    if record.paired() && !record.mate_unmapped() {
        end.read2_reference_index = record.mate_reference_index;
        // `SAMUtils.getMateUnclippedStart` / `getMateUnclippedEnd`, off the `MC` tag: the mate's
        // clipped alignment walked back over its own clips.
        end.read2_coordinate = mate_five_prime(record);
        end.orientation = orientation_byte(record.reverse_strand(), record.mate_reverse_strand());
        end.orientation_for_optical = if record.first_of_pair() {
            orientation_byte(record.reverse_strand(), record.mate_reverse_strand())
        } else {
            orientation_byte(record.mate_reverse_strand(), record.reverse_strand())
        };
    }

    if record.unmapped() || (record.paired() && record.mate_unmapped()) {
        end.has_unmapped = 1;
    }
    if options.parse_read_names {
        end.location = location(&record.name);
        if end.location.known {
            end.read_group = record.read_group;
        }
    }
    end
}

/// The mate's 5' coordinate, from its `MC` cigar and its alignment start.
fn mate_five_prime(record: &Record) -> i32 {
    let Some(cigar) = &record.mate_cigar else {
        // Without the tag there is nothing to read; the caller has already refused or skipped.
        return -1;
    };
    let start = record.mate_alignment_start;
    if record.mate_reverse_strand() {
        let mut end = start + cigar.reference_length() as i32 - 1;
        for element in cigar.elements.iter().rev() {
            match element.op {
                htsjdk_bam::cigar::Op::S | htsjdk_bam::cigar::Op::H => end += element.length as i32,
                _ => break,
            }
        }
        end
    } else {
        let mut begin = start;
        for element in &cigar.elements {
            match element.op {
                htsjdk_bam::cigar::Op::S | htsjdk_bam::cigar::Op::H => {
                    begin -= element.length as i32
                }
                _ => break,
            }
        }
        begin
    }
}

/// `MarkQueue.shouldBeInLocations`: only a mapped pair is tracked for optical duplicates.
fn should_be_in_locations(end: &End) -> bool {
    end.paired && end.has_unmapped == 0
}

/// `MarkQueue.add`, which is where a duplicate is decided the moment a second end arrives at a
/// position rather than after the whole set has been gathered.
#[allow(clippy::too_many_arguments)]
fn add(
    non_duplicate: &mut SortedEnds,
    other_end: &mut SortedEnds,
    location_sets: &mut Vec<Vec<End>>,
    duplicate: &mut [bool],
    better: &dyn Fn(&End, &End) -> Ordering,
    mut other: End,
) {
    let mut location_set: Option<usize> = None;
    let mut add_to_location_set = true;
    let mut found_duplicate: Option<End> = None;

    if let Some(current) = non_duplicate.get(&other, &compare_queue).cloned() {
        // The position is taken. Which of the two ends keeps it is the fine comparator's verdict.
        let comparison = better(&current, &other);
        if current.name == other.name {
            // The two ends of one pair at one position: keep the better, remember the other.
            if comparison == Ordering::Greater {
                non_duplicate.remove(&current, &compare_queue);
                non_duplicate.insert(other.clone(), &compare_queue);
                other_end.insert(current.clone(), &compare_queue);
                if should_be_in_locations(&other) {
                    if let Some(set) = current.location_set {
                        // `replace`: the set stays, its identifier changes.
                        if let Some(entries) = location_sets.get_mut(set) {
                            entries.retain(|entry| entry.index != current.index);
                            entries.push(other.clone());
                        }
                        location_set = Some(set);
                        other.location_set = Some(set);
                        add_to_location_set = false;
                    }
                }
            } else {
                other_end.insert(other.clone(), &compare_queue);
                if should_be_in_locations(&current) {
                    location_set = current.location_set;
                    add_to_location_set = false;
                }
            }
        } else if comparison == Ordering::Greater {
            // "other" is the better end: it takes the position, and "current" is a duplicate.
            let set = if should_be_in_locations(&current) {
                current.location_set
            } else {
                location_sets.push(Vec::new());
                Some(location_sets.len() - 1)
            };
            let set = set.unwrap_or_else(|| {
                location_sets.push(Vec::new());
                location_sets.len() - 1
            });
            other.location_set = Some(set);
            location_set = Some(set);
            non_duplicate.remove(&current, &compare_queue);
            non_duplicate.insert(other.clone(), &compare_queue);

            // The loser's own pair, if it was being held, is not a duplicate.
            if let Some(pair) = other_end.get(&current, &compare_queue).cloned() {
                other_end.remove(&current, &compare_queue);
                duplicate[pair.index] = true;
            }
            found_duplicate = Some(current);
        } else {
            // "current" keeps the position; "other" is the duplicate.
            if should_be_in_locations(&current) {
                location_set = current.location_set;
            }
            found_duplicate = Some(other.clone());
        }
    } else {
        // The first end at this position.
        if should_be_in_locations(&other) {
            location_sets.push(Vec::new());
            let set = location_sets.len() - 1;
            other.location_set = Some(set);
            location_set = Some(set);
        }
        non_duplicate.insert(other.clone(), &compare_queue);
    }

    if other.paired && other.has_unmapped == 0 && add_to_location_set {
        if let Some(set) = location_set {
            location_sets[set].push(other.clone());
        }
    }

    if let Some(duplicated) = found_duplicate {
        duplicate[duplicated.index] = true;
    }
}

/// `MarkDuplicatesWithMateCigarIterator.tryPollingTheToMarkQueue` and `MarkQueue.poll` together:
/// take from the queue every end the window has left behind, and with each one every fragment and
/// unpaired read that shared its position.
#[allow(clippy::too_many_arguments)]
fn drain(
    non_duplicate: &mut SortedEnds,
    other_end: &mut SortedEnds,
    location_sets: &mut [Vec<End>],
    duplicate: &mut [bool],
    optical: &mut [bool],
    options: &Options,
    records: &[Record],
    flush: bool,
    current: Option<&End>,
    minimum_distance: i32,
    reference_index: i32,
) {
    loop {
        if non_duplicate.is_empty() {
            return;
        }
        if !flush {
            let Some(current) = current else { return };
            let head = non_duplicate.first().expect("not empty");
            let far_enough = minimum_distance < current.read1_coordinate - head.read1_coordinate;
            if reference_index == current.read1_reference_index && !far_enough {
                return;
            }
        }

        let polled = non_duplicate.poll_first().expect("not empty");

        if polled.is_paired() {
            // Its own other end, if held, is not a duplicate either.
            other_end.remove(&polled, &compare_queue);

            // "Remove from the set fragments and unpaired, which only have two possible
            // orientations": a mapped pair leaving the queue takes them with it.
            let mut probe = polled.clone();
            probe.read2_reference_index = -1;
            probe.read2_coordinate = -1;
            probe.orientation = if matches!(polled.orientation, FF | FR | F) {
                F
            } else {
                R
            };
            // The probe stands for a position, not a record, so it must not win the score
            // tie-break against what is there.
            if let Some(found) = non_duplicate.get(&probe, &compare_queue).cloned() {
                duplicate[found.index] = true;
                non_duplicate.remove(&probe, &compare_queue);
            }
        }

        // `trackOpticalDuplicates` over this end's location set, for the first end of each pair.
        if should_be_in_locations(&polled) && polled.first_of_pair && options.parse_read_names {
            if let Some(set) = polled.location_set {
                let entries = &location_sets[set];
                if !entries.is_empty() {
                    track_optical(entries, options, records, optical);
                }
            }
        }
    }
}

/// `AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates` over one location set: the
/// ends within the pixel distance of the best-scoring one are optical duplicates of it.
fn track_optical(entries: &[End], options: &Options, records: &[Record], optical: &mut [bool]) {
    let known: Vec<&End> = entries.iter().filter(|end| end.location.known).collect();
    if known.len() < 2 {
        return;
    }
    // The keeper is the highest score, first one winning, as `MarkDuplicates` picks it.
    let mut keeper = 0;
    for (position, end) in known.iter().enumerate() {
        if end.score > known[keeper].score {
            keeper = position;
        }
    }
    for (position, end) in known.iter().enumerate() {
        if position == keeper {
            continue;
        }
        let other = known[keeper];
        if end.read_group == other.read_group
            && end.location.tile == other.location.tile
            && (end.location.x - other.location.x).abs() <= options.optical_duplicate_pixel_distance
            && (end.location.y - other.location.y).abs() <= options.optical_duplicate_pixel_distance
        {
            optical[end.index] = true;
            let _ = records;
        }
    }
}
