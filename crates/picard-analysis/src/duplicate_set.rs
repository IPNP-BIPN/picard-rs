//! htsjdk's `DuplicateSetIterator`: the third way this repository marks duplicates.
//!
//! `MarkDuplicates` sorts every read end and cuts sets on equal keys.
//! [`crate::mate_cigar_iterator`] keeps a window and decides as ends arrive. This one walks the
//! coordinate-sorted file once and cuts a new set wherever the record it is looking at is no
//! longer comparable to the set's REPRESENTATIVE, then marks inside the set by read name.
//!
//! `SimpleMarkDuplicatesWithMateCigar` is `MarkDuplicates`'s subclass driving this, which is why
//! it can differ from both of the others on the same file.
//!
//! Two comparators do the work, and the difference between them is one boolean each:
//!
//!   * `duplicateSetCompare` collapses the orientation, so a fragment on the forward strand is
//!     comparable to a pair whose first end is forward, and it does not look at how many ends are
//!     mapped. That is set MEMBERSHIP;
//!   * `compare` keeps the orientations apart, does look at the ends, and then continues into the
//!     duplicate score (with the mate cigar), the read name and which end of the pair it is. That
//!     is ORDER, and its first record is the representative.
//!
//! Ported from `htsjdk.samtools.DuplicateSetIterator`, `htsjdk.samtools.DuplicateSet` and
//! `htsjdk.samtools.SAMRecordDuplicateComparator` at tag 4.2.0.

use std::cmp::Ordering;

use crate::mark_duplicates::{duplicate_score_with, Options, Record, ScoringStrategy};

/// `SAMRecordDuplicateComparator`'s own orientation constants, which are NOT `ReadEnds`': here a
/// fragment sits BETWEEN the paired orientations rather than before them.
const FF: u8 = 0;
const FR: u8 = 1;
const F: u8 = 2;
const RF: u8 = 3;
const RR: u8 = 4;
const R: u8 = 5;

fn paired_and_both_mapped(record: &Record) -> bool {
    record.paired() && !record.unmapped() && !record.mate_unmapped()
}

fn has_unmapped_end(record: &Record) -> bool {
    record.unmapped() || (record.paired() && record.mate_unmapped())
}

fn has_mapped_end(record: &Record) -> bool {
    !record.unmapped() || (record.paired() && !record.mate_unmapped())
}

/// `Attr.ReadCoordinate`: the unclipped end the strand points away from.
fn read_coordinate(record: &Record) -> i32 {
    if record.reverse_strand() {
        record.unclipped_end()
    } else {
        record.unclipped_start()
    }
}

/// `Attr.MateCoordinate`, off the `MC` tag, and `-1` where there is no mapped mate.
fn mate_coordinate(record: &Record) -> i32 {
    if !paired_and_both_mapped(record) {
        return -1;
    }
    let Some(cigar) = &record.mate_cigar else {
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

fn mate_reference_index(record: &Record) -> i32 {
    if paired_and_both_mapped(record) {
        record.mate_reference_index
    } else {
        -1
    }
}

/// `getPairedOrientation`: the pair's two strands, or the fragment's one.
fn paired_orientation(record: &Record) -> u8 {
    if paired_and_both_mapped(record) {
        match (record.reverse_strand(), record.mate_reverse_strand()) {
            (true, true) => RR,
            (true, false) => RF,
            (false, true) => FR,
            (false, false) => FF,
        }
    } else if record.reverse_strand() {
        R
    } else {
        F
    }
}

/// `compareOrientationByteCollapseOrientation`: a FRAGMENT is comparable to any pair whose own
/// first end points the same way, which is how a fragment joins a pair's set.
fn compare_orientation_collapsed(left: u8, right: u8) -> Ordering {
    if left == F || left == R {
        if left == F && matches!(right, F | FR | FF) {
            return Ordering::Equal;
        }
        if left == R && matches!(right, R | RF | RR) {
            return Ordering::Equal;
        }
    } else if right == F || right == R {
        return compare_orientation_collapsed(right, left).reverse();
    }
    left.cmp(&right)
}

/// `SAMRecordDuplicateComparator.fileOrderCompare(a, b, collapseOrientation, considerEnds)`.
fn file_order_compare(
    left: &Record,
    right: &Record,
    left_library: i32,
    right_library: i32,
    collapse_orientation: bool,
    consider_ends: bool,
) -> Ordering {
    let mut cmp = left_library.cmp(&right_library);

    if cmp == Ordering::Equal {
        // An unmapped record sorts LAST here, as it does everywhere in htsjdk.
        cmp = match (left.reference_index, right.reference_index) {
            (-1, -1) => Ordering::Equal,
            (-1, _) => Ordering::Greater,
            (_, -1) => Ordering::Less,
            (a, b) => a.cmp(&b),
        };
    }
    if cmp == Ordering::Equal {
        cmp = read_coordinate(left).cmp(&read_coordinate(right));
    }
    if cmp == Ordering::Equal {
        let (a, b) = (paired_orientation(left), paired_orientation(right));
        cmp = if collapse_orientation {
            compare_orientation_collapsed(a, b)
        } else {
            a.cmp(&b)
        };
    }
    if paired_and_both_mapped(left) && paired_and_both_mapped(right) {
        if cmp == Ordering::Equal {
            cmp = mate_reference_index(left).cmp(&mate_reference_index(right));
        }
        if cmp == Ordering::Equal {
            cmp = mate_coordinate(left).cmp(&mate_coordinate(right));
        }
    }
    if cmp == Ordering::Equal {
        cmp = (!has_mapped_end(left)).cmp(&!has_mapped_end(right));
    }
    if cmp == Ordering::Equal && consider_ends {
        cmp = if left.paired() == right.paired() {
            has_unmapped_end(left).cmp(&has_unmapped_end(right))
        } else if left.paired() {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    cmp
}

/// `duplicateSetCompare`: membership, which collapses the orientation and ignores the ends.
fn duplicate_set_compare(
    left: &Record,
    right: &Record,
    left_library: i32,
    right_library: i32,
) -> Ordering {
    file_order_compare(left, right, left_library, right_library, true, false)
}

/// `SAMRecordDuplicateComparator.compare`: order within a set, and so which record represents it.
///
/// It is also the comparator a writer uses for `SO:duplicate`, which is why it is public: a tool
/// that writes a file in that order sorts by exactly this.
pub fn compare(
    left: &Record,
    right: &Record,
    left_library: i32,
    right_library: i32,
    scoring: ScoringStrategy,
) -> Ordering {
    let mut cmp = file_order_compare(left, right, left_library, right_library, false, true);
    if cmp == Ordering::Equal {
        cmp = duplicate_score_with(right, scoring, true)
            .cmp(&duplicate_score_with(left, scoring, true));
    }
    if cmp == Ordering::Equal {
        cmp = left.name.cmp(&right.name);
    }
    if cmp == Ordering::Equal && left.paired() && right.paired() {
        cmp = (!left.first_of_pair()).cmp(&!right.first_of_pair());
    }
    cmp
}

/// Which records are duplicates, under htsjdk's duplicate sets.
///
/// A set is closed when the next record is no longer comparable to the representative, or when the
/// representative itself is unmapped or secondary -- such a record is a set of its own. Inside a
/// closed set, every record whose read NAME differs from the representative's is a duplicate, and
/// the first record after sorting never is: the marking is by TEMPLATE, so both ends of the
/// representative's pair survive.
/// One set's records in the order `DuplicateSet.getRecords` returns them, which is the FULL
/// comparator's and not the file's.
fn sorted_set(
    set: &[usize],
    _representative: usize,
    records: &[Record],
    library_of: &[i32],
    options: &Options,
) -> Vec<usize> {
    let mut sorted = set.to_vec();
    sorted.sort_by(|a, b| {
        compare(
            &records[*a],
            &records[*b],
            library_of[*a],
            library_of[*b],
            options.scoring,
        )
    });
    sorted
}

/// The duplicate sets themselves, as indices into `records`, in the order the iterator yields
/// them and with each set's records in the order it returns them.
///
/// `mark_duplicate_sets` is the same walk with a verdict written at the end of each set; a tool
/// that reads the SETS -- `CollectUmiPrevalenceMetrics` counts the distinct barcodes in one --
/// needs the grouping and not the flags.
pub fn duplicate_sets(records: &[Record], options: &Options) -> Vec<Vec<usize>> {
    let mut libraries: Vec<String> = Vec::new();
    let mut library_of: Vec<i32> = Vec::with_capacity(records.len());
    for record in records {
        let id = match libraries.iter().position(|known| *known == record.library) {
            Some(at) => at as i32,
            None => {
                libraries.push(record.library.clone());
                (libraries.len() - 1) as i32
            }
        };
        library_of.push(id);
    }
    // The same whole-file re-sort `mark_duplicate_sets` does, by the FULL comparator: the
    // iterator is built with `preSorted = false`.
    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_by(|a, b| {
        compare(
            &records[*a],
            &records[*b],
            library_of[*a],
            library_of[*b],
            options.scoring,
        )
    });

    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut set: Vec<usize> = Vec::new();
    let mut representative: usize = 0;
    for index in order {
        let record = &records[index];
        if set.is_empty() {
            set.push(index);
            representative = index;
            continue;
        }
        let head = &records[representative];
        let same = !head.unmapped()
            && !head.secondary_or_supplementary()
            && duplicate_set_compare(head, record, library_of[representative], library_of[index])
                == Ordering::Equal;
        if same {
            if compare(
                head,
                record,
                library_of[representative],
                library_of[index],
                options.scoring,
            ) == Ordering::Greater
            {
                representative = index;
            }
            set.push(index);
        } else {
            out.push(sorted_set(
                &set,
                representative,
                records,
                &library_of,
                options,
            ));
            set.clear();
            set.push(index);
            representative = index;
        }
    }
    if !set.is_empty() {
        out.push(sorted_set(
            &set,
            representative,
            records,
            &library_of,
            options,
        ));
    }
    out
}

pub fn mark_duplicate_sets(records: &[Record], options: &Options) -> Vec<bool> {
    let mut duplicate = vec![false; records.len()];
    let mut libraries: Vec<String> = Vec::new();
    let mut library_of: Vec<i32> = Vec::with_capacity(records.len());
    for record in records {
        let id = match libraries.iter().position(|known| *known == record.library) {
            Some(at) => at as i32,
            None => {
                libraries.push(record.library.clone());
                (libraries.len() - 1) as i32
            }
        };
        library_of.push(id);
    }

    // `new DuplicateSetIterator(..., preSorted = false, ...)`: the WHOLE file is re-sorted by the
    // duplicate comparator before a single set is cut. That is not the file's own order -- the
    // comparator keys on the UNCLIPPED coordinate, so a soft-clipped read moves -- and cutting
    // sets in file order marks a different file.
    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_by(|a, b| {
        compare(
            &records[*a],
            &records[*b],
            library_of[*a],
            library_of[*b],
            options.scoring,
        )
    });

    let mut set: Vec<usize> = Vec::new();
    let mut representative: usize = 0;

    let close = |set: &mut Vec<usize>, representative: usize, duplicate: &mut Vec<bool>| {
        if set.is_empty() {
            return;
        }
        let sorted = sorted_set(set, representative, records, &library_of, options);
        let name = records[representative].name.clone();
        for index in &sorted {
            let record = &records[*index];
            if !record.unmapped() && !record.secondary_or_supplementary() && record.name != name {
                duplicate[*index] = true;
            }
        }
        // `records.get(0).setDuplicateReadFlag(false)`, after the sort.
        duplicate[sorted[0]] = false;
        set.clear();
    };

    for index in order {
        let record = &records[index];
        if set.is_empty() {
            set.push(index);
            representative = index;
            continue;
        }
        let head = &records[representative];
        if head.unmapped() || head.secondary_or_supplementary() {
            close(&mut set, representative, &mut duplicate);
            set.push(index);
            representative = index;
            continue;
        }
        let cmp =
            duplicate_set_compare(head, record, library_of[representative], library_of[index]);
        if cmp == Ordering::Equal {
            if compare(
                head,
                record,
                library_of[representative],
                library_of[index],
                options.scoring,
            ) == Ordering::Greater
            {
                representative = index;
            }
            set.push(index);
        } else {
            close(&mut set, representative, &mut duplicate);
            set.push(index);
            representative = index;
        }
    }
    close(&mut set, representative, &mut duplicate);

    duplicate
}
