//! The two mate-cigar duplicate markers: `MarkDuplicatesWithMateCigar` and its simple sibling.
//!
//! [`crate::mark_duplicates`] decides a pair's position from the two ends it has seen, so it holds
//! every unpaired end until its mate arrives and needs two passes over the file. These two read
//! the other end's position off the `MC` tag instead, which makes one pass over a
//! coordinate-sorted file enough.
//!
//! # What the measurement says
//!
//! Seven of the golden's eleven inputs come out identical under all three tools, soft-clipped
//! mates and soft-clipped first ends included: reading the mate's cigar does not change the answer
//! when the two ends are in the file anyway. What changes is what happens when they are not, and
//! that is the whole of this module:
//!
//!  * **a pair with no `MC` is SKIPPED**, and skipping it removes it from its set, so a set of two
//!    pairs where one lacks the tag marks NEITHER. `MarkDuplicates` on the same file marks one;
//!  * **`SKIP_PAIRS_WITH_NO_MATE_CIGAR=false` is a refusal**, not a second algorithm: a
//!    `PicardException` naming the read;
//!  * **the simple one refuses the same file outright**, with htsjdk's own wording rather than
//!    Picard's, because it asks the record for a mate cigar it does not have;
//!  * **and both refuse a queryname-sorted file** that `MarkDuplicates` accepts.
//!
//! # What is not ported
//!
//! `MINIMUM_DISTANCE` is the width of the window `MarkDuplicatesWithMateCigar` buffers records in,
//! and no fixture separates a run that used it from one that did not: the golden's distant-mate
//! case is marked the same at the default window and at three thousand. It is an argument this
//! port accepts and does not act on, which is said here rather than implied by its absence. The
//! same goes for `BLOCK_SIZE`, which is a buffer size.
//!
//! Ported from `picard.sam.markduplicates.MarkDuplicatesWithMateCigar` and
//! `picard.sam.markduplicates.SimpleMarkDuplicatesWithMateCigar` in Picard 3.4.0.

use htsjdk_bam::cigar::Cigar;

use crate::mark_duplicates::{mark, Marking, Options, Record};

/// The sort order the header declares, which both tools check before anything else.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Queryname,
    Unsorted,
}

/// `MarkDuplicatesWithMateCigar`'s own arguments, beside the ones it shares.
#[derive(Debug, Clone)]
pub struct MateCigarOptions {
    pub base: Options,
    /// `SKIP_PAIRS_WITH_NO_MATE_CIGAR`, true by default.
    pub skip_pairs_with_no_mate_cigar: bool,
    /// `MINIMUM_DISTANCE`, accepted and not acted on: see the module's note.
    pub minimum_distance: i32,
}

impl Default for MateCigarOptions {
    fn default() -> Self {
        Self {
            base: Options::default(),
            skip_pairs_with_no_mate_cigar: true,
            minimum_distance: -1,
        }
    }
}

/// What either tool refuses a file for, in the words it refuses it with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// `PicardException`, thrown by both before a record is read.
    NotCoordinateSorted,
    /// `PicardException`, thrown by `MarkDuplicatesWithMateCigar` with the skip turned off.
    NoMateCigar { read: String },
    /// `SAMException`, thrown by the simple one whatever the skip says, because it asks htsjdk for
    /// the mate cigar and htsjdk is the one that refuses.
    MateCigarNotFound { read: String },
}

impl Refusal {
    /// The message, which is the reference's and not this port's.
    pub fn message(&self) -> String {
        match self {
            Refusal::NotCoordinateSorted => {
                "This program requires inputs in coordinate SortOrder".to_string()
            }
            Refusal::NoMateCigar { read } => format!(
                "Read {read} was mapped and had a mapped mate, but no mate cigar (\"MC\") tag."
            ),
            Refusal::MateCigarNotFound { read } => {
                format!("Mate CIGAR (Tag MC) not found: {read}")
            }
        }
    }

    /// The exception class the reference throws, which is not the same for the two tools.
    pub fn exception(&self) -> &'static str {
        match self {
            Refusal::NotCoordinateSorted | Refusal::NoMateCigar { .. } => "picard.PicardException",
            Refusal::MateCigarNotFound { .. } => "htsjdk.samtools.SAMException",
        }
    }
}

/// A record's mate cigar, where it carries one.
pub fn mate_cigar(record: &Record) -> Option<&Cigar> {
    record.mate_cigar.as_ref()
}

/// Whether a record is a mapped read with a mapped mate, which is the only shape that needs `MC`.
pub fn needs_mate_cigar(record: &Record) -> bool {
    record.paired() && !record.unmapped() && !record.mate_unmapped()
}

/// `MarkDuplicatesWithMateCigar.doWork`, over records already in memory.
///
/// The skip is the algorithm's, not a convenience: a pair without the tag leaves the run entirely,
/// which is why a set of two pairs where one lacks it marks neither. The metrics still count it,
/// because they are counted in the writing pass over every record.
pub fn mark_with_mate_cigar(
    records: &[Record],
    order: SortOrder,
    options: &MateCigarOptions,
) -> Result<Marking, Refusal> {
    if order != SortOrder::Coordinate {
        return Err(Refusal::NotCoordinateSorted);
    }
    let mut skipped: Vec<usize> = Vec::new();
    for (index, record) in records.iter().enumerate() {
        if !needs_mate_cigar(record) || mate_cigar(record).is_some() {
            continue;
        }
        if !options.skip_pairs_with_no_mate_cigar {
            return Err(Refusal::NoMateCigar {
                read: record.name.clone(),
            });
        }
        skipped.push(index);
    }
    Ok(marking_without(records, &skipped, &options.base))
}

/// `SimpleMarkDuplicatesWithMateCigar.doWork`, over records already in memory.
///
/// It is a `MarkDuplicates` subclass driven by htsjdk's duplicate-set iterator, and that iterator
/// is what refuses a record with no mate cigar: the refusal is htsjdk's `SAMException` and carries
/// htsjdk's wording, whatever `SKIP_PAIRS_WITH_NO_MATE_CIGAR` says.
pub fn simple_mark_with_mate_cigar(
    records: &[Record],
    order: SortOrder,
    options: &Options,
) -> Result<Marking, Refusal> {
    if order != SortOrder::Coordinate {
        return Err(Refusal::NotCoordinateSorted);
    }
    if let Some(record) = records
        .iter()
        .find(|record| needs_mate_cigar(record) && mate_cigar(record).is_none())
    {
        return Err(Refusal::MateCigarNotFound {
            read: record.name.clone(),
        });
    }
    Ok(mark(records, options))
}

/// The marking of a file with some pairs taken out of it, and the metrics of the whole.
///
/// The skipped records are removed before the sets are cut, so they change what the sets hold, and
/// they are put back unmarked for the writing pass: the reference writes them out untouched.
fn marking_without(records: &[Record], skipped: &[usize], options: &Options) -> Marking {
    if skipped.is_empty() {
        return mark(records, options);
    }
    let kept: Vec<Record> = records
        .iter()
        .enumerate()
        .filter(|(index, _)| !skipped.contains(index))
        .map(|(_, record)| record.clone())
        .collect();
    let inner = mark(&kept, options);
    // The metrics are counted over every record, skipped ones included, because the writing pass
    // walks the whole file.
    let whole = mark(records, options);

    let mut duplicate = vec![false; records.len()];
    let mut optical = vec![false; records.len()];
    let mut duplicate_type = vec![None; records.len()];
    let mut written = vec![true; records.len()];
    let mut position = 0;
    for index in 0..records.len() {
        if skipped.contains(&index) {
            continue;
        }
        duplicate[index] = inner.duplicate[position];
        optical[index] = inner.optical[position];
        duplicate_type[index] = inner.duplicate_type[position].clone();
        written[index] = inner.written[position];
        position += 1;
    }
    let mut metrics = whole.metrics;
    for row in &mut metrics {
        // The duplicate counters are the inner run's; the examined ones are the whole file's.
        let inner_row = inner
            .metrics
            .iter()
            .find(|other| other.library == row.library);
        row.unpaired_read_duplicates = inner_row.map(|r| r.unpaired_read_duplicates).unwrap_or(0);
        row.read_pair_duplicates = inner_row.map(|r| r.read_pair_duplicates).unwrap_or(0);
        row.read_pair_optical_duplicates = inner_row
            .map(|r| r.read_pair_optical_duplicates)
            .unwrap_or(0);
        row.estimated_library_size = crate::mark_duplicates::estimate_library_size(
            row.read_pairs_examined - row.read_pair_optical_duplicates,
            row.read_pairs_examined - row.read_pair_duplicates,
        );
        let examined = row.unpaired_reads_examined + row.read_pairs_examined;
        row.percent_duplication = if examined == 0 {
            0.0
        } else {
            (row.unpaired_read_duplicates + row.read_pair_duplicates * 2) as f64
                / (row.unpaired_reads_examined + row.read_pairs_examined * 2) as f64
        };
    }
    Marking {
        duplicate,
        optical,
        duplicate_type,
        written,
        metrics,
        all_sets: inner.all_sets,
        optical_sets: inner.optical_sets,
        non_optical_sets: inner.non_optical_sets,
    }
}
