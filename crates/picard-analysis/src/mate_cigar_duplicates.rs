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

use crate::mark_duplicates::{Marking, Options, Record};

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
    /// `PicardException`, thrown by `MarkDuplicatesWithMateCigarIterator`'s constructor: this tool
    /// refuses one of the three scoring strategies outright, because the two ends of a pair are
    /// scored separately here and summing base qualities would let them disagree.
    SumOfBaseQualitiesUnsupported,
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
            // htsjdk appends the RECORD and not its name: `SAMRecord.toString` renders
            // `<name> <1|2>/2 <length>b aligned to <contig>:<start>-<end>.`. This port carries the
            // name alone, because a `Record` here has no contig NAME to render -- the difference
            // is visible in the covering array and is tracked rather than hidden.
            Refusal::MateCigarNotFound { read } => {
                format!("Mate CIGAR (Tag MC) not found: {read}")
            }
            Refusal::SumOfBaseQualitiesUnsupported => {
                "SUM_OF_BASE_QUALITIES not supported as this \
                 may cause inconsistencies across ends in a pair.  Please use a different scoring \
                 strategy."
                    .to_string()
            }
        }
    }

    /// The exception class the reference throws, which is not the same for the two tools.
    pub fn exception(&self) -> &'static str {
        match self {
            Refusal::NotCoordinateSorted
            | Refusal::NoMateCigar { .. }
            | Refusal::SumOfBaseQualitiesUnsupported => "picard.PicardException",
            Refusal::MateCigarNotFound { .. } => "htsjdk.samtools.SAMException",
        }
    }
}

/// `SAMRecord.toString`, which is what htsjdk appends to the refusal rather than the read name:
/// `<name> <1|2>/2 <length>b aligned to <contig>:<start>-<end>.`, or `unmapped read.`.
///
/// The contig is the record's own reference name, which this crate's `Record` does not carry, so
/// the caller supplies it through [`describe_with_contig`]; without one the index is printed,
/// which no reference output ever contains and so cannot be mistaken for a match.
/// [`describe`] with the contig name the header gives.
pub fn describe_with_contig(record: &Record, contig: Option<&str>) -> String {
    let mut out = record.name.clone();
    if record.paired() {
        out.push_str(if record.first_of_pair() {
            " 1/2"
        } else {
            " 2/2"
        });
    }
    out.push(' ');
    out.push_str(&record.qualities.len().to_string());
    out.push('b');
    if record.unmapped() {
        out.push_str(" unmapped read.");
    } else {
        let name = contig
            .map(str::to_string)
            .unwrap_or_else(|| record.reference_index.to_string());
        let end = record.alignment_start + record.cigar.reference_length() as i32 - 1;
        out.push_str(&format!(
            " aligned to {name}:{}-{end}.",
            record.alignment_start
        ));
    }
    out
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
    // The iterator's constructor refuses the strategy before it reads a record, and it refuses it
    // whatever the file holds.
    if options.base.scoring == crate::mark_duplicates::ScoringStrategy::SumOfBaseQualities {
        return Err(Refusal::SumOfBaseQualitiesUnsupported);
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
    // The mate-cigar path scores a pair from one end, which the base options do not say.
    let base = Options {
        assume_mate_cigar: true,
        ..options.base.clone()
    };
    // The marking itself is the iterator's, not `MarkDuplicates`'s: a sliding window with a queue
    // that decides a duplicate as the second end arrives, rather than a sorted whole cut into
    // sets. The two disagree wherever three ends share a position and arrival order is not score
    // order (picard-rs #307).
    let decisions = crate::mate_cigar_iterator::mark_with_queue(records, &base, &skipped);
    Ok(crate::mark_duplicates::marking_from(
        records,
        &base,
        &decisions.duplicate,
        &decisions.optical,
    ))
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
    contigs: &[String],
) -> Result<Marking, Refusal> {
    if order != SortOrder::Coordinate {
        return Err(Refusal::NotCoordinateSorted);
    }
    // The refusal is htsjdk's, from `SAMUtils.getMateUnclippedStart`/`End`, and WHICH record it
    // names is decided by the sort that `DuplicateSetIterator` runs before it cuts a single set.
    // Java's sort begins by comparing `a[1]` against `a[0]`, in that order, and the comparator
    // resolves the mate coordinate of its first argument first -- so the record named is the
    // SECOND of the file, not the first. The order below is that examination order; a file whose
    // failure lies deeper than the first comparison would need the whole of TimSort to predict.
    let examination = [1usize, 0]
        .into_iter()
        .chain(2..records.len())
        .filter(|index| *index < records.len());
    for index in examination {
        let record = &records[index];
        if needs_mate_cigar(record) && mate_cigar(record).is_none() {
            return Err(Refusal::MateCigarNotFound {
                read: describe_with_contig(
                    record,
                    contigs
                        .get(record.reference_index.max(0) as usize)
                        .map(String::as_str),
                ),
            });
        }
    }
    // The marking is htsjdk's duplicate-set iterator, not `MarkDuplicates`'s sorted whole: this
    // tool is a `MarkDuplicates` subclass whose iterator does the grouping for it.
    let base = Options {
        assume_mate_cigar: true,
        ..options.clone()
    };
    let duplicate = crate::duplicate_set::mark_duplicate_sets(records, &base);
    let optical = vec![false; records.len()];
    Ok(crate::mark_duplicates::marking_from(
        records, &base, &duplicate, &optical,
    ))
}
