//! `SinglePassSamProgram.makeItSo` for the metrics tools that extend it and need nothing from the
//! reference but the walk: the input is read, the sort order is checked, the reference walker is
//! asked for every mapped record's contig, and the records are handed to the tool until
//! `STOP_AFTER` of them have been seen.
//!
//! The order is the reference's own: a sort order that is not coordinate is refused unless
//! `ASSUME_SORTED` (which defaults to TRUE for these tools), and taking that escape moves the
//! refusal to the walker when the input is not in contig order.
//!
//! Ported from `picard.analysis.SinglePassSamProgram` in Picard 3.4.0.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::metrics_cli::{check_coordinate_sorted, thrown, Args, ReferenceWalker};

/// `makeItSo`, for a tool whose `usesNoRefReads()` is true (it reads past the mapped records into
/// the unmapped tail). An `Err` from `accept` is a Java throwable, class and message, and ends
/// the run the way an uncaught one does.
pub fn run(
    args: &Args,
    input: &str,
    header: &SamHeader,
    records: &[BamRecord],
    mut accept: impl FnMut(&SamHeader, &BamRecord) -> Result<(), String>,
) {
    let assume_sorted = args.bool("ASSUME_SORTED", true);
    let stop_after = args.int("STOP_AFTER", 0);
    let with_reference = args.get("REFERENCE_SEQUENCE").is_some();
    check_coordinate_sorted(input, header, assume_sorted);

    let mut walker = ReferenceWalker::default();
    let mut count = 0i64;
    for record in records {
        if with_reference && record.reference_index != -1 {
            walker.get(record.reference_index);
        }
        if let Err(message) = accept(header, record) {
            thrown(&message);
        }
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
    }
}
