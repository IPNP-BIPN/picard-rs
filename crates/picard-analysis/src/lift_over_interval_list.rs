//! `LiftOverIntervalList`.
//!
//! Ports `picard.util.LiftOverIntervalList.doWork` at tag 3.4.0: adjust the coordinates of an
//! interval list from one reference build to another using a UCSC chain file. Each input interval
//! is lifted with [`htsjdk_bam::liftover::LiftOver`]; those that lift are collected, sorted by the
//! **target** dictionary, and written, and those that do not are dropped (and would go to the
//! optional `REJECT` file, which this port does not emit because it never reaches `OUTPUT`).
//!
//! The output header is a fresh `@HD VN:1.6 SO:coordinate` followed by the `@SQ` lines of the
//! `SEQUENCE_DICTIONARY` verbatim: `toIntervals.sorted().write` re-serializes the same
//! `SAMSequenceDictionary` that the `.dict` was parsed from, an idempotent round trip (the same one
//! [`crate::bed_to_interval_list`] relies on), so `M5`/`UR` and their order survive unchanged. The
//! input list's own header is irrelevant to the output: its intervals are read by contig name and
//! the "from" dictionary is never consulted by the lift.
//!
//! The return code mirrors the tool: `0` if every interval lifted, `1` if any was rejected
//! (`rejects.getIntervals().isEmpty() ? 0 : 1`).

use htsjdk_bam::interval::{Interval, IntervalList, ParseError};
use htsjdk_bam::liftover::{ChainParseError, LiftOver, MissingToSequence};

/// The lifted interval list plus the tool's process return code.
#[derive(Debug, PartialEq, Eq)]
pub struct LiftOverResult {
    /// The output `.interval_list` text.
    pub output: String,
    /// `0` if all intervals lifted, `1` if any were rejected.
    pub return_code: i32,
}

/// Why the liftover could not run.
#[derive(Debug, PartialEq, Eq)]
pub enum LiftOverIntervalListError {
    /// The chain file could not be parsed.
    Chain(ChainParseError),
    /// The input interval list could not be parsed.
    Input(ParseError),
    /// A chain names a "to" sequence absent from `SEQUENCE_DICTIONARY` (`validateToSequences`).
    MissingToSequence(MissingToSequence),
}

/// The `@SQ` `SN:` names in header order (for the coordinate sort) and the verbatim `@SQ` lines
/// (for the output header).
fn target_dictionary(sequence_dictionary: &str) -> (Vec<String>, Vec<String>) {
    let mut names = Vec::new();
    let mut sq_lines = Vec::new();
    for line in sequence_dictionary.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }
        for field in line.split('\t') {
            if let Some(v) = field.strip_prefix("SN:") {
                names.push(v.to_string());
            }
        }
        sq_lines.push(line.to_string());
    }
    (names, sq_lines)
}

/// `LiftOverIntervalList.doWork` over the input interval list, target sequence dictionary, and
/// chain file.
pub fn lift_over_interval_list(
    input_interval_list: &str,
    sequence_dictionary: &str,
    chain: &str,
    min_liftover_pct: f64,
) -> Result<LiftOverResult, LiftOverIntervalListError> {
    let mut lift = LiftOver::load(chain).map_err(LiftOverIntervalListError::Chain)?;
    lift.set_lift_over_min_match(min_liftover_pct);

    // The input intervals are read by contig name; the input's own dictionary plays no part in the
    // lift, so an empty dictionary suffices to parse them.
    let input = IntervalList::parse_body(Vec::new(), input_interval_list)
        .map_err(LiftOverIntervalListError::Input)?;

    let (names, sq_lines) = target_dictionary(sequence_dictionary);
    lift.validate_to_sequences(&names)
        .map_err(LiftOverIntervalListError::MissingToSequence)?;

    let mut lifted: Vec<Interval> = Vec::new();
    let mut rejected = false;
    for interval in &input.intervals {
        match lift.lift_over(interval) {
            Some(to) => lifted.push(to),
            None => rejected = true,
        }
    }

    // Sort the lifted intervals by the target dictionary and render the file.
    let to_list = IntervalList {
        dictionary: names,
        intervals: lifted,
    }
    .sorted();

    let mut output = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for sq in &sq_lines {
        output.push_str(sq);
        output.push('\n');
    }
    output.push_str(&to_list.write_body());

    Ok(LiftOverResult {
        output,
        return_code: if rejected { 1 } else { 0 },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SD: &str = "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:200\n";

    #[test]
    fn a_same_strand_interval_lifts_by_the_block_offset() {
        // chr1[0,100) -> chr2[50,150) '+', so chr1:11-20 -> chr2:61-70.
        let input = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t11\t20\t+\tA\n";
        let chain = "chain 1000 chr1 100 + 0 100 chr2 200 + 50 150 1\n100\n\n";
        let r = lift_over_interval_list(input, SD, chain, 0.95).unwrap();
        assert_eq!(
            r.output,
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:200\nchr2\t61\t70\t+\tA\n"
        );
        assert_eq!(r.return_code, 0);
    }

    #[test]
    fn a_negative_chain_flips_the_strand_and_coordinates() {
        // chr1[0,100) -> chr2[0,100) '-', size 200: chr1:11-20 -> chr2:181-190 on '-'.
        let input = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t11\t20\t+\tA\n";
        let chain = "chain 1000 chr1 100 + 0 100 chr2 200 - 0 100 2\n100\n\n";
        let r = lift_over_interval_list(input, SD, chain, 0.95).unwrap();
        assert_eq!(
            r.output,
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:200\nchr2\t181\t190\t-\tA\n"
        );
        assert_eq!(r.return_code, 0);
    }

    #[test]
    fn an_interval_outside_the_chain_is_rejected_and_sets_rc_1() {
        // chr1[0,50) -> chr2[0,50): chr1:11-20 lifts, chr1:71-80 does not.
        let input = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t11\t20\t+\tA\nchr1\t71\t80\t+\tB\n";
        let chain = "chain 1000 chr1 100 + 0 50 chr2 200 + 0 50 3\n50\n\n";
        let r = lift_over_interval_list(input, SD, chain, 0.95).unwrap();
        assert_eq!(
            r.output,
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:200\nchr2\t11\t20\t+\tA\n"
        );
        assert_eq!(r.return_code, 1);
    }
}
