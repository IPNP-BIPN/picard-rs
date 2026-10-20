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

use std::io::{self, BufRead, Write};

use htsjdk_bam::interval::ParseError;
use htsjdk_bam::liftover::{ChainParseError, LiftOver, MissingToSequence};

use crate::compact_interval::{self as ci, Contigs, Iv};

/// The lifted interval list plus the tool's process return code.
#[derive(Debug, PartialEq, Eq)]
pub struct LiftOverResult {
    /// The output `.interval_list` text.
    pub output: String,
    /// `0` if all intervals lifted, `1` if any were rejected.
    pub return_code: i32,
    /// What `REJECT` receives: the input's own header, re-encoded the way every `IntervalList.write`
    /// encodes one, and the intervals that did not lift, in input order.
    pub rejects: String,
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

/// The lifted list and the rejects, ready to write.
pub struct Lifted {
    /// The output's header: a fresh `@HD` and the target dictionary's `@SQ` lines.
    header: String,
    contigs: Contigs,
    /// The lifted intervals, sorted by the target dictionary.
    intervals: Vec<Iv>,
    /// The input's own header, re-encoded the way every `IntervalList.write` encodes one.
    rejects_header: String,
    input_contigs: Contigs,
    /// The intervals that did not lift, in input order.
    rejected: Vec<Iv>,
}

impl Lifted {
    /// `0` if all intervals lifted, `1` if any were rejected.
    pub fn return_code(&self) -> i32 {
        i32::from(!self.rejected.is_empty())
    }

    pub fn write_output(&self, out: &mut impl Write) -> io::Result<()> {
        out.write_all(self.header.as_bytes())?;
        ci::write_lines(&self.contigs, &self.intervals, out)
    }

    pub fn write_rejects(&self, out: &mut impl Write) -> io::Result<()> {
        out.write_all(self.rejects_header.as_bytes())?;
        ci::write_lines(&self.input_contigs, &self.rejected, out)
    }
}

/// `LiftOverIntervalList.doWork` over the input interval list (read line by line), the target
/// sequence dictionary, and the chain file. A refusal is reported once the input has read
/// cleanly to its end, as a whole-file read of it came first.
pub fn lift(
    input: impl BufRead,
    sequence_dictionary: &str,
    chain: &str,
    min_liftover_pct: f64,
) -> io::Result<Result<Lifted, LiftOverIntervalListError>> {
    let mut lift = match LiftOver::load(chain) {
        Ok(lift) => lift,
        Err(error) => {
            ci::drain(input)?;
            return Ok(Err(LiftOverIntervalListError::Chain(error)));
        }
    };
    lift.set_lift_over_min_match(min_liftover_pct);

    // The input intervals are read by contig name; the input's own dictionary plays no part in the
    // lift.
    let mut input_contigs = Contigs::default();
    let mut intervals = Vec::new();
    let input_header = match ci::read_list(input, &mut input_contigs, &mut intervals)? {
        Ok(header) => header,
        Err(error) => return Ok(Err(LiftOverIntervalListError::Input(error))),
    };

    let (names, sq_lines) = target_dictionary(sequence_dictionary);
    if let Err(missing) = lift.validate_to_sequences(&names) {
        return Ok(Err(LiftOverIntervalListError::MissingToSequence(missing)));
    }

    // Each interval that lifts is replaced in place by its lifted self.
    let mut contigs = Contigs::with_dictionary(&names);
    let mut rejected = Vec::new();
    let mut lifted = intervals;
    lifted.retain_mut(
        |interval| match lift.lift_over(&interval.to_interval(&input_contigs)) {
            Some(to) => {
                *interval = Iv::from_interval(to, &mut contigs);
                true
            }
            None => {
                rejected.push(interval.clone());
                false
            }
        },
    );
    let input_header: String = input_header
        .leading_lines()
        .iter()
        .map(|line| format!("{line}\n"))
        .collect();
    let rejects_header =
        htsjdk_bam::reader::parse_header_text(&input_header).encode_replacing_version();

    // Sort the lifted intervals by the target dictionary.
    ci::sort(&contigs, &mut lifted);

    let mut header = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for sq in &sq_lines {
        header.push_str(sq);
        header.push('\n');
    }
    Ok(Ok(Lifted {
        header,
        contigs,
        intervals: lifted,
        rejects_header,
        input_contigs,
        rejected,
    }))
}

/// [`lift`] over the input interval list's text, rendering both files.
pub fn lift_over_interval_list(
    input_interval_list: &str,
    sequence_dictionary: &str,
    chain: &str,
    min_liftover_pct: f64,
) -> Result<LiftOverResult, LiftOverIntervalListError> {
    let lifted = lift(
        input_interval_list.as_bytes(),
        sequence_dictionary,
        chain,
        min_liftover_pct,
    )
    .expect("a text reads")?;
    let render = |write: &dyn Fn(&mut Vec<u8>) -> io::Result<()>| {
        let mut out = Vec::new();
        write(&mut out).expect("a Vec writes");
        String::from_utf8(out).expect("the inputs were text")
    };
    Ok(LiftOverResult {
        output: render(&|out| lifted.write_output(out)),
        return_code: lifted.return_code(),
        rejects: render(&|out| lifted.write_rejects(out)),
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
