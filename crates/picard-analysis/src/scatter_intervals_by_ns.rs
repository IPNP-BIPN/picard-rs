//! `ScatterIntervalsByNs`.
//!
//! Ports `picard.util.ScatterIntervalsByNs.doWork`/`segregateReference` at tag 3.4.0: walk a reference
//! into alternating runs of no-call (`N`) and called (`ACGT`) bases, labelling each run `Nmer` or
//! `ACGTmer`, then emit them as a coordinate-sorted interval list restricted to the requested
//! `OUTPUT_TYPE` (`N`, `ACGT`, or `BOTH`).
//!
//! Each contig is scanned base by base with `SequenceUtil.isNoCall` (a base is a no-call if, once
//! upper-cased, it is `N` or `.`). A run `[start, i)` (0-based) becomes the interval `start + 1 ..= i`
//! (1-based inclusive), named for the run that just closed. Short `N` runs are then folded away: while
//! the front three intervals are `ACGTmer`, `Nmer`, `ACGTmer`, all mutually abutting, and the middle
//! `Nmer` is at most `MAX_TO_MERGE` bases, they are replaced by one `ACGTmer` spanning all three (and
//! that merged interval is reconsidered against what follows). Everything else is emitted in order, so
//! the result is already in coordinate order.
//!
//! The output header is a fresh `@HD VN:1.6 SO:coordinate` followed by the reference dictionary's
//! `@SQ` lines verbatim (the same idempotent `.dict` round trip [`crate::bed_to_interval_list`] relies
//! on), so `M5`/`UR` survive unchanged. Contigs are taken from the FASTA in file order, which is the
//! dictionary order for a `.dict` built from that FASTA.

use std::collections::VecDeque;

use htsjdk_bam::fasta::{read_fasta, FastaError};
use htsjdk_bam::interval::Interval;

const NMER: &str = "Nmer";
const ACGTMER: &str = "ACGTmer";

/// `OUTPUT_TYPE`: which run labels reach the output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputType {
    /// `N`: only `Nmer` runs.
    N,
    /// `ACGT`: only `ACGTmer` runs.
    Acgt,
    /// `BOTH`: every run (the default).
    Both,
}

impl OutputType {
    /// `OutputType.accepts`.
    fn accepts(&self, name: &str) -> bool {
        match self {
            OutputType::N => name == NMER,
            OutputType::Acgt => name == ACGTMER,
            OutputType::Both => name == NMER || name == ACGTMER,
        }
    }
}

/// `ScatterIntervalsByNs`'s options.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `OUTPUT_TYPE`.
    pub output_type: OutputType,
    /// `MAX_TO_MERGE`: the longest `N` run that keeps its flanking `ACGT` runs joined.
    pub max_to_merge: i32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            output_type: OutputType::Both,
            max_to_merge: 1,
        }
    }
}

/// `SequenceUtil.isNoCall`: a base is a no-call if, upper-cased, it is `N` or `.`.
fn is_no_call(base: u8) -> bool {
    let upper = base.to_ascii_uppercase();
    upper == b'N' || upper == b'.'
}

/// `Interval.abuts`: same contig, one ending exactly where the other begins.
fn abuts(a: &Interval, b: &Interval) -> bool {
    a.contig == b.contig && (a.end + 1 == b.start || b.end + 1 == a.start)
}

/// The reference dictionary's `@SQ` lines, verbatim and in order.
fn sq_lines(sequence_dictionary: &str) -> Vec<String> {
    sequence_dictionary
        .lines()
        .filter(|line| line.starts_with("@SQ"))
        .map(str::to_string)
        .collect()
}

/// `ScatterIntervalsByNs.doWork` over a reference FASTA and its dictionary, returning the interval
/// list text.
pub fn scatter_intervals_by_ns(
    reference_fasta: &str,
    sequence_dictionary: &str,
    opts: &Options,
) -> Result<String, FastaError> {
    let contigs = read_fasta(reference_fasta.as_bytes())?;

    // Build the alternating N / ACGT runs across every contig, in dictionary (file) order.
    let mut preliminary: Vec<Interval> = Vec::new();
    for contig in &contigs {
        let bases = &contig.bases;
        if bases.is_empty() {
            continue;
        }
        let mut n_open = is_no_call(bases[0]);
        let mut start = 0usize;
        for (i, &base) in bases.iter().enumerate() {
            let current_is_n = is_no_call(base);
            if n_open != current_is_n {
                preliminary.push(run(&contig.name, start + 1, i, n_open));
                start = i;
                n_open = !n_open;
            }
        }
        preliminary.push(run(&contig.name, start + 1, bases.len(), n_open));
    }

    // Fold `ACGT, short-N, ACGT` trios back into one `ACGT`, reconsidering the merged run each time.
    let mut queue: VecDeque<Interval> = preliminary.into();
    let mut segregated: Vec<Interval> = Vec::new();
    while !queue.is_empty() {
        if queue.len() >= 3
            && queue[0].name.as_deref() == Some(ACGTMER)
            && queue[1].name.as_deref() == Some(NMER)
            && queue[2].name.as_deref() == Some(ACGTMER)
            && abuts(&queue[0], &queue[1])
            && abuts(&queue[1], &queue[2])
            && queue[1].length() <= opts.max_to_merge
        {
            let contig = queue[0].contig.clone();
            let start = queue[0].start;
            let end = queue[2].end;
            queue.pop_front();
            queue.pop_front();
            queue.pop_front();
            queue.push_front(Interval::with_strand_and_name(
                &contig,
                start,
                end,
                false,
                Some(ACGTMER),
            ));
        } else {
            segregated.push(queue.pop_front().unwrap());
        }
    }

    let mut out = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for sq in sq_lines(sequence_dictionary) {
        out.push_str(&sq);
        out.push('\n');
    }
    for interval in &segregated {
        if opts
            .output_type
            .accepts(interval.name.as_deref().unwrap_or(""))
        {
            out.push_str(&interval.to_file_line());
            out.push('\n');
        }
    }
    Ok(out)
}

/// One run `[start, end]` (1-based inclusive) labelled for the block that just closed.
fn run(contig: &str, start: usize, end: usize, n_open: bool) -> Interval {
    Interval::with_strand_and_name(
        contig,
        start as i32,
        end as i32,
        false,
        Some(if n_open { NMER } else { ACGTMER }),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const DICT: &str = "@HD\tVN:1.6\n\
        @SQ\tSN:chr1\tLN:20\tM5:aaa\tUR:file:///ref.fasta\n";

    fn header() -> String {
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\tM5:aaa\tUR:file:///ref.fasta\n"
            .to_string()
    }

    #[test]
    fn alternating_blocks_become_labelled_intervals() {
        // 4 ACGT, 4 N, 4 ACGT: with MAX_TO_MERGE=1 the 4-base N run is too long to fold.
        let fasta = ">chr1\nACGTNNNNACGT\n";
        let opts = Options {
            max_to_merge: 1,
            ..Options::default()
        };
        let out = scatter_intervals_by_ns(fasta, DICT, &opts).unwrap();
        assert_eq!(
            out,
            format!(
                "{}chr1\t1\t4\t+\tACGTmer\nchr1\t5\t8\t+\tNmer\nchr1\t9\t12\t+\tACGTmer\n",
                header()
            )
        );
    }

    #[test]
    fn a_short_n_run_folds_its_flanking_acgt_runs_together() {
        // 4 ACGT, 1 N, 4 ACGT: the single N is <= MAX_TO_MERGE, so all nine bases become one ACGT.
        let fasta = ">chr1\nACGTNACGT\n";
        let out = scatter_intervals_by_ns(fasta, DICT, &Options::default()).unwrap();
        assert_eq!(out, format!("{}chr1\t1\t9\t+\tACGTmer\n", header()));
    }

    #[test]
    fn output_type_n_keeps_only_n_runs() {
        let fasta = ">chr1\nACGTNNNNACGT\n";
        let opts = Options {
            output_type: OutputType::N,
            max_to_merge: 1,
        };
        let out = scatter_intervals_by_ns(fasta, DICT, &opts).unwrap();
        assert_eq!(out, format!("{}chr1\t5\t8\t+\tNmer\n", header()));
    }

    #[test]
    fn a_dot_is_a_no_call_like_n() {
        let fasta = ">chr1\nACGT....ACGT\n";
        let opts = Options {
            output_type: OutputType::N,
            max_to_merge: 1,
        };
        let out = scatter_intervals_by_ns(fasta, DICT, &opts).unwrap();
        assert_eq!(out, format!("{}chr1\t5\t8\t+\tNmer\n", header()));
    }
}
