//! `SamToFastq` as a runnable binary: the covering array's port side.
//!
//! The whole tool is `picard_analysis::sam_to_fastq_full`: which reads survive
//! (`INCLUDE_NON_PRIMARY_ALIGNMENTS`, `INCLUDE_NON_PF_READS`), where each goes (one set of files, or
//! one per read group named by `RG_TAG`, interleaved or split, unpaired reads apart), and what is
//! written for it (clipped at `CLIPPING_ATTRIBUTE`, put back in sequencing orientation by
//! `RE_REVERSE`, trimmed, quality-trimmed and capped). The array compares every file a run leaves.

fn main() {
    picard_analysis::output::install();
    picard_analysis::sam_to_fastq_full::main_with(false)
}
