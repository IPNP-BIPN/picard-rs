//! `NonNFastaSize`.
//!
//! Ports `picard.reference.NonNFastaSize.doWork` at tag 3.4.0 for the default (whole-genome) case:
//! count the bases of a FASTA that are not `N`, and write that count followed by a newline.
//!
//! `doWork` uppercases each contig's bases (`StringUtil.toUpperCase`) and counts every base that is
//! not `N` (`SequenceUtil.N`), summed over every contig, with a `WholeGenomeReferenceSequenceMask`
//! that admits every position. So the count is simply the number of non-`N` bases across the file,
//! case-insensitively (`n` counts as `N`). The port reads the FASTA with [`htsjdk_bam::fasta`], which
//! preserves case, so it uppercases as it counts.
//!
//! Scope of this slice: no `INTERVALS` (the default whole-genome mask). The interval-restricted count
//! is a separate surface.

use htsjdk_bam::fasta::{read_fasta, FastaError};

/// `NonNFastaSize.doWork` over FASTA text: the non-`N` base count, followed by a newline.
pub fn non_n_fasta_size(fasta: &str) -> Result<String, FastaError> {
    let sequences = read_fasta(fasta.as_bytes())?;
    let mut non_n: u64 = 0;
    for sequence in &sequences {
        for &base in &sequence.bases {
            if !base.eq_ignore_ascii_case(&b'N') {
                non_n += 1;
            }
        }
    }
    Ok(format!("{non_n}\n"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_non_n_bases_case_insensitively() {
        // 20 + 8 bases across chr1 with four N; 25 across chr2 with none.
        let fasta = ">chr1\nACGTacgtACGTNNNNacgt\nACGTACGT\n>chr2\nTTTTTTTTTTGGGGGGGGGGCCCCC\n";
        // chr1: 28 bases, 4 are N -> 24; chr2: 25 -> total 49.
        assert_eq!(non_n_fasta_size(fasta).unwrap(), "49\n");
    }

    #[test]
    fn lowercase_n_also_counts_as_n() {
        assert_eq!(non_n_fasta_size(">c\nAnCnG\n").unwrap(), "3\n");
    }
}
