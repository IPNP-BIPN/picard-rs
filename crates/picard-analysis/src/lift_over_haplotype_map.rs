//! `LiftOverHaplotypeMap`: a haplotype database moved onto another reference through a UCSC chain
//! file.
//!
//! The chain itself is [`htsjdk_bam::liftover::LiftOver`]; the table is
//! [`crate::haplotype_map`]. What is here is the loop between them, which is the whole tool.
//!
//! Ported from `picard.fingerprint.LiftOverHaplotypeMap` in Picard 3.4.0.

use htsjdk_bam::interval::Interval;
use htsjdk_bam::liftover::{LiftOver, MissingToSequence};

use crate::haplotype_map::{rows, HaplotypeBlock, Row, Snp};

/// `LiftOverHaplotypeMap.LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS`, which is NOT the usual 1.
pub const LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS: i32 = 101;

/// What one run produced: the rows of the lifted table and the tool's exit code.
#[derive(Debug, Clone, PartialEq)]
pub struct LiftOverHaplotypeMapResult {
    pub rows: Vec<Row>,
    pub return_code: i32,
}

/// `doWork`: every block's SNPs lifted one at a time, then written out.
///
/// A SNP that does not lift is DROPPED and the run carries on, so the output is a shorter table
/// rather than none at all and a database whose every SNP fails leaves a file holding its header
/// and its column line alone. The exit code for any failure at all is 101.
///
/// The alleles are carried over unchanged. The chain may put the SNP on the negative strand and
/// the tool still writes the bases it read, never their complements, and the frequency stays the
/// minor allele's whatever the strand.
///
/// A block whose every SNP fails is still added to the map, as an empty one, which is why it
/// contributes no row rather than an empty one.
pub fn lift_over_haplotype_map(
    blocks: &[HaplotypeBlock],
    lift: &LiftOver,
    sequence_order: &[String],
) -> Result<LiftOverHaplotypeMapResult, MissingToSequence> {
    lift.validate_to_sequences(sequence_order)?;
    let mut lifted: Vec<HaplotypeBlock> = Vec::with_capacity(blocks.len());
    let mut any_failed = false;
    for block in blocks {
        let mut to = HaplotypeBlock::default();
        for snp in &block.snps {
            let interval = Interval::new(&snp.chromosome, snp.position, snp.position);
            match lift.lift_over(&interval) {
                Some(lifted_interval) => {
                    // The addition cannot fail: a block's SNPs land on one contig or the chain
                    // would have had to name two, which the reference does not check either.
                    let _ = to.add_snp(Snp {
                        name: snp.name.clone(),
                        chromosome: lifted_interval.contig.clone(),
                        position: lifted_interval.start,
                        major_allele: snp.major_allele,
                        minor_allele: snp.minor_allele,
                        minor_allele_frequency: snp.minor_allele_frequency,
                        panels: snp.panels.clone(),
                    });
                }
                None => any_failed = true,
            }
        }
        lifted.push(to);
    }
    Ok(LiftOverHaplotypeMapResult {
        rows: rows(&lifted, sequence_order),
        return_code: if any_failed {
            LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS
        } else {
            0
        },
    })
}
