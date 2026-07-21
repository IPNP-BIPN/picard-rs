//! Adapter-sequence detection.
//!
//! Ported from `picard.analysis.AdapterUtility` and the `picard.util.IlluminaUtil`
//! `IlluminaAdapterPair` constants, tag 3.4.0.
//!
//! The shape is unusual enough to be worth stating: the adapter *sequences* are never matched
//! whole. `prepareAdapterSequences` slices each one into every 16-mer, throws away any 16-mer
//! with more than one `N`, adds the reverse complement of each survivor, and matches a read's
//! **first 16 bases** against that set with at most one mismatch. So a read is called an adapter
//! read on the strength of 16 bases, and the 58- to 67-base sequences in the defaults are only
//! ever a source of 16-mers.
//!
//! Two things follow that a reimplementation is likely to get wrong:
//!
//!   * The kmer set is a `HashSet<String>` and the match array is filled by iterating it, so the
//!     array's order is Java's hash order. It does not reach the output, because `isAdapter` is
//!     an existential over the array and returns the same answer in any order - but that has to
//!     be *checked*, not assumed, and it is why this port can use a sorted set without a
//!     divergence.
//!   * The `N` filter counts `N` in the *kmer*, not in the read, and the threshold it compares
//!     against is `MAX_ADAPTER_ERRORS`, the same constant used for read mismatches. A kmer with
//!     two `N`s is dropped entirely rather than being allowed two free mismatches.

use std::collections::BTreeSet;

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::bases_equal;

/// `ADAPTER_MATCH_LENGTH`: the number of read bases checked.
const ADAPTER_MATCH_LENGTH: usize = 16;
/// `MAX_ADAPTER_ERRORS`: mismatches allowed, and also the maximum `N` count in a usable kmer.
const MAX_ADAPTER_ERRORS: usize = 1;

const READ_UNMAPPED: u16 = 0x4;
const READ_REVERSE: u16 = 0x10;

/// `IlluminaAdapterPair.SINGLE_END`, `PAIRED_END` and `INDEXED`, in the order
/// `AdapterUtility.DEFAULT_ADAPTER_SEQUENCE` lists them.
///
/// The 5' sequence is identical in all three pairs, so the list contains it three times. That
/// duplication is harmless here - the kmers go into a set - but it is reproduced because the
/// list is also what Picard echoes into its own log line.
pub const DEFAULT_ADAPTER_SEQUENCE: &[&str] = &[
    // SINGLE_END
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    // PAIRED_END
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    // INDEXED, whose 3' adapter carries eight Ns for the barcode
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
];

/// `SequenceUtil.complement`.
fn complement(base: u8) -> u8 {
    match base {
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

pub struct AdapterUtility {
    /// The kmers, held in a sorted set rather than a hash set. See the module note: the match is
    /// existential, so the order cannot reach the output, and a deterministic order makes a
    /// divergence traceable.
    kmers: Vec<Vec<u8>>,
}

impl AdapterUtility {
    /// `new AdapterUtility(adapterSequence)`.
    pub fn new(adapter_sequences: &[&str]) -> Self {
        let mut set: BTreeSet<Vec<u8>> = BTreeSet::new();
        for seq in adapter_sequences {
            let seq = seq.as_bytes();
            if seq.len() < ADAPTER_MATCH_LENGTH {
                continue;
            }
            for i in 0..=(seq.len() - ADAPTER_MATCH_LENGTH) {
                let kmer = seq[i..i + ADAPTER_MATCH_LENGTH].to_ascii_uppercase();
                // The N count is on the kmer and the threshold is MAX_ADAPTER_ERRORS, so a kmer
                // with two Ns is dropped rather than tolerated.
                let ns = kmer.iter().filter(|&&c| c == b'N').count();
                if ns <= MAX_ADAPTER_ERRORS {
                    set.insert(reverse_complement(&kmer));
                    set.insert(kmer);
                }
            }
        }
        AdapterUtility {
            kmers: set.into_iter().collect(),
        }
    }

    /// `AdapterUtility` with Picard's default adapter list.
    pub fn with_defaults() -> Self {
        Self::new(DEFAULT_ADAPTER_SEQUENCE)
    }

    /// `isAdapter(record)`.
    ///
    /// A read that mapped with a non-zero mapping quality is never an adapter read, whatever its
    /// bases say. So the same bases can be an adapter read in one BAM and not in another,
    /// depending only on how the aligner scored them.
    pub fn is_adapter(&self, rec: &BamRecord) -> bool {
        let unmapped = rec.flags & READ_UNMAPPED != 0;
        if !unmapped && rec.mapping_quality != 0 {
            return false;
        }
        // A mapped read on the reverse strand is stored in reference orientation, so it has to
        // be read back the other way to recover the order the machine saw.
        let rev_comp = !unmapped && (rec.flags & READ_REVERSE != 0);
        self.is_adapter_sequence(&rec.read_bases, rev_comp)
    }

    /// `isAdapterSequence(read, revCompRead)`.
    pub fn is_adapter_sequence(&self, read: &[u8], rev_comp_read: bool) -> bool {
        if read.len() < ADAPTER_MATCH_LENGTH {
            return false;
        }
        for adapter in &self.kmers {
            let mut errors = 0;
            for (i, &a) in adapter.iter().enumerate() {
                let base = if rev_comp_read {
                    complement(read[read.len() - i - 1])
                } else {
                    read[i]
                };
                if !bases_equal(base, a) {
                    errors += 1;
                    if errors > MAX_ADAPTER_ERRORS {
                        break;
                    }
                }
            }
            if errors <= MAX_ADAPTER_ERRORS {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn util() -> AdapterUtility {
        AdapterUtility::with_defaults()
    }

    #[test]
    fn the_first_sixteen_bases_of_a_default_adapter_match() {
        let read = b"AATGATACGGCGACCACCGAGATCT";
        assert!(util().is_adapter_sequence(read, false));
    }

    /// One mismatch is tolerated, two are not. Built against a single made-up adapter so the
    /// distances are exact; the default list is dense enough in 16-mers that a "two off" read
    /// can land within one error of some *other* kmer.
    #[test]
    fn one_mismatch_is_allowed_and_two_are_not() {
        let u = AdapterUtility::new(&["ACGTACGTACGTACGT"]);
        assert!(u.is_adapter_sequence(b"ACGTACGTACGTACGT", false), "exact");
        assert!(u.is_adapter_sequence(b"TCGTACGTACGTACGT", false), "one off");
        assert!(
            !u.is_adapter_sequence(b"TTGTACGTACGTACGT", false),
            "two off"
        );
    }

    #[test]
    fn a_read_shorter_than_the_match_length_is_never_an_adapter() {
        assert!(
            !util().is_adapter_sequence(b"AATGATACGGCGACC", false),
            "15 bases"
        );
    }

    /// Only kmers with at most one `N` survive, so the eight-N stretch of the indexed adapter
    /// contributes nothing, while the kmers that clip its edge with a single `N` do.
    #[test]
    fn kmers_with_more_than_one_n_are_dropped() {
        let u = AdapterUtility::new(&["AAAAAAAAAAAAAAANNAAAAAAAAAAAAAAA"]);
        assert!(
            !u.kmers
                .iter()
                .any(|k| k.iter().filter(|&&c| c == b'N').count() > 1),
            "no surviving kmer has two Ns"
        );
        assert!(
            u.kmers.iter().any(|k| k.contains(&b'N')),
            "the single-N kmers at the edges do survive"
        );
    }

    #[test]
    fn every_kmer_has_its_reverse_complement() {
        let u = AdapterUtility::new(&["ACGTACGTACGTACGTAC"]);
        for k in &u.kmers {
            assert!(
                u.kmers.contains(&reverse_complement(k)),
                "{:?} has no reverse complement",
                String::from_utf8_lossy(k)
            );
        }
    }

    /// A mapped read with a non-zero mapping quality is exempt, whatever its bases are.
    #[test]
    fn mapping_quality_overrides_the_bases() {
        let u = util();
        let mut rec = BamRecord {
            read_name: "r".to_string(),
            flags: 0,
            reference_index: 0,
            alignment_start: 1,
            mapping_quality: 60,
            cigar: Default::default(),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: b"AATGATACGGCGACCACCGAGATCT".to_vec(),
            base_qualities: vec![30; 25],
            tags: Default::default(),
        };
        assert!(!u.is_adapter(&rec), "mapped with MAPQ 60");
        rec.mapping_quality = 0;
        assert!(u.is_adapter(&rec), "the same bases at MAPQ 0");
        rec.mapping_quality = 60;
        rec.flags = READ_UNMAPPED;
        assert!(u.is_adapter(&rec), "unmapped, so MAPQ is not consulted");
    }

    /// A mapped reverse-strand read is stored in reference orientation and is read back the
    /// other way, so the adapter has to be at the *end* of the stored bases.
    #[test]
    fn a_reverse_strand_read_is_matched_from_the_other_end() {
        let u = util();
        let forward = b"AATGATACGGCGACCACCGAGATCT".to_vec();
        let stored = reverse_complement(&forward);
        let rec = BamRecord {
            read_name: "r".to_string(),
            flags: READ_REVERSE,
            reference_index: 0,
            alignment_start: 1,
            mapping_quality: 0,
            cigar: Default::default(),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: stored,
            base_qualities: vec![30; 25],
            tags: Default::default(),
        };
        assert!(u.is_adapter(&rec));
    }
}
