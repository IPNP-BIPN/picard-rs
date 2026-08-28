//! `SortGff`: a GFF3 file sorted by contig and then by start.
//!
//! Reading and writing GFF3 is not ported. What is ported is the order the comparator produces
//! and what a sequence dictionary changes about it.
//!
//! Ported from `picard.annotation.SortGff` in Picard 3.4.0.

/// One feature, reduced to what the comparator reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    /// The order it was read in, which is what a stable sort falls back on.
    pub index: usize,
}

/// `FeatureComparator`, which sorts by contig and then by start.
///
/// Without a dictionary the contigs are compared as STRINGS, so `chr10` sorts before `chr2`. With
/// one they are compared by their index in it, and a contig the dictionary does not name has
/// index -1 and therefore sorts BEFORE every contig it does name.
///
/// The comparison SUBTRACTS two ints in both halves, so it is a difference and not a sign: a
/// start past two billion would overflow it, which no real GFF reaches.
pub fn compare(a: &Feature, b: &Feature, dictionary: Option<&[String]>) -> std::cmp::Ordering {
    let first = match dictionary {
        None => a.contig.cmp(&b.contig),
        Some(dictionary) => {
            sequence_index(dictionary, &a.contig).cmp(&sequence_index(dictionary, &b.contig))
        }
    };
    if first != std::cmp::Ordering::Equal {
        return first;
    }
    a.start.cmp(&b.start)
}

/// `SAMSequenceDictionary.getSequenceIndex`, which answers -1 for a name it does not hold.
pub fn sequence_index(dictionary: &[String], contig: &str) -> i32 {
    dictionary
        .iter()
        .position(|name| name == contig)
        .map_or(-1, |index| index as i32)
}

/// The features in order.
///
/// The sort is STABLE, so two features that start together keep the order they were read in and
/// their ends are never compared.
pub fn sort(features: &[Feature], dictionary: Option<&[String]>) -> Vec<Feature> {
    let mut sorted = features.to_vec();
    sorted.sort_by(|a, b| compare(a, b, dictionary));
    sorted
}

/// The version directive the codec writes, which is ITS OWN and not the input's.
///
/// A file that opens with `##gff-version 3.1.26` comes back opening with this.
pub const GFF_VERSION_DIRECTIVE: &str = "##gff-version 3.1.25";

/// The refusal a file the codec cannot decode produces.
///
/// A file with NO FEATURE gets the same one as a file that is not GFF at all: `canDecode` wants a
/// feature and not only a directive, so an empty file is refused rather than sorted into an empty
/// file.
pub fn cannot_decode_message(path: &str) -> String {
    format!("Input file {path} cannot be read by Gff3Codec")
}

/// `nRecordsInMemory`, which decides only where the sort holds its records.
pub const DEFAULT_RECORDS_IN_MEMORY: usize = 50000;
