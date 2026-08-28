//! `CollectIndependentReplicateMetrics`: which duplicate sets are not duplicates at all.
//!
//! A duplicate set whose reads disagree on the allele at a heterozygous site cannot have come from
//! one molecule, so it is an independent replicate rather than a duplicate. Everything the tool
//! does is deciding which sets get to answer that question.
//!
//! The filters sit at four levels and only the first is a property of the VCF:
//!
//!  * the SITE, by its genotype and its `--MINIMUM_GQ`;
//!  * the READ, by its `--MINIMUM_MQ` and, unless `--FILTER_UNPAIRED_READS` is turned off, by
//!    being paired at all;
//!  * the BASE, by its `--MINIMUM_BQ`, which leaves the read in the set and takes its allele away;
//!  * and the BARCODE, by every one of its own qualities against `--MINIMUM_BARCODE_BQ`, which
//!    decides whether the UMI counters are touched and nothing else.
//!
//! Ported from `picard.analysis.replicates.CollectIndependentReplicateMetrics`.

/// `MINIMUM_GQ`, on the site's genotype quality.
pub const DEFAULT_MINIMUM_GQ: i32 = 90;
/// `MINIMUM_MQ`, on the read.
pub const DEFAULT_MINIMUM_MQ: i32 = 40;
/// `MINIMUM_BQ`, on the base at the site.
pub const DEFAULT_MINIMUM_BQ: i32 = 17;
/// `MINIMUM_BARCODE_BQ`, on every base of the barcode.
pub const DEFAULT_MINIMUM_BARCODE_BQ: i32 = 30;
/// `BARCODE_TAG` and `BARCODE_BQ`, the two tags a UMI is read from.
pub const DEFAULT_BARCODE_TAG: &str = "RX";
pub const DEFAULT_BARCODE_QUALITY_TAG: &str = "QX";

/// The two set sizes the tool has counters for.
pub const DOUBLETON_SIZE: usize = 2;
pub const TRIPLETON_SIZE: usize = 3;

/// What a set of a given size is counted as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SetSize {
    /// One read, counted in `nTotalReads` and nowhere else.
    Singleton,
    Doubleton,
    Tripleton,
    /// More than three, whose reads are counted in `nReadsInBigSets`.
    Big,
}

/// `setSize` as the counters see it.
pub fn set_size(reads: usize) -> SetSize {
    match reads {
        0 | 1 => SetSize::Singleton,
        DOUBLETON_SIZE => SetSize::Doubleton,
        TRIPLETON_SIZE => SetSize::Tripleton,
        _ => SetSize::Big,
    }
}

/// `SetClassification`, from the three counts a set produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SetClassification {
    /// A read carried a base that is neither of the site's two alleles.
    MismatchingAllele,
    /// Both alleles are in the set, which is what an independent replicate looks like.
    DifferentAlleles,
    ReferenceAllele,
    AlternateAllele,
}

/// `classifySet`, in the reference's own order.
///
/// The order matters: a set that carries a third allele is MISMATCHING even when it also carries
/// both of the site's own, so the third-allele test comes first and the heterogeneity test second.
/// A set whose reads all lost their base to the quality floor counts as a reference set, both
/// counts being zero.
pub fn classify_set(reference: usize, alternate: usize, other: usize) -> SetClassification {
    if other != 0 {
        return SetClassification::MismatchingAllele;
    }
    if alternate > 0 && reference > 0 {
        return SetClassification::DifferentAlleles;
    }
    if reference == 0 {
        return SetClassification::AlternateAllele;
    }
    SetClassification::ReferenceAllele
}

/// `calculateEditDistance`: a Hamming distance over two barcodes of equal length.
///
/// The reference validates the lengths rather than padding, so two barcodes of different lengths
/// are an argument error and not a distance.
pub fn edit_distance(left: &str, right: &str) -> Option<u8> {
    if left.len() != right.len() {
        return None;
    }
    Some(
        left.bytes()
            .zip(right.bytes())
            .filter(|(a, b)| a != b)
            .count() as u8,
    )
}

/// Whether a doubleton's barcodes are used at all.
///
/// Every base of every barcode has to reach the floor. A read with no quality tag contributes an
/// EMPTY string, which has no base under the floor, so a set whose reads carry no barcode at all
/// is a set with good barcodes.
pub fn barcodes_are_usable(qualities: &[&str], minimum_barcode_bq: i32) -> bool {
    !qualities.iter().any(|written| {
        written
            .bytes()
            .any(|quality| i32::from(quality as i8 - 33) < minimum_barcode_bq)
    })
}

/// Whether a site is looked at at all: heterozygous, and over the genotype-quality floor.
pub fn site_is_used(is_heterozygous: bool, genotype_quality: i32, minimum_gq: i32) -> bool {
    is_heterozygous && genotype_quality >= minimum_gq
}

/// Whether a read reaches the sets, which is decided before any set is built.
pub fn read_is_used(
    mapping_quality: i32,
    paired: bool,
    minimum_mq: i32,
    filter_unpaired_reads: bool,
) -> bool {
    mapping_quality >= minimum_mq && (paired || !filter_unpaired_reads)
}

/// Whether the base at the site is counted, which leaves the read in its set either way.
///
/// The test is strictly greater: `read.getBaseQualities()[offset] <= MINIMUM_BQ` skips the base, so
/// a base exactly at the floor is skipped and not kept.
pub fn base_is_counted(base_quality: i32, minimum_bq: i32) -> bool {
    base_quality > minimum_bq
}

/// The `nThreeAllelesSites` the tail of `doWork` adds, which is one for a run that examined
/// nothing at all.
///
/// The counter is incremented whenever the loop ended with no locus pending a merge. A run whose
/// site list is empty is in that state from the start, so a homozygous site, a site under the
/// quality floor and a file of unpaired reads each report one three-allele site having looked at
/// nothing. This is the reference's own arithmetic and not a rounding of it.
pub fn three_allele_sites_from_the_tail(a_locus_is_pending: bool) -> i32 {
    if a_locus_is_pending {
        0
    } else {
        1
    }
}

/// The extension the output carries, which is none: `OUTPUT` is a file name and not a basename.
pub const OUTPUT_IS_A_BASENAME: bool = false;
