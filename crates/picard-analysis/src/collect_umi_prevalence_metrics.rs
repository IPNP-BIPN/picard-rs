//! `CollectUmiPrevalenceMetrics`: how many distinct UMIs each duplicate set holds.
//!
//! Grouping the reads into duplicate sets is not ported; the six filters that decide which reads
//! reach a set are, along with the counting.
//!
//! Ported from `picard.analysis.CollectUmiPrevalenceMetrics` in Picard 3.4.0.

use std::collections::BTreeMap;
use std::collections::HashSet;

/// `CollectUmiPrevalenceMetrics.MINIMUM_MQ`.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 30;
/// `CollectUmiPrevalenceMetrics.MINIMUM_BARCODE_BQ`.
pub const DEFAULT_MINIMUM_BARCODE_BASE_QUALITY: i32 = 30;
/// `CollectUmiPrevalenceMetrics.BARCODE_TAG`.
pub const DEFAULT_BARCODE_TAG: &str = "RX";
/// `CollectUmiPrevalenceMetrics.BARCODE_BQ`.
pub const DEFAULT_BARCODE_QUALITY_TAG: &str = "BQ";

/// One read, reduced to what the six filters read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub unmapped: bool,
    pub mapping_quality: i32,
    pub secondary_or_supplementary: bool,
    pub paired: bool,
    pub barcode: Option<String>,
    /// The barcode's per-base qualities, already decoded from the tag's FASTQ characters.
    pub barcode_qualities: Option<Vec<i32>>,
}

/// The arguments the filters consult.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arguments {
    pub minimum_mapping_quality: i32,
    pub minimum_barcode_base_quality: i32,
    pub filter_unpaired_reads: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            minimum_mapping_quality: DEFAULT_MINIMUM_MAPPING_QUALITY,
            minimum_barcode_base_quality: DEFAULT_MINIMUM_BARCODE_BASE_QUALITY,
            filter_unpaired_reads: true,
        }
    }
}

/// `BarcodeQualityFilter.filterOut`, transcribed with its predicate INVERTED as it stands.
///
/// It computes `badQuality`, whether any base of the barcode is under the floor, and then returns
/// its NEGATION. So a barcode entirely above the floor is filtered OUT and one carrying a bad base
/// is kept, which is the opposite of what the name and the argument's documentation say.
///
/// Two consequences follow, and both are visible in the golden. A file of well-formed barcodes
/// reports nothing at all. And LOWERING the floor drops more reads rather than fewer, because a
/// lower floor is one fewer base under it.
///
/// A read with no quality tag returns false before any of that, so an absent tag is the only other
/// way past the filter.
pub fn barcode_quality_filters_out(read: &Read, minimum: i32) -> bool {
    let Some(qualities) = &read.barcode_qualities else {
        return false;
    };
    let bad_quality = qualities.iter().any(|quality| *quality < minimum);
    !bad_quality
}

/// `UMITagPresentFilter.filterOut`: a read with no barcode tag is dropped.
pub fn barcode_tag_filters_out(read: &Read) -> bool {
    read.barcode.is_none()
}

/// The whole `AggregateFilter`, in the order the tool assembles it.
///
/// The sixth is added only when `--FILTER_UNPAIRED_READS` is set, which it is by default.
pub fn filters_out(read: &Read, arguments: &Arguments) -> bool {
    if read.unmapped {
        return true;
    }
    if read.mapping_quality < arguments.minimum_mapping_quality {
        return true;
    }
    if read.secondary_or_supplementary {
        return true;
    }
    if barcode_tag_filters_out(read) {
        return true;
    }
    if barcode_quality_filters_out(read, arguments.minimum_barcode_base_quality) {
        return true;
    }
    arguments.filter_unpaired_reads && !read.paired
}

/// `SAMUtils.fastqToPhred`: a FASTQ quality character is its code point less thirty-three.
///
/// The tool strips spaces from the tag first, which is what lets a `BQ` written with separators
/// be read at all.
pub fn decode_barcode_qualities(tag: &str) -> Vec<i32> {
    tag.chars()
        .filter(|c| *c != ' ')
        .map(|c| c as i32 - 33)
        .collect()
}

/// The histogram: how many duplicate sets hold each number of distinct UMIs.
///
/// The key is the UMI COUNT and the value is the number of sets, so three sets of one UMI each is
/// `1 -> 3` and not `3 -> 1`. The UMIs of a set are a SET, so two reads carrying the same tag are
/// one UMI. A set every one of whose reads was filtered contributes nothing, not a zero.
pub fn histogram(sets: &[Vec<Read>], arguments: &Arguments) -> BTreeMap<usize, i64> {
    let mut counts: BTreeMap<usize, i64> = BTreeMap::new();
    for set in sets {
        let kept: Vec<&Read> = set
            .iter()
            .filter(|read| !filters_out(read, arguments))
            .collect();
        if kept.is_empty() {
            continue;
        }
        let barcodes: HashSet<&String> = kept
            .iter()
            .filter_map(|read| read.barcode.as_ref())
            .collect();
        *counts.entry(barcodes.len()).or_default() += 1;
    }
    counts
}
