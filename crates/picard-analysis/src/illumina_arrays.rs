//! Illumina's genotyping-array files, and the two smallest tools that read them.
//!
//! Three formats, none of them self-describing, and one of them not even a stream:
//!
//!  * a `.bpm` (bead pool manifest) opens with `BPM`, a version byte of one, an int version between
//!    three and five, two strings, the locus count, an index block the parser SKIPS, the names, one
//!    normalization id apiece, and then a locus entry each;
//!  * an `.egt` (cluster file) is a header and then, per locus, three counts and four triples of
//!    floats and fifteen floats nobody reads;
//!  * a `.gtc` (genotype calls) is a TABLE OF CONTENTS: an id and an ABSOLUTE offset per kind of
//!    data, and the reader seeks. The number of SNPs is not a payload at all, its OFFSET is the
//!    value.
//!
//! Ported from `picard.arrays.illumina.IlluminaBPMFile`, `InfiniumEGTFile`, `InfiniumGTCFile`,
//! `InfiniumDataFile`, `BpmToNormalizationManifestCsv` and `CompareGtcFiles` in Picard 3.4.0.

/// `InfiniumDataFile.parseString`: a varint length and then that many bytes.
///
/// The length is seven bits a byte, the high bit saying another follows, so a string of under a
/// hundred and twenty-eight bytes costs one byte of length and a longer one costs two.
pub fn parse_string(bytes: &[u8], at: usize) -> Option<(String, usize)> {
    let mut length = 0usize;
    let mut shift = 0;
    let mut cursor = at;
    loop {
        let byte = *bytes.get(cursor)?;
        cursor += 1;
        length += ((byte & 0x7F) as usize) << shift;
        if byte & 0x80 == 0 {
            break;
        }
        shift += 7;
    }
    let text = String::from_utf8(bytes.get(cursor..cursor + length)?.to_vec()).ok()?;
    Some((text, cursor + length))
}

/// One locus of a bead pool manifest, reduced to what the tools read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub name: String,
    pub index: usize,
    pub chrom: String,
    pub position: i32,
    pub address_a: i32,
    pub address_b: i32,
    pub assay_type: u8,
    /// The id as the file carries it, before the assay type is folded in.
    pub raw_normalization_id: u8,
}

impl Locus {
    /// `locusEntry.normalizationId = normId + 100 * assayType`.
    ///
    /// Which is why the CSV reports 101 for a locus whose file says 1: the number in the output is
    /// not the number in the file.
    pub fn normalization_id(&self) -> i32 {
        i32::from(self.raw_normalization_id) + 100 * i32::from(self.assay_type)
    }
}

/// What a manifest can be wrong about, in the reference's own words.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ManifestRefusal {
    /// `Invalid normalization ID: <id> for name: <name>`.
    NormalizationId { id: u8, name: String },
    /// `Invalid assay_type '<type>' for address B '<address>' in BPM file`.
    AssayType { assay_type: u8, address_b: i32 },
}

/// The two cross-checks the parser makes of every locus, in its own order.
///
/// The normalization id is checked BEFORE the assay type is folded in, which is what makes an id
/// of a hundred and one an error rather than a locus of assay type one; and an assay type of zero
/// must have no B address while any other must have one.
pub fn validate(locus: &Locus) -> Result<(), ManifestRefusal> {
    if locus.raw_normalization_id > 100 {
        return Err(ManifestRefusal::NormalizationId {
            id: locus.raw_normalization_id,
            name: locus.name.clone(),
        });
    }
    if (locus.assay_type != 0 && locus.address_b == 0)
        || (locus.assay_type == 0 && locus.address_b != 0)
    {
        return Err(ManifestRefusal::AssayType {
            assay_type: locus.assay_type,
            address_b: locus.address_b,
        });
    }
    Ok(())
}

/// One row of the normalization manifest CSV.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizationRow {
    pub index: usize,
    pub name: String,
    pub chromosome: String,
    pub position: i32,
    pub gentrain_score: f32,
    pub snp: String,
    pub illumina_strand: String,
    pub customer_strand: String,
    pub normalization_id: i32,
}

/// The CSV's header, which the tool writes verbatim.
pub const NORMALIZATION_HEADER: &str =
    "Index,Name,Chromosome,Position,GenTrain Score,SNP,ILMN Strand,Customer Strand,NormID";

/// One row as the file writes it: the score to four decimal places and nothing quoted.
pub fn normalization_line(row: &NormalizationRow) -> String {
    format!(
        "{},{},{},{},{:.4},{},{},{},{}",
        row.index,
        row.name,
        row.chromosome,
        row.position,
        row.gentrain_score,
        row.snp,
        row.illumina_strand,
        row.customer_strand,
        row.normalization_id
    )
}

/// What `CompareGtcFiles` answers a caller: a status rather than a list.
///
/// The differences themselves go through log4j, which holds the stream it was initialised with, so
/// they cannot be captured in-process; the status is what a caller reads and what the golden
/// records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comparison {
    /// Every compared getter agreed.
    Same,
    /// At least one did not.
    Different,
}

/// The fields the comparison EXCLUDES, which is what makes two runs of one chip comparable.
///
/// A different sample name is not a difference: it is one of the fields expected to differ between
/// two runs, so it is left out of the comparison rather than compared and forgiven.
pub const EXCLUDED_FROM_COMPARISON: [&str; 4] = [
    "getSampleName",
    "getSamplePlate",
    "getSampleWell",
    "getAutoCallDate",
];

/// `compareGTCFiles`, over whatever a caller can put side by side.
///
/// Every pair of values is compared but the excluded ones, and an array of a different length is a
/// difference like any other.
pub fn compare(fields: &[(&str, Vec<String>, Vec<String>)]) -> Comparison {
    for (name, left, right) in fields {
        if EXCLUDED_FROM_COMPARISON.contains(name) {
            continue;
        }
        if left != right {
            return Comparison::Different;
        }
    }
    Comparison::Same
}
