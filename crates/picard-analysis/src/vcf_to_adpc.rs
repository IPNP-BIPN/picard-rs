//! `VcfToAdpc`: the `adpc.bin` an array VCF becomes.
//!
//! A sixteen-byte header and then EIGHTEEN bytes per sample per locus: two unsigned shorts of
//! raw intensity, three floats, one more unsigned short for the genotype. Everything is
//! little-endian and nothing is compressed, so the file is byte-comparable.
//!
//! Ported from `picard.arrays.VcfToAdpc`, `picard.arrays.illumina.IlluminaAdpcFileWriter`,
//! `picard.arrays.illumina.IlluminaGenotype` and `picard.arrays.illumina.InfiniumDataFile`.

/// The header, which is a literal in the writer rather than a magic number.
pub const HEADER: &[u8; 16] = b"1234567890123456";

/// `InfiniumDataFile.MAX_UNSIGNED_SHORT`, which an intensity is truncated to rather than refused.
pub const MAX_UNSIGNED_SHORT: i32 = 65535;

/// The four codes an Illumina genotype takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IlluminaGenotype {
    Aa = 0,
    Ab = 1,
    Bb = 2,
    Nn = 3,
}

impl IlluminaGenotype {
    pub fn code(self) -> u16 {
        self as u16
    }
}

/// Which genotype a call maps to, decided by the array's OWN alleles.
///
/// `ALLELE_A` and `ALLELE_B` are header fields, so the same `0/0` is an `AA` or a `BB` depending
/// on which allele the array called A. An allele that is the reference carries a trailing `*` in
/// those fields, which is stripped before the match, and an uncalled genotype is `NN` without
/// either field being read.
pub fn illumina_genotype(
    called: Option<(&str, &str)>,
    allele_a: &str,
    allele_b: &str,
) -> Option<IlluminaGenotype> {
    let Some((first, second)) = called else {
        return Some(IlluminaGenotype::Nn);
    };
    let a = allele_a.trim_end_matches('*');
    let b = allele_b.trim_end_matches('*');
    match (first == a, first == b) {
        (true, _) if second == a => Some(IlluminaGenotype::Aa),
        (true, _) if second == b => Some(IlluminaGenotype::Ab),
        (_, true) if second == a => Some(IlluminaGenotype::Ab),
        (_, true) if second == b => Some(IlluminaGenotype::Bb),
        _ => None,
    }
}

/// `getUnsignedShortAttributeAsInt`: over the limit is TRUNCATED, under zero is refused.
pub fn raw_intensity(value: i32) -> Option<u16> {
    if value < 0 {
        return None;
    }
    Some(value.min(MAX_UNSIGNED_SHORT) as u16)
}

/// One record, as the writer lays it out.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Record {
    pub a_intensity: u16,
    pub b_intensity: u16,
    /// The normalized intensities are OPTIONAL where the raw ones are not, and an absent one is
    /// written as a NaN rather than skipped: the record's width does not move.
    pub a_normalized: f32,
    pub b_normalized: f32,
    pub gc_score: f32,
    pub genotype: IlluminaGenotype,
}

/// The eighteen bytes of one record, little-endian throughout: 2 + 2 + 4 + 4 + 4 + 2.
pub fn write_record(record: &Record) -> Vec<u8> {
    let mut out = Vec::with_capacity(18);
    out.extend(record.a_intensity.to_le_bytes());
    out.extend(record.b_intensity.to_le_bytes());
    out.extend(record.a_normalized.to_le_bytes());
    out.extend(record.b_normalized.to_le_bytes());
    out.extend(record.gc_score.to_le_bytes());
    out.extend(record.genotype.code().to_le_bytes());
    out
}

/// The whole file: the header, then every record in the order they were written.
///
/// The order is sample-major: the tool walks the VCF once per sample, so two samples of one locus
/// are two records in a row and not one interleaved pair.
pub fn write_file(records: &[Record]) -> Vec<u8> {
    let mut out = HEADER.to_vec();
    for record in records {
        out.extend(write_record(record));
    }
    out
}

/// The samples file: one name per line, with no trailing newline.
pub fn samples_file(samples: &[&str]) -> String {
    samples.join("\n")
}

/// The marker file: the number of loci, as a bare number.
pub fn markers_file(loci: usize) -> String {
    loci.to_string()
}

/// What a refusal leaves behind, which is an exit code and nothing else.
///
/// The tool catches its own exception, logs it and returns one, and the log reaches no stream a
/// caller can capture, so these messages are carried here for a reader rather than compared: the
/// golden holds the code.
pub const NO_RECORDS_MESSAGE: &str = "Found no records in VCF";
pub const DIFFERING_LOCI_MESSAGE: &str = "VCFs have differing number of loci";
/// The exit code a refusal returns.
pub const REFUSAL_EXIT_CODE: i32 = 1;
