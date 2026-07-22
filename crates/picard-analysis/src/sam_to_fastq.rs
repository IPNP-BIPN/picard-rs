//! `SamToFastq`, the unpaired default path.
//!
//! Ports `picard.sam.SamToFastq.handleRecord` and `writeRecord` at tag 3.4.0, for **unpaired reads
//! with default options** written to a single FASTQ. This is the opener of the read-data-manipulation
//! archetype: a record transform whose byte-identity is the whole output file, here FASTQ text.
//!
//! The paired handling (`SECOND_END_FASTQ`, `INTERLEAVE`, `OUTPUT_PER_RG`), the clipping
//! (`CLIPPING_ATTRIBUTE`), and the trimming (`READ1_TRIM`, `QUALITY`, `*_MAX_BASES_TO_WRITE`) are
//! separate surfaces and are **not** claimed here; the collector asserts it was given only unpaired
//! reads so it cannot silently emit a half-right paired file.
//!
//! The default per-read transform, from `writeRecord`:
//!   - the sequence is `getReadString`, the qualities `getBaseQualityString`;
//!   - with `RE_REVERSE` on (the default) a negative-strand read is reverse-complemented and its
//!     qualities reversed, so the FASTQ is in original sequencing orientation;
//!   - the record is `@name` / sequence / `+` / qualities.
//!
//! And the filtering, from `handleRecord`: a secondary or supplementary read is dropped unless
//! `INCLUDE_NON_PRIMARY_ALIGNMENTS`, and a vendor-failed read unless `INCLUDE_NON_PF_READS`, both
//! off by default.

use htsjdk_bam::fastq::{phred_to_fastq, write_record, FastqRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::reverse_complement;

const READ_PAIRED: u16 = 0x1;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;
const FIRST_OF_PAIR: u16 = 0x40;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

const NULL_STRING: &str = "*";

/// The options that gate the unpaired path, all defaulting to Picard's defaults.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `RE_REVERSE`, default true: put negative-strand reads back in sequencing orientation.
    pub re_reverse: bool,
    /// `INCLUDE_NON_PRIMARY_ALIGNMENTS`, default false.
    pub include_non_primary: bool,
    /// `INCLUDE_NON_PF_READS`, default false.
    pub include_non_pf: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            re_reverse: true,
            include_non_primary: false,
            include_non_pf: false,
        }
    }
}

/// `SAMRecord.getReadString()`.
fn read_string(rec: &BamRecord) -> String {
    if rec.read_bases.is_empty() {
        NULL_STRING.to_string()
    } else {
        String::from_utf8_lossy(&rec.read_bases).into_owned()
    }
}

/// `SAMRecord.getBaseQualityString()`.
fn base_quality_string(rec: &BamRecord) -> String {
    if rec.base_qualities.is_empty() {
        NULL_STRING.to_string()
    } else {
        phred_to_fastq(&rec.base_qualities)
    }
}

/// `writeRecord(read, mateNumber, writer, 0, null)`: the FASTQ text for one read.
///
/// `mate_number` is `None` for an unpaired read (`@name`) and `Some(1)`/`Some(2)` for a paired
/// read (`@name/1`, `@name/2`).
fn write_one(rec: &BamRecord, mate_number: Option<u32>, opts: &Options) -> String {
    let mut seq = read_string(rec);
    let mut quals = base_quality_string(rec);

    if opts.re_reverse && rec.flags & READ_NEGATIVE_STRAND != 0 {
        // SequenceUtil.reverseComplement(String) on the bases, StringUtil.reverseString on the quals.
        let mut bytes = seq.into_bytes();
        reverse_complement(&mut bytes);
        seq = String::from_utf8(bytes).expect("bases stay ASCII");
        quals = quals.chars().rev().collect();
    }

    let name = match mate_number {
        None => rec.read_name.clone(),
        Some(n) => format!("{}/{}", rec.read_name, n),
    };
    write_record(&FastqRecord {
        read_name: Some(name),
        read_string: Some(seq),
        quality_header: Some(String::new()),
        quality_string: Some(quals),
    })
}

/// Whether `handleRecord` drops this record before writing: a secondary or supplementary read
/// without `INCLUDE_NON_PRIMARY_ALIGNMENTS`, or a vendor-failed read without `INCLUDE_NON_PF_READS`.
fn is_dropped(rec: &BamRecord, opts: &Options) -> bool {
    let secondary_or_supplementary =
        rec.flags & (SECONDARY_ALIGNMENT | SUPPLEMENTARY_ALIGNMENT) != 0;
    (secondary_or_supplementary && !opts.include_non_primary)
        || (rec.flags & READ_FAILS_VENDOR_QUALITY != 0 && !opts.include_non_pf)
}

/// Runs the unpaired default path over the records in file order, returning the FASTQ file text.
///
/// Panics on a paired read, since the paired writers are not part of this port and emitting a
/// fragment file for paired input would be a silent, wrong result rather than a missing feature.
pub fn sam_to_fastq_unpaired(records: &[BamRecord], opts: &Options) -> String {
    let mut out = String::new();
    for rec in records {
        if is_dropped(rec, opts) {
            continue;
        }
        assert!(
            rec.flags & READ_PAIRED == 0,
            "sam_to_fastq_unpaired given a paired read; the paired path is not ported"
        );
        out.push_str(&write_one(rec, None, opts));
    }
    out
}

/// The paired default path with `SECOND_END_FASTQ` (two separate files), returning
/// `(first_of_pair_fastq, second_of_pair_fastq)`.
///
/// Ports `handleRecord`'s paired branch: reads are matched by name through a first-seen map, and a
/// pair is emitted only when its second mate arrives, so the output order in each file follows the
/// order pairs *complete* in the input. The first-of-pair mate is written to the first file with a
/// `/1` suffix and the second-of-pair to the second file with `/2`, regardless of which physically
/// arrived second. A read filtered by [`is_dropped`] never enters the map, so a pair with one
/// dropped mate is left incomplete and unwritten, matching htsjdk's MATE_NOT_FOUND handling.
pub fn sam_to_fastq_paired(records: &[BamRecord], opts: &Options) -> (String, String) {
    use std::collections::HashMap;

    let mut first_seen: HashMap<&str, &BamRecord> = HashMap::new();
    let mut first_end = String::new();
    let mut second_end = String::new();

    for rec in records {
        if is_dropped(rec, opts) {
            continue;
        }
        assert!(
            rec.flags & READ_PAIRED != 0,
            "sam_to_fastq_paired given an unpaired read"
        );
        match first_seen.remove(rec.read_name.as_str()) {
            None => {
                first_seen.insert(rec.read_name.as_str(), rec);
            }
            Some(first) => {
                let (read1, read2) = if rec.flags & FIRST_OF_PAIR != 0 {
                    (rec, first)
                } else {
                    (first, rec)
                };
                first_end.push_str(&write_one(read1, Some(1), opts));
                second_end.push_str(&write_one(read2, Some(2), opts));
            }
        }
    }
    (first_end, second_end)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, flags: u16, bases: &[u8], quals: &[u8]) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            read_bases: bases.to_vec(),
            base_qualities: quals.to_vec(),
            ..Default::default()
        }
    }

    #[test]
    fn a_forward_read_is_emitted_verbatim() {
        let r = rec("r1", 0, b"ACGT", &[40, 40, 30, 20]);
        assert_eq!(
            sam_to_fastq_unpaired(&[r], &Options::default()),
            "@r1\nACGT\n+\nII?5\n"
        );
    }

    #[test]
    fn a_negative_strand_read_is_reverse_complemented_and_its_quals_reversed() {
        // Stored ACGT with quals 40,30,20,10; revcomp = ACGT, quals reversed = 10,20,30,40.
        let r = rec("r1", READ_NEGATIVE_STRAND, b"ACGT", &[40, 30, 20, 10]);
        // revcomp(ACGT) = ACGT; quals "I?5+" reversed = "+5?I".
        assert_eq!(
            sam_to_fastq_unpaired(&[r], &Options::default()),
            "@r1\nACGT\n+\n+5?I\n"
        );
    }

    #[test]
    fn secondary_supplementary_and_vendor_fail_are_dropped_by_default() {
        let recs = [
            rec("keep", 0, b"AC", &[40, 40]),
            rec("sec", SECONDARY_ALIGNMENT, b"AC", &[40, 40]),
            rec("sup", SUPPLEMENTARY_ALIGNMENT, b"AC", &[40, 40]),
            rec("qc", READ_FAILS_VENDOR_QUALITY, b"AC", &[40, 40]),
        ];
        assert_eq!(
            sam_to_fastq_unpaired(&recs, &Options::default()),
            "@keep\nAC\n+\nII\n"
        );
    }

    #[test]
    #[should_panic(expected = "paired read")]
    fn a_paired_read_is_rejected() {
        let r = rec("r1", READ_PAIRED, b"AC", &[40, 40]);
        sam_to_fastq_unpaired(&[r], &Options::default());
    }

    #[test]
    fn a_pair_is_split_across_the_two_files_by_first_second_of_pair() {
        // The second-of-pair mate physically arrives first, but read1 is still the first-of-pair.
        let second = rec("p1", READ_PAIRED | 0x80, b"TT", &[30, 30]);
        let first = rec("p1", READ_PAIRED | FIRST_OF_PAIR, b"AA", &[40, 40]);
        let (r1, r2) = sam_to_fastq_paired(&[second, first], &Options::default());
        assert_eq!(r1, "@p1/1\nAA\n+\nII\n");
        assert_eq!(r2, "@p1/2\nTT\n+\n??\n");
    }

    #[test]
    fn a_pair_with_a_dropped_mate_is_not_written() {
        // The first-of-pair is a supplementary read (dropped), so the pair never completes.
        let first = rec(
            "p1",
            READ_PAIRED | FIRST_OF_PAIR | SUPPLEMENTARY_ALIGNMENT,
            b"AA",
            &[40, 40],
        );
        let second = rec("p1", READ_PAIRED | 0x80, b"TT", &[30, 30]);
        let (r1, r2) = sam_to_fastq_paired(&[first, second], &Options::default());
        assert_eq!(r1, "");
        assert_eq!(r2, "");
    }
}
