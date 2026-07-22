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

const READ_NEGATIVE_STRAND: u16 = 0x10;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;
const READ_PAIRED: u16 = 0x1;
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

/// `writeRecord(read, null, writer, 0, null)`: the FASTQ text for one unpaired read.
fn write_one(rec: &BamRecord, opts: &Options) -> String {
    let mut seq = read_string(rec);
    let mut quals = base_quality_string(rec);

    if opts.re_reverse && rec.flags & READ_NEGATIVE_STRAND != 0 {
        // SequenceUtil.reverseComplement(String) on the bases, StringUtil.reverseString on the quals.
        let mut bytes = seq.into_bytes();
        reverse_complement(&mut bytes);
        seq = String::from_utf8(bytes).expect("bases stay ASCII");
        quals = quals.chars().rev().collect();
    }

    write_record(&FastqRecord {
        read_name: Some(rec.read_name.clone()),
        read_string: Some(seq),
        quality_header: Some(String::new()),
        quality_string: Some(quals),
    })
}

/// Runs the unpaired default path over the records in file order, returning the FASTQ file text.
///
/// Panics on a paired read, since the paired writers are not part of this port and emitting a
/// fragment file for paired input would be a silent, wrong result rather than a missing feature.
pub fn sam_to_fastq_unpaired(records: &[BamRecord], opts: &Options) -> String {
    let mut out = String::new();
    for rec in records {
        // isSecondaryOrSupplementary && !INCLUDE_NON_PRIMARY_ALIGNMENTS
        let secondary_or_supplementary =
            rec.flags & (SECONDARY_ALIGNMENT | SUPPLEMENTARY_ALIGNMENT) != 0;
        if secondary_or_supplementary && !opts.include_non_primary {
            continue;
        }
        // vendor fail && !INCLUDE_NON_PF_READS
        if rec.flags & READ_FAILS_VENDOR_QUALITY != 0 && !opts.include_non_pf {
            continue;
        }
        assert!(
            rec.flags & READ_PAIRED == 0,
            "sam_to_fastq_unpaired given a paired read; the paired path is not ported"
        );
        out.push_str(&write_one(rec, opts));
    }
    out
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
}
