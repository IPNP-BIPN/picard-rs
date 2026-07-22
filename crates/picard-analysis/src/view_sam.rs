//! `ViewSam`.
//!
//! Ports `picard.sam.ViewSam.writeSamText` at tag 3.4.0, for the **default I/O path** (no
//! `INTERVAL_LIST`): print a SAM/BAM as SAM text to stdout, optionally filtering records by whether
//! they are aligned and whether they passed vendor quality. The header is written **verbatim** (from
//! `getTextHeader`) and no `@PG` and no timestamp are added, so the whole output is comparable raw.
//!
//! The filter is a pure per-record predicate reading only the record's flags, so it runs on all cores
//! and stays byte-identical: rayon's `collect` keeps the surviving records in their original order
//! (decision 0006).
//!
//! The `INTERVAL_LIST` path (which walks an index) is a separate surface and is not ported here.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text::write_alignment;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use rayon::prelude::*;

const READ_UNMAPPED: u16 = 0x4;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;

/// `ViewSam.AlignmentStatus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlignmentStatus {
    #[default]
    All,
    Aligned,
    Unaligned,
}

/// `ViewSam.PfStatus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PfStatus {
    #[default]
    All,
    /// Passed vendor quality (the vendor-fail flag is clear).
    Pf,
    /// Failed vendor quality (the vendor-fail flag is set).
    NonPf,
}

/// `ViewSam`'s options, defaulting to Picard's defaults (print everything).
#[derive(Debug, Clone, Copy, Default)]
pub struct Options {
    pub alignment_status: AlignmentStatus,
    pub pf_status: PfStatus,
    /// `HEADER_ONLY`: print only the header.
    pub header_only: bool,
    /// `RECORDS_ONLY`: print only the records.
    pub records_only: bool,
}

/// Whether `writeSamText` prints this record given the two status filters.
fn keep(rec: &BamRecord, opts: &Options) -> bool {
    let unmapped = rec.flags & READ_UNMAPPED != 0;
    match opts.alignment_status {
        AlignmentStatus::Aligned if unmapped => return false,
        AlignmentStatus::Unaligned if !unmapped => return false,
        _ => {}
    }

    let vendor_fail = rec.flags & READ_FAILS_VENDOR_QUALITY != 0;
    match opts.pf_status {
        PfStatus::Pf if vendor_fail => return false,
        PfStatus::NonPf if !vendor_fail => return false,
        _ => {}
    }
    true
}

/// `ViewSam.writeSamText` for the default (no-interval) path, returning the SAM text.
pub fn view_sam(input_sam: &str, opts: &Options) -> Result<String, ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    // RECORDS_ONLY suppresses the header; otherwise print it verbatim (getTextHeader).
    if opts.records_only {
        // No header: encode the records only, keeping the writer's per-record format.
        if opts.header_only {
            return Ok(String::new());
        }
        let body: String = records
            .par_iter()
            .filter(|rec| keep(rec, opts))
            .map(|rec| render(&header, rec))
            .collect();
        return Ok(body);
    }

    let mut out = header.encode();
    // HEADER_ONLY stops after the header.
    if opts.header_only {
        return Ok(out);
    }
    // Each surviving record renders independently; rayon's collect preserves order (decision 0006).
    let body: String = records
        .par_iter()
        .filter(|rec| keep(rec, opts))
        .map(|rec| render(&header, rec))
        .collect();
    out.push_str(&body);
    Ok(out)
}

/// One record as `getSAMString` renders it: the alignment line and its trailing newline, with
/// reference names resolved against the header.
fn render(header: &htsjdk_bam::header::SamHeader, rec: &BamRecord) -> String {
    let name_of = |index: i32| -> &str {
        if index < 0 {
            "*"
        } else {
            header
                .sequences
                .get(index as usize)
                .map(|s| s.name.as_str())
                .unwrap_or("*")
        }
    };
    let mut line = write_alignment(
        rec,
        name_of(rec.reference_index),
        name_of(rec.mate_reference_index),
    )
    .expect("records that parsed re-encode as SAM text");
    line.push('\n');
    line
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        m1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        qc\t512\tchr1\t150\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        m2\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

    fn names(sam: &str) -> Vec<&str> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap())
            .collect()
    }

    #[test]
    fn the_default_prints_the_header_and_every_record() {
        let out = view_sam(INPUT, &Options::default()).unwrap();
        assert!(out.starts_with("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"));
        assert_eq!(names(&out), ["m1", "qc", "m2", "u1"]);
    }

    #[test]
    fn aligned_drops_the_unmapped_read_only() {
        let opts = Options {
            alignment_status: AlignmentStatus::Aligned,
            ..Default::default()
        };
        // The vendor-fail read qc is still aligned, so it stays; only u1 goes.
        assert_eq!(names(&view_sam(INPUT, &opts).unwrap()), ["m1", "qc", "m2"]);
    }

    #[test]
    fn unaligned_keeps_only_the_unmapped_read() {
        let opts = Options {
            alignment_status: AlignmentStatus::Unaligned,
            ..Default::default()
        };
        assert_eq!(names(&view_sam(INPUT, &opts).unwrap()), ["u1"]);
    }

    #[test]
    fn pf_drops_the_vendor_fail_read() {
        let opts = Options {
            pf_status: PfStatus::Pf,
            ..Default::default()
        };
        assert_eq!(names(&view_sam(INPUT, &opts).unwrap()), ["m1", "m2", "u1"]);
    }

    #[test]
    fn nonpf_keeps_only_the_vendor_fail_read() {
        let opts = Options {
            pf_status: PfStatus::NonPf,
            ..Default::default()
        };
        assert_eq!(names(&view_sam(INPUT, &opts).unwrap()), ["qc"]);
    }

    #[test]
    fn header_only_stops_after_the_header() {
        let opts = Options {
            header_only: true,
            ..Default::default()
        };
        let out = view_sam(INPUT, &opts).unwrap();
        assert_eq!(out, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    }
}
