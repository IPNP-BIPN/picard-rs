//! `CleanSam`.
//!
//! Ports `picard.sam.CleanSam.doWork` at tag 3.4.0, with the per-record work in
//! `AbstractAlignmentMerger.createNewCigarsIfMapsOffEndOfReference` and
//! `CigarUtil.softClipEndOfRead`. Two edits per record: soft-clip an alignment that hangs off the
//! end of its reference, and set `MAPQ` to 0 on an unmapped read that still carries a nonzero one.
//!
//! `doWork` keeps the input header unchanged, adds no `@PG` and no timestamp, and does not re-sort,
//! so the whole output is comparable raw.
//!
//! This ports the **read's** off-the-end clip and the MAPQ fix. The mate's clip (via the `MC` tag)
//! is a separate surface: a record carrying an `MC` tag would need it too, so the port asserts there
//! is none rather than silently leaving a mate cigar unclipped.

use htsjdk_bam::cigar::{soft_clip_end_of_read, Cigar, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::tag::Tag;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

const READ_UNMAPPED: u16 = 0x4;

/// `createNewCigarIfMapsOffEndOfReference` for the read, then the unmapped-MAPQ fix.
fn clean_record(rec: &mut BamRecord, header: &SamHeader) {
    assert!(
        rec.tags.get(Tag::new(b"MC")).is_none(),
        "CleanSam: a record carries an MC tag; the mate clip is not ported"
    );

    if rec.flags & READ_UNMAPPED == 0 {
        let ref_len = header.sequences[rec.reference_index as usize].length;
        // getAlignmentEnd = alignmentStart + referenceLength - 1.
        let alignment_end = rec.alignment_start + rec.cigar.reference_length() as i32 - 1;
        let overhang = alignment_end - ref_len;
        if overhang > 0 {
            // 1-based first read base to clip.
            let mut clip_from = rec.read_length() as i32 - overhang + 1;
            // If the last element is already a soft clip, subtract it from clipFrom.
            if let Some(last) = rec.cigar.elements.last() {
                if last.op == Op::S {
                    clip_from -= last.length as i32;
                }
            }
            let new_elements = soft_clip_end_of_read(clip_from, &rec.cigar.elements);
            rec.cigar = Cigar::new(new_elements);
        }
    }

    if rec.flags & READ_UNMAPPED != 0 && rec.mapping_quality != 0 {
        rec.mapping_quality = 0;
    }
}

/// `CleanSam.doWork` for SAM I/O: clean every record in input order and rewrite.
///
/// Reads at LENIENT stringency, as CleanSam does (it downgrades STRICT to LENIENT), so it can
/// accept the very records it exists to fix, such as an unmapped read with a nonzero MAPQ.
pub fn clean_sam(input_sam: &str) -> Result<String, ParseError> {
    let (header, records) = clean(input_sam)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// `CleanSam.doWork` up to the write: the header and the cleaned records.
fn clean(input_sam: &str) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), ParseError> {
    use rayon::prelude::*;
    let (header, mut records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    // Each record is cleaned independently from the (immutable, shared) header, and `par_iter_mut`
    // mutates them in place without reordering, so the result is byte-identical to a serial loop
    // regardless of how the work is split across cores. See decision 0006.
    records
        .par_iter_mut()
        .for_each(|rec| clean_record(rec, &header));
    Ok((header, records))
}

/// `CleanSam.doWork` for SAM input and **BAM** output. The same cleaning, written through htsjdk-rs's
/// byte-identical `BamWriter`; CleanSam adds no `@PG`. Byte-identity to Picard's
/// `USE_JDK_DEFLATER=true` output follows transitively: the records are the ones `clean_sam` already
/// reproduces (the CleanSam oracle), and the `BamWriter` is proven byte-identical over arbitrary
/// records (the SamFormatConverter oracle and htsjdk-rs's `every_file_is_byte_identical_to_htsjdks`).
pub fn clean_sam_to_bam(input_sam: &str) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, records) = clean(input_sam)?;
    let mut writer = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        writer
            .write(rec)
            .expect("records that parsed re-encode as BAM");
    }
    Ok(writer
        .finish()
        .expect("in-memory BAM writer never fails to finish"))
}

#[cfg(test)]
mod tests {
    use super::*;

    // chr1 is 150 long, so a 36M read at 130 ends at 165, hanging 15 off the end.
    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:150\n\
        offend\t0\tchr1\t130\t60\t36M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGT\t\
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        inside\t0\tchr1\t10\t60\t36M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGT\t\
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        unmapped\t4\t*\t0\t30\t*\t*\t0\t0\tACGT\tIIII\n";

    fn cigars(sam: &str) -> Vec<(String, String)> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0].to_string(), f[5].to_string())
            })
            .collect()
    }

    #[test]
    fn a_read_hanging_off_the_end_is_soft_clipped() {
        let out = clean_sam(INPUT).unwrap();
        let c = cigars(&out);
        // 130..165 over a 150bp contig: overhang 15, clipFrom = 36-15+1 = 22 -> 21M15S.
        assert_eq!(c.iter().find(|(n, _)| n == "offend").unwrap().1, "21M15S");
        // A read comfortably inside is untouched.
        assert_eq!(c.iter().find(|(n, _)| n == "inside").unwrap().1, "36M");
    }

    /// The parallel per-record cleaning is byte-identical to a serial loop, on an input large
    /// enough that rayon actually splits the work. This is the invariant that lets the transform go
    /// multicore without weakening the byte-identity claim (decision 0006).
    #[test]
    fn parallel_and_serial_cleaning_produce_the_same_bytes() {
        use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
        use htsjdk_bam::header::{SamHeader, SequenceRecord};
        use htsjdk_bam::sam_file::write_sam;
        use rayon::prelude::*;

        let mut header = SamHeader::new();
        header.set_sort_order("coordinate");
        header.sequences.push(SequenceRecord::new("chr1", 150));

        let make = || {
            (0..5000)
                .map(|i| BamRecord {
                    read_name: format!("r{i}"),
                    reference_index: 0,
                    alignment_start: 100 + (i % 45), // most reads hang off the 150bp contig
                    mapping_quality: 60,
                    cigar: Cigar::new(vec![CigarElement {
                        length: 36,
                        op: Op::M,
                    }]),
                    read_bases: vec![b'A'; 36],
                    base_qualities: vec![40; 36],
                    ..Default::default()
                })
                .collect::<Vec<_>>()
        };

        let mut serial = make();
        for r in &mut serial {
            clean_record(r, &header);
        }

        let mut parallel = make();
        parallel
            .par_iter_mut()
            .for_each(|r| clean_record(r, &header));

        assert_eq!(write_sam(&header, &serial), write_sam(&header, &parallel));
    }

    #[test]
    fn an_unmapped_read_has_its_mapping_quality_zeroed() {
        let out = clean_sam(INPUT).unwrap();
        let row = out.lines().find(|l| l.starts_with("unmapped")).unwrap();
        let mapq = row.split('\t').nth(4).unwrap();
        assert_eq!(mapq, "0", "unmapped read kept a nonzero MAPQ: {row}");
    }

    /// The BAM output decodes back to exactly the SAM output; the writer's byte-identity to htsjdk is
    /// proven elsewhere.
    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        let sam = clean_sam(INPUT).unwrap();
        let bam = clean_sam_to_bam(INPUT).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bam).expect("bam decompresses");
        let reader = htsjdk_bam::reader::BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(
            htsjdk_bam::sam_file::write_sam(&header, &records).unwrap(),
            sam
        );
    }
}
