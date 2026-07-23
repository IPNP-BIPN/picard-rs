//! `BuildBamIndex`.
//!
//! Ports `picard.sam.BuildBamIndex.doWork` at tag 3.4.0: read a coordinate-sorted BAM and write its
//! `.bai`. The tool opens the BAM with source-in-records and calls `BAMIndexer.createIndex`; that
//! index build is [`htsjdk_bam::build_bam_index`], proven byte-identical to htsjdk's read-side
//! `BAMIndexer` (distinct from the write-side `setCreateIndex` index). This tool is therefore a
//! pass-through, and byte-identity against Picard with `USE_JDK_DEFLATER=true` holds transitively;
//! unit tests only, no new Picard oracle. CRAM inputs and the sort-order assertion (the primitive
//! already produces the correct index for a coordinate-sorted BAM) are out of scope here.

use htsjdk_bam::build_index::{build_bam_index as build_index, BuildIndexError};

/// `BuildBamIndex.doWork` for a BAM given as its raw bytes: the `.bai` index.
pub fn build_bam_index(bam: &[u8]) -> Result<Vec<u8>, BuildIndexError> {
    build_index(bam)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::{SamHeader, SequenceRecord};
    use htsjdk_bam::record::BamRecord;
    use htsjdk_bam::writer::BamWriter;

    fn coord_bam() -> Vec<u8> {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 100_000));
        let mut w = BamWriter::new(Vec::new(), &h).unwrap();
        for (name, start) in [("a", 10), ("b", 500)] {
            let rec = BamRecord {
                read_name: name.to_string(),
                reference_index: 0,
                alignment_start: start,
                mapping_quality: 60,
                cigar: htsjdk_bam::cigar::Cigar::new(vec![htsjdk_bam::cigar::CigarElement {
                    length: 10,
                    op: htsjdk_bam::cigar::Op::M,
                }]),
                ..BamRecord::default()
            };
            w.write(&rec).unwrap();
        }
        w.finish().unwrap()
    }

    #[test]
    fn writes_a_bai_for_a_coordinate_sorted_bam() {
        let bai = build_bam_index(&coord_bam()).unwrap();
        assert_eq!(&bai[..4], b"BAI\x01");
    }

    #[test]
    fn a_non_bam_input_is_rejected() {
        assert!(build_bam_index(b"not a bam").is_err());
    }
}
