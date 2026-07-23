//! `BamIndexStats`.
//!
//! Ports `picard.sam.BamIndexStats.doWork` (which delegates to `BAMIndexMetaData.printIndexStats`) at
//! tag 3.4.0: print the per-reference aligned/unaligned record counts of a BAM's index, then the
//! total no-coordinate record count. The counts come from the `.bai`; the port builds the index with
//! [`htsjdk_bam::build_bam_index`] (whose metadata is identical to the write-side index it would read,
//! since the counts do not depend on the virtual offsets) and reads it with
//! [`htsjdk_bam::parse_bai_metadata`], both oracle-proven, so this tool is a thin composition.
//!
//! A reference with no reads has no metadata chunk in the index, and htsjdk prints `Aligned= 0
//! Unaligned= 0` for it (its `getMetaData` returns a zeroed record, not null), so `None` maps to
//! `(0, 0)`. CRAM/CSI are out of scope.

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::{parse_bai_metadata, BuildIndexError};

/// Why `BamIndexStats` could not run.
#[derive(Debug)]
pub enum BamIndexStatsError {
    Build(BuildIndexError),
    Bgzf(String),
    Bai(String),
}

impl From<BuildIndexError> for BamIndexStatsError {
    fn from(e: BuildIndexError) -> Self {
        BamIndexStatsError::Build(e)
    }
}

/// `printIndexStats` for a BAM given as its raw bytes: the exact text the tool prints.
pub fn bam_index_stats(bam: &[u8]) -> Result<String, BamIndexStatsError> {
    let decoded =
        htsjdk_bgzf::decompress_all(bam).map_err(|e| BamIndexStatsError::Bgzf(format!("{e:?}")))?;
    let reader =
        BamReader::new(&decoded).map_err(|e| BamIndexStatsError::Bgzf(format!("{e:?}")))?;
    let sequences = reader.header.text.sequences.clone();

    let bai = build_bam_index(bam)?;
    let stats = parse_bai_metadata(&bai).map_err(|e| BamIndexStatsError::Bai(format!("{e:?}")))?;

    let mut out = String::new();
    for (i, seq) in sequences.iter().enumerate() {
        let (aligned, unaligned) = stats
            .references
            .get(i)
            .and_then(|m| m.as_ref())
            .map(|m| (m.aligned, m.unaligned))
            .unwrap_or((0, 0));
        out.push_str(&format!(
            "{} length=\t{}\tAligned= {}\tUnaligned= {}\n",
            seq.name, seq.length, aligned, unaligned
        ));
    }
    out.push_str(&format!(
        "NoCoordinateCount= {}\n",
        stats.no_coordinate_records
    ));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
    use htsjdk_bam::header::{SamHeader, SequenceRecord};
    use htsjdk_bam::record::{BamRecord, READ_UNMAPPED_FLAG};
    use htsjdk_bam::writer::BamWriter;

    fn m10() -> Cigar {
        Cigar::new(vec![CigarElement {
            length: 10,
            op: Op::M,
        }])
    }

    #[test]
    fn prints_aligned_unaligned_and_no_coordinate_counts() {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 100_000));
        h.sequences.push(SequenceRecord::new("chr2", 50_000));
        let mut w = BamWriter::new(Vec::new(), &h).unwrap();
        for start in [10, 500] {
            w.write(&BamRecord {
                read_name: "m".into(),
                reference_index: 0,
                alignment_start: start,
                mapping_quality: 60,
                cigar: m10(),
                ..BamRecord::default()
            })
            .unwrap();
        }
        for _ in 0..2 {
            w.write(&BamRecord {
                read_name: "u".into(),
                reference_index: -1,
                flags: READ_UNMAPPED_FLAG,
                ..BamRecord::default()
            })
            .unwrap();
        }
        let bam = w.finish().unwrap();

        let text = bam_index_stats(&bam).unwrap();
        assert_eq!(
            text,
            "chr1 length=\t100000\tAligned= 2\tUnaligned= 0\n\
             chr2 length=\t50000\tAligned= 0\tUnaligned= 0\n\
             NoCoordinateCount= 2\n"
        );
    }
}
