//! `BuildBamIndex`.
//!
//! Ports `picard.sam.BuildBamIndex.doWork` at tag 3.4.0: read a coordinate-sorted BAM and write its
//! `.bai`. The tool opens the BAM with source-in-records and calls `BAMIndexer.createIndex`; that
//! index build is [`htsjdk_bam::build_bam_index`], proven byte-identical to htsjdk's read-side
//! `BAMIndexer` (distinct from the write-side `setCreateIndex` index). This tool is therefore a
//! pass-through, and byte-identity against Picard with `USE_JDK_DEFLATER=true` holds transitively;
//! unit tests only, no new Picard oracle. CRAM inputs and the sort-order assertion (the primitive
//! already produces the correct index for a coordinate-sorted BAM) are out of scope here.

use htsjdk_bam::bin::{compute_indexing_bin, BIN_GENOMIC_SPAN, NO_ALIGNMENT_START};
use htsjdk_bam::build_index::BuildIndexError;
use htsjdk_bam::index::{BamIndexer, Chunk};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BAM_MAGIC;
use htsjdk_bgzf::BgzfReader;
use std::io::Read;

/// `BuildBamIndex.doWork` for a BAM given as its raw bytes: the `.bai` index.
pub fn build_bam_index(bam: &[u8]) -> Result<Vec<u8>, BuildIndexError> {
    build_bam_index_from(bam)
}

/// [`htsjdk_bam::build_index::build_bam_index`] over a reader rather than the whole file in memory:
/// the same loop, the same virtual file pointers around each record's bytes, the same errors. The
/// record's bytes go through one buffer reused from record to record.
pub fn build_bam_index_from(bam: impl Read) -> Result<Vec<u8>, BuildIndexError> {
    let mut reader = BgzfReader::new(bam);
    let reference_lengths = consume_header(&mut reader)?;
    let mut indexer = BamIndexer::new(reference_lengths.clone());
    let mut full = Vec::new();
    loop {
        let start_offset = reader.virtual_pos();
        let mut size_bytes = [0u8; 4];
        if !read_exact_or_eof(&mut reader, &mut size_bytes)? {
            break;
        }
        let block_size = i32::from_le_bytes(size_bytes);
        full.clear();
        full.extend_from_slice(&size_bytes);
        full.resize(4 + block_size as usize, 0);
        reader
            .read_exact(&mut full[4..])
            .map_err(|_| BuildIndexError::Truncated)?;
        let stop_offset = reader.virtual_pos();
        let (record, _) = BamRecord::decode(&full)?.ok_or(BuildIndexError::Truncated)?;
        let index_bin = index_bin_for(&record, &reference_lengths);
        indexer.process(
            record.reference_index,
            record.alignment_start,
            record.alignment_end(),
            index_bin,
            record.read_unmapped(),
            Chunk {
                start: start_offset,
                end: stop_offset,
            },
        );
    }
    Ok(indexer.finish())
}

fn read_exact_or_eof<R: Read>(r: &mut R, buf: &mut [u8]) -> Result<bool, BuildIndexError> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..]) {
            Ok(0) if filled == 0 => return Ok(false),
            Ok(0) | Err(_) => return Err(BuildIndexError::Truncated),
            Ok(n) => filled += n,
        }
    }
    Ok(true)
}

fn read_i32<R: Read>(r: &mut R) -> Result<i32, BuildIndexError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)
        .map_err(|_| BuildIndexError::Truncated)?;
    Ok(i32::from_le_bytes(b))
}

/// The header's binary part skipped, its reference lengths kept, as htsjdk-rs's `consume_header`.
fn consume_header<R: Read>(reader: &mut BgzfReader<R>) -> Result<Vec<i32>, BuildIndexError> {
    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .map_err(|_| BuildIndexError::Truncated)?;
    if magic != BAM_MAGIC {
        return Err(BuildIndexError::NotABam);
    }
    let skip = |reader: &mut BgzfReader<R>, n: i32| {
        std::io::copy(&mut reader.by_ref().take(n as u64), &mut std::io::sink())
            .map_err(|_| BuildIndexError::Truncated)
    };
    let l_text = read_i32(reader)?;
    skip(reader, l_text)?;
    let n_ref = read_i32(reader)?;
    let mut lengths = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref {
        let l_name = read_i32(reader)?;
        skip(reader, l_name)?;
        lengths.push(read_i32(reader)?);
    }
    Ok(lengths)
}

/// The indexing bin `BAMFileWriter` would record, as htsjdk-rs's `index_bin_for`.
fn index_bin_for(record: &BamRecord, reference_lengths: &[i32]) -> i32 {
    let too_large = record.reference_index >= 0
        && reference_lengths
            .get(record.reference_index as usize)
            .is_some_and(|&len| len > BIN_GENOMIC_SPAN);
    if too_large || record.alignment_start == NO_ALIGNMENT_START {
        0
    } else {
        compute_indexing_bin(record.alignment_start, record.alignment_end()).unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::{SamHeader, SequenceRecord};
    use htsjdk_bam::record::BamRecord;

    fn coord_bam() -> Vec<u8> {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 100_000));
        let mut w = crate::output::bam_writer(Vec::new(), &h).unwrap();
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
    fn the_streamed_index_is_htsjdk_rs_own() {
        let bam = coord_bam();
        assert_eq!(
            build_bam_index(&bam).unwrap(),
            htsjdk_bam::build_index::build_bam_index(&bam).unwrap()
        );
    }

    #[test]
    fn a_non_bam_input_is_rejected() {
        assert!(build_bam_index(b"not a bam").is_err());
    }
}
