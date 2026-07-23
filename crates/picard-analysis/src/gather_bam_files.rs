//! `GatherBamFiles`.
//!
//! Ports `picard.sam.GatherBamFiles.doWork` at tag 3.4.0 for its block-copying fast path. When every
//! input is a BAM file, `determineBlockCopyingStatus` returns true and the tool calls
//! `BamFileIoUtils.gatherWithBlockCopying`, concatenating the inputs by copying their compressed
//! blocks (first file's header kept, the rest dropped) rather than re-encoding. That block copy is
//! [`htsjdk_bam::gather_bam_files`], proven byte-identical to htsjdk's `gatherWithBlockCopying`, so
//! this tool is a pass-through and byte-identity holds transitively against Picard with
//! `USE_JDK_DEFLATER=true`.
//!
//! The `gatherNormally` path (any input a SAM or CRAM, i.e. `determineBlockCopyingStatus == false`)
//! re-encodes through a `SAMFileWriter` and is a separate slice, deferred here. `CREATE_INDEX` /
//! `CREATE_MD5_FILE` default false. `IOUtil.unrollFiles` (expanding `.list` inputs) is a CLI concern;
//! this function takes the already-resolved BAM byte streams in order.

use htsjdk_bam::gather_bam_files as gather_block_copy;
use htsjdk_bam::reheader::ReheaderError;

/// Why `GatherBamFiles` could not run.
#[derive(Debug)]
pub enum GatherError {
    /// A block copy failed (an input was not a valid BAM, or had a defective tail).
    BlockCopy(ReheaderError),
}

impl From<ReheaderError> for GatherError {
    fn from(e: ReheaderError) -> Self {
        GatherError::BlockCopy(e)
    }
}

/// `GatherBamFiles.doWork` (block-copy fast path): the inputs concatenated into one BAM. Each input
/// is a whole raw BAM file (BGZF-framed, terminator included); they must already share a header, as
/// the fast path requires.
pub fn gather_bam_files(inputs: &[&[u8]]) -> Result<Vec<u8>, GatherError> {
    Ok(gather_block_copy(inputs)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::reader::BamReader;
    use htsjdk_bam::sam_file::read_sam_with;
    use htsjdk_bam::text_parse::ValidationStringency;
    use htsjdk_bam::writer::BamWriter;

    const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";

    fn build_bam(reads: &str) -> Vec<u8> {
        let sam = format!("{HEADER}{reads}");
        let (header, records) = read_sam_with(&sam, ValidationStringency::Lenient).unwrap();
        let mut w = BamWriter::new(Vec::new(), &header).unwrap();
        for r in &records {
            w.write(r).unwrap();
        }
        w.finish().unwrap()
    }

    #[test]
    fn gathers_records_from_all_inputs_under_the_first_header() {
        let a = build_bam("a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let b = build_bam("b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tTTTT\tIIII\n");
        let out = gather_bam_files(&[&a, &b]).unwrap();
        let decoded = htsjdk_bgzf::decompress_all(&out).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        assert!(reader.header.raw_text.contains("@SQ\tSN:chr1\tLN:1000"));
        let names: Vec<_> = reader.map(|r| r.unwrap().read_name).collect();
        assert_eq!(names, vec!["a".to_string(), "b".to_string()]);
    }

    #[test]
    fn a_non_bam_input_is_rejected() {
        let a = build_bam("a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        assert!(gather_bam_files(&[&a, b"not a bam"]).is_err());
    }
}
