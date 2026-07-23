//! `AddCommentsToBam`.
//!
//! Ports `picard.sam.AddCommentsToBam.doWork` at tag 3.4.0: copy a BAM file, adding `@CO` comment
//! lines to its header. The tool reads the header, rejects any comment containing a newline
//! (`PicardException`), appends each with `SAMFileHeader.addComment`, and rewrites the file with
//! `BamFileIoUtils.reheaderBamFile` (a block copy of the record data, only the header re-encoded).
//! That block copy is [`htsjdk_bam::reheader_bam`], proven byte-identical to htsjdk's
//! `reheaderBamFile`; this tool only parses the header, validates and adds the comments, and calls
//! it. Byte-identity therefore holds against Picard with `USE_JDK_DEFLATER=true`, transitively.
//!
//! `CREATE_MD5_FILE` and `CREATE_INDEX` default to false, so the sidecar `.md5` / `.bai` outputs are
//! out of scope here. The tool's `.sam` rejection is a filename check in `doWork`; this function
//! takes BAM bytes, and a non-BAM input is rejected by `reheader_bam` (`ReheaderError::NotABam`).

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::reheader::{reheader_bam, ReheaderError};
use htsjdk_bgzf::decompress_all;

/// Why `AddCommentsToBam` could not run.
#[derive(Debug)]
pub enum AddCommentsError {
    /// A comment contained a newline (`"Comments can not contain a new line"`).
    NewlineInComment,
    /// The input could not be read as a BAM, or the reheader failed.
    Reheader(ReheaderError),
    /// The input BGZF could not be decoded.
    Bgzf(String),
}

impl From<ReheaderError> for AddCommentsError {
    fn from(e: ReheaderError) -> Self {
        AddCommentsError::Reheader(e)
    }
}

/// `AddCommentsToBam.doWork` for a BAM given as its raw bytes: the copied BAM with the `@CO`
/// comments added to its header.
pub fn add_comments_to_bam(
    input_bam: &[u8],
    comments: &[&str],
) -> Result<Vec<u8>, AddCommentsError> {
    let decoded =
        decompress_all(input_bam).map_err(|e| AddCommentsError::Bgzf(format!("{e:?}")))?;
    let reader = BamReader::new(&decoded).map_err(|e| AddCommentsError::Bgzf(format!("{e:?}")))?;
    let mut header = reader.header.text.clone();

    for comment in comments {
        if comment.contains('\n') {
            return Err(AddCommentsError::NewlineInComment);
        }
        header.add_comment(comment);
    }

    Ok(reheader_bam(&header, input_bam)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::sam_file::read_sam_with;
    use htsjdk_bam::text_parse::ValidationStringency;
    use htsjdk_bam::writer::BamWriter;

    const SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
        a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    fn build_bam(sam: &str) -> Vec<u8> {
        let (header, records) = read_sam_with(sam, ValidationStringency::Lenient).unwrap();
        let mut w = BamWriter::new(Vec::new(), &header).unwrap();
        for r in &records {
            w.write(r).unwrap();
        }
        w.finish().unwrap()
    }

    #[test]
    fn comments_are_added_and_records_preserved() {
        let input = build_bam(SAM);
        let out = add_comments_to_bam(&input, &["first", "second"]).unwrap();
        let decoded = decompress_all(&out).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        assert!(reader.header.raw_text.contains("@CO\tfirst"));
        assert!(reader.header.raw_text.contains("@CO\tsecond"));
        let names: Vec<_> = reader.map(|r| r.unwrap().read_name).collect();
        assert_eq!(names, vec!["a".to_string()]);
    }

    #[test]
    fn a_newline_in_a_comment_is_rejected() {
        let input = build_bam(SAM);
        assert!(matches!(
            add_comments_to_bam(&input, &["bad\ncomment"]),
            Err(AddCommentsError::NewlineInComment)
        ));
    }

    #[test]
    fn a_non_bam_input_is_rejected() {
        assert!(add_comments_to_bam(b"not a bam", &["c"]).is_err());
    }
}
