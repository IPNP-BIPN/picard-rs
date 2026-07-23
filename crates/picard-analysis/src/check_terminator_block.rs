//! `CheckTerminatorBlock`.
//!
//! Ports `picard.sam.CheckTerminatorBlock.doWork` at tag 3.4.0: inspect a BAM file's BGZF tail and
//! report its [`FileTermination`], returning exit code `100` when the file is defective and `0`
//! otherwise. The classification itself is [`htsjdk_bgzf::check_termination`], which is proven
//! against htsjdk's `BlockCompressedInputStream.checkTermination`; this tool only maps the result to
//! the enum name Picard prints and to the exit code.

use htsjdk_bgzf::{check_termination, FileTermination};

/// The exit code Picard returns for a defective file (`FileTermination.DEFECTIVE`).
pub const DEFECTIVE_EXIT_CODE: i32 = 100;

/// The `FileTermination` name Picard writes, exactly as `Enum.name()` renders it.
fn termination_name(t: FileTermination) -> &'static str {
    match t {
        FileTermination::HasTerminatorBlock => "HAS_TERMINATOR_BLOCK",
        FileTermination::HasHealthyLastBlock => "HAS_HEALTHY_LAST_BLOCK",
        FileTermination::Defective => "DEFECTIVE",
    }
}

/// `CheckTerminatorBlock.doWork` for BAM input given as its raw bytes: the termination name Picard
/// prints to stderr, and the exit code (`100` if defective, else `0`).
pub fn check_terminator_block(data: &[u8]) -> (String, i32) {
    let term = check_termination(data);
    let code = if term == FileTermination::Defective {
        DEFECTIVE_EXIT_CODE
    } else {
        0
    };
    (termination_name(term).to_string(), code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::SamHeader;
    use htsjdk_bam::writer::BamWriter;

    /// A BAM produced by the byte-identical writer ends in the terminator block and exits 0.
    fn a_bam() -> Vec<u8> {
        let header = SamHeader::new();
        let writer = BamWriter::new(Vec::new(), &header).unwrap();
        writer.finish().unwrap()
    }

    #[test]
    fn a_well_formed_bam_has_a_terminator_block() {
        let bam = a_bam();
        assert_eq!(
            check_terminator_block(&bam),
            ("HAS_TERMINATOR_BLOCK".to_string(), 0)
        );
    }

    #[test]
    fn a_bam_missing_its_terminator_has_a_healthy_last_block() {
        let bam = a_bam();
        let trimmed = &bam[..bam.len() - 28];
        assert_eq!(
            check_terminator_block(trimmed),
            ("HAS_HEALTHY_LAST_BLOCK".to_string(), 0)
        );
    }

    #[test]
    fn a_truncated_bam_is_defective_and_exits_100() {
        let bam = a_bam();
        let trimmed = &bam[..bam.len() - 33];
        assert_eq!(
            check_terminator_block(trimmed),
            ("DEFECTIVE".to_string(), DEFECTIVE_EXIT_CODE)
        );
    }
}
