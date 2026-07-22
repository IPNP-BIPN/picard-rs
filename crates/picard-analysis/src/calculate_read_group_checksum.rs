//! `CalculateReadGroupChecksum`.
//!
//! Ports `picard.sam.CalculateReadGroupChecksum.doWork` at tag 3.4.0. The whole tool is a thin
//! wrapper: read the header, compute an MD5 over its read groups
//! ([`calculate_read_group_record_checksum`](htsjdk_bam::read_group_checksum::calculate_read_group_record_checksum),
//! ported in htsjdk-rs), and write the hash text. The output file (`<input>.read_group_md5`) holds
//! **only** the 32-character hex digest, with no trailing newline, so it is compared raw.
//!
//! Changing, adding, or removing a read group changes the digest, which is the point: it is a quick
//! fingerprint of a file's read-group metadata.

use htsjdk_bam::read_group_checksum::calculate_read_group_record_checksum;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

/// `CalculateReadGroupChecksum.doWork` for a SAM input: the MD5 hex digest of the read groups.
pub fn calculate_read_group_checksum(input_sam: &str) -> Result<String, ParseError> {
    let (header, _) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    Ok(calculate_read_group_record_checksum(&header.read_groups))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_digest_is_a_bare_32_char_hex_string() {
        let input = "@HD\tVN:1.6\n@RG\tID:rg1\tSM:s\n\
            r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
        let digest = calculate_read_group_checksum(input).unwrap();
        assert_eq!(digest.len(), 32);
        assert!(digest.chars().all(|c| c.is_ascii_hexdigit()));
        assert!(!digest.ends_with('\n'));
    }

    #[test]
    fn adding_a_read_group_changes_the_digest() {
        let one = "@HD\tVN:1.6\n@RG\tID:rg1\tSM:s\n";
        let two = "@HD\tVN:1.6\n@RG\tID:rg1\tSM:s\n@RG\tID:rg2\tSM:t\n";
        assert_ne!(
            calculate_read_group_checksum(one).unwrap(),
            calculate_read_group_checksum(two).unwrap()
        );
    }
}
