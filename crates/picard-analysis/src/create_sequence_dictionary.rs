//! `CreateSequenceDictionary`.
//!
//! Ports `picard.sam.CreateSequenceDictionary.makeSequenceDictionary` at tag 3.4.0, for the default
//! path: read a reference FASTA and write a `.dict` (a headerless SAM header) with one `@SQ` line per
//! contig. Opens the **reference** archetype.
//!
//! Each `@SQ` line carries the contig name (`SN`, the FASTA header up to the first whitespace), its
//! length (`LN`), an `M5` checksum, and the reference URI (`UR`). The `M5` is
//! `SequenceUtil.calculateMD5`: the MD5 of the contig's bases **uppercased**, as lower-case hex. The
//! `UR` is the reference's `file:` URI, which is path-dependent and therefore canonicalized away in
//! comparison; everything else (`@HD VN:1.6`, and each `SN`/`LN`/`M5`) is exact. No `@PG` and no
//! timestamp are written.
//!
//! # The optional attributes, and the order they are written in
//!
//! `makeSequenceRecord` sets `M5`, then `AS` when `GENOME_ASSEMBLY` was given, then `UR`
//! unconditionally, then `SP` when `SPECIES` was given. A `SAMSequenceRecord`'s attributes keep
//! insertion order (htsjdk-rs decision 0009), so that IS the order they appear in, and a port that
//! sorted them or that wrote `AS` after `UR` would produce a different file for the same inputs.
//!
//! `URI` overrides the reference's own `file:` URI when it is given, which is why the `UR` value is
//! a parameter here rather than derived: the tool has an argument for it.
//!
//! `NUM_SEQUENCES` truncates the dictionary to its first N contigs.
//!
//! `TRUNCATE_NAMES_AT_WHITESPACE=false` is NOT ported, and the reason is one layer down:
//! `ReferenceSequenceFileFactory.getReferenceSequenceFile` takes the flag and the reader keeps the
//! whole header when it is false, while `htsjdk_bam::fasta::ReferenceSequence` truncates and drops
//! the description. Supporting it here would mean inventing the untruncated name from something
//! this port no longer has. The array holds the argument at Picard's default until the reader
//! carries the description. `ALT_NAMES` and the alias handling are a separate surface.

use htsjdk_bam::fasta::{read_fasta, FastaError};
use md5::{Digest, Md5};

/// `SequenceUtil.calculateMD5`: the MD5 of the bases uppercased, rendered as lower-case hex.
fn calculate_md5(bases: &[u8]) -> String {
    let mut digest = Md5::new();
    // Uppercase byte by byte, as htsjdk does, so a soft-masked (lower-case) reference hashes the same
    // as its upper-case form.
    let upper: Vec<u8> = bases.iter().map(u8::to_ascii_uppercase).collect();
    digest.update(&upper);
    let mut hex = String::with_capacity(32);
    for b in digest.finalize() {
        hex.push_str(&format!("{b:02x}"));
    }
    hex
}

/// `CreateSequenceDictionary.makeSequenceDictionary` for a FASTA input. `fasta` is the reference
/// bytes; `reference_uri` is the `UR` value (its `file:` URI).
pub fn create_sequence_dictionary(fasta: &[u8], reference_uri: &str) -> Result<String, FastaError> {
    create_sequence_dictionary_with(fasta, reference_uri, &Options::default())
}

/// `CreateSequenceDictionary`'s optional arguments, defaulting to Picard's defaults.
#[derive(Debug, Clone, Default)]
pub struct Options {
    /// `GENOME_ASSEMBLY`, written as `AS` when present.
    pub genome_assembly: Option<String>,
    /// `SPECIES`, written as `SP` when present.
    pub species: Option<String>,
    /// `NUM_SEQUENCES`: keep only the first N contigs. Picard's default is `Integer.MAX_VALUE`.
    pub num_sequences: Option<usize>,
}

/// `makeSequenceDictionary` with the optional attributes the tool accepts.
pub fn create_sequence_dictionary_with(
    fasta: &[u8],
    reference_uri: &str,
    options: &Options,
) -> Result<String, FastaError> {
    let contigs = read_fasta(fasta)?;
    let mut out = String::from("@HD\tVN:1.6\n");
    let keep = options.num_sequences.unwrap_or(usize::MAX);
    for contig in contigs.iter().take(keep) {
        // `SamSequenceRecordsIterator.makeSequenceRecord`'s order: M5, then AS when given, then
        // UR, then SP when given. The attributes keep insertion order, so this is the file's.
        out.push_str(&format!(
            "@SQ\tSN:{}\tLN:{}\tM5:{}",
            contig.name,
            contig.bases.len(),
            calculate_md5(&contig.bases),
        ));
        if let Some(assembly) = &options.genome_assembly {
            out.push_str(&format!("\tAS:{assembly}"));
        }
        out.push_str(&format!("\tUR:{reference_uri}"));
        if let Some(species) = &options.species {
            out.push_str(&format!("\tSP:{species}"));
        }
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA: &[u8] = b">chr1 desc\nACGTacgtNNNN\nACGT\n>chr2\nTTTTGGGGCCCC\n";

    fn sq(dict: &str, sn: &str) -> String {
        dict.lines()
            .find(|l| l.contains(&format!("SN:{sn}\t")) || l.contains(&format!("SN:{sn}")))
            .unwrap()
            .to_string()
    }

    #[test]
    fn the_header_leads_with_hd_and_no_pg() {
        let dict = create_sequence_dictionary(FASTA, "file:///ref.fasta").unwrap();
        assert!(dict.starts_with("@HD\tVN:1.6\n"));
        assert!(!dict.contains("@PG"));
    }

    #[test]
    fn a_contig_reports_its_name_length_and_uppercased_md5() {
        let dict = create_sequence_dictionary(FASTA, "file:///ref.fasta").unwrap();
        // chr1 = ACGTacgtNNNN + ACGT = 16 bases; the description after the space is dropped from SN.
        let line = sq(&dict, "chr1");
        assert!(line.contains("SN:chr1\t"), "got {line}");
        assert!(line.contains("LN:16\t"), "got {line}");
        // M5 is the MD5 of the uppercased bases ACGTACGTNNNNACGT.
        assert!(
            line.contains("M5:6fe448fd331944d9de0964282ec020dd\t"),
            "got {line}"
        );
        assert!(line.ends_with("UR:file:///ref.fasta"), "got {line}");
    }

    #[test]
    fn the_md5_uppercases_so_soft_masking_does_not_change_it() {
        assert_eq!(calculate_md5(b"acgt"), calculate_md5(b"ACGT"));
    }
}
