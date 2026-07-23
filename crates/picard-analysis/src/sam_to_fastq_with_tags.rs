//! `SamToFastqWithTags`, the per-tag-group FASTQ output.
//!
//! `SamToFastqWithTags` extends `SamToFastq` (Picard 3.4.0): it writes the normal read FASTQ exactly
//! as [`crate::sam_to_fastq`] does, and **in addition** writes one FASTQ per `SEQUENCE_TAG_GROUP`,
//! whose reads are built from the record's tag values rather than its bases. This module ports that
//! additional output (`writeTagRecords` and the tag-writer file naming); the base read FASTQ is the
//! already-ported `sam_to_fastq` path and is unchanged.
//!
//! For each group and each surviving read, `writeTagRecords` emits a FASTQ record:
//!   - the name is the read name, or `name/1` / `name/2` for a paired read's two ends;
//!   - the sequence is the group's sequence tags looked up on the record and joined by the group's
//!     separator (default empty);
//!   - the quality is, when no `QUALITY_TAG_GROUP` was given, `~` repeated to the sequence length;
//!     otherwise the group's quality tags joined by `~` repeated to the separator's length.
//!
//! A referenced tag that the record lacks is an error (`assertTagExists` throws).
//!
//! The output file for a group is `IOUtil.makeFileNameSafe(spec.replace(',', '_')) + ".fastq"`, e.g.
//! `SEQUENCE_TAG_GROUP=CB,UR` writes `CB_UR.fastq`. Unpaired input writes one record per read; paired
//! input (`SECOND_END_FASTQ`) writes both ends (`name/1`, `name/2`) into the **same** group file, per
//! completed pair.
//!
//! Scope of this slice: unpaired and paired reads with default options (no `OUTPUT_PER_RG`, no
//! `COMPRESS_OUTPUTS_PER_TAG_GROUP`). The per-read-group file fan-out is a separate surface; each
//! entry point asserts it was given input of the matching pairing rather than emitting a half-right
//! file.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const READ_PAIRED: u16 = 0x1;
const FIRST_OF_PAIR: u16 = 0x40;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

/// The quality character used to fill in missing qualities (`TAG_SPLIT_QUAL`).
const TAG_SPLIT_QUAL: char = '~';

/// One `SEQUENCE_TAG_GROUP` (with its matching `QUALITY_TAG_GROUP` and `TAG_GROUP_SEPERATOR`).
#[derive(Debug, Clone)]
pub struct TagGroup {
    /// The raw group spec, e.g. `"CB,UR"`; only used to name the output file.
    pub spec: String,
    /// The sequence tags, in order, whose values are concatenated into the read sequence.
    pub sequence_tags: Vec<String>,
    /// The quality tags, when `QUALITY_TAG_GROUP` was given; `None` fills quality with `~`.
    pub quality_tags: Option<Vec<String>>,
    /// The separator placed between sequence tag values (`TAG_GROUP_SEPERATOR`, default empty).
    pub separator: String,
}

impl TagGroup {
    /// A group of comma-joined `sequence` tags with default (empty) separator and `~` qualities.
    pub fn new(sequence: &str) -> Self {
        TagGroup {
            spec: sequence.to_string(),
            sequence_tags: sequence.split(',').map(|s| s.trim().to_string()).collect(),
            quality_tags: None,
            separator: String::new(),
        }
    }

    /// Set the matching `QUALITY_TAG_GROUP` (comma-joined quality tags).
    pub fn with_quality(mut self, quality: &str) -> Self {
        self.quality_tags = Some(quality.split(',').map(|s| s.trim().to_string()).collect());
        self
    }

    /// Set the `TAG_GROUP_SEPERATOR`.
    pub fn with_separator(mut self, separator: &str) -> Self {
        self.separator = separator.to_string();
        self
    }

    /// `makeFileNameSafe(spec.replace(',', '_')) + ".fastq"`.
    fn file_name(&self) -> String {
        format!(
            "{}.fastq",
            make_file_name_safe(&self.spec.replace(',', "_"))
        )
    }
}

/// A read is missing a tag that a group references (`assertTagExists`).
#[derive(Debug)]
pub struct MissingTag {
    pub read_name: String,
    pub tag: String,
}

/// `IOUtil.makeFileNameSafe`: trim, then replace each whitespace or reserved punctuation character
/// with `_`.
fn make_file_name_safe(s: &str) -> String {
    const RESERVED: &str = "!\"#$%&'()*/:;<=>?@[]\\^`{|}~";
    s.trim()
        .chars()
        .map(|c| {
            if c.is_whitespace() || RESERVED.contains(c) {
                '_'
            } else {
                c
            }
        })
        .collect()
}

/// `SAMRecord.getStringAttribute`, erroring if absent (`assertTagExists`).
fn tag_value<'a>(rec: &'a BamRecord, tag: &str) -> Result<&'a str, MissingTag> {
    let code = tag.as_bytes();
    let key = Tag::new(&[code[0], code[1]]);
    match rec.tags.get(key) {
        Some(TagValue::Str(v)) => Ok(v),
        _ => Err(MissingTag {
            read_name: rec.read_name.clone(),
            tag: tag.to_string(),
        }),
    }
}

/// A record is dropped from the FASTQ unless secondary/supplementary and vendor-fail are included;
/// this mirrors `SamToFastq.handleRecord`'s filter (defaults off), so the tag FASTQ covers exactly the
/// reads the base FASTQ does.
fn is_dropped(rec: &BamRecord) -> bool {
    rec.flags & (SECONDARY_ALIGNMENT | SUPPLEMENTARY_ALIGNMENT) != 0
        || rec.flags & READ_FAILS_VENDOR_QUALITY != 0
}

/// `writeTagRecords` for one read and one group: the FASTQ record built from the group's tags.
fn write_tag_group(
    rec: &BamRecord,
    mate_number: Option<u32>,
    group: &TagGroup,
) -> Result<String, MissingTag> {
    let header = match mate_number {
        None => rec.read_name.clone(),
        Some(m) => format!("{}/{}", rec.read_name, m),
    };

    let seq_parts: Vec<&str> = group
        .sequence_tags
        .iter()
        .map(|t| tag_value(rec, t))
        .collect::<Result<_, _>>()?;
    let sequence = seq_parts.join(&group.separator);

    let quality = match &group.quality_tags {
        None => TAG_SPLIT_QUAL.to_string().repeat(sequence.len()),
        Some(quality_tags) => {
            let qual_sep = TAG_SPLIT_QUAL.to_string().repeat(group.separator.len());
            let qual_parts: Vec<&str> = quality_tags
                .iter()
                .map(|t| tag_value(rec, t))
                .collect::<Result<_, _>>()?;
            qual_parts.join(&qual_sep)
        }
    };

    Ok(format!("@{header}\n{sequence}\n+\n{quality}\n"))
}

/// `SamToFastqWithTags`' additional output for **unpaired** reads with default options: one
/// `(file name, FASTQ text)` per `SEQUENCE_TAG_GROUP`, over the reads in file order.
///
/// Panics on a paired read, since the paired split is not part of this slice and a fragment tag file
/// for paired input would be a silent, wrong result rather than a missing feature.
pub fn sam_to_fastq_with_tags_unpaired(
    records: &[BamRecord],
    groups: &[TagGroup],
) -> Result<Vec<(String, String)>, MissingTag> {
    let mut files: Vec<(String, String)> = groups
        .iter()
        .map(|g| (g.file_name(), String::new()))
        .collect();

    for rec in records.iter().filter(|r| !is_dropped(r)) {
        assert!(
            rec.flags & READ_PAIRED == 0,
            "sam_to_fastq_with_tags_unpaired given a paired read; the paired path is not ported"
        );
        for (i, group) in groups.iter().enumerate() {
            files[i].1.push_str(&write_tag_group(rec, None, group)?);
        }
    }

    Ok(files)
}

/// `SamToFastqWithTags`' additional output for **paired** reads (with `SECOND_END_FASTQ`) with default
/// options: one `(file name, FASTQ text)` per `SEQUENCE_TAG_GROUP`. Unlike the base read FASTQ, which
/// splits the two ends into two files, each group's tag FASTQ is a **single** file carrying both ends,
/// the first-of-pair (`name/1`) then the second-of-pair (`name/2`), per completed pair.
///
/// Mirrors [`crate::sam_to_fastq::sam_to_fastq_paired`]'s pairing: reads are matched by name through a
/// first-seen map, a pair is emitted only when its second mate arrives (so ordering follows pair
/// completion), and a pair with a dropped mate is left unwritten. Panics on an unpaired read.
pub fn sam_to_fastq_with_tags_paired(
    records: &[BamRecord],
    groups: &[TagGroup],
) -> Result<Vec<(String, String)>, MissingTag> {
    use std::collections::HashMap;

    let mut files: Vec<(String, String)> = groups
        .iter()
        .map(|g| (g.file_name(), String::new()))
        .collect();
    let mut first_seen: HashMap<&str, &BamRecord> = HashMap::new();

    for rec in records.iter().filter(|r| !is_dropped(r)) {
        assert!(
            rec.flags & READ_PAIRED != 0,
            "sam_to_fastq_with_tags_paired given an unpaired read"
        );
        match first_seen.remove(rec.read_name.as_str()) {
            None => {
                first_seen.insert(rec.read_name.as_str(), rec);
            }
            Some(first) => {
                let (read1, read2) = if rec.flags & FIRST_OF_PAIR != 0 {
                    (rec, first)
                } else {
                    (first, rec)
                };
                for (i, group) in groups.iter().enumerate() {
                    files[i]
                        .1
                        .push_str(&write_tag_group(read1, Some(1), group)?);
                    files[i]
                        .1
                        .push_str(&write_tag_group(read2, Some(2), group)?);
                }
            }
        }
    }

    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::sam_file::read_sam;

    // Two unpaired reads carrying cell-barcode-style tags, matching the SamToFastqWithTags oracle.
    const INPUT: &str = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n\
        r1\t4\t*\t0\t0\t*\t*\t0\t0\tAACCGGTT\tIIIIIIII\tCR:Z:ACGT\tCY:Z:FFFF\tCB:Z:TTGG\tUR:Z:CC\tUY:Z:!!\n\
        r2\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGCCCC\tJJJJJJJJ\tCR:Z:TGCA\tCY:Z:####\tCB:Z:AAAA\tUR:Z:GG\tUY:Z:@@\n";

    fn records() -> Vec<BamRecord> {
        read_sam(INPUT).unwrap().1
    }

    #[test]
    fn a_single_tag_group_uses_the_sequence_and_quality_tags() {
        let groups = [TagGroup::new("CR").with_quality("CY")];
        let out = sam_to_fastq_with_tags_unpaired(&records(), &groups).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].0, "CR.fastq");
        assert_eq!(out[0].1, "@r1\nACGT\n+\nFFFF\n@r2\nTGCA\n+\n####\n");
    }

    #[test]
    fn a_multi_tag_group_concatenates_the_tags_and_names_the_file_with_underscores() {
        let groups = [TagGroup::new("CB,UR").with_quality("CY,UY")];
        let out = sam_to_fastq_with_tags_unpaired(&records(), &groups).unwrap();
        assert_eq!(out[0].0, "CB_UR.fastq");
        assert_eq!(out[0].1, "@r1\nTTGGCC\n+\nFFFF!!\n@r2\nAAAAGG\n+\n####@@\n");
    }

    #[test]
    fn without_a_quality_group_the_quality_is_filled_with_tildes() {
        let groups = [TagGroup::new("CB,UR")];
        let out = sam_to_fastq_with_tags_unpaired(&records(), &groups).unwrap();
        // TTGG+CC = 6 bases -> six '~'.
        assert_eq!(out[0].1, "@r1\nTTGGCC\n+\n~~~~~~\n@r2\nAAAAGG\n+\n~~~~~~\n");
    }

    #[test]
    fn multiple_groups_each_produce_a_file() {
        let groups = [
            TagGroup::new("CR").with_quality("CY"),
            TagGroup::new("CB,UR").with_quality("CY,UY"),
        ];
        let out = sam_to_fastq_with_tags_unpaired(&records(), &groups).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].0, "CR.fastq");
        assert_eq!(out[1].0, "CB_UR.fastq");
    }

    #[test]
    fn a_missing_tag_is_an_error() {
        let groups = [TagGroup::new("ZZ")];
        let err = sam_to_fastq_with_tags_unpaired(&records(), &groups).unwrap_err();
        assert_eq!(err.tag, "ZZ");
        assert_eq!(err.read_name, "r1");
    }

    #[test]
    fn a_pair_writes_both_ends_into_one_file_per_group() {
        // A first/second-of-pair template; the tag FASTQ carries both ends interleaved (name/1, name/2).
        let paired = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n\
            p1\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\tCR:Z:ACGT\tCY:Z:FFFF\n\
            p1\t141\t*\t0\t0\t*\t*\t0\t0\tCCCC\tJJJJ\tCR:Z:TTTT\tCY:Z:####\n";
        let recs = read_sam(paired).unwrap().1;
        let groups = [TagGroup::new("CR").with_quality("CY")];
        let out = sam_to_fastq_with_tags_paired(&recs, &groups).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].0, "CR.fastq");
        assert_eq!(out[0].1, "@p1/1\nACGT\n+\nFFFF\n@p1/2\nTTTT\n+\n####\n");
    }
}
