//! `AddOrReplaceReadGroups`.
//!
//! Ports `picard.sam.AddOrReplaceReadGroups.doWork` at tag 3.4.0 for SAM I/O with the default sort
//! order (keep the input's): clone the input header, **replace all its read groups with one new
//! `@RG`**, stamp the `RG` tag on every record, and rewrite them in input order.
//!
//! `doWork` adds no `@PG` and no timestamp, and with `SORT_ORDER` unset the records are not
//! re-ordered (`makeWriter` is told the output is presorted), so the whole SAM is comparable raw.
//! The one new read group carries the four required fields in setter order, `LB PL SM PU`, after
//! its `ID`; the optional fields (CN, DS, DT, PI, PG, PM, KS, FO) and an explicit `SORT_ORDER` are
//! separate surfaces.
//!
//! Replacing the `RG` tag is a sorted-insert that overwrites the existing tag in place, so a record
//! that already carried an `RG` keeps its tag position and only its value changes.

use htsjdk_bam::header::{ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::ParseError;

/// The read-group fields: `RGID` (default `"1"`) and the four required attributes.
#[derive(Debug, Clone)]
pub struct Options {
    pub rgid: String,
    pub rglb: String,
    pub rgpl: String,
    pub rgpu: String,
    pub rgsm: String,
}

impl Options {
    pub fn new(rglb: &str, rgpl: &str, rgpu: &str, rgsm: &str) -> Self {
        Options {
            rgid: "1".to_string(),
            rglb: rglb.to_string(),
            rgpl: rgpl.to_string(),
            rgpu: rgpu.to_string(),
            rgsm: rgsm.to_string(),
        }
    }
}

/// Builds the single replacement `@RG`. The attribute order is the setter order htsjdk uses:
/// `setLibrary`, `setPlatform`, `setSample`, `setPlatformUnit` -> `LB PL SM PU`.
fn build_read_group(opts: &Options) -> ReadGroup {
    let mut rg = ReadGroup::new(&opts.rgid);
    rg.attributes.set("LB", &opts.rglb);
    rg.attributes.set("PL", &opts.rgpl);
    rg.attributes.set("SM", &opts.rgsm);
    rg.attributes.set("PU", &opts.rgpu);
    rg
}

fn replace_read_groups(header: &mut SamHeader, opts: &Options) {
    header.read_groups = vec![build_read_group(opts)];
}

/// `AddOrReplaceReadGroups.doWork` up to the write: the reheadered records with the new `RG` tag.
fn apply(
    input_sam: &str,
    opts: &Options,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), ParseError> {
    use rayon::prelude::*;
    let (mut header, mut records) = read_sam(input_sam)?;
    replace_read_groups(&mut header, opts);
    // Stamping the RG tag is per-record and independent; par_iter_mut keeps the order, so the bytes
    // match the serial loop (decision 0006).
    records.par_iter_mut().for_each(|rec| {
        rec.tags
            .insert(Tag::new(b"RG"), TagValue::Str(opts.rgid.clone()));
    });
    Ok((header, records))
}

/// `AddOrReplaceReadGroups.doWork` for SAM I/O, default sort order.
pub fn add_or_replace_read_groups(input_sam: &str, opts: &Options) -> Result<String, ParseError> {
    let (header, records) = apply(input_sam, opts)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// `AddOrReplaceReadGroups.doWork` for SAM input and **BAM** output. Same reheader + RG stamp, written
/// through htsjdk-rs's byte-identical `BamWriter`; the tool adds no `@PG`. Byte-identity to Picard's
/// `USE_JDK_DEFLATER=true` output follows transitively: the records are the ones
/// `add_or_replace_read_groups` already reproduces (its oracle), and the `BamWriter` is proven
/// byte-identical over arbitrary records (the SamFormatConverter oracle and htsjdk-rs's
/// `every_file_is_byte_identical_to_htsjdks`).
pub fn add_or_replace_read_groups_to_bam(
    input_sam: &str,
    opts: &Options,
) -> Result<Vec<u8>, ParseError> {
    use htsjdk_bam::writer::BamWriter;
    let (header, records) = apply(input_sam, opts)?;
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

    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @RG\tID:old\tSM:otherSample\n\
        r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n\
        r2\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n";

    fn opts() -> Options {
        Options {
            rgid: "2".to_string(),
            ..Options::new("lib1", "ILLUMINA", "unit1", "sample1")
        }
    }

    #[test]
    fn the_single_new_read_group_replaces_the_old_one() {
        let out = add_or_replace_read_groups(INPUT, &opts()).unwrap();
        // Exactly one @RG, in setter order, and the old one is gone.
        assert!(out.contains("@RG\tID:2\tLB:lib1\tPL:ILLUMINA\tSM:sample1\tPU:unit1"));
        assert!(!out.contains("ID:old"));
        assert_eq!(out.matches("@RG").count(), 1);
    }

    #[test]
    fn every_record_gets_the_new_rg_tag() {
        let out = add_or_replace_read_groups(INPUT, &opts()).unwrap();
        for row in out.lines().filter(|l| !l.starts_with('@')) {
            assert!(row.ends_with("RG:Z:2"), "row lost its RG tag: {row}");
        }
        assert!(!out.contains("RG:Z:old"));
    }

    #[test]
    fn the_records_keep_their_input_order_and_the_header_is_otherwise_preserved() {
        let out = add_or_replace_read_groups(INPUT, &opts()).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:coordinate")); // sort order unchanged
        assert!(out.contains("@SQ\tSN:chr1\tLN:1000"));
        let names: Vec<&str> = out
            .lines()
            .filter(|l| !l.starts_with('@'))
            .filter_map(|l| l.split('\t').next())
            .collect();
        assert_eq!(names, ["r1", "r2"]);
    }

    /// The BAM output decodes back to exactly the SAM output; the writer's byte-identity to htsjdk is
    /// proven elsewhere.
    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        let sam = add_or_replace_read_groups(INPUT, &opts()).unwrap();
        let bam = add_or_replace_read_groups_to_bam(INPUT, &opts()).unwrap();
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
