//! `AddOrReplaceReadGroups`.
//!
//! Ports `picard.sam.AddOrReplaceReadGroups.doWork` at tag 3.4.0 for SAM I/O with the default sort
//! order (keep the input's): clone the input header, **replace all its read groups with one new
//! `@RG`**, stamp the `RG` tag on every record, and rewrite them in input order.
//!
//! `doWork` adds no `@PG` and no timestamp, and with `SORT_ORDER` unset the records are not
//! re-ordered (`makeWriter` is told the output is presorted), so the whole SAM is comparable raw.
//! The one new read group carries the four required fields in setter order, `LB PL SM PU`, after
//! its `ID`, then whichever optional fields were given, in `doWork`'s setter order: `CN DS DT PI PG
//! PM KS FO`. An explicit `SORT_ORDER` is a separate surface.
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
    /// The optional fields that were given, as `(tag, value)` in setter order, values already in
    /// the form htsjdk writes them (`DT` reformatted by [`iso8601`]).
    pub optional: Vec<(&'static str, String)>,
}

impl Options {
    pub fn new(rglb: &str, rgpl: &str, rgpu: &str, rgsm: &str) -> Self {
        Options {
            rgid: "1".to_string(),
            rglb: rglb.to_string(),
            rgpl: rgpl.to_string(),
            rgpu: rgpu.to_string(),
            rgsm: rgsm.to_string(),
            optional: Vec::new(),
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
    for (tag, value) in &opts.optional {
        rg.attributes.set(tag, value);
    }
    rg
}

/// `new Iso8601Date(String).toString()` in a JVM whose zone is UTC, which is the oracle's: htsjdk's
/// `DateParser` reads `YYYY-MM-DD`, optionally followed by `THH:MM[:SS]` and a zone (`Z`, `+HH:MM`
/// or `+HHMM`; none means UTC), and the date is written back as `yyyy-MM-dd'T'HH:mm:ssZ` in UTC.
/// `None` for a shape this does not read.
pub fn iso8601(text: &str) -> Option<String> {
    let num = |s: &str| -> Option<i64> {
        (!s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())).then(|| s.parse().ok())?
    };
    let (date, rest) = text.split_at(text.len().min(10));
    let mut parts = date.split('-');
    let (y, mo, d) = (
        num(parts.next()?)?,
        num(parts.next()?)?,
        num(parts.next()?)?,
    );
    let (mut secs, mut offset) = (0i64, 0i64);
    if let Some(time) = rest.strip_prefix('T') {
        let zone_at = time.find(['Z', '+', '-']).unwrap_or(time.len());
        let (clock, zone) = time.split_at(zone_at);
        let fields: Vec<i64> = clock.split(':').map(num).collect::<Option<_>>()?;
        if !(2..=3).contains(&fields.len()) {
            return None;
        }
        secs = fields[0] * 3600 + fields[1] * 60 + fields.get(2).copied().unwrap_or(0);
        if let Some(sign) = zone.chars().next().filter(|c| *c != 'Z') {
            let digits: String = zone[1..].chars().filter(|c| *c != ':').collect();
            if digits.len() != 4 {
                return None;
            }
            let minutes = num(&digits[..2])? * 60 + num(&digits[2..])?;
            offset = if sign == '-' {
                -minutes * 60
            } else {
                minutes * 60
            };
        }
    } else if !rest.is_empty() {
        return None;
    }
    // Days from the civil date (Howard Hinnant's algorithm), then back after the shift to UTC.
    let (yy, mm) = if mo <= 2 {
        (y - 1, mo + 9)
    } else {
        (y, mo - 3)
    };
    let era = yy.div_euclid(400);
    let yoe = yy - era * 400;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + (153 * mm + 2) / 5 + d - 1;
    let total = (era * 146097 + doe) * 86400 + secs - offset;
    let (days, sod) = (total.div_euclid(86400), total.rem_euclid(86400));
    let era = days.div_euclid(146097);
    let doe = days - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = era * 400 + yoe + i64::from(month <= 2);
    Some(format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}+0000",
        sod / 3600,
        sod % 3600 / 60,
        sod % 60
    ))
}

fn replace_read_groups(header: &mut SamHeader, opts: &Options) {
    header.read_groups = vec![build_read_group(opts)];
}

/// The reheader and the `RG` stamp, over records that are already in memory.
///
/// The SAM entry points below read their records first and then call this; the runnable binary
/// calls it directly, because its input is a BAM as often as it is a SAM.
///
/// Stamping the `RG` tag is per-record and independent; `par_iter_mut` keeps the order, so the
/// bytes match the serial loop (decision 0006).
pub fn replace_and_stamp(header: &mut SamHeader, records: &mut [BamRecord], opts: &Options) {
    use rayon::prelude::*;
    replace_read_groups(header, opts);
    records.par_iter_mut().for_each(|rec| {
        rec.tags
            .insert(Tag::new(b"RG"), TagValue::Str(opts.rgid.clone()));
    });
}

/// `AddOrReplaceReadGroups.doWork` up to the write: the reheadered records with the new `RG` tag.
fn apply(
    input_sam: &str,
    opts: &Options,
) -> Result<(htsjdk_bam::header::SamHeader, Vec<BamRecord>), ParseError> {
    let (mut header, mut records) = read_sam(input_sam)?;
    replace_and_stamp(&mut header, &mut records, opts);
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
