//! `FastqToSam`, the unpaired default path.
//!
//! Ports `picard.sam.FastqToSam.createSamFileHeader`, `createSamRecord`, and `doUnpaired` at tag
//! 3.4.0, for **a single FASTQ input with default options** written as SAM. The counterpart to
//! `SamToFastq`, and the other half of the FASTQ round trip at the tool level.
//!
//! Two things make the whole SAM output comparable raw, with no canonicalization. First,
//! `createSamFileHeader` builds the header by hand with only an `@HD` and one `@RG`, and adds **no
//! `@PG`**, so there is no command line to strip. Second, the header carries no timestamp. What is
//! left is exact: `@HD VN:1.6 SO:queryname`, then `@RG ID:<rg> SM:<sample>`, then the reads.
//!
//! Each read becomes an **unmapped** record whose name is the FASTQ header cleaned by
//! `getSamReadNameFromFastqHeader`, whose bases and Standard-format qualities are decoded from the
//! record, and which carries the `RG` tag. The reads are then **queryname-sorted** (the default
//! `SORT_ORDER`), which is why this waited on `SAMRecordQueryNameComparator`.
//!
//! Only the unpaired, Standard-quality path is claimed. Paired input (two FASTQs), the Solexa and
//! Illumina quality formats, and the other read-group and sort-order options are separate surfaces.

use htsjdk_bam::fastq::{as_sam_record, parse_fastq, FastqError};
use htsjdk_bam::header::{ReadGroup, SamHeader};
use htsjdk_bam::query_name;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use htsjdk_bam::tag::{Tag, TagValue};

/// The options that gate the unpaired path.
#[derive(Debug, Clone)]
pub struct Options {
    /// `READ_GROUP_NAME`, default `"A"`.
    pub read_group_name: String,
    /// `SAMPLE_NAME` (required, no default).
    pub sample_name: String,
    /// `SORT_ORDER`, `queryname` by default. It is written into the header AND it decides the
    /// order: the writer is given the header's order, so asking for coordinate on a file of
    /// unmapped reads keeps them in the order the FASTQ had them.
    pub sort_order: String,
}

impl Options {
    pub fn new(sample_name: &str) -> Self {
        Options {
            read_group_name: "A".to_string(),
            sample_name: sample_name.to_string(),
            sort_order: "queryname".to_string(),
        }
    }
}

/// `createSamFileHeader`: an `@HD` with `SO:queryname` and a single `@RG`, nothing else.
fn build_header(read_group_name: &str, sample_name: &str) -> SamHeader {
    let mut header = SamHeader::new();
    header.set_sort_order("queryname");
    let mut rg = ReadGroup::new(read_group_name);
    rg.attributes.set("SM", sample_name);
    header.read_groups.push(rg);
    header
}

/// `doUnpaired`: convert every FASTQ read to an unmapped record, queryname-sort, and write SAM.
///
/// Blank-line skipping is off (`ALLOW_AND_IGNORE_EMPTY_LINES` defaults false), matching the reader
/// FastqToSam constructs.
pub fn fastq_to_sam_unpaired(fastq_text: &str, opts: &Options) -> Result<String, FastqError> {
    let (header, records) = fastq_to_records_unpaired(fastq_text, opts)?;
    Ok(write_sam(&header, &records).expect("unmapped records always encode as SAM text"))
}

/// `doUnpaired` up to the write: the header and the queryname-sorted unmapped records. Shared by the
/// SAM and BAM renderers so they cannot drift.
fn fastq_to_records_unpaired(
    fastq_text: &str,
    opts: &Options,
) -> Result<(SamHeader, Vec<BamRecord>), FastqError> {
    let mut header = build_header(&opts.read_group_name, &opts.sample_name);
    header.set_sort_order(&opts.sort_order);

    use rayon::prelude::*;
    // Converting each FASTQ record to an unmapped SAM record is independent; the parallel map's
    // `collect` preserves order, so the pre-sort record list is identical to the serial one
    // (decision 0006). The queryname sort below then imposes the final order.
    let mut records: Vec<BamRecord> = parse_fastq(fastq_text, false)?
        .par_iter()
        .map(|frec| {
            // createSamRecord: asSAMRecord gives the unmapped read (name cleaned, bases, quals);
            // FastqToSam then tags it with the read group.
            let mut rec = as_sam_record(frec);
            rec.tags
                .insert(Tag::new(b"RG"), TagValue::Str(opts.read_group_name.clone()));
            rec
        })
        .collect();

    // The writer sorts by whatever the header says. Every record here is unmapped, so a
    // coordinate sort has nothing to order by and the FASTQ's own order survives.
    if opts.sort_order == "queryname" {
        records.sort_by(query_name::compare);
    }
    Ok((header, records))
}

/// `FastqToSam` `doUnpaired` for **BAM** output: the same unmapped, queryname-sorted records written
/// through htsjdk-rs's byte-identical `BamWriter`. FastqToSam adds no `@PG`, so byte-identity to
/// Picard with `USE_JDK_DEFLATER=true` follows transitively (the records are those the SAM path
/// already reproduces, and `BamWriter` is oracle-gated over arbitrary records).
pub fn fastq_to_sam_unpaired_to_bam(
    fastq_text: &str,
    opts: &Options,
) -> Result<Vec<u8>, FastqError> {
    let (header, records) = fastq_to_records_unpaired(fastq_text, opts)?;
    let mut w =
        crate::output::bam_writer(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        w.write(rec).expect("unmapped record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

/// Every `FastqToSam` argument that decides an output byte, for [`fastq_to_bam_full`].
#[derive(Debug, Clone)]
pub struct FullOptions {
    pub read_group_name: String,
    pub sample_name: String,
    pub sort_order: String,
    pub library_name: Option<String>,
    pub platform: Option<String>,
    pub platform_unit: Option<String>,
    pub sequencing_center: Option<String>,
    pub predicted_insert_size: Option<i32>,
    pub description: Option<String>,
    /// `RUN_DATE`, already in the form `Iso8601Date.toString` writes.
    pub run_date: Option<String>,
    pub platform_model: Option<String>,
    pub program_group: Option<String>,
    pub comments: Vec<String>,
    pub min_q: i32,
    pub max_q: i32,
    /// `Standard`, `Illumina` or `Solexa`: the encoding the qualities are read in.
    pub quality_format: String,
    pub allow_and_ignore_empty_lines: bool,
}

/// `createSamFileHeader`: the read group's fields in setter order, then the comments.
fn full_header(opts: &FullOptions) -> SamHeader {
    let mut header = build_header(&opts.read_group_name, &opts.sample_name);
    let group = &mut header.read_groups[0];
    let optional = [
        ("LB", opts.library_name.clone()),
        ("PL", opts.platform.clone()),
        ("PU", opts.platform_unit.clone()),
        ("CN", opts.sequencing_center.clone()),
        ("PI", opts.predicted_insert_size.map(|v| v.to_string())),
        ("DS", opts.description.clone()),
        ("DT", opts.run_date.clone()),
        ("PM", opts.platform_model.clone()),
        ("PG", opts.program_group.clone()),
    ];
    for (tag, value) in optional {
        if let Some(value) = value {
            group.attributes.set(tag, &value);
        }
    }
    for comment in &opts.comments {
        header.comments.push(if comment.starts_with("@CO\t") {
            comment.clone()
        } else {
            format!("@CO\t{comment}")
        });
    }
    header.set_sort_order(&opts.sort_order);
    header
}

/// `convertQuality`: the quality characters as Phred values, in the encoding they were read in.
fn convert_quality(text: &str, format: &str) -> Vec<u8> {
    text.bytes()
        .map(|q| match format {
            "Illumina" => q.wrapping_sub(64),
            "Solexa" => {
                // `SolexaQualityConverter`: below the addend a quality is zero, above it the
                // log-odds conversion, rounded.
                if q < 64 {
                    0
                } else {
                    let solexa = (q as f64 - 64.0) / 10.0;
                    (10.0 * (1.0 + 10f64.powf(solexa)).log10()).round() as u8
                }
            }
            _ => q.wrapping_sub(33),
        })
        .collect()
}

/// `getReadNameTokens`: the name and, when it ends in `/1` or `/2`, that number.
fn read_name_tokens(name: &str) -> (String, Option<String>) {
    match name.rfind('/') {
        Some(idx) if matches!(&name[idx + 1..], "1" | "2") => {
            (name[..idx].to_string(), Some(name[idx + 1..].to_string()))
        }
        _ => (name.to_string(), None),
    }
}

/// `doUnpaired` or `doPaired`, from the FASTQ text to the BAM the writer produces. An `Err` is a
/// Java throwable, class and message.
pub fn fastq_to_bam_full(
    fastq1: &str,
    fastq2: Option<&str>,
    opts: &FullOptions,
) -> Result<Vec<u8>, String> {
    use htsjdk_bam::fastq::get_sam_read_name_from_fastq_header;

    let header = full_header(opts);
    let parse = |text: &str| {
        parse_fastq(text, opts.allow_and_ignore_empty_lines)
            .map_err(|e| format!("htsjdk.samtools.SAMException: {e}"))
    };
    let reads1 = parse(fastq1)?;
    let reads2 = match fastq2 {
        Some(text) => Some(parse(text)?),
        None => None,
    };
    let make = |frec: &htsjdk_bam::fastq::FastqRecord, name: String| -> Result<BamRecord, String> {
        let quals = convert_quality(
            frec.quality_string.as_deref().unwrap_or(""),
            &opts.quality_format,
        );
        for &q in &quals {
            let q = q as i32;
            if q < opts.min_q || q > opts.max_q {
                return Err(format!(
                    "picard.PicardException: Base quality {q} is not in the range {}..{} for read {}",
                    opts.min_q,
                    opts.max_q,
                    frec.read_name.as_deref().unwrap_or("")
                ));
            }
        }
        let mut rec = as_sam_record(frec);
        rec.read_name = name;
        rec.base_qualities = quals;
        rec.tags
            .insert(Tag::new(b"RG"), TagValue::Str(opts.read_group_name.clone()));
        Ok(rec)
    };
    let mut records = Vec::new();
    match reads2 {
        None => {
            for frec in &reads1 {
                let name =
                    get_sam_read_name_from_fastq_header(frec.read_name.as_deref().unwrap_or(""));
                records.push(make(frec, name)?);
            }
        }
        Some(reads2) => {
            for (frec1, frec2) in reads1.iter().zip(reads2.iter()) {
                let name1 =
                    get_sam_read_name_from_fastq_header(frec1.read_name.as_deref().unwrap_or(""));
                let name2 =
                    get_sam_read_name_from_fastq_header(frec2.read_name.as_deref().unwrap_or(""));
                if name1.is_empty() || name2.is_empty() {
                    return Err(
                        "picard.PicardException: Pair read name cannot be empty".to_string()
                    );
                }
                let (base1, num1) = read_name_tokens(&name1);
                let (base2, num2) = read_name_tokens(&name2);
                if base1 != base2 {
                    return Err(format!(
                        "picard.PicardException: In paired mode, read name 1 ({base1}) does not match read name 2 ({base2})"
                    ));
                }
                if num1.is_some() != num2.is_some()
                    || num1.as_deref().is_some_and(|n| n != "1")
                    || num2.as_deref().is_some_and(|n| n != "2")
                {
                    return Err(format!(
                        "picard.PicardException: the pair numbers of {name1} and {name2} do not read 1 and 2"
                    ));
                }
                let mut first = make(frec1, base1.clone())?;
                first.flags |= 0x1 | 0x8 | 0x40;
                records.push(first);
                let mut second = make(frec2, base1)?;
                second.flags |= 0x1 | 0x8 | 0x80;
                records.push(second);
            }
            if reads1.len() != reads2.len() {
                return Err(
                    "picard.PicardException: Input paired fastq files must be the same length"
                        .to_string(),
                );
            }
        }
    }
    if opts.sort_order == "queryname" {
        records.sort_by(query_name::compare);
    }
    let mut w =
        crate::output::bam_writer(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        w.write(rec).expect("unmapped record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_header_has_no_program_record_and_is_queryname_sorted() {
        let text = fastq_to_sam_unpaired("@r1\nACGT\n+\nIIII\n", &Options::new("s1")).unwrap();
        assert!(text.contains("@HD\tVN:1.6\tSO:queryname"), "got {text}");
        assert!(text.contains("@RG\tID:A\tSM:s1"), "got {text}");
        assert!(!text.contains("@PG"), "FastqToSam writes no @PG: {text}");
    }

    #[test]
    fn the_bam_output_round_trips_to_the_sam_output() {
        use htsjdk_bam::reader::BamReader;
        let input = "@r2\nAA\n+\nII\n@r1\nCC\n+\nII\n@r10\nGG\n+\nII\n";
        let opts = Options::new("s1");
        let sam = fastq_to_sam_unpaired(input, &opts).unwrap();
        let bam = fastq_to_sam_unpaired_to_bam(input, &opts).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bam).unwrap();
        let reader = BamReader::new(&plain).unwrap();
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(write_sam(&header, &records).unwrap(), sam);
    }

    #[test]
    fn reads_are_emitted_in_queryname_order_not_input_order() {
        // Input order r2, r1, r10; queryname (String.compareTo) order is r1, r10, r2.
        let input = "@r2\nAA\n+\nII\n@r1\nCC\n+\nII\n@r10\nGG\n+\nII\n";
        let text = fastq_to_sam_unpaired(input, &Options::new("s1")).unwrap();
        let names: Vec<&str> = text
            .lines()
            .filter(|l| !l.starts_with('@'))
            .filter_map(|l| l.split('\t').next())
            .collect();
        assert_eq!(names, ["r1", "r10", "r2"]);
    }

    #[test]
    fn a_read_is_unmapped_with_its_rg_tag_and_decoded_quals() {
        let text = fastq_to_sam_unpaired("@r/1\nACGT\n+\nII?5\n", &Options::new("s1")).unwrap();
        let row = text.lines().find(|l| !l.starts_with('@')).unwrap();
        let f: Vec<&str> = row.split('\t').collect();
        assert_eq!(f[0], "r"); // /1 suffix stripped by the read-name cleanup
        assert_eq!(f[1], "4"); // unmapped flag, not paired
        assert_eq!(f[9], "ACGT"); // SEQ
        assert_eq!(f[10], "II?5"); // QUAL round-trips
        assert!(row.contains("RG:Z:A"));
    }
}
