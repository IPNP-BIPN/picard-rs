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
}

impl Options {
    pub fn new(sample_name: &str) -> Self {
        Options {
            read_group_name: "A".to_string(),
            sample_name: sample_name.to_string(),
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
    let header = build_header(&opts.read_group_name, &opts.sample_name);

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

    // The writer sorts by queryname because the header is SO:queryname.
    records.sort_by(query_name::compare);
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
    use htsjdk_bam::writer::BamWriter;
    let (header, records) = fastq_to_records_unpaired(fastq_text, opts)?;
    let mut w = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
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
