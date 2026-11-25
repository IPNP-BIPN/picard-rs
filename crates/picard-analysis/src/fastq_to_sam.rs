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
pub fn full_header(opts: &FullOptions) -> SamHeader {
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
    let header = full_header(opts);
    let mut records = fastq_records(fastq1, fastq2, opts)?.collect::<Result<Vec<_>, _>>()?;
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

/// The exception `FastqReader.error` throws for `error`, raised while reading the record that
/// starts at `line` (htsjdk counts four lines a record, skipped blank lines aside). Picard opens
/// its readers with no `File`, so the message ends `in fastq ` and a space.
pub fn fastq_exception(error: &FastqError, line: usize) -> String {
    // shortcut: htsjdk-rs's checkLine names the line kinds by their enum constants, not by
    // `LineType.toString`; drop the renaming once htsjdk-rs prints the printable names.
    let message = error
        .0
        .replace("SequenceLine", "Sequence Line")
        .replace("QualityHeader", "Quality Header")
        .replace("QualityLine", "Quality Line");
    format!("htsjdk.samtools.SAMException: {message} at line {line} in fastq ")
}

/// `htsjdk.samtools.fastq.FastqReader` as Picard drives it: the constructor and every `next()`
/// read one record ahead, so a malformed record throws when the record before it is taken.
pub struct ReadAhead<'a> {
    reader: htsjdk_bam::fastq::FastqReader<'a>,
    read: usize,
    next: Option<htsjdk_bam::fastq::FastqRecord>,
}

impl<'a> ReadAhead<'a> {
    /// The constructor, which reads the first record.
    pub fn open(text: &'a str, skip_blank_lines: bool) -> Result<Self, String> {
        let mut reader = ReadAhead {
            reader: htsjdk_bam::fastq::FastqReader::new(text, skip_blank_lines),
            read: 0,
            next: None,
        };
        reader.advance()?;
        Ok(reader)
    }

    fn advance(&mut self) -> Result<(), String> {
        let line = 1 + 4 * self.read;
        self.next = self
            .reader
            .next_record()
            .map_err(|e| fastq_exception(&e, line))?;
        self.read += usize::from(self.next.is_some());
        Ok(())
    }

    pub fn has_next(&self) -> bool {
        self.next.is_some()
    }

    /// `next()`: the record read ahead, once the one after it is read.
    pub fn next_record(&mut self) -> Result<htsjdk_bam::fastq::FastqRecord, String> {
        let record = self.next.take().expect("next() called when !hasNext()");
        self.advance()?;
        Ok(record)
    }
}

/// `QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE`.
const MAX_RECORDS_TO_ITERATE: usize = 10_000;

/// `FastqToSam.determineQualityFormat`: the quality characters of the first 10000 records, the
/// readers taken in turn, and `generateBestGuess(FASTQ, expected)` over them. A record the
/// readers throw on in that window is thrown from here.
pub fn determine_quality_format(
    readers: Vec<ReadAhead<'_>>,
    expected: Option<&str>,
) -> Result<&'static str, String> {
    let mut queue: std::collections::VecDeque<_> = readers.into();
    let mut seen = [false; 256];
    let mut count = 0;
    loop {
        while queue.front().is_some_and(|r| !r.has_next()) {
            queue.pop_front();
        }
        let Some(mut reader) = queue.pop_front() else {
            break;
        };
        if count == MAX_RECORDS_TO_ITERATE {
            break;
        }
        count += 1;
        let record = reader.next_record()?;
        for &b in record.quality_string.as_deref().unwrap_or("").as_bytes() {
            seen[usize::from(b)] = true;
        }
        queue.push_back(reader);
    }
    best_guess(&seen, expected)
}

/// `generateCandidateQualities` and `generateBestGuess` for a FASTQ, over the set of quality
/// characters seen.
fn best_guess(seen: &[bool; 256], expected: Option<&str>) -> Result<&'static str, String> {
    let exception = |m: &str| format!("htsjdk.samtools.SAMException: {m}");
    if !seen.iter().any(|&s| s) {
        return Err(exception(
            "Cannot determine candidate qualities: no qualities found.",
        ));
    }
    // `QualityScheme`, in `FastqQualityFormat`'s order: the ASCII range every quality must fall
    // in, and the range at least one must fall in when that is checked.
    let schemes: [(
        &'static str,
        std::ops::RangeInclusive<usize>,
        Option<std::ops::RangeInclusive<usize>>,
    ); 3] = [
        ("Solexa", 59..=126, None),
        ("Illumina", 64..=126, None),
        ("Standard", 33..=126, Some(33..=58)),
    ];
    let candidates = |check_expected: bool| -> Vec<&'static str> {
        schemes
            .iter()
            .filter(|(_, range, wanted)| {
                (0..256).all(|q| !seen[q] || range.contains(&q))
                    && !(check_expected
                        && wanted.as_ref().is_some_and(|w| !w.clone().any(|q| seen[q])))
            })
            .map(|(name, ..)| *name)
            .collect()
    };
    if let Some(expected) = expected {
        return candidates(false)
            .into_iter()
            .find(|c| *c == expected)
            .ok_or_else(|| {
                exception(&format!(
                    "The quality values do not fall in the range appropriate for the expected \
                     quality of {expected}."
                ))
            });
    }
    match candidates(true).as_slice() {
        [one] => Ok(one),
        // {Solexa, Illumina} and, for a FASTQ, {Illumina, Standard} are Illumina.
        ["Solexa", "Illumina"] | ["Illumina", "Standard"] => Ok("Illumina"),
        ["Solexa", "Standard"] => Ok("Standard"),
        [] => Err(exception(
            "The quality format cannot be determined: all formats were excluded.",
        )),
        _ => Err(exception(
            "The quality format cannot be determined: no formats were excluded.",
        )),
    }
}

/// `IOUtil.openFileForBufferedReading`: a `.gz`, `.gzip` or `.bfq` name is gunzipped.
pub fn read_fastq(path: &str) -> std::io::Result<String> {
    use std::io::Read;
    let file = std::fs::File::open(path)?;
    let mut text = String::new();
    if [".gz", ".gzip", ".bfq"].iter().any(|e| path.ends_with(e)) {
        flate2::read::MultiGzDecoder::new(file).read_to_string(&mut text)?;
    } else {
        std::io::BufReader::new(file).read_to_string(&mut text)?;
    }
    Ok(text)
}

/// [`fastq_to_bam_full`]'s records in the order the FASTQs give them, built one at a time. Each
/// file is parsed through first, so a parse failure anywhere outranks a record's, as reading them
/// whole did.
pub fn fastq_records<'a>(
    fastq1: &'a str,
    fastq2: Option<&'a str>,
    opts: &'a FullOptions,
) -> Result<impl Iterator<Item = Result<BamRecord, String>> + 'a, String> {
    use htsjdk_bam::fastq::{get_sam_read_name_from_fastq_header, FastqReader};

    let reader = move |text: &'a str| FastqReader::new(text, opts.allow_and_ignore_empty_lines);
    let count = |text: &'a str| -> Result<usize, String> {
        let mut reader = reader(text);
        let mut count = 0;
        while reader
            .next_record()
            .map_err(|e| fastq_exception(&e, 1 + 4 * count))?
            .is_some()
        {
            count += 1;
        }
        Ok(count)
    };
    let count1 = count(fastq1)?;
    let count2 = fastq2.map(count).transpose()?;
    let records = move |text: &'a str| {
        let mut reader = reader(text);
        std::iter::from_fn(move || reader.next_record().expect("parsed once already"))
    };
    let make = move |frec: &htsjdk_bam::fastq::FastqRecord,
                     name: String|
          -> Result<BamRecord, String> {
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
    let all: Box<dyn Iterator<Item = Result<BamRecord, String>> + 'a> = match fastq2 {
        None => Box::new(records(fastq1).map(move |frec| {
            let name = get_sam_read_name_from_fastq_header(frec.read_name.as_deref().unwrap_or(""));
            make(&frec, name)
        })),
        Some(text2) => {
            let pairs = records(fastq1).zip(records(text2)).flat_map(move |(frec1, frec2)| {
                let pair = (|| {
                    let name1 = get_sam_read_name_from_fastq_header(
                        frec1.read_name.as_deref().unwrap_or(""),
                    );
                    let name2 = get_sam_read_name_from_fastq_header(
                        frec2.read_name.as_deref().unwrap_or(""),
                    );
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
                    let mut first = make(&frec1, base1.clone())?;
                    first.flags |= 0x1 | 0x8 | 0x40;
                    let mut second = make(&frec2, base1)?;
                    second.flags |= 0x1 | 0x8 | 0x80;
                    Ok([first, second])
                })();
                match pair {
                    Ok(both) => both.map(Ok).to_vec(),
                    Err(e) => vec![Err(e)],
                }
            });
            // The lengths are compared once every pair the shorter file makes has been built.
            let unequal = (Some(count1) != count2).then(|| {
                Err(
                    "picard.PicardException: Input paired fastq files must be the same length"
                        .to_string(),
                )
            });
            Box::new(pairs.chain(unequal))
        }
    };
    Ok(all)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn guess(qualities: &str, expected: Option<&str>) -> Result<&'static str, String> {
        let text = format!("@r\n{}\n+\n{qualities}\n", "A".repeat(qualities.len()));
        determine_quality_format(vec![ReadAhead::open(&text, false).unwrap()], expected)
    }

    #[test]
    fn the_detector_guesses_as_generate_best_guess_does() {
        // '#' (35) is in Standard's expected range 33..=58 and below Solexa's and Illumina's.
        assert_eq!(guess("#I", None), Ok("Standard"));
        // 'h' (104) fits all three, and nothing falls in 33..=58: Solexa and Illumina, so Illumina.
        assert_eq!(guess("h", None), Ok("Illumina"));
        // ';' (59) excludes Illumina and is not in 33..=58: Solexa alone.
        assert_eq!(guess(";h", None), Ok("Solexa"));
        // A claim is checked against the ranges only.
        assert_eq!(guess("h", Some("Standard")), Ok("Standard"));
        assert!(guess("#", Some("Illumina"))
            .unwrap_err()
            .ends_with("the expected quality of Illumina."));
        assert!(guess("\u{7f}", None)
            .unwrap_err()
            .ends_with("all formats were excluded."));
    }

    #[test]
    fn a_malformed_record_throws_with_its_line_when_the_one_before_is_taken() {
        let text = "@a\nA\n+\nI\n@b\nA\n+\n";
        let mut reader = ReadAhead::open(text, false).unwrap();
        assert_eq!(
            reader.next_record().unwrap_err(),
            "htsjdk.samtools.SAMException: File is too short - missing Quality Line at line 5 \
             in fastq "
        );
    }

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
