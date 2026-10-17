//! `SamToFastq.doWork` and `SamToFastqWithTags`, whole.
//!
//! Ports `picard.sam.SamToFastq` and `picard.sam.SamToFastqWithTags` at tag 3.4.0: the writers
//! `generateWriters` builds (one set for every read group, or one set per read group under
//! `OUTPUT_PER_RG`), `handleRecord`'s filters and mate matching, `writeRecord`'s clipping,
//! re-reversal, trimming, quality trimming and length cap, and the tag-group files the subclass
//! adds beside them. [`crate::sam_to_fastq`] keeps the default paths it started with; this is the
//! one the runnable binaries drive.
//!
//! A writer is a file the moment it is built, so a FASTQ the run never writes a read to is still
//! left behind, empty, exactly as the reference's is. The one writer built lazily is a read group's
//! second-of-pair file, which exists only once a pair reaches it.
//!
//! An `Err` is a Java throwable, class and message, as the reference's uncaught one prints it.

use htsjdk_bam::fastq::{phred_to_fastq, write_record, FastqRecord};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::reverse_complement;
use htsjdk_bam::tag::{Tag, TagValue};

const READ_PAIRED: u16 = 0x1;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const READ_FAILS_VENDOR_QUALITY: u16 = 0x200;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

/// `SamToFastqWithTags`' arguments: one entry per `SEQUENCE_TAG_GROUP`.
#[derive(Debug, Clone, Default)]
pub struct TagGroups {
    pub sequence: Vec<String>,
    pub quality: Vec<String>,
    pub separator: Vec<String>,
    pub compress: bool,
}

/// Every `SamToFastq` argument that decides an output byte.
#[derive(Debug, Clone)]
pub struct Options {
    pub input: String,
    pub fastq: Option<String>,
    pub second_end_fastq: Option<String>,
    pub unpaired_fastq: Option<String>,
    pub output_per_rg: bool,
    pub compress_outputs_per_rg: bool,
    pub rg_tag: String,
    pub output_dir: Option<String>,
    pub re_reverse: bool,
    pub interleave: bool,
    pub include_non_pf_reads: bool,
    pub include_non_primary_alignments: bool,
    pub clipping_attribute: Option<String>,
    pub clipping_action: Option<String>,
    pub clipping_min_length: i32,
    pub read1_trim: i32,
    pub read1_max_bases_to_write: Option<i32>,
    pub read2_trim: i32,
    pub read2_max_bases_to_write: Option<i32>,
    pub quality: Option<i32>,
    /// `None` for `SamToFastq`, the subclass's groups for `SamToFastqWithTags`.
    pub tags: Option<TagGroups>,
}

impl Options {
    /// Picard's defaults, with the input named.
    pub fn new(input: &str) -> Self {
        Options {
            input: input.to_string(),
            fastq: None,
            second_end_fastq: None,
            unpaired_fastq: None,
            output_per_rg: false,
            compress_outputs_per_rg: false,
            rg_tag: "PU".to_string(),
            output_dir: None,
            re_reverse: true,
            interleave: false,
            include_non_pf_reads: false,
            include_non_primary_alignments: false,
            clipping_attribute: None,
            clipping_action: None,
            clipping_min_length: 0,
            read1_trim: 0,
            read1_max_bases_to_write: None,
            read2_trim: 0,
            read2_max_bases_to_write: None,
            quality: None,
            tags: None,
        }
    }

    /// `customCommandLineValidation`: the subclass's own failures when it has any, else the base
    /// tool's, every one of them, in the order the reference collects them.
    pub fn validate(&self) -> Vec<String> {
        if let Some(tags) = &self.tags {
            let mut errors = Vec::new();
            if !tags.quality.is_empty() && tags.sequence.len() != tags.quality.len() {
                errors.push(
                    "QUALITY_TAG_GROUP size must be equal to SEQUENCE_TAG_GROUP or not be specified at all."
                        .to_string(),
                );
            }
            if !tags.separator.is_empty() && tags.sequence.len() != tags.separator.len() {
                errors.push(
                    "TAG_GROUP_SEPERATOR size must be equal to SEQUENCE_TAG_GROUP or not be specified at all."
                        .to_string(),
                );
            }
            if !errors.is_empty() {
                return errors;
            }
        }
        let mut errors = Vec::new();
        if self.interleave && self.second_end_fastq.is_some() {
            errors.push("Cannot set INTERLEAVE to true and pass in a SECOND_END_FASTQ".to_string());
        }
        if self.unpaired_fastq.is_some() && self.second_end_fastq.is_none() {
            errors.push("UNPAIRED_FASTQ may only be set when also emitting read1 and read2 fastqs (so SECOND_END_FASTQ must also be set).".to_string());
        }
        if self.clipping_attribute.is_some() != self.clipping_action.is_some() {
            errors.push(
                "Both or neither of CLIPPING_ATTRIBUTE and CLIPPING_ACTION should be set."
                    .to_string(),
            );
        }
        if let Some(action) = &self.clipping_action {
            if action != "N" && action != "X" && action.parse::<i32>().is_err() {
                errors.push("CLIPPING ACTION must be one of: N, X, or an integer".to_string());
            }
        }
        if self.output_per_rg != self.output_dir.is_some() {
            errors.push("If OUTPUT_PER_RG is true, then OUTPUT_DIR should be set. If ".to_string());
        }
        if self.output_per_rg
            && !(self.rg_tag.eq_ignore_ascii_case("PU") || self.rg_tag.eq_ignore_ascii_case("ID"))
        {
            errors.push("RG_TAG must be: PU or ID".to_string());
        }
        errors
    }
}

/// The files a run leaves, in the order their writers were built, each with what was written.
#[derive(Debug, Default)]
pub struct Outputs {
    pub files: Vec<(String, String)>,
}

impl Outputs {
    fn open(&mut self, path: String) -> usize {
        self.files.push((path, String::new()));
        self.files.len() - 1
    }
}

/// One read group's writers: `FastqWriters` in the reference.
#[derive(Debug, Clone, Copy)]
struct Writers {
    first: usize,
    /// `None` with no second-of-pair writer at all; `Some(None)` for a lazy one not yet built.
    second: Option<Option<usize>>,
    unpaired: usize,
}

/// `IOUtil.makeFileNameSafe`: trimmed, then every whitespace character and every one of
/// ``!"#$%&'()*/:;<=>?@[]\^`{|}~`` replaced by `_`, one for one.
fn file_name_safe(name: &str) -> String {
    name.trim()
        .chars()
        .map(|c| {
            if c.is_whitespace() || "!\"#$%&'()*/:;<=>?@[]\\^`{|}~".contains(c) {
                '_'
            } else {
                c
            }
        })
        .collect()
}

/// `new File(OUTPUT_DIR, name)` or `new File(name)`.
fn in_directory(dir: Option<&str>, name: &str) -> String {
    match dir {
        Some(dir) => format!("{}/{name}", dir.trim_end_matches('/')),
        None => name.to_string(),
    }
}

/// A two-character tag name as a key; any other length names no tag a record can carry.
fn tag_named(name: &str) -> Option<Tag> {
    match name.as_bytes() {
        [a, b] => Some(Tag::new(&[*a, *b])),
        _ => None,
    }
}

/// `SAMRecord.getReadGroup()`: the header's record for the read's `RG`, or `None`.
fn read_group_of(header: &SamHeader, rec: &BamRecord) -> Option<usize> {
    match rec.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => header.read_groups.iter().position(|g| g.id == *id),
        _ => None,
    }
}

/// `clip(src, point, replacement, posStrand)`, including the out-of-range substring the reference
/// throws on.
fn clip(
    src: &str,
    point: i32,
    replacement: Option<char>,
    pos_strand: bool,
) -> Result<String, String> {
    let len = src.chars().count() as i32;
    let chars: Vec<char> = src.chars().collect();
    let (begin, end) = if pos_strand {
        (0, point - 1)
    } else {
        (len - point + 1, len)
    };
    if begin < 0 || end > len || begin > end {
        return Err(format!(
            "java.lang.StringIndexOutOfBoundsException: begin {begin}, end {end}, length {len}"
        ));
    }
    let mut out: String = chars[begin as usize..end as usize].iter().collect();
    if let Some(r) = replacement {
        let count = (len - point + 1).max(0) as usize;
        let pad: String = std::iter::repeat_n(r, count).collect();
        out = if pos_strand { out + &pad } else { pad + &out };
    }
    Ok(out)
}

/// `String.substring(begin)` and `substring(0, end)`, with Java's exception for a bad bound.
fn substring(s: &str, begin: i32, end: Option<i32>) -> Result<String, String> {
    let len = s.chars().count() as i32;
    let end = end.unwrap_or(len);
    if begin < 0 || end > len || begin > end {
        return Err(format!(
            "java.lang.StringIndexOutOfBoundsException: begin {begin}, end {end}, length {len}"
        ));
    }
    Ok(s.chars()
        .skip(begin as usize)
        .take((end - begin) as usize)
        .collect())
}

/// `TrimmingUtil.findQualityTrimPoint`.
fn quality_trim_point(quals: &[i32], trim_qual: i32) -> i32 {
    let length = quals.len() as i32;
    if trim_qual < 1 || length == 0 {
        return 0;
    }
    let (mut score, mut max_score, mut trim_point) = (0, 0, length);
    for i in (0..quals.len()).rev() {
        score += trim_qual - quals[i];
        if score < 0 {
            break;
        }
        if score > max_score {
            max_score = score;
            trim_point = i as i32;
        }
    }
    trim_point
}

/// `writeRecord`: one read's FASTQ block.
fn fastq_block(
    rec: &BamRecord,
    mate: Option<u32>,
    opts: &Options,
    bases_to_trim: i32,
    max_bases: Option<i32>,
) -> Result<String, String> {
    let header = match mate {
        None => rec.read_name.clone(),
        Some(n) => format!("{}/{n}", rec.read_name),
    };
    let mut read = if rec.read_bases.is_empty() {
        "*".to_string()
    } else {
        String::from_utf8_lossy(&rec.read_bases).into_owned()
    };
    let mut quals = if rec.base_qualities.is_empty() {
        "*".to_string()
    } else {
        phred_to_fastq(&rec.base_qualities)
    };
    let negative = rec.flags & READ_NEGATIVE_STRAND != 0;
    if let (Some(attribute), Some(action)) = (&opts.clipping_attribute, &opts.clipping_action) {
        let point = match tag_named(attribute).and_then(|t| rec.tags.get(t)) {
            None => None,
            Some(TagValue::Int(v)) => Some(*v as i32),
            Some(_) => {
                return Err(format!(
                    "java.lang.ClassCastException: the {attribute} attribute is not an Integer"
                ))
            }
        };
        if let Some(mut point) = point {
            if point < opts.clipping_min_length {
                point = (read.chars().count() as i32).min(opts.clipping_min_length);
            }
            if action.eq_ignore_ascii_case("X") {
                read = clip(&read, point, None, !negative)?;
                quals = clip(&quals, point, None, !negative)?;
            } else if action.eq_ignore_ascii_case("N") {
                read = clip(&read, point, Some('N'), !negative)?;
            } else {
                let q = action.parse::<i32>().map_err(|e| format!("{e}"))?;
                let new_qual = phred_to_fastq(&[q as u8]).chars().next().unwrap_or('!');
                quals = clip(&quals, point, Some(new_qual), !negative)?;
            }
        }
    }
    if opts.re_reverse && negative {
        let mut bytes = read.into_bytes();
        reverse_complement(&mut bytes);
        read = String::from_utf8(bytes).expect("bases stay ASCII");
        quals = quals.chars().rev().collect();
    }
    if bases_to_trim > 0 {
        read = substring(&read, bases_to_trim, None)?;
        quals = substring(&quals, bases_to_trim, None)?;
    }
    if let Some(quality) = opts.quality {
        let phred: Vec<i32> = quals.bytes().map(|b| b as i32 - 33).collect();
        let index = quality_trim_point(&phred, quality).max(1);
        if (index as usize) < phred.len() {
            read = substring(&read, 0, Some(index))?;
            quals = substring(&quals, 0, Some(index))?;
        }
    }
    if let Some(max) = max_bases {
        if (max as usize) < read.chars().count() {
            read = substring(&read, 0, Some(max))?;
            quals = substring(&quals, 0, Some(max))?;
        }
    }
    Ok(write_record(&FastqRecord {
        read_name: Some(header),
        read_string: Some(read),
        quality_header: Some(String::new()),
        quality_string: Some(quals),
    }))
}

/// `writeTagRecords`: one block per tag group, for the group's file.
fn tag_blocks(rec: &BamRecord, mate: Option<u32>, tags: &TagGroups) -> Result<Vec<String>, String> {
    let header = match mate {
        None => rec.read_name.clone(),
        Some(n) => format!("{}/{n}", rec.read_name),
    };
    let value = |tag: &str| -> Result<String, String> {
        match tag_named(tag).and_then(|t| rec.tags.get(t)) {
            Some(TagValue::Str(s)) => Ok(s.clone()),
            None => Err(format!(
                "picard.PicardException: Record: {} does have a value for tag: {tag}",
                rec.read_name
            )),
            Some(_) => Err(format!(
                "java.lang.ClassCastException: the {tag} attribute is not a String"
            )),
        }
    };
    let mut blocks = Vec::new();
    for (i, group) in tags.sequence.iter().enumerate() {
        let separator = tags.separator.get(i).map(String::as_str).unwrap_or("");
        let sequence = group
            .trim()
            .split(',')
            .map(value)
            .collect::<Result<Vec<_>, _>>()?
            .join(separator);
        let quality = if tags.quality.is_empty() {
            "~".repeat(sequence.chars().count())
        } else {
            tags.quality[i]
                .trim()
                .split(',')
                .map(value)
                .collect::<Result<Vec<_>, _>>()?
                .join(&"~".repeat(separator.chars().count()))
        };
        blocks.push(write_record(&FastqRecord {
            read_name: Some(header.clone()),
            read_string: Some(sequence),
            quality_header: Some(String::new()),
            quality_string: Some(quality),
        }));
    }
    Ok(blocks)
}

/// `doWork`, from the writers to the `MATE_NOT_FOUND` check. `strict` is whether
/// `VALIDATION_STRINGENCY` is `STRICT`, which is the only stringency that turns reads left without
/// a mate into a failure, and only once everything has been written.
pub fn run(
    header: &SamHeader,
    records: &[BamRecord],
    opts: &Options,
    strict: bool,
) -> (Outputs, Result<(), String>) {
    let mut out = Outputs::default();
    let result = run_into(header, records, opts, strict, &mut out);
    (out, result)
}

fn run_into(
    header: &SamHeader,
    records: &[BamRecord],
    opts: &Options,
    strict: bool,
    out: &mut Outputs,
) -> Result<(), String> {
    let groups = &header.read_groups;
    // `generateWriters`: the key `None` is the reads without a read group.
    let mut writers: Vec<(Option<usize>, Writers)> = Vec::new();
    let rg_file = |index: usize, suffix: &str| -> Result<String, String> {
        let group = &groups[index];
        let name = if opts.rg_tag.eq_ignore_ascii_case("PU") {
            group.attributes.get("PU").map(str::to_string)
        } else if opts.rg_tag.eq_ignore_ascii_case("ID") {
            Some(group.id.clone())
        } else {
            None
        };
        let name = name.ok_or_else(|| {
            format!(
                "picard.PicardException: The selected RG_TAG: {} is not present in the header.",
                opts.rg_tag
            )
        })?;
        let extension = if opts.compress_outputs_per_rg {
            ".fastq.gz"
        } else {
            ".fastq"
        };
        Ok(in_directory(
            opts.output_dir.as_deref(),
            &format!("{}{suffix}{extension}", file_name_safe(&name)),
        ))
    };
    if !opts.output_per_rg {
        // `IOUtil.assertFileIsWritable(FASTQ)`, whose null check is worded for reading.
        let Some(fastq) = opts.fastq.clone() else {
            return Err(
                "java.lang.IllegalArgumentException: Cannot check readability of null file."
                    .to_string(),
            );
        };
        let first = out.open(fastq);
        let second = if opts.interleave {
            Some(Some(first))
        } else {
            opts.second_end_fastq
                .clone()
                .map(|path| Some(out.open(path)))
        };
        let unpaired = match &opts.unpaired_fastq {
            Some(path) => out.open(path.clone()),
            None => first,
        };
        let set = Writers {
            first,
            second,
            unpaired,
        };
        writers.push((None, set));
        for index in 0..groups.len() {
            writers.push((Some(index), set));
        }
    } else {
        for index in 0..groups.len() {
            let first = out.open(rg_file(index, "_1")?);
            let second = if opts.interleave {
                Some(Some(first))
            } else {
                Some(None)
            };
            writers.push((
                Some(index),
                Writers {
                    first,
                    second,
                    unpaired: first,
                },
            ));
        }
    }
    // `generateAdditionalWriters`: the tag files, after the read files.
    let mut tag_writers: Vec<(Option<usize>, Vec<usize>)> = Vec::new();
    if let Some(tags) = &opts.tags {
        let parent = opts.fastq.as_deref().and_then(|f| {
            std::path::Path::new(f)
                .parent()
                .map(|p| p.display().to_string())
        });
        let make = |prefix: &str, out: &mut Outputs| -> Vec<usize> {
            tags.sequence
                .iter()
                .map(|group| {
                    let mut name = file_name_safe(&format!("{prefix}{}", group.replace(',', "_")));
                    name += if tags.compress { ".fastq.gz" } else { ".fastq" };
                    let dir = opts.output_dir.as_deref().or(parent.as_deref());
                    out.open(in_directory(dir, &name))
                })
                .collect()
        };
        if !opts.output_per_rg {
            let set = make("", out);
            tag_writers.push((None, set.clone()));
            for index in 0..groups.len() {
                tag_writers.push((Some(index), set.clone()));
            }
        } else {
            for (index, group) in groups.iter().enumerate() {
                // `getPlatformUnit() + "_"`: a missing PU is the string "null".
                let prefix = if opts.rg_tag.eq_ignore_ascii_case("PU") {
                    format!("{}_", group.attributes.get("PU").unwrap_or("null"))
                } else {
                    format!("{}_", group.id)
                };
                let set = make(&prefix, out);
                tag_writers.push((Some(index), set));
            }
        }
    }
    if writers.is_empty() {
        let base = format!("{} does not contain Read Groups", opts.input);
        return Err(format!(
            "picard.PicardException: {}",
            if opts.output_per_rg {
                format!("{base}, consider not using the OUTPUT_PER_RG option")
            } else {
                base
            }
        ));
    }

    let mut first_seen: std::collections::HashMap<&str, &BamRecord> =
        std::collections::HashMap::new();
    for rec in records {
        if rec.flags & (SECONDARY_ALIGNMENT | SUPPLEMENTARY_ALIGNMENT) != 0
            && !opts.include_non_primary_alignments
        {
            continue;
        }
        if rec.flags & READ_FAILS_VENDOR_QUALITY != 0 && !opts.include_non_pf_reads {
            continue;
        }
        let group = read_group_of(header, rec);
        let slot = writers.iter().position(|(g, _)| *g == group);
        let (mut read1, mut read2) = (None, None);
        if rec.flags & READ_PAIRED != 0 {
            match first_seen.remove(rec.read_name.as_str()) {
                None => {
                    first_seen.insert(rec.read_name.as_str(), rec);
                }
                Some(first) => {
                    let mates = (first.flags & FIRST_OF_PAIR != 0
                        && rec.flags & SECOND_OF_PAIR != 0)
                        || (rec.flags & FIRST_OF_PAIR != 0 && first.flags & SECOND_OF_PAIR != 0);
                    if !mates {
                        return Err(format!(
                            "picard.PicardException: Illegal mate state: {}",
                            first.read_name
                        ));
                    }
                    let (r1, r2) = if rec.flags & FIRST_OF_PAIR != 0 {
                        (rec, first)
                    } else {
                        (first, rec)
                    };
                    let Some(slot) = slot else {
                        return Err("java.lang.NullPointerException: Cannot invoke \"picard.sam.SamToFastq$FastqWriters.getFirstOfPair()\" because \"fq\" is null".to_string());
                    };
                    let set = writers[slot].1;
                    let block = fastq_block(
                        r1,
                        Some(1),
                        opts,
                        opts.read1_trim,
                        opts.read1_max_bases_to_write,
                    )?;
                    out.files[set.first].1.push_str(&block);
                    let second = match set.second {
                        None => {
                            return Err("picard.PicardException: Input contains paired reads but no SECOND_END_FASTQ specified.".to_string())
                        }
                        Some(Some(index)) => index,
                        Some(None) => {
                            let index = out.open(rg_file(group.expect("per-group writers"), "_2")?);
                            for entry in writers.iter_mut() {
                                if entry.0 == group {
                                    entry.1.second = Some(Some(index));
                                }
                            }
                            index
                        }
                    };
                    let block = fastq_block(
                        r2,
                        Some(2),
                        opts,
                        opts.read2_trim,
                        opts.read2_max_bases_to_write,
                    )?;
                    out.files[second].1.push_str(&block);
                    read1 = Some(r1);
                    read2 = Some(r2);
                }
            }
        } else {
            let Some(slot) = slot else {
                return Err("java.lang.NullPointerException: Cannot invoke \"picard.sam.SamToFastq$FastqWriters.getUnpaired()\" because \"fq\" is null".to_string());
            };
            let block = fastq_block(
                rec,
                None,
                opts,
                opts.read1_trim,
                opts.read1_max_bases_to_write,
            )?;
            out.files[writers[slot].1.unpaired].1.push_str(&block);
        }
        // `handleAdditionalRecords`.
        if let Some(tags) = &opts.tags {
            let Some(files) = tag_writers
                .iter()
                .find(|(g, _)| *g == group)
                .map(|(_, f)| f.clone())
            else {
                continue;
            };
            let blocks: Vec<(usize, String)> = if rec.flags & READ_PAIRED != 0 {
                match (read1, read2) {
                    (Some(r1), Some(r2)) => {
                        let mut both = Vec::new();
                        if !tags.sequence.is_empty() {
                            for (i, b) in tag_blocks(r1, Some(1), tags)?.into_iter().enumerate() {
                                both.push((files[i], b));
                            }
                            for (i, b) in tag_blocks(r2, Some(2), tags)?.into_iter().enumerate() {
                                both.push((files[i], b));
                            }
                        }
                        both
                    }
                    _ => Vec::new(),
                }
            } else if tags.sequence.is_empty() {
                Vec::new()
            } else {
                tag_blocks(rec, None, tags)?
                    .into_iter()
                    .enumerate()
                    .map(|(i, b)| (files[i], b))
                    .collect()
            };
            for (file, block) in blocks {
                out.files[file].1.push_str(&block);
            }
        }
    }
    if !first_seen.is_empty() && strict {
        return Err(format!(
            "htsjdk.samtools.SAMFormatException: SAM validation error: ERROR::MATE_NOT_FOUND:Found {} unpaired mates",
            first_seen.len()
        ));
    }
    Ok(())
}

/// The runnable binary, for both tools: `NAME=value` arguments in, the files and the exit code
/// out. `with_tags` is `SamToFastqWithTags`.
pub fn main_with(with_tags: bool) -> ! {
    let args: Vec<String> = std::env::args().collect();
    let one = |key: &str| -> Option<String> {
        args.iter()
            .find_map(|a| a.strip_prefix(&format!("{key}=")).map(str::to_string))
    };
    let all = |key: &str| -> Vec<String> {
        args.iter()
            .filter_map(|a| a.strip_prefix(&format!("{key}=")).map(str::to_string))
            .collect()
    };
    let flag = |key: &str, default: bool| one(key).map_or(default, |v| crate::args::is_true(&v));
    let int = |key: &str| one(key).and_then(|v| v.parse::<i32>().ok());
    let Some(input) = one("INPUT").or_else(|| one("I")) else {
        eprintln!("INPUT= is required");
        std::process::exit(1);
    };
    let mut opts = Options::new(&input);
    opts.fastq = one("FASTQ").or_else(|| one("F"));
    opts.second_end_fastq = one("SECOND_END_FASTQ").or_else(|| one("F2"));
    opts.unpaired_fastq = one("UNPAIRED_FASTQ").or_else(|| one("FU"));
    opts.output_per_rg = flag("OUTPUT_PER_RG", false);
    opts.compress_outputs_per_rg = flag("COMPRESS_OUTPUTS_PER_RG", false);
    if let Some(tag) = one("RG_TAG") {
        opts.rg_tag = tag;
    }
    opts.output_dir = one("OUTPUT_DIR");
    opts.re_reverse = flag("RE_REVERSE", true);
    opts.interleave = flag("INTERLEAVE", false);
    opts.include_non_pf_reads = flag("INCLUDE_NON_PF_READS", false);
    opts.include_non_primary_alignments = flag("INCLUDE_NON_PRIMARY_ALIGNMENTS", false);
    opts.clipping_attribute = one("CLIPPING_ATTRIBUTE");
    opts.clipping_action = one("CLIPPING_ACTION");
    opts.clipping_min_length = int("CLIPPING_MIN_LENGTH").unwrap_or(0);
    opts.read1_trim = int("READ1_TRIM").unwrap_or(0);
    opts.read1_max_bases_to_write = int("READ1_MAX_BASES_TO_WRITE");
    opts.read2_trim = int("READ2_TRIM").unwrap_or(0);
    opts.read2_max_bases_to_write = int("READ2_MAX_BASES_TO_WRITE");
    opts.quality = int("QUALITY");
    if with_tags {
        opts.tags = Some(TagGroups {
            sequence: all("SEQUENCE_TAG_GROUP"),
            quality: all("QUALITY_TAG_GROUP"),
            separator: all("TAG_GROUP_SEPERATOR"),
            compress: flag("COMPRESS_OUTPUTS_PER_TAG_GROUP", false),
        });
    }
    // Barclay's required-argument check, before the tool's own validation: FASTQ may be left out
    // only when one of the arguments it is mutually exclusive with is given.
    let given = |key: &str| one(key).is_some();
    if opts.fastq.is_none()
        && !given("OUTPUT_PER_RG")
        && !given("COMPRESS_OUTPUTS_PER_RG")
        && !given("OUTPUT_DIR")
    {
        eprintln!("Argument FASTQ was missing: Argument 'FASTQ' is required unless one of {{[OUTPUT_PER_RG, COMPRESS_OUTPUTS_PER_RG, OUTPUT_DIR]}} are provided");
        std::process::exit(1);
    }
    let stringency = one("VALIDATION_STRINGENCY").unwrap_or_else(|| "STRICT".to_string());
    if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
        eprintln!("unknown VALIDATION_STRINGENCY: {stringency}");
        std::process::exit(1);
    }
    // Barclay prints the usage and then every failure, one per line.
    let errors = opts.validate();
    if !errors.is_empty() {
        for error in errors {
            eprintln!("{error}");
        }
        std::process::exit(1);
    }

    let read = || -> Result<(SamHeader, Vec<BamRecord>), String> {
        let raw = std::fs::read(&input).map_err(|e| e.to_string())?;
        let raw = crate::cram_input::uncram(raw);
        if raw.starts_with(&[0x1f, 0x8b]) {
            let plain = crate::bgzf::decompress_owned(raw).map_err(|e| format!("{e:?}"))?;
            let reader = crate::record_validation::reader(&plain).map_err(|e| format!("{e:?}"))?;
            let header = reader.header.text.clone();
            let records = reader
                .map(|r| r.map_err(|e| format!("{e:?}")))
                .collect::<Result<_, _>>()?;
            Ok((header, records))
        } else {
            let text = String::from_utf8(raw).map_err(|e| e.to_string())?;
            htsjdk_bam::sam_file::read_sam(&text).map_err(|e| format!("{e:?}"))
        }
    };
    let (header, records) = match read() {
        Ok(read) => read,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };
    let (outputs, result) = run(&header, &records, &opts, stringency == "STRICT");
    for (path, text) in &outputs.files {
        let written = if path.ends_with(".gz") {
            let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
            std::io::Write::write_all(&mut writer, text.as_bytes())
                .and_then(|()| writer.into_inner())
                .and_then(|bytes| std::fs::write(path, bytes))
        } else {
            std::fs::write(path, text)
        };
        if let Err(e) = written {
            eprintln!("{path}: {e}");
            std::process::exit(1);
        }
    }
    if let Err(thrown) = result {
        eprintln!("Exception in thread \"main\" {thrown}");
        std::process::exit(1);
    }
    std::process::exit(0);
}
