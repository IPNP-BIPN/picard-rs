//! `BamToBfq` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fastq.BamToBfq` and `picard.fastq.SamToBfqWriter` at tag 3.4.0, with the record
//! encoding in `picard_analysis::bam_to_bfq`.
//!
//! What decides the bytes:
//!
//! * Barclay validates the arguments in field order, and it counts an argument given as `null` as
//!   SET: `FLOWCELL_BARCODE=null` beside `OUTPUT_FILE_PREFIX` is the mutex refusal, and so is
//!   `RUN_BARCODE=null` beside `READ_NAME_PREFIX`. `customCommandLineValidation` then fills the
//!   prefix in as `FLOWCELL_BARCODE + "." + LANE` and the name prefix as `RUN_BARCODE + ":"`, with
//!   Java's `null` spelled out.
//! * `READS_TO_ALIGN` first counts the writable records -- which wants a queryname-sorted header --
//!   and keeps every `floor(count / READS_TO_ALIGN)`th. The count is not the writer's filter: it
//!   takes any `XN` as noise where the writer takes only `XN=1`, it knows nothing of the
//!   whole-read clip, and it runs in the constructor before `includeNonPfReads` is assigned, so
//!   it always drops the reads that fail vendor QC.
//! * Each read is written as a name (with `/1` or `/2`), a length and one byte per base; a
//!   paired run sends a read to the file its first-of-pair flag names. Every `READ_CHUNK_SIZE`
//!   records the files roll over to the next index, so a run that ends on a multiple leaves an
//!   empty pair of files behind. A `.bfq` is gzip (`IOUtil.hasGzipFileExtension`), and a gzip
//!   stream is what this writes: the harness compares what it inflates to.

use std::fs::File;
use std::io::Write;

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bgzf::BgzfWriter;
use picard_analysis::bam_to_bfq::{encode_base_and_quality, MAX_SEED_REGION_NOCALL_FIXES};
use picard_analysis::metrics_cli::{absolute, read_input, refuse_validation, sort_order, thrown};

const SEED_REGION_LENGTH: usize = 28;

/// The command line as Barclay sees it: every argument given, `null` included.
struct CommandLine {
    pairs: Vec<(String, String)>,
}

impl CommandLine {
    fn from_env() -> Self {
        const ALIASES: [(&str, &str); 9] = [
            ("I", "INPUT"),
            ("F", "FLOWCELL_BARCODE"),
            ("L", "LANE"),
            ("NUM", "READS_TO_ALIGN"),
            ("CHUNK", "READ_CHUNK_SIZE"),
            ("PE", "PAIRED_RUN"),
            ("RB", "RUN_BARCODE"),
            ("NONPF", "INCLUDE_NON_PF_READS"),
            ("R", "REFERENCE_SEQUENCE"),
        ];
        let mut pairs = Vec::new();
        for raw in std::env::args().skip(1) {
            let raw = raw.trim_start_matches('-');
            if let Some((name, value)) = raw.split_once('=') {
                let long = ALIASES
                    .iter()
                    .find(|(short, _)| *short == name)
                    .map_or(name, |(_, long)| long);
                pairs.push((long.to_string(), value.to_string()));
            }
        }
        CommandLine { pairs }
    }

    fn has_been_set(&self, name: &str) -> bool {
        self.pairs.iter().any(|(n, _)| n == name)
    }

    /// The value, `None` for absent or `null`.
    fn get(&self, name: &str) -> Option<&str> {
        self.pairs
            .iter()
            .rev()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_str())
            .filter(|v| *v != "null")
    }

    fn integer(&self, name: &str) -> Option<i32> {
        self.get(name).map(|v| {
            v.parse().unwrap_or_else(|_| {
                refuse_validation(
                    "BamToBfq",
                    &[format!(
                        "Argument '{name}' cannot be set to '{v}': it is not a number"
                    )],
                )
            })
        })
    }

    fn boolean(&self, name: &str) -> Option<bool> {
        self.get(name).map(|v| v.eq_ignore_ascii_case("true"))
    }
}

/// Barclay's `validateValues`, one argument at a time in field order; `mutex` is each argument's
/// mutex list after `validateArgumentDefinitions` has made the relation symmetric.
fn validate(args: &CommandLine) {
    let fields: [(&str, bool, &[&str]); 13] = [
        ("INPUT", false, &[]),
        ("ANALYSIS_DIR", false, &[]),
        ("FLOWCELL_BARCODE", false, &["OUTPUT_FILE_PREFIX"]),
        ("LANE", true, &["OUTPUT_FILE_PREFIX"]),
        ("OUTPUT_FILE_PREFIX", false, &["FLOWCELL_BARCODE", "LANE"]),
        ("READS_TO_ALIGN", true, &[]),
        ("READ_CHUNK_SIZE", true, &[]),
        ("PAIRED_RUN", false, &[]),
        ("RUN_BARCODE", true, &["READ_NAME_PREFIX"]),
        ("READ_NAME_PREFIX", true, &["RUN_BARCODE"]),
        ("INCLUDE_NON_PF_READS", true, &[]),
        ("CLIP_ADAPTERS", true, &[]),
        ("BASES_TO_WRITE", true, &[]),
    ];
    for (name, optional, mutex) in fields {
        let provided: Vec<&str> = mutex
            .iter()
            .copied()
            .filter(|m| args.has_been_set(m))
            .collect();
        if args.has_been_set(name) && !provided.is_empty() {
            refuse_validation(
                "BamToBfq",
                &[format!(
                    "Argument '{name}' cannot be used in conjunction with argument(s) {}",
                    provided.join(" ")
                )],
            );
        }
        if !optional && !args.has_been_set(name) && provided.is_empty() {
            let why = if mutex.is_empty() {
                format!("Argument '{name}' is required")
            } else {
                format!(
                    "Argument '{name}' is required unless one of {{[{}]}} are provided",
                    mutex.join(", ")
                )
            };
            refuse_validation("BamToBfq", &[format!("Argument {name} was missing: {why}")]);
        }
    }
}

fn int_tag(record: &BamRecord, tag: &[u8; 2]) -> Option<i64> {
    match record.tags.get(Tag::new(tag)) {
        Some(TagValue::Int(v)) => Some(*v),
        _ => None,
    }
}

fn has_tag(record: &BamRecord, tag: &[u8; 2]) -> bool {
    record.tags.get(Tag::new(tag)).is_some()
}

/// `TagFilter(XN, 1)`.
fn is_noise(record: &BamRecord) -> bool {
    int_tag(record, b"XN") == Some(1)
}

/// `FailsVendorReadQualityFilter`.
fn fails_vendor_quality(record: &BamRecord) -> bool {
    record.flags & 0x200 != 0
}

/// `WholeReadClippedFilter`.
fn whole_read_clipped(record: &BamRecord) -> bool {
    int_tag(record, b"XT") == Some(1)
}

fn first_of_pair(record: &BamRecord) -> bool {
    record.flags & 0x40 != 0
}

/// One `.bfq` being written: a gzip stream, which `BinaryCodec` writes little-endian into.
struct Codec {
    writer: BgzfWriter<File>,
}

impl Codec {
    fn open(path: &str) -> Codec {
        let file = File::create(path).unwrap_or_else(|e| {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Error opening file for writing: file://{path}: {e}"
            ))
        });
        Codec {
            writer: BgzfWriter::new(file),
        }
    }

    fn write(&mut self, bytes: &[u8]) {
        self.writer
            .write_all(bytes)
            .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}")));
    }

    fn close(mut self) {
        self.writer
            .finish()
            .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}")));
    }
}

struct BfqWriter {
    input: String,
    output_prefix: String,
    paired: bool,
    increment: i32,
    chunk: i32,
    name_prefix: String,
    include_non_pf: bool,
    clip_adapters: bool,
    bases_to_write: Option<i32>,
    wrote: i32,
    codec1: Option<Codec>,
    codec2: Option<Codec>,
}

impl BfqWriter {
    /// `initializeNextBfqFiles`.
    fn next_files(&mut self, index: i32) {
        if let Some(codec) = self.codec1.take() {
            codec.close();
        }
        if let Some(codec) = self.codec2.take() {
            codec.close();
        }
        self.codec1 = Some(Codec::open(&format!("{}{index}.1.bfq", self.output_prefix)));
        if self.paired {
            self.codec2 = Some(Codec::open(&format!("{}{index}.2.bfq", self.output_prefix)));
        }
    }

    /// `writeFastqRecord` into codec 1 or 2.
    fn write_record(&mut self, to_first: bool, record: &BamRecord, name: &str) {
        let name = name.strip_prefix(self.name_prefix.as_str()).unwrap_or(name);
        let mut out: Vec<u8> = Vec::new();
        let name_length = name.encode_utf16().count() as i32 + 1;
        out.extend_from_slice(&name_length.to_le_bytes());
        out.extend(name.chars().map(|c| c as u32 as u8));
        out.push(0);

        let seqs = &record.read_bases;
        // `getBaseQualityString`: FASTQ-encoded.
        let quals: Vec<i32> = record
            .base_qualities
            .iter()
            .map(|q| i32::from(*q) + 33)
            .collect();
        let mut retained = seqs.len();
        if self.clip_adapters {
            if let Some(trim) = int_tag(record, b"XT") {
                retained = seqs
                    .len()
                    .min((SEED_REGION_LENGTH as i64).max(trim - 1).max(0) as usize);
            }
        }
        let length = self.bases_to_write.unwrap_or(seqs.len() as i32);
        out.extend_from_slice(&length.to_le_bytes());
        if length < 0 {
            thrown(&format!("java.lang.NegativeArraySizeException: {length}"));
        }
        let mut encoded = vec![0u8; length as usize];
        let mut fixes = 0;
        for i in 0..retained.min(encoded.len()) {
            let mut quality = (quals[i] - 33).min(63);
            let base = match seqs[i] {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                b'N' | b'n' | b'.' => {
                    if i < SEED_REGION_LENGTH {
                        if fixes < MAX_SEED_REGION_NOCALL_FIXES {
                            quality = 1;
                            fixes += 1;
                        } else {
                            quality = 0;
                        }
                    } else {
                        quality = 1;
                    }
                    0
                }
                other => thrown(&format!(
                    "picard.PicardException: Unknown base when writing bfq file: {}",
                    other as char
                )),
            };
            encoded[i] = encode_base_and_quality(base, quality as u8);
        }
        for byte in encoded.iter_mut().skip(retained) {
            *byte = encode_base_and_quality(0, 1);
        }
        out.extend_from_slice(&encoded);
        let codec = if to_first {
            self.codec1.as_mut()
        } else {
            self.codec2.as_mut()
        };
        codec.expect("the codec is open").write(&out);
    }

    /// After a record is written: count it, and roll the files over at a chunk boundary.
    fn wrote_one(&mut self, file_index: &mut i32) {
        self.wrote += 1;
        if self.chunk > 0 && self.wrote % self.chunk == 0 {
            self.next_files(*file_index);
            *file_index += 1;
        }
    }

    fn write_single_end(&mut self, records: &[BamRecord]) {
        let mut file_index = 0;
        self.next_files(file_index);
        file_index += 1;
        let mut count: i32 = 0;
        for record in records {
            if is_noise(record)
                || whole_read_clipped(record)
                || (!self.include_non_pf && fails_vendor_quality(record))
            {
                continue;
            }
            count = count.wrapping_add(1);
            if count % self.increment == 0 {
                let name = format!("{}/1", record.read_name);
                self.write_record(true, record, &name);
                self.wrote_one(&mut file_index);
            }
        }
    }

    fn write_paired_end(&mut self, records: &[BamRecord]) {
        let mut file_index = 0;
        self.next_files(file_index);
        file_index += 1;
        let mut count: i32 = 0;
        let mut iter = records.iter();
        while let Some(first) = iter.next() {
            let Some(second) = iter.next() else {
                thrown(&format!(
                    "picard.PicardException: Mismatched number of records in {}",
                    absolute(&self.input)
                ));
            };
            if second.read_name != first.read_name || first_of_pair(first) == first_of_pair(second)
            {
                thrown(&format!(
                    "picard.PicardException: Unmatched read pairs in {}: {}, {}.",
                    absolute(&self.input),
                    first.read_name,
                    second.read_name
                ));
            }
            if is_noise(first) && is_noise(second) {
                continue;
            }
            if !self.include_non_pf && (fails_vendor_quality(first) || fails_vendor_quality(second))
            {
                continue;
            }
            if whole_read_clipped(first) || whole_read_clipped(second) {
                continue;
            }
            count = count.wrapping_add(1);
            if count % self.increment == 0 {
                let name = format!("{}/1", first.read_name);
                self.write_record(first_of_pair(first), first, &name);
                let name = format!("{}/2", second.read_name);
                self.write_record(first_of_pair(second), second, &name);
                self.wrote_one(&mut file_index);
            }
        }
    }
}

/// `countWritableRecords`, which `READS_TO_ALIGN` divides.
fn count_writable_records(
    input: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    paired: bool,
    include_non_pf: bool,
) -> i32 {
    if sort_order(header) != "queryname" {
        thrown(&format!(
            "picard.PicardException: Input file ({}) needs to be sorted by queryname.",
            absolute(input)
        ));
    }
    let mut count = 0;
    if !paired {
        for record in records {
            if is_noise(record) || (!include_non_pf && fails_vendor_quality(record)) {
                continue;
            }
            count += 1;
        }
    } else {
        let mut iter = records.iter();
        while let Some(first) = iter.next() {
            let Some(second) = iter.next() else {
                thrown("java.util.NoSuchElementException");
            };
            if has_tag(first, b"XN") && has_tag(second, b"XN") {
                continue;
            }
            if !include_non_pf && (fails_vendor_quality(first) || fails_vendor_quality(second)) {
                continue;
            }
            count += 1;
        }
    }
    count
}

/// `(int) Math.floor(double)`: NaN is 0 and the infinities saturate.
fn java_floor_to_int(value: f64) -> i32 {
    if value.is_nan() {
        0
    } else {
        value.floor().clamp(i32::MIN as f64, i32::MAX as f64) as i32
    }
}

fn main() {
    let args = CommandLine::from_env();
    validate(&args);

    let input = args.get("INPUT").unwrap_or("null").to_string();
    let analysis_dir = args.get("ANALYSIS_DIR").unwrap_or("null").to_string();
    let lane = args.integer("LANE");
    let reads_to_align = args.integer("READS_TO_ALIGN");
    let chunk = args.integer("READ_CHUNK_SIZE");
    let chunk = if args.has_been_set("READ_CHUNK_SIZE") {
        chunk
    } else {
        Some(2_000_000)
    };
    let paired = args.boolean("PAIRED_RUN").unwrap_or_else(|| {
        thrown("java.lang.NullPointerException");
    });
    let include_non_pf = args.boolean("INCLUDE_NON_PF_READS").unwrap_or(false);
    let clip_adapters = args.boolean("CLIP_ADAPTERS").unwrap_or(true);
    let bases_to_write = args.integer("BASES_TO_WRITE");

    // `customCommandLineValidation`.
    let output_file_prefix = match args.get("OUTPUT_FILE_PREFIX") {
        Some(prefix) => prefix.to_string(),
        None => format!(
            "{}.{}",
            args.get("FLOWCELL_BARCODE").unwrap_or("null"),
            lane.map_or("null".to_string(), |l| l.to_string())
        ),
    };
    let name_prefix = match args.get("READ_NAME_PREFIX") {
        Some(prefix) => prefix.to_string(),
        None => format!("{}:", args.get("RUN_BARCODE").unwrap_or("null")),
    };

    // `doWork`.
    let mut output_prefix = absolute(&analysis_dir);
    if !output_prefix.ends_with('/') {
        output_prefix.push('/');
    }
    output_prefix.push_str(&output_file_prefix);
    output_prefix.push('.');

    if !std::path::Path::new(&input).is_file() {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
            absolute(&input)
        ));
    }
    let (header, records) = read_input(&input);

    let mut increment = 1;
    if let Some(total) = reads_to_align {
        // The constructor counts BEFORE it assigns `includeNonPfReads`, so the count always
        // sees the field's initial `false`, whatever was asked for.
        let writable = count_writable_records(&input, &header, &records, paired, false);
        increment = java_floor_to_int(f64::from(writable) / f64::from(total));
        if increment == 0 {
            increment = 1;
        }
    }

    let mut writer = BfqWriter {
        input: input.clone(),
        output_prefix,
        paired,
        increment,
        chunk: chunk.unwrap_or(0),
        name_prefix,
        include_non_pf,
        clip_adapters,
        bases_to_write,
        wrote: 0,
        codec1: None,
        codec2: None,
    };
    if paired {
        writer.write_paired_end(&records);
    } else {
        writer.write_single_end(&records);
    }
    if let Some(codec) = writer.codec1.take() {
        codec.close();
    }
    if let Some(codec) = writer.codec2.take() {
        codec.close();
    }
}
