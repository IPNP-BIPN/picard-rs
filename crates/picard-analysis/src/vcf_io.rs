//! Reading a VCF the way `VCFFileReader` hands it to a tool that writes it back out, and writing it.
//!
//! The five `picard.vcf` manipulation tools (`MakeSitesOnlyVcf`, `RenameSampleInVcf`, `SortVcf`,
//! `UpdateVcfSequenceDictionary`, and the reading half of `VcfToIntervalList`) share one
//! round trip, and one fact about it decides the bytes of every genotype column they write.
//!
//! # A genotype block is copied, not re-encoded, when the samples were already sorted
//!
//! `AbstractVCFCodec.parseVCFLine` wraps `parts[8]` -- the FORMAT column and every sample column,
//! joined back with tabs -- in a `LazyGenotypesContext`, and decodes it at once only when
//!
//! ```java
//! if ( !header.samplesWereAlreadySorted() )
//!     lazy.decode();
//! ```
//!
//! `VCFEncoder.write` then has a fast path for a context that is still lazy:
//!
//! ```java
//! if (gc.isLazyWithData() && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String) {
//!     vcfOutput.append(VCFConstants.FIELD_SEPARATOR);
//!     vcfOutput.append(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
//! ```
//!
//! So the same record comes out two different ways. In a file whose sample names are in sorted
//! order the block is copied verbatim -- `GT:GQ:DP` stays in that order and `0/0:.:.` keeps its
//! trailing dots -- and in a file whose names are not, it is decoded and re-encoded: the FORMAT keys
//! sorted after `GT`, the trailing missing fields trimmed, and the columns written in the order of
//! the header the writer was given. Measured on the oracle with the same records under the two
//! sample orders.
//!
//! It also explains how `RenameSampleInVcf` works at all: the header it writes names a sample
//! nobody's genotypes are keyed by, and the copy never looks the name up.
//!
//! Anything that touches the genotypes decodes them (`MakeSitesOnlyVcf`'s `subsetToSamples`), and a
//! decoded record is re-encoded whatever the sample order was.

use htsjdk_vcf::encoder::{EncodeError, MissingFields, VcfEncoder};
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::reader::{read_vcf, ReadFailure, VcfFile};
use htsjdk_vcf::record_parse::{split_condensed, NUM_STANDARD_FIELDS};
use htsjdk_vcf::variant::VariantContext;

/// One record as the codec hands it over: decoded, plus the block it would copy if nothing decodes
/// it first.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub variant: VariantContext,
    /// `parts[8]` while the context is still lazy; `None` once it has been decoded, or when the
    /// line had no genotype columns at all.
    pub lazy_genotypes: Option<String>,
}

impl Record {
    /// `LazyGenotypesContext.decode()`, or any access that forces it: the copy is gone and the
    /// genotypes will be encoded from what was parsed.
    pub fn decode(&mut self) {
        self.lazy_genotypes = None;
    }
}

/// A whole file, read by `VCFFileReader(file, false)` and iterated.
#[derive(Debug, Clone, PartialEq)]
pub struct LazyVcf {
    pub file: VcfFile,
    pub records: Vec<Record>,
}

/// `ParsingUtils.isSorted`: each name no greater than the next, by `String.compareTo`, which is
/// UTF-16 code unit order.
pub fn is_sorted(names: &[String]) -> bool {
    names
        .windows(2)
        .all(|pair| java_compare(&pair[0], &pair[1]) != std::cmp::Ordering::Greater)
}

/// `String.compareTo`.
pub fn java_compare(a: &str, b: &str) -> std::cmp::Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

/// `VCFHeader.getSampleNamesInOrder`: the genotype samples, sorted.
pub fn sample_names_in_order(header: &VcfHeader) -> Vec<String> {
    let mut names = header.samples.clone();
    names.sort_by(|a, b| java_compare(a, b));
    names
}

/// Read a VCF, keeping each record's genotype block while the codec would have kept it lazy.
pub fn read_lazy(text: &str) -> Result<LazyVcf, ReadFailure> {
    let file = read_vcf(text)?;
    let lazy = is_sorted(&file.header.samples) && !file.header.samples.is_empty();

    // The body starts after the `#CHROM` line, exactly where `read_vcf` started decoding, and the
    // records are its non-`#` lines in order.
    let body = text
        .lines()
        .skip_while(|line| !line.starts_with("#CHROM"))
        .skip(1)
        .filter(|line| !line.starts_with('#'));
    let records = file
        .records
        .iter()
        .zip(body)
        .map(|(variant, line)| {
            let lazy_genotypes = if lazy {
                split_condensed(line, '\t', NUM_STANDARD_FIELDS + 1, true)
                    .into_iter()
                    .nth(NUM_STANDARD_FIELDS)
            } else {
                None
            };
            Record {
                variant: variant.clone(),
                lazy_genotypes,
            }
        })
        .collect();
    Ok(LazyVcf { file, records })
}

/// `VCFWriter.writeHeader` and then `add` per record: the whole file.
pub fn write_records(header: &VcfHeader, records: &[Record]) -> Result<String, EncodeError> {
    write_records_with(header, records, MissingFields::default())
}

/// [`write_records`] under the writer's `ALLOW_MISSING_FIELDS_IN_HEADER` option.
pub fn write_records_with(
    header: &VcfHeader,
    records: &[Record],
    missing: MissingFields,
) -> Result<String, EncodeError> {
    let encoders = Encoders::new(header, missing);
    let mut out = header.write();
    for record in records {
        encoders.encode(record, &mut out)?;
    }
    Ok(out)
}

/// The two encoders one record line goes through: the full one, and one over the same lines with
/// no samples for the site columns of a lazy record, which writes no FORMAT of its own because the
/// copied block supplies it.
struct Encoders {
    header: VcfHeader,
    sites_header: VcfHeader,
    missing: MissingFields,
}

impl Encoders {
    fn new(header: &VcfHeader, missing: MissingFields) -> Self {
        Self {
            missing,
            header: header.clone(),
            sites_header: VcfHeader {
                lines: header.lines.clone(),
                samples: Vec::new(),
            },
        }
    }

    /// One record's line, newline included, appended to `out`.
    fn encode(&self, record: &Record, out: &mut String) -> Result<(), EncodeError> {
        match &record.lazy_genotypes {
            Some(block) => {
                let mut site = record.variant.clone();
                site.genotypes.clear();
                VcfEncoder::new(&self.sites_header)
                    .with_missing_fields(self.missing)
                    .encode_into(&site, out)?;
                out.push('\t');
                out.push_str(block);
            }
            // The reader keeps its own copy of the block whatever the sample order, so a record
            // this module decoded is marked decoded there too, or the encoder would copy it.
            None => {
                record.variant.genotypes.decode();
                VcfEncoder::new(&self.header)
                    .with_missing_fields(self.missing)
                    .encode_into(&record.variant, out)?
            }
        }
        out.push('\n');
        Ok(())
    }
}

/// `VCFWriter` with `INDEX_ON_THE_FLY`, and the `.idx` it leaves beside the file at `close()`.
///
/// `IndexingVariantContextWriter.add` hands the indexer each record's position BEFORE the line is
/// written, so a record's block starts where its line does and the header is counted in every
/// position. Nothing here checks the order: the on-the-fly creator takes what it is given, and a
/// file out of coordinate order is written and indexed without complaint -- measured on the
/// oracle, which is why MakeSitesOnlyVcf can index an unsorted input that SortVcf exists to fix.
pub struct Indexed<'a> {
    /// `setReferenceDictionary`: one `DICT:` property per sequence, in dictionary order.
    pub dictionary: &'a [Sequence],
    /// The output path; the index records it as the absolute `file:` URI.
    pub path: &'a std::path::Path,
    /// The written file's modification time in milliseconds, which the header records beside its
    /// size. The size is the text's length; the time is the file system's, so no two runs agree
    /// on it.
    pub timestamp: i64,
}

/// The whole text, and the index when one was asked for.
pub fn write_indexed(
    header: &VcfHeader,
    records: &[Record],
    indexed: Option<Indexed>,
) -> Result<(String, Option<Vec<u8>>), EncodeError> {
    let text = write_records(header, records)?;
    let index = indexed.and_then(|indexed| tribble_index(header, records, &text, indexed));
    Ok((text, index))
}

/// The `.idx` of a written text.
fn tribble_index(
    header: &VcfHeader,
    records: &[Record],
    text: &str,
    indexed: Indexed,
) -> Option<Vec<u8>> {
    use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR, VERSION};
    use htsjdk_tribble::index_write::{BalanceApproach, BuiltIndex, DynamicIndexCreator, Feature};

    let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
    let header_length = header.write().len();
    let mut at = header_length as i64;
    for (record, line) in records
        .iter()
        .zip(text[header_length..].split_terminator('\n'))
    {
        creator.add_feature(
            &Feature {
                contig: record.variant.contig.clone(),
                start: record.variant.start as i32,
                end: record.variant.stop as i32,
            },
            at,
        );
        at += line.len() as i64 + 1;
    }
    let mut properties: Vec<(String, String)> = indexed
        .dictionary
        .iter()
        .map(|sequence| {
            (
                format!("DICT:{}", sequence.name),
                sequence.length.to_string(),
            )
        })
        .collect();
    properties.extend(creator.properties());
    // The creator starts a new per-contig index whenever the contig changes, so a contig that
    // comes back after another is two of them; the index then puts each into a `LinkedHashMap` by
    // name, where the second REPLACES the first and keeps the first's place.
    let (index_type, contigs, interval_contigs) = match creator.finalize(text.len() as i64) {
        Ok(BuiltIndex::Linear(contigs)) => (LINEAR, by_name(contigs, |c| &c.name), Vec::new()),
        Ok(BuiltIndex::IntervalTree(intervals)) => {
            (INTERVAL_TREE, Vec::new(), by_name(intervals, |c| &c.name))
        }
        // Only an empty linear index refuses, and a VCF always has a header.
        Err(_) => return None,
    };
    let absolute = std::path::absolute(indexed.path).unwrap_or_else(|_| indexed.path.into());
    let index = TribbleIndex {
        index_type,
        version: VERSION,
        indexed_path: format!("file://{}", absolute.display()),
        indexed_file_size: text.len() as i64,
        indexed_file_timestamp: indexed.timestamp,
        indexed_file_md5: String::new(),
        flags: 0,
        properties,
        contigs,
        interval_contigs,
    };
    index.write().ok()
}

/// `LinkedHashMap.put` over a list: a repeated key keeps its first position and its last value.
fn by_name<T>(items: Vec<T>, name: impl Fn(&T) -> &String) -> Vec<T> {
    let mut kept: Vec<T> = Vec::with_capacity(items.len());
    for item in items {
        match kept.iter().position(|k| name(k) == name(&item)) {
            Some(at) => kept[at] = item,
            None => kept.push(item),
        }
    }
    kept
}

/// `VariantContextWriterBuilder.determineOutputTypeFromFile`: a name in
/// `FileExtensions.COMPRESSED_VCF_LIST` is a `BLOCK_COMPRESSED_VCF`.
pub fn is_compressed_vcf(path: &str) -> bool {
    path.ends_with(".vcf.gz") || path.ends_with(".vcf.bgz")
}

/// Write the file, and its index beside it when there is one: a `.tbi` for a block-compressed
/// name, a `.idx` for any other.
pub fn write_output(
    path: &str,
    header: &VcfHeader,
    records: &[Record],
    index_dictionary: Option<&[Sequence]>,
) -> Result<(), Box<dyn std::error::Error>> {
    Ok(write_output_with(
        path,
        header,
        records,
        index_dictionary,
        MissingFields::default(),
    )?)
}

/// [`write_output`] under the writer's `ALLOW_MISSING_FIELDS_IN_HEADER` option.
pub fn write_output_with(
    path: &str,
    header: &VcfHeader,
    records: &[Record],
    index_dictionary: Option<&[Sequence]>,
    missing: MissingFields,
) -> Result<(), WriteError> {
    if is_compressed_vcf(path) {
        return write_block_compressed(path, header, records, index_dictionary, missing);
    }
    let output = std::path::Path::new(path);
    let text = write_records_with(header, records, missing)?;
    std::fs::write(output, &text)?;
    if let Some(dictionary) = index_dictionary {
        let timestamp = std::fs::metadata(output)?
            .modified()?
            .duration_since(std::time::UNIX_EPOCH)
            .map(|elapsed| elapsed.as_millis() as i64)
            .unwrap_or(0);
        let indexed = Indexed {
            dictionary,
            path: output,
            timestamp,
        };
        if let Some(index) = tribble_index(header, records, &text, indexed) {
            std::fs::write(format!("{path}.idx"), index)?;
        }
    }
    Ok(())
}

/// `Defaults.BUFFER_SIZE`, the `BufferedOutputStream` the builder puts between the BGZF stream and
/// the file.
const BUFFER_SIZE: usize = 128 * 1024;

/// A `BLOCK_COMPRESSED_VCF`: `VCFWriter` over a `BlockCompressedOutputStream` at the installed level
/// and deflater, and with `INDEX_ON_THE_FLY` a `TabixIndexCreator(refDict, TabixFormat.VCF)` whose
/// `.tbi` is written at `close()` through the same deflater.
///
/// `IndexingVariantContextWriter.add` hands the creator the BGZF virtual offset BEFORE the line is
/// written. Unlike the `.idx` creator, the tabix one checks the order: a record that goes back on
/// its contig, or a contig that comes back, is an uncaught `IllegalArgumentException`, and the JVM
/// exits with whatever the buffered stream had already passed to the file.
pub fn write_block_compressed(
    path: &str,
    header: &VcfHeader,
    records: &[Record],
    index_dictionary: Option<&[Sequence]>,
    missing: MissingFields,
) -> Result<(), WriteError> {
    use htsjdk_tribble::tabix::{FeatureRef, TabixFormat, TabixIndexCreator};
    use std::io::Write;

    let file = std::io::BufWriter::with_capacity(BUFFER_SIZE, std::fs::File::create(path)?);
    let mut bgzf = crate::output::bgzf_writer(file);
    bgzf.write_all(header.write().as_bytes())?;
    let mut creator = index_dictionary.map(|_| TabixIndexCreator::new(TabixFormat::VCF));
    let lengths: std::collections::HashMap<&str, i32> = index_dictionary
        .unwrap_or_default()
        .iter()
        .map(|sequence| (sequence.name.as_str(), sequence.length as i32))
        .collect();
    let encoders = Encoders::new(header, missing);
    let mut line = String::new();
    for record in records {
        line.clear();
        encoders.encode(record, &mut line)?;
        if let Some(creator) = creator.as_mut() {
            let variant = &record.variant;
            // `sequenceDictionary.getSequence(name).getSequenceLength()` when a new contig opens,
            // which throws on a contig the dictionary lacks.
            let sequence_length = match (index_dictionary, lengths.get(variant.contig.as_str())) {
                (Some([_, ..]), None) => die(MISSING_SEQUENCE),
                (_, length) => length.copied().unwrap_or(0),
            };
            let feature = FeatureRef {
                contig: &variant.contig,
                start: variant.start as i32,
                end: variant.stop as i32,
                // `VariantContext.toString()`, which only the contig-out-of-order refusal quotes
                // and which this port does not reproduce.
                description: &variant.contig,
                sequence_length,
            };
            if let Err(error) = creator.add_feature(feature, bgzf.file_pointer() as i64) {
                die(&format!("{}: {}", error.java_class(), error.message()));
            }
        }
        bgzf.write_all(line.as_bytes())?;
    }
    // `VCFWriter.close()` flushes before `IndexingVariantContextWriter.close()` asks for the final
    // position, and `BlockCompressedOutputStream.flush()` deflates the pending block: the last
    // chunk ends at the next block's address with offset zero, not inside the last block.
    bgzf.flush()?;
    let end = bgzf.file_pointer() as i64;
    bgzf.into_inner()?
        .into_inner()
        .map_err(|e| e.into_error())?;
    if let Some(creator) = creator {
        let index = creator
            .finish(end)
            .unwrap_or_else(|error| die(&format!("{}: {}", error.java_class(), error.message())));
        let (level, deflater) = crate::output::current();
        std::fs::write(
            format!("{path}.tbi"),
            index.write_with_level(level, deflater),
        )?;
    }
    Ok(())
}

/// The helpful `NullPointerException` the tabix creator throws for a contig its dictionary lacks.
const MISSING_SEQUENCE: &str = "java.lang.NullPointerException: Cannot invoke \
     \"htsjdk.samtools.SAMSequenceRecord.getSequenceLength()\" because the return value of \
     \"htsjdk.samtools.SAMSequenceDictionary.getSequence(String)\" is null";

/// Why a block-compressed VCF could not be written: a record the header cannot carry, or the file.
#[derive(Debug)]
pub enum WriteError {
    Encode(EncodeError),
    Io(std::io::Error),
}

impl From<EncodeError> for WriteError {
    fn from(error: EncodeError) -> Self {
        WriteError::Encode(error)
    }
}

impl From<std::io::Error> for WriteError {
    fn from(error: std::io::Error) -> Self {
        WriteError::Io(error)
    }
}

impl std::fmt::Display for WriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WriteError::Encode(error) => write!(f, "{error:?}"),
            WriteError::Io(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for WriteError {}

/// `std::fs::read_to_string` for a VCF, which inflates a gzip file whatever its name as
/// `VCFFileReader` does, and fails as `read_to_string` would on anything else.
pub fn read_to_string(path: &str) -> std::io::Result<String> {
    let bytes = std::fs::read(path)?;
    let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{e:?}")))?
    } else {
        bytes
    };
    String::from_utf8(bytes).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "stream did not contain valid UTF-8",
        )
    })
}

/// A VCF's text, inflated when the file is gzip, which `VCFFileReader` reads whatever its name.
pub fn read_text(path: &str) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("java.io.IOException: {e}"))?;
    let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .map_err(|e| format!("htsjdk.samtools.SAMException: {e:?}"))?
    } else {
        bytes
    };
    String::from_utf8(bytes).map_err(|e| format!("java.io.IOException: {e}"))
}

/// `VCFFileReader(file, false)` on a path, with a read failure as the exception it was upstream.
///
/// A gzip file is a `.vcf.gz`, which the reader inflates whatever its name.
pub fn read_path(path: &str) -> Result<LazyVcf, String> {
    let text = read_text(path)?;
    read_lazy(&text)
        .map_err(|failure| format!("{}: {}", failure.error.class(), failure.error.message()))
}

/// `FileExtensions.VCF_LIST`: what `IOUtil.unrollPaths` takes as a variant file rather than a list.
pub const VCF_LIST: [&str; 4] = [".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"];

/// `IOUtil.unrollPaths(paths, VCF_LIST)`: a path whose file name ends in a VCF extension is itself,
/// and any other is read as a list of paths, one per non-blank trimmed line, recursively.
///
/// The Java walks a stack and reverses what it collected, which is the input order, lists
/// expanded in place.
pub fn unroll_paths(inputs: &[String]) -> Result<Vec<String>, String> {
    let mut stack: Vec<String> = inputs.to_vec();
    let mut output: Vec<String> = Vec::new();
    while let Some(path) = stack.pop() {
        let name = std::path::Path::new(&path)
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        if VCF_LIST.iter().any(|ext| name.ends_with(ext)) {
            output.push(path);
            continue;
        }
        let text = std::fs::read_to_string(&path).map_err(|_| {
            let absolute = std::path::absolute(&path).unwrap_or_else(|_| path.clone().into());
            format!(
                "java.lang.IllegalArgumentException: had trouble reading from file://{}",
                absolute.display()
            )
        })?;
        stack.extend(
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty())
                .map(str::to_string),
        );
    }
    output.reverse();
    Ok(output)
}

/// `VCFHeader.addMetaDataLine` for an unstructured line: it goes in only when no line of the
/// "other" kind (anything but INFO, FORMAT, FILTER and contig) already has its key, so of two
/// comments with one key the first is the one that stays.
pub fn add_other_meta_data_line(header: &mut VcfHeader, key: &str, value: &str) {
    let taken = header.lines.iter().any(|line| match line {
        HeaderLine::Unstructured { key: k, .. } | HeaderLine::Structured { key: k, .. } => k == key,
        _ => false,
    });
    if !taken {
        header.lines.push(HeaderLine::Unstructured {
            key: key.to_string(),
            value: value.to_string(),
        });
    }
}

/// The line htsjdk's `Log.error` prints: level, time, the tool's class name and the message parts
/// run together, a throwable among them by its `toString()`.
///
/// The time is `yyyy-MM-dd HH:mm:ss` in the JVM's zone, which the oracle image leaves at UTC.
pub fn log_error(tool: &str, message: &str) {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (days, seconds) = (now.div_euclid(86_400), now.rem_euclid(86_400));
    // Howard Hinnant's `civil_from_days`.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    eprintln!(
        "ERROR\t{year:04}-{month:02}-{day:02} {:02}:{:02}:{:02}\t{tool}\t{message}",
        seconds / 3_600,
        seconds % 3_600 / 60,
        seconds % 60
    );
}

/// Print the line the JVM prints for an uncaught exception, and exit the way it does.
pub fn die(exception: &str) -> ! {
    eprintln!("Exception in thread \"main\" {exception}");
    std::process::exit(1);
}

/// One `SAMSequenceRecord`, with the one attribute a VCF contig line can carry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    /// `AS`, which `VCFContigHeaderLine` writes as `assembly`.
    pub assembly: Option<String>,
    /// `M5`, which only matters to `isSameSequence`.
    pub md5: Option<String>,
}

impl Sequence {
    /// `SAMSequenceRecord.toString`, with the index the dictionary gives it.
    pub fn describe(&self, index: usize) -> String {
        format!(
            "SAMSequenceRecord(name={},length={},dict_index={},assembly={},alternate_names=[])",
            self.name,
            self.length,
            index,
            self.assembly.as_deref().unwrap_or("null")
        )
    }
}

/// `VCFHeader.getSequenceDictionary`: one record per contig line, or `None` when there are none.
///
/// A contig line without a length reads as `UNKNOWN_SEQUENCE_LENGTH`, zero.
pub fn header_dictionary(header: &VcfHeader) -> Option<Vec<Sequence>> {
    let sequences: Vec<Sequence> = header
        .lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Contig { fields, .. } => {
                let field = |key: &str| {
                    fields
                        .iter()
                        .find(|(k, _)| k == key)
                        .map(|(_, v)| v.clone())
                };
                Some(Sequence {
                    name: field("ID").unwrap_or_default(),
                    length: field("length").and_then(|l| l.parse().ok()).unwrap_or(0),
                    assembly: field("assembly"),
                    md5: None,
                })
            }
            _ => None,
        })
        .collect();
    if sequences.is_empty() {
        None
    } else {
        Some(sequences)
    }
}

/// The `@SQ` lines of a SAM text header, which is what a `.dict` and an interval list's header are.
pub fn parse_sam_dictionary(text: &str) -> Vec<Sequence> {
    text.lines()
        .filter(|line| line.starts_with("@SQ\t"))
        .map(|line| {
            let mut sequence = Sequence {
                name: String::new(),
                length: 0,
                assembly: None,
                md5: None,
            };
            for field in line.split('\t').skip(1) {
                let Some((tag, value)) = field.split_once(':') else {
                    continue;
                };
                match tag {
                    "SN" => sequence.name = value.to_string(),
                    "LN" => sequence.length = value.parse().unwrap_or(0),
                    "AS" => sequence.assembly = Some(value.to_string()),
                    "M5" => sequence.md5 = Some(value.to_string()),
                    _ => {}
                }
            }
            sequence
        })
        .collect()
}

/// `VCFHeader.setSequenceDictionary`: every contig line goes and one per sequence comes back, at the
/// sequence's dictionary index, with `ID`, `length` and -- only when the record has one --
/// `assembly`, in that order.
pub fn set_sequence_dictionary(header: &mut VcfHeader, dictionary: &[Sequence]) {
    header
        .lines
        .retain(|line| !matches!(line, HeaderLine::Contig { .. }));
    for (index, sequence) in dictionary.iter().enumerate() {
        let mut fields = vec![
            ("ID".to_string(), sequence.name.clone()),
            ("length".to_string(), sequence.length.to_string()),
        ];
        if let Some(assembly) = &sequence.assembly {
            fields.push(("assembly".to_string(), assembly.clone()));
        }
        header.lines.push(HeaderLine::Contig {
            index: index as i32,
            fields,
        });
    }
}

/// `SAMSequenceRecord.isSameSequence`, for two records at the same position in their dictionaries.
///
/// A zero length is unknown and matches any. The MD5s decide when both carry one, and the names do
/// otherwise; the assembly is not compared at all.
pub fn is_same_sequence(this: &Sequence, that: &Sequence) -> bool {
    if this.length != 0 && that.length != 0 && this.length != that.length {
        return false;
    }
    match (&this.md5, &that.md5) {
        (Some(a), Some(b)) => {
            a.trim_start_matches('0').to_lowercase() == b.trim_start_matches('0').to_lowercase()
        }
        _ => this.name == that.name,
    }
}

/// `SAMSequenceDictionary.assertSameDictionary`, as the `AssertionError` message it throws.
pub fn assert_same_dictionary(this: &[Sequence], that: &[Sequence]) -> Result<(), String> {
    let template = |detail: String| format!("SAM dictionaries are not the same: {detail}.");
    let mut those = that.iter().enumerate();
    for (index, this_sequence) in this.iter().enumerate() {
        match those.next() {
            None => {
                return Err(template(format!(
                    "{} is present in only one dictionary",
                    this_sequence.describe(index)
                )))
            }
            Some((that_index, that_sequence)) => {
                if !is_same_sequence(that_sequence, this_sequence) {
                    return Err(template(format!(
                        "{} was found when {} was expected",
                        that_sequence.describe(that_index),
                        this_sequence.describe(index)
                    )));
                }
            }
        }
    }
    if let Some((index, sequence)) = those.next() {
        return Err(template(format!(
            "{} is present in only one dictionary",
            sequence.describe(index)
        )));
    }
    Ok(())
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`, for the shapes this repository's corpus has.
///
/// A FASTA is not read: `ReferenceSequenceFileFactory` looks for the `.dict` beside it. A VCF
/// contributes its contig lines -- and `null` when it has none, which is the one shape that comes
/// back empty-handed -- and a `.dict` or an interval list its `@SQ` lines.
pub fn extract_dictionary(path: &str) -> std::io::Result<Option<Vec<Sequence>>> {
    let candidate = std::path::Path::new(path);
    let extension = candidate
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or_default();
    match extension {
        "fasta" | "fa" | "fna" => {
            let text = std::fs::read_to_string(candidate.with_extension("dict"))?;
            Ok(Some(parse_sam_dictionary(&text)))
        }
        "vcf" => {
            let text = std::fs::read_to_string(path)?;
            let file = read_vcf(&text).map_err(|failure| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, failure.error.message())
            })?;
            Ok(header_dictionary(&file.header))
        }
        _ => Ok(Some(parse_sam_dictionary(&std::fs::read_to_string(path)?))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
        ##contig=<ID=chr1,length=2000>\n";

    fn file(samples: &str) -> String {
        format!(
            "{HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n\
             chr1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ:DP\t0/1:30:8\t0/0:.:.\n"
        )
    }

    #[test]
    fn sorted_samples_copy_the_genotype_block() {
        let vcf = read_lazy(&file("a\tb")).unwrap();
        let text = write_records(&vcf.file.header, &vcf.records).unwrap();
        assert!(text.ends_with("GT:GQ:DP\t0/1:30:8\t0/0:.:.\n"), "{text}");
    }

    #[test]
    fn unsorted_samples_are_decoded_and_re_encoded() {
        let vcf = read_lazy(&file("b\ta")).unwrap();
        assert!(vcf.records[0].lazy_genotypes.is_none());
        let text = write_records(&vcf.file.header, &vcf.records).unwrap();
        assert!(text.ends_with("GT:DP:GQ\t0/1:8:30\t0/0\n"), "{text}");
    }

    #[test]
    fn a_dictionary_mismatch_names_the_found_record_first() {
        let expected = parse_sam_dictionary("@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n");
        let found = parse_sam_dictionary("@SQ\tSN:chr1\tLN:2000\tAS:x\n@SQ\tSN:chr2\tLN:1000\n");
        assert_eq!(
            assert_same_dictionary(&expected, &found).unwrap_err(),
            "SAM dictionaries are not the same: SAMSequenceRecord(name=chr2,length=1000,\
             dict_index=1,assembly=null,alternate_names=[]) was found when \
             SAMSequenceRecord(name=chr2,length=1500,dict_index=1,assembly=null,\
             alternate_names=[]) was expected."
        );
        assert!(assert_same_dictionary(&expected[..1], &found[..1]).is_ok());
    }

    #[test]
    fn a_compressed_name_is_bgzf_with_a_tabix_index_beside_it() {
        let dir = std::env::temp_dir().join(format!("vcf_io_gz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("out.vcf.gz");
        let path = path.to_str().unwrap();
        let vcf = read_lazy(&file("a\tb")).unwrap();
        let dictionary = header_dictionary(&vcf.file.header).unwrap();
        write_output(path, &vcf.file.header, &vcf.records, Some(&dictionary)).unwrap();

        let bytes = std::fs::read(path).unwrap();
        assert!(bytes.ends_with(&htsjdk_bgzf::EMPTY_GZIP_BLOCK));
        let text = write_records(&vcf.file.header, &vcf.records).unwrap();
        assert_eq!(read_to_string(path).unwrap(), text);
        assert_eq!(read_path(path).unwrap().records, vcf.records);
        let index = std::fs::read(format!("{path}.tbi")).unwrap();
        assert!(htsjdk_bgzf::decompress_all(&index)
            .unwrap()
            .starts_with(b"TBI\x01"));
        assert!(!std::path::Path::new(&format!("{path}.idx")).exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn only_the_compressed_vcf_list_is_block_compressed() {
        assert!(is_compressed_vcf("a.vcf.gz"));
        assert!(is_compressed_vcf("a.vcf.bgz"));
        assert!(!is_compressed_vcf("a.vcf"));
        assert!(!is_compressed_vcf("a.gz"));
    }
}
