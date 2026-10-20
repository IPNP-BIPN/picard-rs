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
    /// `isFullyDecoded()`: true for a record the BCF2 codec made, whose values are already typed
    /// and whose attributes are already in the order the record iterates them. A BCF writer runs
    /// `fullyDecode` on any other record first, and that decides the INFO order it writes.
    pub fully_decoded: bool,
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
                fully_decoded: false,
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
    match written(header, records, missing) {
        (text, None) => Ok(text),
        (_, Some(error)) => Err(error),
    }
}

/// [`write_records_with`] up to the first refusal: the header and every whole line before it, and
/// the refusal. A refused line never leaves `VCFWriter`'s line buffer.
fn written(
    header: &VcfHeader,
    records: &[Record],
    missing: MissingFields,
) -> (String, Option<EncodeError>) {
    let encoders = Encoders::new(header, missing);
    let mut out = header.write();
    for record in records {
        let start = out.len();
        if let Err(error) = encoders.encode(record, &mut out) {
            out.truncate(start);
            return (out, Some(error));
        }
    }
    (out, None)
}

/// The exception `VCFEncoder` throws for a refusal.
pub fn encode_exception(error: &EncodeError) -> String {
    match error {
        EncodeError::MissingFromHeader {
            key,
            field,
            contig,
            start,
        } => format!(
            "java.lang.IllegalStateException: Key {key} found in VariantContext field {field} at \
             {contig}:{start} but this key isn't defined in the VCFHeader.  We require all VCFs \
             to have complete VCF headers by default."
        ),
        other => format!("{other:?}"),
    }
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
    let index = indexed.and_then(|indexed| text_index(header, &text, records, &indexed));
    Ok((text, index))
}

/// The `.idx` of `text`, already written: each line's feature at the position its bytes start at.
fn text_index(
    header: &VcfHeader,
    text: &str,
    records: &[Record],
    indexed: &Indexed,
) -> Option<Vec<u8>> {
    let header_length = header.write().len();
    let mut at = header_length as i64;
    let mut positions = Vec::with_capacity(records.len());
    for line in text[header_length..].split_terminator('\n') {
        positions.push(at);
        at += line.len() as i64 + 1;
    }
    tribble_index(features_of(records), &positions, text.len() as i64, indexed)
}

/// The `.idx` `IndexingVariantContextWriter` writes at `close()`: each record's feature at the
/// position its bytes start at, and the file's length as the end.
fn tribble_index(
    features: impl IntoIterator<Item = (String, i64, i64)>,
    positions: &[i64],
    file_length: i64,
    indexed: &Indexed,
) -> Option<Vec<u8>> {
    use htsjdk_tribble::index_write::{BalanceApproach, DynamicIndexCreator, Feature};

    let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
    for ((contig, start, end), &at) in features.into_iter().zip(positions) {
        creator.add_feature(
            &Feature {
                contig,
                start: start as i32,
                end: end as i32,
            },
            at,
        );
    }
    finish_tribble(creator, file_length, indexed)
}

/// The `.idx` of the features `creator` was given, the file being `file_length` bytes long.
fn finish_tribble(
    creator: htsjdk_tribble::index_write::DynamicIndexCreator,
    file_length: i64,
    indexed: &Indexed,
) -> Option<Vec<u8>> {
    use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR, VERSION};
    use htsjdk_tribble::index_write::BuiltIndex;

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
    let (index_type, contigs, interval_contigs) = match creator.finalize(file_length) {
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
        indexed_file_size: file_length,
        indexed_file_timestamp: indexed.timestamp,
        indexed_file_md5: String::new(),
        flags: 0,
        properties,
        contigs,
        interval_contigs,
    };
    index.write().ok()
}

/// Write a VCF the tool rendered as text, and its `.idx` beside it the way
/// `VariantContextWriterBuilder`'s default `INDEX_ON_THE_FLY` writes one for a `.vcf`: whatever
/// `CREATE_INDEX` says, given the dictionary the writer was built with. A block-compressed or BCF
/// output is written as it is, without an index here.
///
/// The features are read back off the record lines: the contig, `POS`, and the end, which is
/// `INFO`'s `END` when there is one and the last base of `REF` otherwise, as `getEnd` has it.
pub fn write_text_indexed(
    path: &str,
    text: &str,
    dictionary: Option<&[Sequence]>,
) -> std::io::Result<()> {
    std::fs::write(path, text)?;
    let lower = path.to_ascii_lowercase();
    let Some(dictionary) = dictionary else {
        return Ok(());
    };
    if [".gz", ".bgz", ".bcf"]
        .iter()
        .any(|ext| lower.ends_with(ext))
    {
        return Ok(());
    }
    let header_length = match text.find("\n#CHROM") {
        Some(at) => text[at + 1..]
            .find('\n')
            .map_or(text.len(), |end| at + 1 + end + 1),
        None => return Ok(()),
    };
    let mut positions = Vec::new();
    let mut at = header_length as i64;
    let features = text[header_length..].split_terminator('\n').map(|line| {
        positions.push(at);
        at += line.len() as i64 + 1;
        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        let start: i64 = fields.get(1).and_then(|p| p.parse().ok()).unwrap_or(0);
        let reference = fields.get(3).map_or(1, |r| r.len() as i64);
        let end = fields
            .get(7)
            .and_then(|info| {
                info.split(';')
                    .find_map(|entry| entry.strip_prefix("END="))
                    .and_then(|end| end.parse().ok())
            })
            .unwrap_or(start + reference - 1);
        (fields[0].to_string(), start, end)
    });
    let output = std::path::Path::new(path);
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
    let features: Vec<_> = features.collect();
    if let Some(index) = tribble_index(features, &positions, text.len() as i64, &indexed) {
        std::fs::write(format!("{path}.idx"), index)?;
    }
    Ok(())
}

/// Each record's feature: its contig, start and end.
fn features_of(records: &[Record]) -> impl Iterator<Item = (String, i64, i64)> + '_ {
    records.iter().map(|record| {
        (
            record.variant.contig.clone(),
            record.variant.start,
            record.variant.stop,
        )
    })
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

/// [`write_output`] under the writer's `ALLOW_MISSING_FIELDS_IN_HEADER` option: every record through
/// one [`VcfWriter`]. What the writer refuses is the JVM's exception and exit, with what the
/// buffered stream had flushed left on disk; an I/O failure is returned.
pub fn write_output_with(
    path: &str,
    header: &VcfHeader,
    records: &[Record],
    index_dictionary: Option<&[Sequence]>,
    missing: MissingFields,
) -> Result<(), WriteError> {
    let mut writer = VcfWriter::create(path, header, index_dictionary, missing).or_else(thrown)?;
    for record in records {
        writer.add(record).or_else(thrown)?;
    }
    writer.finish().or_else(thrown)
}

/// An I/O failure as it is, and anything the JVM would have thrown, thrown.
pub fn thrown<T>(error: WriteError) -> Result<T, WriteError> {
    match error {
        WriteError::Io(_) | WriteError::Other(_) => Err(error),
        other => die(&other.exception()),
    }
}

/// `Defaults.BUFFER_SIZE`, the `BufferedOutputStream` the builder puts between the BGZF stream and
/// the file.
const BUFFER_SIZE: usize = 128 * 1024;

/// What `VariantContextWriterBuilder` builds for a name, one record at a time: `BCF2Writer` for a
/// `.bcf`, `VCFWriter` over BGZF for a block-compressed name, `VCFWriter` over text for any other,
/// each with `INDEX_ON_THE_FLY` when it is given a dictionary.
///
/// `IndexingVariantContextWriter.add` hands the indexer each record's position BEFORE the record is
/// encoded and written, so a record's block starts where its line does and the header is counted
/// in every position. The `.idx` creator takes features in any order: a file out of coordinate
/// order is written and indexed without complaint -- measured on the oracle, which is why
/// MakeSitesOnlyVcf can index an unsorted input that SortVcf exists to fix. The tabix creator a
/// block-compressed name gets checks the order: a record that goes back on its contig, or a contig
/// that comes back, is an uncaught `IllegalArgumentException`.
///
/// A refusal is returned rather than thrown, so that a caller still reading its input can let a
/// read failure further on win, as it does when the whole input is read first. Nothing is flushed
/// after one: the file holds what the stream had already passed on.
pub struct VcfWriter {
    path: String,
    sink: Sink,
    index: Index,
    encoders: Encoders,
    line: String,
}

enum Sink {
    /// The byte count, which is a text line's position.
    Text { out: JavaBuffered, at: u64 },
    /// `VCFWriter` over a `BlockCompressedOutputStream` at the installed level and deflater.
    Bgzf(htsjdk_bgzf::BgzfWriter<std::io::BufWriter<std::fs::File>>),
    Bcf {
        out: JavaBuffered,
        writer: Box<htsjdk_vcf::bcf2::writer::Bcf2Writer>,
        at: u64,
    },
}

enum Index {
    None,
    /// The `.idx` of a `.vcf` or `.bcf`, and the dictionary its `DICT:` properties come from.
    Tribble {
        creator: htsjdk_tribble::index_write::DynamicIndexCreator,
        dictionary: Vec<Sequence>,
    },
    /// `TabixIndexCreator(refDict, TabixFormat.VCF)`, written as `.tbi` at `close()` through the
    /// same deflater.
    Tabix {
        creator: htsjdk_tribble::tabix::TabixIndexCreator,
        dictionary: Vec<Sequence>,
        lengths: std::collections::HashMap<String, i32>,
    },
}

impl VcfWriter {
    /// The writer, its header written. The builder opens the file before the writer exists, so a
    /// header `BCF2Writer` refuses leaves an empty file behind.
    pub fn create(
        path: &str,
        header: &VcfHeader,
        index_dictionary: Option<&[Sequence]>,
        missing: MissingFields,
    ) -> Result<Self, WriteError> {
        use htsjdk_tribble::index_write::{BalanceApproach, DynamicIndexCreator};
        use htsjdk_tribble::tabix::{TabixFormat, TabixIndexCreator};

        let bcf = is_bcf(path);
        let index = match index_dictionary {
            None => Index::None,
            Some(dictionary) if !bcf && is_compressed_vcf(path) => Index::Tabix {
                creator: TabixIndexCreator::new(TabixFormat::VCF),
                lengths: dictionary
                    .iter()
                    .map(|sequence| (sequence.name.clone(), sequence.length as i32))
                    .collect(),
                dictionary: dictionary.to_vec(),
            },
            Some(dictionary) => Index::Tribble {
                creator: DynamicIndexCreator::new(BalanceApproach::ForSeekTime),
                dictionary: dictionary.to_vec(),
            },
        };
        let file = std::fs::File::create(path)?;
        let sink = if bcf {
            let writer = htsjdk_vcf::bcf2::writer::Bcf2Writer::new(header, false)
                .map_err(|e| WriteError::Exception(e.to_string()))?;
            let bytes = writer.header_bytes();
            let mut out = JavaBuffered::new(file);
            // `BCFVersion.write` is the magic as one call and each version byte as another, and
            // `BCF2Type.INT32.write` is four single bytes.
            let mut at = 0;
            for len in [3, 1, 1, 1, 1, 1, 1, bytes.len() - 9] {
                out.write(&bytes[at..at + len])?;
                at += len;
            }
            Sink::Bcf {
                out,
                writer: Box::new(writer),
                at: at as u64,
            }
        } else if is_compressed_vcf(path) {
            let mut bgzf =
                crate::output::bgzf_writer(std::io::BufWriter::with_capacity(BUFFER_SIZE, file));
            std::io::Write::write_all(&mut bgzf, header.write().as_bytes())?;
            Sink::Bgzf(bgzf)
        } else {
            let text = header.write();
            let mut out = JavaBuffered::new(file);
            out.write(text.as_bytes())?;
            Sink::Text {
                out,
                at: text.len() as u64,
            }
        };
        Ok(VcfWriter {
            path: path.to_string(),
            sink,
            index,
            encoders: Encoders::new(header, missing),
            line: String::new(),
        })
    }

    /// `IndexingVariantContextWriter.add`: the indexer, then the writer.
    pub fn add(&mut self, record: &Record) -> Result<(), WriteError> {
        let at = self.position()?;
        self.index(record, at)?;
        self.write(record)
    }

    /// Where the next record starts: a byte count, or the BGZF virtual offset.
    pub fn position(&mut self) -> Result<u64, WriteError> {
        Ok(match &mut self.sink {
            Sink::Text { at, .. } | Sink::Bcf { at, .. } => *at,
            Sink::Bgzf(bgzf) => bgzf.file_pointer()?,
        })
    }

    /// This writer's indexer given `record` at `position`, which may be another file's: a
    /// builder can hand one tabix creator to two writers.
    pub fn index(&mut self, record: &Record, position: u64) -> Result<(), WriteError> {
        let variant = &record.variant;
        match &mut self.index {
            Index::None => {}
            Index::Tribble { creator, .. } => creator.add_feature(
                &htsjdk_tribble::index_write::Feature {
                    contig: variant.contig.clone(),
                    start: variant.start as i32,
                    end: variant.stop as i32,
                },
                position as i64,
            ),
            Index::Tabix {
                creator,
                dictionary,
                lengths,
            } => {
                // `sequenceDictionary.getSequence(name).getSequenceLength()` when a new contig
                // opens, which throws on a contig the dictionary lacks.
                let sequence_length = match lengths.get(&variant.contig) {
                    None if !dictionary.is_empty() => {
                        return Err(WriteError::Exception(MISSING_SEQUENCE.to_string()))
                    }
                    length => length.copied().unwrap_or(0),
                };
                let feature = htsjdk_tribble::tabix::FeatureRef {
                    contig: &variant.contig,
                    start: variant.start as i32,
                    end: variant.stop as i32,
                    // `VariantContext.toString()`, which only the contig-out-of-order refusal
                    // quotes and which this port does not reproduce.
                    description: &variant.contig,
                    sequence_length,
                };
                creator
                    .add_feature(feature, position as i64)
                    .map_err(|error| {
                        WriteError::Exception(format!(
                            "{}: {}",
                            error.java_class(),
                            error.message()
                        ))
                    })?;
            }
        }
        Ok(())
    }

    /// `VCFWriter.add` or `BCF2Writer.add` once the indexer has the position: the record encoded,
    /// or the encoder's refusal, and its bytes written.
    pub fn write(&mut self, record: &Record) -> Result<(), WriteError> {
        match &mut self.sink {
            Sink::Text { out, at } => {
                self.line.clear();
                self.encoders.encode(record, &mut self.line)?;
                out.write(self.line.as_bytes())?;
                *at += self.line.len() as u64;
            }
            Sink::Bgzf(bgzf) => {
                self.line.clear();
                self.encoders.encode(record, &mut self.line)?;
                std::io::Write::write_all(bgzf, self.line.as_bytes())?;
            }
            Sink::Bcf { out, writer, at } => {
                use htsjdk_vcf::bcf2::writer::Origin;
                let origin = if record.fully_decoded {
                    Origin::Bcf2Codec
                } else {
                    Origin::VcfCodec
                };
                let encoded = writer
                    .encode(&record.variant, origin)
                    .map_err(|e| WriteError::Exception(e.to_string()))?;
                // `writeBlock`: two INT32 lengths, then the site block and the genotype block.
                let sites = u32::from_le_bytes(encoded[..4].try_into().expect("a length")) as usize;
                let mut from = 0;
                for len in [1, 1, 1, 1, 1, 1, 1, 1, sites, encoded.len() - 8 - sites] {
                    out.write(&encoded[from..from + len])?;
                    from += len;
                }
                *at += encoded.len() as u64;
            }
        }
        Ok(())
    }

    /// `close()`: the file, then its index beside it.
    pub fn finish(self) -> Result<(), WriteError> {
        let output = std::path::Path::new(&self.path);
        let length = match self.sink {
            Sink::Text { out, at } | Sink::Bcf { out, at, .. } => {
                out.finish()?;
                at as i64
            }
            Sink::Bgzf(mut bgzf) => {
                // `VCFWriter.close()` flushes before `IndexingVariantContextWriter.close()` asks
                // for the final position, and `BlockCompressedOutputStream.flush()` deflates the
                // pending block: the last chunk ends at the next block's address with offset
                // zero, not inside the last block.
                std::io::Write::flush(&mut bgzf)?;
                let end = bgzf.file_pointer()? as i64;
                bgzf.into_inner()?
                    .into_inner()
                    .map_err(|e| e.into_error())?;
                end
            }
        };
        match self.index {
            Index::None => {}
            Index::Tribble {
                creator,
                dictionary,
            } => {
                let indexed = Indexed {
                    dictionary: &dictionary,
                    path: output,
                    timestamp: modified_millis(output)?,
                };
                if let Some(index) = finish_tribble(creator, length, &indexed) {
                    std::fs::write(format!("{}.idx", self.path), index)?;
                }
            }
            Index::Tabix { creator, .. } => {
                let index = creator.finish(length).map_err(|error| {
                    WriteError::Exception(format!("{}: {}", error.java_class(), error.message()))
                })?;
                let (level, deflater) = crate::output::current();
                std::fs::write(
                    format!("{}.tbi", self.path),
                    index.write_with_level(level, deflater),
                )?;
            }
        }
        Ok(())
    }
}

/// What reaches the file through `IOUtil.maybeBufferOutputStream(out, Defaults.BUFFER_SIZE)`,
/// counted per `write` call as the writer makes them: a write that does not fit flushes what is
/// buffered first, and one at least as large as the buffer goes straight through. Nothing else
/// flushes before `close()`, which a run that throws never reaches, so a refusal leaves exactly
/// `flushed` bytes on disk -- none on a small file, whole buffers' worth on a large one.
struct JavaBuffered {
    file: std::fs::File,
    buffer: Vec<u8>,
}

impl JavaBuffered {
    /// `samjdk.buffer_size`'s default.
    const SIZE: usize = 128 * 1024;

    fn new(file: std::fs::File) -> Self {
        JavaBuffered {
            file,
            buffer: Vec::with_capacity(Self::SIZE),
        }
    }

    fn write(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        use std::io::Write;
        if bytes.len() > Self::SIZE - self.buffer.len() {
            self.file.write_all(&self.buffer)?;
            self.buffer.clear();
        }
        if bytes.len() >= Self::SIZE {
            self.file.write_all(bytes)
        } else {
            self.buffer.extend_from_slice(bytes);
            Ok(())
        }
    }

    fn finish(mut self) -> std::io::Result<()> {
        std::io::Write::write_all(&mut self.file, &self.buffer)
    }
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
    /// What the BCF2 writer refused with.
    Other(Box<dyn std::error::Error>),
    /// Anything else the writer or its indexer throws, as the JVM prints it.
    Exception(String),
}

impl WriteError {
    /// The exception the JVM prints for it.
    pub fn exception(&self) -> String {
        match self {
            WriteError::Encode(error) => encode_exception(error),
            WriteError::Exception(exception) => exception.clone(),
            other => other.to_string(),
        }
    }
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
            WriteError::Other(error) => write!(f, "{error}"),
            WriteError::Exception(exception) => write!(f, "{exception}"),
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
        crate::bgzf::decompress_owned(bytes)
            .map_err(|e| format!("htsjdk.samtools.SAMException: {e:?}"))?
    } else {
        bytes
    };
    String::from_utf8(bytes).map_err(|e| format!("java.io.IOException: {e}"))
}

/// `VCFFileReader.isBCF` and `VariantContextWriterBuilder.isBCF`: the name decides, not the bytes.
pub fn is_bcf(path: &str) -> bool {
    std::path::Path::new(path)
        .file_name()
        .is_some_and(|name| name.to_string_lossy().ends_with(".bcf"))
}

/// The modification time `IndexingVariantContextWriter` records, in milliseconds.
fn modified_millis(path: &std::path::Path) -> std::io::Result<i64> {
    Ok(std::fs::metadata(path)?
        .modified()?
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_millis() as i64)
        .unwrap_or(0))
}

/// `VCFFileReader` on a `.bcf`: the BCF2 codec, every record fully decoded.
pub fn read_bcf_path(path: &str) -> Result<LazyVcf, String> {
    let file = read_bcf_file(path)?;
    let records = file
        .records
        .iter()
        .map(|variant| Record {
            variant: variant.clone(),
            lazy_genotypes: None,
            fully_decoded: true,
        })
        .collect();
    Ok(LazyVcf { file, records })
}

/// The BCF2 codec's whole file, with a refusal as `VCFFileReader` surfaces it.
fn read_bcf_file(path: &str) -> Result<VcfFile, String> {
    use htsjdk_vcf::bcf2::codec::read_bcf;

    let bytes = std::fs::read(path).map_err(|e| format!("java.io.IOException: {e}"))?;
    let absolute = std::path::absolute(path).unwrap_or_else(|_| path.into());
    let source = format!("file://{}", absolute.display());
    let file = read_bcf(&bytes).map_err(|failure| failure.reader_exception(&source).to_string())?;
    Ok(VcfFile {
        header: file.header,
        codec_version: file.codec_version,
        header_version: file.header_version,
        records: file.records,
        skipped: Vec::new(),
    })
}

/// The whole file as [`read_vcf`] gives it, for a tool that reads the records and writes none:
/// a `.bcf` through the BCF2 codec, anything else as VCF text.
pub fn read_vcf_file(path: &str) -> Result<VcfFile, Box<dyn std::error::Error>> {
    if is_bcf(path) {
        return Ok(read_bcf_file(path)?);
    }
    Ok(read_vcf(&read_to_string(path)?).map_err(|e| format!("{e:?}"))?)
}

/// `VCFFileReader(file, false)` on a path, with a read failure as the exception it was upstream.
///
/// A gzip file is a `.vcf.gz`, which the reader inflates whatever its name, and a `.bcf` is read
/// by the BCF2 codec whatever its contents.
pub fn read_path(path: &str) -> Result<LazyVcf, String> {
    if is_bcf(path) {
        return read_bcf_path(path);
    }
    let text = read_text(path)?;
    read_lazy(&text)
        .map_err(|failure| format!("{}: {}", failure.error.class(), failure.error.message()))
}

/// [`read_path`] a record at a time: the header when it is opened, then each record as it is read,
/// with `lazy_genotypes` as [`read_lazy`] keeps it.
///
/// A tool that streams still answers a bad input with the exception the whole-file read gave,
/// which is the first of: the gzip stream's failure anywhere in the file, the first byte that is
/// not UTF-8, the header's refusal, the first record's refusal. So any failure, the reader's own
/// or one the tool raises on the way through ([`VcfReader::fail`]), first reads the rest of the
/// input for a failure that outranks it.
///
/// shortcut: a `.bcf` input is still read whole by the BCF2 codec, then handed over a record at a
/// time; stream it when a BCF input that size turns up.
pub struct VcfReader {
    /// The header and versions; `records` stays empty.
    pub file: VcfFile,
    source: Source,
    /// `WFIterator.currentRecord`: the record after the one last handed out, read as that one was,
    /// once the first has been asked for. `Some(None)` once the input is done.
    ahead: Option<Option<Result<Record, String>>>,
    /// The line of the record last handed out, which [`VcfReader::keep`] keeps; the source's own
    /// line is the one read ahead.
    handed: Vec<u8>,
}

enum Source {
    Bcf(std::vec::IntoIter<Record>),
    Text(Box<TextSource>),
}

struct TextSource {
    lines: Lines,
    /// Bytes read so far, which is where a UTF-8 failure is counted from.
    offset: u64,
    /// `lineNo` as the codec counts it.
    line_number: usize,
    samples: std::sync::Arc<[String]>,
    /// `samplesWereAlreadySorted`, which decides whether the codec decodes each record.
    sorted: bool,
    /// Whether [`read_lazy`] keeps the genotype block.
    lazy: bool,
    line: Vec<u8>,
    /// Set once a failure was reported: nothing further is read.
    done: bool,
}

/// The input's lines, `\n` included, inflated when the file is gzip.
enum Lines {
    Plain(std::io::BufReader<std::fs::File>),
    Bgzf {
        reader: htsjdk_bgzf::BgzfReader<std::io::BufReader<std::fs::File>>,
        block: Vec<u8>,
        at: usize,
    },
}

/// A failure, by what outranks what in a whole-file read.
enum Failure {
    /// The file or its gzip stream: nothing is read after it.
    Bytes(String),
    Utf8(String),
    Parse(String),
}

impl Lines {
    /// The next line appended to `line`, or false at the end of the input.
    fn read(&mut self, line: &mut Vec<u8>) -> Result<bool, String> {
        use std::io::BufRead;
        match self {
            Lines::Plain(reader) => reader
                .read_until(b'\n', line)
                .map(|n| n > 0)
                .map_err(|e| format!("java.io.IOException: {e}")),
            Lines::Bgzf { reader, block, at } => loop {
                let rest = &block[*at..];
                if let Some(end) = rest.iter().position(|&b| b == b'\n') {
                    line.extend_from_slice(&rest[..=end]);
                    *at += end + 1;
                    return Ok(true);
                }
                line.extend_from_slice(rest);
                match reader
                    .next_block()
                    .map_err(|e| format!("htsjdk.samtools.SAMException: {e:?}"))?
                {
                    Some(next) => {
                        *block = next.data;
                        *at = 0;
                    }
                    None => {
                        block.clear();
                        *at = 0;
                        return Ok(!line.is_empty());
                    }
                }
            },
        }
    }
}

/// `String.from_utf8` over the whole file, for one line `offset` bytes in.
fn utf8_at(line: &[u8], offset: u64) -> Result<&str, String> {
    std::str::from_utf8(line).map_err(|e| {
        let at = offset + e.valid_up_to() as u64;
        match e.error_len() {
            Some(len) => {
                format!(
                    "java.io.IOException: invalid utf-8 sequence of {len} bytes from index {at}"
                )
            }
            None => format!("java.io.IOException: incomplete utf-8 byte sequence from index {at}"),
        }
    })
}

/// `str::lines`' line: without its `\n`, and without a `\r` before that.
fn without_terminator(line: &str) -> &str {
    match line.strip_suffix('\n') {
        Some(line) => line.strip_suffix('\r').unwrap_or(line),
        None => line,
    }
}

impl TextSource {
    /// The exception the whole-file read would have thrown, given the first failure met here:
    /// the rest of the input is read for one that outranks it.
    fn resolve(&mut self, first: Failure) -> String {
        self.done = true;
        let (mut utf8, parse) = match first {
            Failure::Bytes(exception) => return exception,
            Failure::Utf8(exception) => (Some(exception), None),
            Failure::Parse(exception) => (None, Some(exception)),
        };
        loop {
            self.line.clear();
            match self.lines.read(&mut self.line) {
                Err(exception) => return exception,
                Ok(false) => break,
                Ok(true) => {}
            }
            if utf8.is_none() {
                utf8 = utf8_at(&self.line, self.offset).err();
            }
            self.offset += self.line.len() as u64;
        }
        utf8.or(parse).expect("a failure")
    }

    /// One line of the input `offset` bytes in, as the codec decodes it: `None` for a `#` line.
    fn decode(
        &mut self,
        file: &VcfFile,
        buffer: &[u8],
        offset: u64,
    ) -> Result<Option<Record>, Failure> {
        use htsjdk_vcf::genotype_parse::{parse_genotypes, GenotypeContext};
        use htsjdk_vcf::genotypes_context::GenotypesContext;
        use htsjdk_vcf::record_parse::{decode_line, split_condensed_borrowed};

        let text = utf8_at(buffer, offset).map_err(Failure::Utf8)?;
        let line = without_terminator(text);
        let decoded = decode_line(line, &file.header, self.line_number, file.codec_version)
            .map_err(|e| Failure::Parse(format!("{}: {}", e.class(), e.message())))?;
        let Some(decoded) = decoded else {
            // A `#` line: no record, no refusal, and no increment of the counter.
            return Ok(None);
        };
        // `parseVCFLine` increments on entry.
        self.line_number += 1;
        let mut variant = decoded.variant;
        if let Some(block) = decoded.genotype_block {
            let context = GenotypeContext {
                line,
                header: &file.header,
                version: file.codec_version,
                contig: &variant.contig,
                pos: variant.start,
                line_number: self.line_number,
            };
            let genotypes = parse_genotypes(&block, &variant.alleles, &context)
                .map_err(|e| Failure::Parse(format!("{}: {}", e.class(), e.message())))?;
            variant.genotypes = GenotypesContext::lazy(genotypes, block, self.samples.clone());
            if !self.sorted {
                variant.genotypes.decode();
            }
        }
        let lazy_genotypes = self
            .lazy
            .then(|| {
                split_condensed_borrowed(line, '\t', NUM_STANDARD_FIELDS + 1, true)
                    .into_iter()
                    .nth(NUM_STANDARD_FIELDS)
                    .map(str::to_string)
            })
            .flatten();
        Ok(Some(Record {
            variant,
            lazy_genotypes,
            fully_decoded: false,
        }))
    }

    fn next(&mut self, file: &VcfFile) -> Option<Result<Record, String>> {
        if self.done {
            return None;
        }
        loop {
            let mut buffer = std::mem::take(&mut self.line);
            buffer.clear();
            match self.lines.read(&mut buffer) {
                Err(exception) => return Some(Err(self.resolve(Failure::Bytes(exception)))),
                Ok(false) => {
                    self.done = true;
                    return None;
                }
                Ok(true) => {}
            }
            let offset = self.offset;
            self.offset += buffer.len() as u64;
            let decoded = self.decode(file, &buffer, offset);
            self.line = buffer;
            match decoded {
                Ok(Some(record)) => return Some(Ok(record)),
                Ok(None) => continue,
                Err(failure) => return Some(Err(self.resolve(failure))),
            }
        }
    }
}

impl VcfReader {
    /// `VCFFileReader(file, false)` up to its header: the file opened and the header read, or the
    /// exception the whole-file read gave.
    pub fn open(path: &str) -> Result<VcfReader, String> {
        use std::io::BufRead;

        if is_bcf(path) {
            let LazyVcf { mut file, records } = read_bcf_path(path)?;
            file.records = Vec::new();
            return Ok(VcfReader {
                file,
                source: Source::Bcf(records.into_iter()),
                ahead: None,
                handed: Vec::new(),
            });
        }
        let mut reader = std::io::BufReader::with_capacity(
            BUFFER_SIZE,
            std::fs::File::open(path).map_err(|e| format!("java.io.IOException: {e}"))?,
        );
        let gzip = reader
            .fill_buf()
            .map_err(|e| format!("java.io.IOException: {e}"))?
            .starts_with(&[0x1f, 0x8b]);
        let lines = if gzip {
            Lines::Bgzf {
                reader: htsjdk_bgzf::BgzfReader::new(reader),
                block: Vec::new(),
                at: 0,
            }
        } else {
            Lines::Plain(reader)
        };
        let mut source = TextSource {
            lines,
            offset: 0,
            line_number: 0,
            samples: std::sync::Arc::from(Vec::new()),
            sorted: true,
            lazy: false,
            line: Vec::new(),
            done: false,
        };

        // The header is every `##` line and the one after them, which `read_header_frame` takes
        // as `#CHROM` or refuses.
        let mut header = Vec::new();
        loop {
            let start = header.len();
            match source.lines.read(&mut header) {
                Err(exception) => return Err(source.resolve(Failure::Bytes(exception))),
                Ok(false) => break,
                Ok(true) => {}
            }
            source.line_number += 1;
            if !header[start..].starts_with(b"##") {
                break;
            }
        }
        source.offset = header.len() as u64;
        let text = match utf8_at(&header, 0) {
            Ok(text) => text,
            Err(exception) => return Err(source.resolve(Failure::Utf8(exception))),
        };
        let file = match read_vcf(text) {
            Ok(file) => file,
            Err(failure) => {
                return Err(source.resolve(Failure::Parse(format!(
                    "{}: {}",
                    failure.error.class(),
                    failure.error.message()
                ))))
            }
        };
        let samples = &file.header.samples;
        source.samples = std::sync::Arc::from(samples.clone());
        source.sorted = samples.windows(2).all(|pair| {
            htsjdk_vcf::genotype_type::compare_sample_names(&pair[0], &pair[1])
                != std::cmp::Ordering::Greater
        });
        source.lazy = is_sorted(samples) && !samples.is_empty();
        Ok(VcfReader {
            file,
            source: Source::Text(Box::new(source)),
            ahead: None,
            handed: Vec::new(),
        })
    }

    /// The next record, or the exception the whole-file read gave.
    ///
    /// `WFIterator` reads the first record when the iterator is made and, in each `next()`, the
    /// one after before it returns: a record that fails to decode throws as the one before it is
    /// handed out, which is then lost.
    pub fn try_next(&mut self) -> Option<Result<Record, String>> {
        match &mut self.source {
            Source::Bcf(records) => records.next().map(Ok),
            Source::Text(source) => {
                let current = match self.ahead.take() {
                    Some(current) => current,
                    None => source.next(&self.file),
                };
                let record = match current {
                    Some(Ok(record)) => record,
                    other => {
                        self.ahead = Some(None);
                        return other;
                    }
                };
                self.handed = source.line.clone();
                match source.next(&self.file) {
                    Some(Err(e)) => {
                        self.ahead = Some(None);
                        Some(Err(e))
                    }
                    next => {
                        self.ahead = Some(next);
                        Some(Ok(record))
                    }
                }
            }
        }
    }

    /// The rest of the input read, and the exception the whole-file read gave if there is one.
    pub fn check(&mut self) -> Result<(), String> {
        while let Some(record) = self.try_next() {
            record?;
        }
        Ok(())
    }

    /// Throw `exception` as a tool that had read its whole input first would: unless the rest of
    /// the input fails to read, which it would have thrown instead.
    pub fn fail(&mut self, exception: &str) -> ! {
        if let Err(read) = self.check() {
            die(&read);
        }
        die(exception)
    }

    /// [`VcfReader::fail`] for what a writer refused: its exception, or an I/O failure as a
    /// `main` returning it prints it.
    pub fn fail_write(&mut self, error: WriteError) -> ! {
        if let Err(read) = self.check() {
            die(&read);
        }
        match error {
            WriteError::Io(_) | WriteError::Other(_) => {
                eprintln!("Error: {error:?}");
                std::process::exit(1)
            }
            refused => die(&refused.exception()),
        }
    }

    /// The record [`VcfReader::try_next`] just returned, kept to be decoded again by
    /// [`VcfReader::restore`]: a text record as its line, a BCF record as it is.
    pub fn keep(&self, record: Record) -> Kept {
        match &self.source {
            Source::Bcf(_) => Kept::Record(Box::new(record)),
            Source::Text(_) => Kept::Line(self.handed.clone().into_boxed_slice()),
        }
    }

    /// A record [`VcfReader::keep`] kept, decoded again as it was the first time.
    pub fn restore(&mut self, kept: Kept) -> Record {
        match (kept, &mut self.source) {
            (Kept::Record(record), _) => *record,
            (Kept::Line(line), Source::Text(source)) => {
                // The counter only names a failing line, and this line has decoded once already.
                let line_number = source.line_number;
                let decoded = source.decode(&self.file, &line, 0);
                source.line_number = line_number;
                match decoded {
                    Ok(Some(record)) => record,
                    _ => unreachable!("a kept line decodes as it did"),
                }
            }
            (Kept::Line(_), Source::Bcf(_)) => unreachable!("a BCF record is kept whole"),
        }
    }

    /// Every record, through `map`, into a writer at `path` with `header`: one record in memory at
    /// a time. What `map` or the writer refuses is thrown through [`VcfReader::fail`].
    pub fn copy_to(
        &mut self,
        path: &str,
        header: &VcfHeader,
        index_dictionary: Option<&[Sequence]>,
        mut map: impl FnMut(Record) -> Result<Record, String>,
    ) {
        let mut writer =
            VcfWriter::create(path, header, index_dictionary, MissingFields::default())
                .unwrap_or_else(|error| self.fail_write(error));
        while let Some(record) = self.next() {
            let record = map(record).unwrap_or_else(|exception| self.fail(&exception));
            writer
                .add(&record)
                .unwrap_or_else(|error| self.fail_write(error));
        }
        writer
            .finish()
            .unwrap_or_else(|error| self.fail_write(error));
    }
}

/// A record held for later at a fraction of its decoded size: see [`VcfReader::keep`].
pub enum Kept {
    Line(Box<[u8]>),
    Record(Box<Record>),
}

/// The records, dying with the whole-file read's exception at the first that fails.
impl Iterator for VcfReader {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        self.try_next()
            .map(|record| record.unwrap_or_else(|exception| die(&exception)))
    }
}

/// [`VcfReader::fail`] for a tool that read several inputs whole, in order, before `exception`:
/// the first of them that fails to read throws instead. A reader's own failure is `exception`
/// with `files` the inputs before it.
pub fn fail_all(files: &mut [VcfReader], exception: &str) -> ! {
    for file in files {
        if let Err(read) = file.check() {
            die(&read);
        }
    }
    die(exception)
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

    /// What [`VcfReader`] hands over, record by record and then the failure, against [`read_path`].
    fn streamed(path: &str) -> Result<(VcfFile, Vec<Record>), String> {
        let mut reader = VcfReader::open(path)?;
        let mut records = Vec::new();
        while let Some(record) = reader.try_next() {
            records.push(record?);
        }
        Ok((reader.file, records))
    }

    #[test]
    fn the_stream_reads_what_the_whole_file_read_reads_and_fails_as_it_fails() {
        let dir = std::env::temp_dir().join(format!("vcf_io_stream_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let record = "chr1\t200\t.\tC\tT\t.\t.\t.\tGT:GQ:DP\t0/0:.:.\t1/1:9:3\n";
        let cases = [
            ("sorted", format!("{}{record}", file("a\tb")).into_bytes()),
            (
                "unsorted",
                format!("{}# x\r\n{record}", file("b\ta")).into_bytes(),
            ),
            ("no-newline", file("a\tb").trim_end().as_bytes().to_vec()),
            (
                "bad-record",
                format!("{}chr1\t5\n{record}", file("a\tb")).into_bytes(),
            ),
            ("no-chrom", HEADER.as_bytes().to_vec()),
            // A record refused early, a byte that is not UTF-8 later: the UTF-8 failure wins.
            (
                "utf8-after-record",
                [
                    format!("{}chr1\t5\n", file("a\tb")).as_bytes(),
                    b"chr1\t9\t\xff\n",
                ]
                .concat(),
            ),
            (
                "utf8-in-header",
                [
                    b"##fileformat=VCFv4.2\n##x=\xe2\n".as_slice(),
                    file("a\tb").as_bytes(),
                ]
                .concat(),
            ),
        ];
        for (name, bytes) in cases {
            for gz in [false, true] {
                let path = dir.join(format!("{name}.vcf{}", if gz { ".gz" } else { "" }));
                let bytes = if gz {
                    let mut bgzf = htsjdk_bgzf::BgzfWriter::new(Vec::new());
                    std::io::Write::write_all(&mut bgzf, &bytes).unwrap();
                    bgzf.into_inner().unwrap()
                } else {
                    bytes.clone()
                };
                std::fs::write(&path, bytes).unwrap();
                let path = path.to_str().unwrap();
                let whole = read_path(path).map(|vcf| {
                    let mut file = vcf.file;
                    file.records.clear();
                    file.skipped.clear();
                    (file, vcf.records)
                });
                assert_eq!(streamed(path), whole, "{name} gz={gz}");
            }
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_refused_record_leaves_what_the_buffered_stream_had_flushed() {
        let dir = std::env::temp_dir().join(format!("vcf_io_refused_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("out.vcf");
        let path = path.to_str().unwrap();
        let vcf = read_lazy(&file("a\tb")).unwrap();
        let mut writer =
            VcfWriter::create(path, &vcf.file.header, None, MissingFields::default()).unwrap();
        // Whole lines past one buffer's worth: the first buffer reaches the file, the rest waits.
        let line = write_records(&vcf.file.header, &vcf.records).unwrap()
            [vcf.file.header.write().len()..]
            .len();
        let lines = JavaBuffered::SIZE / line + 2;
        for _ in 0..lines {
            writer.add(&vcf.records[0]).unwrap();
        }
        let flushed = std::fs::metadata(path).unwrap().len() as usize;
        let header = vcf.file.header.write().len();
        assert_eq!(
            flushed,
            header + (JavaBuffered::SIZE - header) / line * line
        );
        writer.finish().unwrap();
        assert_eq!(
            std::fs::metadata(path).unwrap().len() as usize,
            header + lines * line
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
