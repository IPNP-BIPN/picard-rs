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

use htsjdk_vcf::encoder::{EncodeError, VcfEncoder};
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
    let encoder = VcfEncoder::new(header);
    // The site columns of a lazy record go through the same encoder with no samples in its header,
    // so it writes no FORMAT of its own; the copied block supplies it.
    let sites_header = VcfHeader {
        lines: header.lines.clone(),
        samples: Vec::new(),
    };
    let sites_encoder = VcfEncoder::new(&sites_header);

    let mut out = header.write();
    for record in records {
        match &record.lazy_genotypes {
            Some(block) => {
                let mut site = record.variant.clone();
                site.genotypes.clear();
                sites_encoder.encode_into(&site, &mut out)?;
                out.push('\t');
                out.push_str(block);
            }
            None => encoder.encode_into(&record.variant, &mut out)?,
        }
        out.push('\n');
    }
    Ok(out)
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
    use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR, VERSION};
    use htsjdk_tribble::index_write::{BalanceApproach, BuiltIndex, DynamicIndexCreator, Feature};

    let text = write_records(header, records)?;
    let Some(indexed) = indexed else {
        return Ok((text, None));
    };
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
        Err(_) => return Ok((text, None)),
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
    Ok((text, index.write().ok()))
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

/// Write the file, and its `.idx` beside it when there is one.
pub fn write_output(
    path: &str,
    header: &VcfHeader,
    records: &[Record],
    index_dictionary: Option<&[Sequence]>,
) -> Result<(), Box<dyn std::error::Error>> {
    let output = std::path::Path::new(path);
    let text = write_records(header, records).map_err(|e| format!("{e:?}"))?;
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
        let (_, index) =
            write_indexed(header, records, Some(indexed)).map_err(|e| format!("{e:?}"))?;
        if let Some(index) = index {
            std::fs::write(format!("{path}.idx"), index)?;
        }
    }
    Ok(())
}

/// `VCFFileReader(file, false)` on a path, with a read failure as the exception it was upstream.
///
/// A gzip file is a `.vcf.gz`, which the reader inflates whatever its name.
pub fn read_path(path: &str) -> Result<LazyVcf, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("java.io.IOException: {e}"))?;
    let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .map_err(|e| format!("htsjdk.samtools.SAMException: {e:?}"))?
    } else {
        bytes
    };
    let text = String::from_utf8(bytes).map_err(|e| format!("java.io.IOException: {e}"))?;
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
}
