//! `VcfFormatConverter` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.VcfFormatConverter.doWork` at tag 3.4.0 for VCF and block-compressed VCF,
//! read and written. BCF is not ported: htsjdk-rs has no BCF codec, so a `.bcf` on either side is
//! refused here with a message of the port's own.
//!
//! The records are copied under a copy of the input's header, so the conversion is the shared
//! round trip (`picard_analysis::vcf_io`) and what varies is the container. The output's format
//! is the writer builder's guess from the name: `.vcf.gz` and `.vcf.bgz` are BGZF, `.vcf` is text,
//! and any other name is the builder's `IllegalArgumentException`.
//!
//! # `REQUIRE_INDEX` looks for the index the reader would use
//!
//! `AbstractFeatureReader.getFeatureReader` takes a block-compressed file with a `.tbi` beside it
//! as tabix; anything else goes to the Tribble reader, which with `requireIndex` looks for a
//! `.idx` and refuses without one -- so the message names `.idx` even for a `.vcf.gz`.
//!
//! # A compressed output is indexed with tabix, which checks the order
//!
//! The plain writer's Tribble index takes records in any order. The tabix creator does not: a
//! record that starts before the one written before it on the same contig is an uncaught
//! `IllegalArgumentException` naming both, their file positions BGZF virtual offsets taken before
//! each line is written. That check is reproduced here; the `.tbi` itself is not written, because
//! the htsjdk-rs revision this repository pins predates its tabix module. The index is a side file
//! the covering array does not compare, and its refusals are the part of it a row can observe.

use std::io::Write;

use picard_analysis::vcf_io::{
    die, header_dictionary, read_path, write_output, write_records, Record,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// `IOUtil.hasBlockCompressedExtension`.
fn has_block_compressed_extension(path: &str) -> bool {
    [".gz", ".gzip", ".bgz", ".bgzf"]
        .iter()
        .any(|ext| path.to_lowercase().ends_with(ext))
}

/// `TabixIndexCreator.addFeature`'s two refusals, without the index it builds.
#[derive(Default)]
struct TabixOrder {
    sequence_names: Vec<String>,
    /// `previousFeature`: reference index, start, end, and the position it starts at.
    previous: Option<(usize, i64, i64, u64)>,
}

impl TabixOrder {
    fn add_feature(
        &mut self,
        contig: &str,
        start: i64,
        end: i64,
        position: u64,
    ) -> Result<(), String> {
        let current = self.sequence_names.last().map(String::as_str);
        let reference_index = if current == Some(contig) {
            self.sequence_names.len() - 1
        } else {
            if current.is_some() && self.sequence_names.iter().any(|n| n == contig) {
                // The Java names the feature by `VariantContext.toString()`, which this port does
                // not reproduce; no record in the corpus reaches it.
                return Err(format!(
                    "java.lang.IllegalArgumentException: Sequence {contig} added out sequence of \
                     order"
                ));
            }
            self.sequence_names.len()
        };
        let describe = |(index, start, end, position): (usize, i64, i64, u64)| {
            format!(
                "TabixFeature{{referenceIndex={index}, start={start}, end={end}, \
                 featureStartFilePosition={position}, featureEndFilePosition=-1}}"
            )
        };
        let this = (reference_index, start, end, position);
        if let Some(previous) = self.previous {
            // `TabixFeature.compareTo`: the reference index, then the start.
            if (previous.0, previous.1) > (this.0, this.1) {
                return Err(format!(
                    "java.lang.IllegalArgumentException: Features added out of order: previous \
                     ({}) > next ({})",
                    describe(previous),
                    describe(this)
                ));
            }
        }
        self.previous = Some(this);
        if reference_index == self.sequence_names.len() {
            self.sequence_names.push(contig.to_string());
        }
        Ok(())
    }
}

/// The BGZF writer, with the tabix creator's order checks beside it when indexing, as
/// `VariantContextWriterBuilder` builds it for a `BLOCK_COMPRESSED_VCF`.
fn write_block_compressed(
    path: &str,
    header: &htsjdk_vcf::header::VcfHeader,
    records: &[Record],
    indexing: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let text = write_records(header, records).map_err(|e| format!("{e:?}"))?;
    let header_length = header.write().len();
    let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
    writer.write_all(&text.as_bytes()[..header_length])?;
    let mut order = TabixOrder::default();
    for (record, line) in records
        .iter()
        .zip(text[header_length..].split_inclusive('\n'))
    {
        if indexing {
            let variant = &record.variant;
            if let Err(exception) = order.add_feature(
                &variant.contig,
                variant.start,
                variant.stop,
                writer.file_pointer(),
            ) {
                die(&exception);
            }
        }
        writer.write_all(line.as_bytes())?;
    }
    std::fs::write(path, writer.into_inner()?)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let require_index = arg(&args, "REQUIRE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(true);
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if input.ends_with(".bcf") || output.ends_with(".bcf") {
        return Err("BCF is not ported: htsjdk-rs has no BCF codec".into());
    }

    // `new VCFFileReader(INPUT, REQUIRE_INDEX)`: the header is read before the index is looked for.
    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    if require_index {
        let tabix = has_block_compressed_extension(&input)
            && std::path::Path::new(&format!("{input}.tbi")).exists();
        if !tabix && !std::path::Path::new(&format!("{input}.idx")).exists() {
            let absolute = std::path::absolute(&input)?;
            die(&format!(
                "htsjdk.tribble.TribbleException: An index is required, but none found with file \
                 ending .idx, for input source: file://{}",
                absolute.display()
            ));
        }
    }

    // `new VCFHeader(reader.getFileHeader())`: the same lines and samples.
    let header = vcf.file.header.clone();
    let dictionary = header_dictionary(&header);
    if create_index && dictionary.is_none() {
        die(
            "picard.PicardException: A sequence dictionary must be available in the input file \
             when creating indexed output.",
        );
    }
    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };

    if output.ends_with(".vcf.gz") || output.ends_with(".vcf.bgz") {
        write_block_compressed(&output, &header, &vcf.records, create_index)?;
    } else if output.ends_with(".vcf") {
        write_output(&output, &header, &vcf.records, index_dictionary)?;
    } else {
        die(
            "java.lang.IllegalArgumentException: Output format type is not set, or could not be \
             inferred from the output path. If a path was used, does it have a valid VCF \
             extension (.vcf, .vcf.gz, .vcf.bgz, .bcf)?",
        );
    }
    Ok(())
}
