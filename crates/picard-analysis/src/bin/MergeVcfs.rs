//! `MergeVcfs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.MergeVcfs.doWork` at tag 3.4.0. The read and write round trip is
//! `picard_analysis::vcf_io`, the header merge is htsjdk-rs's `smartMergeHeaders`, and the record
//! merge is `picard_analysis::merge_vcfs`, whose heap decides the order of records that tie.
//!
//! Each input is checked as it is opened, in order: a file without contig lines takes
//! `SEQUENCE_DICTIONARY`'s (and is refused without one), the first file's contig lines become the
//! comparator every later file must be compatible with, and every file's SORTED sample names must
//! equal the first's -- so the same samples in another column order are accepted, and the output's
//! columns are in sorted order.
//!
//! `COMMENT` goes into the first header through `addMetaDataLine`, which keeps only the first line
//! under one key, so of several comments only the first given is written.

use htsjdk_vcf::comparator::VariantContextComparator;
use htsjdk_vcf::encoder::MissingFields;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use picard_analysis::args::arg;
use picard_analysis::merge_vcfs::{merge_streams, MergeFailure};
use picard_analysis::vcf_io::{
    add_other_meta_data_line, die, fail_all, header_dictionary, is_compressed_vcf,
    parse_sam_dictionary, sample_names_in_order, set_sequence_dictionary, unroll_paths, Record,
    Sequence, VcfReader, VcfWriter,
};

/// A collection argument: every occurrence appends, and the value `null` empties it (Barclay).
fn collection(args: &[String], keys: &[&str]) -> Vec<String> {
    let mut values = Vec::new();
    for a in args {
        if let Some(value) = keys.iter().find_map(|key| a.strip_prefix(key)) {
            if value == "null" {
                values.clear();
            } else {
                values.push(value.to_string());
            }
        }
    }
    values
}

const SEQ_DICT_REQUIRED: &str = "A sequence dictionary must be available (either through the \
                                 input file or by setting it explicitly).";

const NPE_UNKNOWN_CONTIG: &str = "java.lang.NullPointerException: Cannot invoke \
     \"java.lang.Integer.intValue()\" because the return value of \"java.util.Map.get(Object)\" \
     is null";

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

fn contig_lines(header: &VcfHeader) -> Vec<HeaderLine> {
    header
        .lines
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .cloned()
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("MergeVcfs");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let inputs = collection(&args, &["INPUT=", "I="]);
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let sequence_dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "D="))
        .filter(|value| value != "null");
    let comments = collection(&args, &["COMMENT=", "CO="]);
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if !output.ends_with(".vcf")
        && !is_compressed_vcf(&output)
        && !picard_analysis::vcf_io::is_bcf(&output)
    {
        return Err("only a .vcf, .vcf.gz or .bcf OUTPUT is ported".into());
    }

    // Read through `SamReaderFactory`, as a SAM header: a `.dict`'s `@SQ` lines.
    let mut dictionary: Option<Vec<Sequence>> = match &sequence_dictionary {
        Some(path) => Some(parse_sam_dictionary(&std::fs::read_to_string(path)?)),
        None => None,
    };

    let paths = unroll_paths(&inputs).unwrap_or_else(|exception| die(&exception));
    let mut files: Vec<VcfReader> = Vec::with_capacity(paths.len());
    let mut comparator: Option<VariantContextComparator> = None;
    let mut sample_list: Vec<String> = Vec::new();
    for path in &paths {
        let vcf =
            VcfReader::open(path).unwrap_or_else(|exception| fail_all(&mut files, &exception));
        files.push(vcf);
        let first = files.len() == 1;
        let header = &mut files.last_mut().expect("just pushed").file.header;
        let admitted = admit(header, path, &dictionary, &mut comparator, &mut sample_list);
        if let Err(exception) = admitted {
            fail_all(&mut files, &exception);
        }
        if dictionary.is_none() {
            dictionary = header_dictionary(header);
        }
        if first {
            for comment in &comments {
                add_other_meta_data_line(header, "MergeVcfs.comment", comment);
            }
        }
    }

    if create_index && dictionary.is_none() {
        fail_all(
            &mut files,
            &format!("picard.PicardException: Index creation failed. {SEQ_DICT_REQUIRED}"),
        );
    }

    // `new VCFHeader(VCFUtils.smartMergeHeaders(headers, false), sampleList)`.
    let versions: Vec<Option<&str>> = files
        .iter()
        .map(|f| f.file.header_version.map(|v| v.version_string()))
        .collect();
    let sources: Vec<Source> = files
        .iter()
        .zip(&versions)
        .map(|(f, version)| Source {
            header: &f.file.header,
            version: *version,
        })
        .collect();
    let merged = smart_merge_headers(&sources, false)
        .map(|(merged, _)| merged)
        .map_err(|e| format!("{}: {}", e.class(), e.message()));
    let merged = merged.unwrap_or_else(|exception| fail_all(&mut files, &exception));
    let header = VcfHeader {
        lines: merged,
        samples: sample_list,
    };

    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    let mut writer =
        VcfWriter::create(&output, &header, index_dictionary, MissingFields::default())
            .unwrap_or_else(|error| fail_all(&mut files, &error.exception()));
    let comparator = comparator.expect("at least one input");
    let inputs = files.len();
    let files = std::cell::RefCell::new(files);
    let failure = merge_streams(
        inputs,
        |input| {
            let mut files = files.borrow_mut();
            // An input is read whole before the next one is opened, in the run this streams.
            match files[input].try_next()? {
                Ok(record) => Some(record),
                Err(exception) => fail_all(&mut files[..input], &exception),
            }
        },
        "htsjdk.variant.variantcontext.VariantContextComparator",
        |a: &Record, b: &Record| {
            comparator
                .compare(&a.variant, &b.variant)
                .map(|c| c.cmp(&0))
        },
        |record| {
            if let Err(error) = writer.add(record) {
                fail_all(&mut files.borrow_mut(), &error.exception());
            }
        },
    );
    let mut files = files.into_inner();
    if let Some(failure) = failure {
        match failure {
            MergeFailure::Compare(_) => fail_all(&mut files, NPE_UNKNOWN_CONTIG),
            MergeFailure::NotSorted(message) => fail_all(&mut files, &message),
        }
    }
    writer
        .finish()
        .unwrap_or_else(|error| fail_all(&mut files, &error.exception()));
    Ok(())
}

/// One input's header checked against the ones before it, as `doWork` checks each file it opens,
/// and its dictionary filled in from `SEQUENCE_DICTIONARY` when it has no contig lines.
fn admit(
    header: &mut VcfHeader,
    path: &str,
    dictionary: &Option<Vec<Sequence>>,
    comparator: &mut Option<VariantContextComparator>,
    sample_list: &mut Vec<String>,
) -> Result<(), String> {
    if contig_lines(header).is_empty() {
        match dictionary {
            None => {
                return Err(format!(
                    "java.lang.IllegalArgumentException: {SEQ_DICT_REQUIRED}"
                ))
            }
            Some(given) => set_sequence_dictionary(header, given),
        }
    }

    let lines = contig_lines(header);
    match comparator {
        None => {
            *comparator = Some(
                VariantContextComparator::from_header_lines(&lines)
                    .map_err(|e| format!("{}: {}", e.class(), e.message()))?,
            );
        }
        Some(existing) => {
            if !existing.is_compatible(&lines) {
                return Err(format!(
                    "java.lang.IllegalArgumentException: The contig entries in input path {} \
                     are not compatible with the others.",
                    absolute(path)
                ));
            }
        }
    }

    let samples = sample_names_in_order(header);
    if sample_list.is_empty() {
        *sample_list = samples;
    } else if *sample_list != samples {
        return Err(format!(
            "java.lang.IllegalArgumentException: Input path {} has sample entries that don't \
             match the other files.",
            absolute(path)
        ));
    }
    Ok(())
}
