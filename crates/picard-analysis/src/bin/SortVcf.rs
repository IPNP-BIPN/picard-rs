//! `SortVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.SortVcf.doWork` at tag 3.4.0. The header merge of one input and the sort live
//! in `picard_analysis::sort_vcf`; several inputs merge their headers as `MergeVcfs` does, through
//! htsjdk-rs's `smartMergeHeaders`. The read and write round trip is `picard_analysis::vcf_io`.
//!
//! The dictionary decides two things. A VCF with no contig lines is refused unless
//! `SEQUENCE_DICTIONARY` gives one, and then that one becomes the file's contig lines. A VCF WITH
//! contig lines keeps them, and a `SEQUENCE_DICTIONARY` beside it is only checked against them:
//! `assertSameDictionary` throws an `AssertionError`, which the tool wraps in an
//! `IllegalArgumentException`, so the message is the error's `toString()` -- class name and all.
//!
//! The dictionary file is read through `SamReaderFactory`, as a SAM header: a `.dict`, whose
//! `@SQ` lines are the sequences.

use htsjdk_vcf::header::VcfHeader;
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use picard_analysis::args::arg;
use picard_analysis::liftover_vcf::LineSpool;
use picard_analysis::sort_sam::DEFAULT_MAX_RECORDS_IN_RAM;
use picard_analysis::sort_vcf::{smart_merge_header_lines, UNKNOWN_CONTIG};
use picard_analysis::vcf_io::{
    assert_same_dictionary, die, fail_all, header_dictionary, is_bcf, parse_sam_dictionary,
    sample_names_in_order, set_sequence_dictionary, thrown, Kept, VcfReader, VcfWriter,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SortVcf");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let inputs = picard_analysis::args::inputs(&args)?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    // Barclay reads the value `null` as no value, which is this argument's default.
    let sequence_dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "SD="))
        .filter(|value| value != "null");
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(true);
    let max_records_in_ram = arg(&args, "MAX_RECORDS_IN_RAM=")
        .map(|value| value.parse::<i64>())
        .transpose()?
        .unwrap_or(DEFAULT_MAX_RECORDS_IN_RAM as i64)
        .max(1) as usize;
    // shortcut: the first TMP_DIR only, where htsjdk picks the one with the most free space.
    let tmp_dir = arg(&args, "TMP_DIR=")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut dictionary = match &sequence_dictionary {
        Some(path) => Some(parse_sam_dictionary(&std::fs::read_to_string(path)?)),
        None => None,
    };

    // `collectFileReadersAndHeaders`: each input's dictionary is the given one, or the first
    // input's, and each input's samples are the first input's.
    let mut files: Vec<VcfReader> = Vec::with_capacity(inputs.len());
    let mut sample_list: Vec<String> = Vec::new();
    for input in &inputs {
        let vcf =
            VcfReader::open(input).unwrap_or_else(|exception| fail_all(&mut files, &exception));
        files.push(vcf);
        let absolute = std::path::absolute(input)?;
        let header = &mut files.last_mut().expect("just pushed").file.header;
        let mut failure = match header_dictionary(header) {
            None => match &dictionary {
                None => Some(format!(
                    "java.lang.IllegalArgumentException: Sequence dictionary was missing or empty \
                     for the VCF: {} Please add a sequence dictionary to this VCF or specify \
                     SEQUENCE_DICTIONARY.",
                    absolute.display()
                )),
                Some(given) => {
                    set_sequence_dictionary(header, given);
                    None
                }
            },
            Some(own) => match &dictionary {
                None => {
                    dictionary = Some(own);
                    None
                }
                Some(given) => assert_same_dictionary(given, &own).err().map(|message| {
                    format!(
                        "java.lang.IllegalArgumentException: java.lang.AssertionError: {message}"
                    )
                }),
            },
        };
        let samples = sample_names_in_order(header);
        if failure.is_none() {
            if sample_list.is_empty() {
                sample_list = samples;
            } else if sample_list != samples {
                failure = Some(format!(
                    "java.lang.IllegalArgumentException: Input file {} has sample names that \
                     don't match the other files.",
                    absolute.display()
                ));
            }
        }
        if let Some(exception) = failure {
            fail_all(&mut files, &exception);
        }
    }

    // `new VCFHeader(VCFUtils.smartMergeHeaders(inputHeaders, false), sampleList)`.
    let mut header = if let [vcf] = &mut files[..] {
        let mut header = vcf.file.header.clone();
        header.lines = smart_merge_header_lines(&header.lines).unwrap_or_else(|message| {
            vcf.fail(&format!("java.lang.IllegalStateException: {message}"));
        });
        header
    } else {
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
        drop(sources);
        let lines = merged.unwrap_or_else(|exception| fail_all(&mut files, &exception));
        VcfHeader {
            lines,
            samples: Vec::new(),
        }
    };
    header.samples = sample_list;

    let output_dictionary = header_dictionary(&header).unwrap_or_default();
    let contigs: Vec<String> = output_dictionary.iter().map(|s| s.name.clone()).collect();
    // `sortInputs`: every input's records, in input order, into one `SortingCollection` that
    // spills past MAX_RECORDS_IN_RAM. A text record is held as its line, decoded again on its way
    // out; a BCF record, which the reader holds whole already, is held whole and spooled by index.
    let mut spool = LineSpool::new(true, max_records_in_ram, tmp_dir);
    let mut whole: Vec<Option<Kept>> = Vec::new();
    let (mut count, mut unknown_contig) = (0usize, false);
    for (input, vcf) in files.iter_mut().enumerate() {
        while let Some(record) = vcf.try_next() {
            let record = record.unwrap_or_else(|exception| die(&exception));
            let variant = &record.variant;
            let contig = contigs.iter().position(|c| *c == variant.contig);
            unknown_contig |= contig.is_none();
            count += 1;
            let start = variant.start;
            let mut held = (input as u32).to_le_bytes().to_vec();
            match vcf.keep(record) {
                Kept::Line(line) => held.extend_from_slice(&line),
                record => {
                    held.extend_from_slice(&whole.len().to_le_bytes());
                    whole.push(Some(record));
                }
            }
            spool.add(contig.unwrap_or(0) as u32, start, &held)?;
        }
    }
    // A record on a contig the dictionary lacks fails any comparison, as in `sort_records`.
    if count > 1 && unknown_contig {
        die(UNKNOWN_CONTIG);
    }

    let index_dictionary = if create_index {
        Some(output_dictionary.as_slice())
    } else {
        None
    };
    let mut writer = VcfWriter::create(&output, &header, index_dictionary, Default::default())
        .or_else(thrown)?;
    for held in spool.into_lines()? {
        let held = held?;
        let input = u32::from_le_bytes(held[..4].try_into().expect("four bytes")) as usize;
        let kept = if !is_bcf(&inputs[input]) {
            Kept::Line(held[4..].into())
        } else {
            let at = usize::from_le_bytes(held[4..].try_into().expect("an index"));
            whole[at].take().expect("each held once")
        };
        writer.add(&files[input].restore(kept)).or_else(thrown)?;
    }
    writer.finish().or_else(thrown)?;
    Ok(())
}
