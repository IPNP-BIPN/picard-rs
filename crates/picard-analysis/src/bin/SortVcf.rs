//! `SortVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.SortVcf.doWork` at tag 3.4.0 for one input. The header merge and the sort
//! live in `picard_analysis::sort_vcf`; the read and write round trip is `picard_analysis::vcf_io`.
//!
//! The dictionary decides two things. A VCF with no contig lines is refused unless
//! `SEQUENCE_DICTIONARY` gives one, and then that one becomes the file's contig lines. A VCF WITH
//! contig lines keeps them, and a `SEQUENCE_DICTIONARY` beside it is only checked against them:
//! `assertSameDictionary` throws an `AssertionError`, which the tool wraps in an
//! `IllegalArgumentException`, so the message is the error's `toString()` -- class name and all.
//!
//! The dictionary file is read through `SamReaderFactory`, as a SAM header: a `.dict`, whose
//! `@SQ` lines are the sequences.

use picard_analysis::args::arg;
use picard_analysis::sort_vcf::{smart_merge_header_lines, sort_records};
use picard_analysis::vcf_io::{
    assert_same_dictionary, die, header_dictionary, parse_sam_dictionary, sample_names_in_order,
    set_sequence_dictionary, thrown, Kept, VcfReader, VcfWriter,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SortVcf");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
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
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let dictionary = match &sequence_dictionary {
        Some(path) => Some(parse_sam_dictionary(&std::fs::read_to_string(path)?)),
        None => None,
    };

    // `collectFileReadersAndHeaders`.
    let mut vcf = VcfReader::open(&input).unwrap_or_else(|exception| die(&exception));
    match header_dictionary(&vcf.file.header) {
        None => match &dictionary {
            None => {
                let absolute = std::path::absolute(&input)?;
                vcf.fail(&format!(
                    "java.lang.IllegalArgumentException: Sequence dictionary was missing or empty \
                     for the VCF: {} Please add a sequence dictionary to this VCF or specify \
                     SEQUENCE_DICTIONARY.",
                    absolute.display()
                ));
            }
            Some(given) => set_sequence_dictionary(&mut vcf.file.header, given),
        },
        Some(own) => match &dictionary {
            // The file's own dictionary becomes the one later inputs are checked against, and
            // with one input there are none.
            None => {}
            Some(given) => {
                if let Err(message) = assert_same_dictionary(given, &own) {
                    vcf.fail(&format!(
                        "java.lang.IllegalArgumentException: java.lang.AssertionError: {message}"
                    ));
                }
            }
        },
    }

    // `new VCFHeader(VCFUtils.smartMergeHeaders(inputHeaders, false), sampleList)`.
    let mut header = vcf.file.header.clone();
    header.lines = smart_merge_header_lines(&header.lines).unwrap_or_else(|message| {
        vcf.fail(&format!("java.lang.IllegalStateException: {message}"));
    });
    header.samples = sample_names_in_order(&vcf.file.header);

    let output_dictionary = header_dictionary(&header).unwrap_or_default();
    let contigs: Vec<String> = output_dictionary.iter().map(|s| s.name.clone()).collect();
    // Each record is held as its line and sorted by its contig and start, then decoded again on
    // its way out: the input's size in memory rather than its decoded records'.
    let mut records: Vec<(usize, i64, Kept)> = Vec::new();
    while let Some(record) = vcf.try_next() {
        let record = record.unwrap_or_else(|exception| die(&exception));
        let variant = &record.variant;
        let contig = contigs.iter().position(|c| *c == variant.contig);
        let start = variant.start;
        records.push((contig.unwrap_or(usize::MAX), start, vcf.keep(record)));
    }
    if let Err(message) = sort_records(&mut records, &contigs, |record| {
        (contigs.get(record.0).map_or("", String::as_str), record.1)
    }) {
        die(&message);
    }

    let index_dictionary = if create_index {
        Some(output_dictionary.as_slice())
    } else {
        None
    };
    let mut writer = VcfWriter::create(&output, &header, index_dictionary, Default::default())
        .or_else(thrown)?;
    for (_, _, kept) in records {
        writer.add(&vcf.restore(kept)).or_else(thrown)?;
    }
    writer.finish().or_else(thrown)?;
    Ok(())
}
