//! `FastqToSam` as a runnable binary: the covering array's port side.
//!
//! The whole tool is `picard_analysis::fastq_to_sam::fastq_to_bam_full`: one FASTQ or a pair, the
//! read group built from every field the command line gives, the comments, the quality range
//! `MIN_Q..MAX_Q` checked read by read, and the sort order. Before any of that, two refusals:
//!
//! * `QUALITY_FORMAT`, when given, is a CLAIM about the file, and `QualityEncodingDetector`
//!   checks it against the first 10000 records; left out, the detector picks the format, as it
//!   does for a regular file only (a pipe must name one);
//! * `USE_SEQUENTIAL_FASTQS` reads the FASTQ's NAME rather than its content: the file must end in
//!   `_001` and one of five extensions so that the next file in the series can be found. The
//!   corpus's `reads_1.fastq` does not, and the refusal quotes the whole list it expected.

use htsjdk_bam::query_name;
use htsjdk_bam::record::BamRecord;
use picard_analysis::add_or_replace_read_groups::iso8601;
use picard_analysis::args::arg;
use picard_analysis::fastq_to_sam::{
    determine_quality_format, fastq_records, full_header, read_fastq, FullOptions, ReadAhead,
};
use picard_analysis::sort_sam::SortingCollection;
use picard_analysis::stream_output::Output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("FastqToSam");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let fastq = arg(&args, "FASTQ=")
        .or_else(|| arg(&args, "F1="))
        .ok_or("FASTQ= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let sample = arg(&args, "SAMPLE_NAME=").ok_or("SAMPLE_NAME= is required")?;
    let flag = |key: &str, default: bool| {
        arg(&args, key)
            .map(|v| picard_analysis::args::is_true(&v))
            .unwrap_or(default)
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let int = |key: &str, default: i32| -> Result<i32, Box<dyn std::error::Error>> {
        Ok(arg(&args, key)
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(default))
    };
    let min_q = int("MIN_Q=", 0)?;
    let max_q = int("MAX_Q=", 93)?;
    // `customCommandLineValidation`, one message and the first that applies.
    if min_q < 0 {
        eprintln!("MIN_Q must be >= 0");
        std::process::exit(1);
    }
    if max_q > 93 {
        eprintln!("MAX_Q must be <= 93");
        std::process::exit(1);
    }
    let second = arg(&args, "FASTQ2=").or_else(|| arg(&args, "F2="));
    let format = arg(&args, "QUALITY_FORMAT=");
    if let Some(format) = &format {
        if !matches!(format.as_str(), "Standard" | "Illumina" | "Solexa") {
            return Err(format!("unknown QUALITY_FORMAT: {format}").into());
        }
    }
    // `PicardHtsPath.isOther`: neither a regular file nor a directory, such as a pipe.
    let other =
        |path: &str| std::fs::metadata(path).is_ok_and(|m| !m.file_type().is_file() && !m.is_dir());
    let regular = !other(&fastq) && !second.as_deref().is_some_and(other);
    if format.is_none() && !regular {
        eprintln!(
            "QUALITY_FORMAT must be specified when either of FASTQ or FASTQ2 is not a regular file"
        );
        std::process::exit(1);
    }
    let text = read_fastq(&fastq)?;
    let second_text = match &second {
        Some(path) => Some(read_fastq(path)?),
        None => None,
    };
    let skip_blank_lines = flag("ALLOW_AND_IGNORE_EMPTY_LINES=", false);
    let refuse = |thrown: String| -> ! {
        eprintln!("Exception in thread \"main\" {thrown}");
        std::process::exit(1);
    };
    // `doWork`: the first reader is opened, and a non-empty regular input's quality format is
    // determined (or checked) before a record is converted.
    let reader1 = ReadAhead::open(&text, skip_blank_lines).unwrap_or_else(|t| refuse(t));
    let format = if !reader1.has_next() {
        if !flag("ALLOW_EMPTY_FASTQ=", false) {
            refuse(
                "picard.PicardException: Input fastq is empty. Set ALLOW_EMPTY_FASTQ if you \
                 still want to write out a file with no reads."
                    .to_string(),
            );
        }
        format
    } else if regular {
        let mut readers = vec![reader1];
        if let Some(text2) = &second_text {
            readers.push(ReadAhead::open(text2, skip_blank_lines).unwrap_or_else(|t| refuse(t)));
        }
        Some(
            determine_quality_format(readers, format.as_deref())
                .unwrap_or_else(|t| refuse(t))
                .to_string(),
        )
    } else {
        format
    };
    // An empty input never reads a quality, so the format it is converted with does not matter.
    let format = format.unwrap_or_else(|| "Standard".to_string());

    // `getSequentialFileList`, which reads the NAME rather than the content. It is checked
    // AFTER the quality format, which is why a row that gets both wrong is refused for the
    // quality and not for the name.
    if flag("USE_SEQUENTIAL_FASTQS=", false) {
        let name = fastq
            .rsplit('/')
            .next()
            .unwrap_or(&fastq)
            .to_ascii_uppercase();
        let sequential = [
            "_001.FASTQ",
            "_001.FASTQ.GZ",
            "_001.FQ",
            "_001.FQ.GZ",
            "_001.BFQ",
        ]
        .iter()
        .any(|suffix| name.ends_with(suffix));
        if !sequential {
            eprintln!(
                "Exception in thread \"main\" picard.PicardException: Could not parse the FASTQ \
                 extension (expected '_001' + '[FASTQ, FASTQ_GZ, FQ, FQ_GZ, BFQ]'): {fastq}"
            );
            std::process::exit(1);
        }
    }

    let sort_order = arg(&args, "SORT_ORDER=")
        .or_else(|| arg(&args, "SO="))
        .unwrap_or_else(|| "queryname".to_string());
    if !matches!(sort_order.as_str(), "queryname" | "coordinate" | "unsorted") {
        return Err(format!("unknown SORT_ORDER: {sort_order}").into());
    }
    let run_date = match arg(&args, "RUN_DATE=").or_else(|| arg(&args, "DT=")) {
        Some(value) => Some(iso8601(&value).ok_or(format!("unreadable RUN_DATE: {value}"))?),
        None => None,
    };
    let options = FullOptions {
        read_group_name: arg(&args, "READ_GROUP_NAME=")
            .or_else(|| arg(&args, "RG="))
            .unwrap_or_else(|| "A".to_string()),
        sample_name: sample,
        sort_order,
        library_name: arg(&args, "LIBRARY_NAME=").or_else(|| arg(&args, "LB=")),
        platform: arg(&args, "PLATFORM=").or_else(|| arg(&args, "PL=")),
        platform_unit: arg(&args, "PLATFORM_UNIT=").or_else(|| arg(&args, "PU=")),
        sequencing_center: arg(&args, "SEQUENCING_CENTER=").or_else(|| arg(&args, "CN=")),
        predicted_insert_size: arg(&args, "PREDICTED_INSERT_SIZE=")
            .or_else(|| arg(&args, "PI="))
            .map(|v| v.parse())
            .transpose()?,
        description: arg(&args, "DESCRIPTION=").or_else(|| arg(&args, "DS=")),
        run_date,
        platform_model: arg(&args, "PLATFORM_MODEL=").or_else(|| arg(&args, "PM=")),
        program_group: arg(&args, "PROGRAM_GROUP=").or_else(|| arg(&args, "PG=")),
        comments: args
            .iter()
            .filter_map(|a| a.strip_prefix("COMMENT=").or_else(|| a.strip_prefix("CO=")))
            .map(str::to_string)
            .collect(),
        min_q,
        max_q,
        quality_format: format,
        allow_and_ignore_empty_lines: skip_blank_lines,
    };
    let records = fastq_records(&text, second_text.as_deref(), &options)
        .unwrap_or_else(|thrown| refuse(thrown));
    // A queryname output is sorted through a SortingCollection, and any other is written in the
    // FASTQ's order; either way a refused record leaves no output.
    let sorted: Box<dyn Iterator<Item = Result<BamRecord, String>>> =
        if options.sort_order == "queryname" {
            let mut sorting = SortingCollection::from_args(query_name::compare, &args)?;
            for record in records {
                sorting.add(record.unwrap_or_else(|thrown| refuse(thrown)))?;
            }
            Box::new(sorting.into_sorted()?.map(|r| r.map_err(|e| e.to_string())))
        } else {
            Box::new(records)
        };
    let mut out = Output::bam(&output, &full_header(&options), false)?;
    for record in sorted {
        match record {
            Ok(record) => out.write(&record)?,
            Err(thrown) => {
                drop(out);
                refuse(thrown);
            }
        }
    }
    out.finish_create_index()?;
    Ok(())
}
