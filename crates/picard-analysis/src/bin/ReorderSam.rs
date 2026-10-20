//! `ReorderSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.ReorderSam.doWork` at tag 3.4.0 down to what the array varies. The reordering
//! and its two refusals live in `picard_analysis::reorder_sam`.
//!
//! Two arguments decide whether a mismatch between the input's dictionary and the new one is an
//! abort or a rewrite, and they are independent:
//!
//! * `ALLOW_INCOMPLETE_DICT_CONCORDANCE` lets a read contig that the new dictionary does not name
//!   become unmapped, instead of aborting the run;
//! * `ALLOW_CONTIG_LENGTH_DISCORDANCE` lets a name that matches with a different LENGTH warn
//!   instead of aborting.
//!
//! `SEQUENCE_DICTIONARY` is read the way `SAMSequenceDictionaryExtractor` reads it: a FASTA does
//! not carry a dictionary, so the `.dict` beside it is what is opened. Passing `ref.fasta` and
//! reading `ref.dict` is not a shortcut here, it is the reference's own rule.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::text_parse::ValidationStringency;
use picard_analysis::args::arg;
use picard_analysis::bam_stream::{Input, Records};
use picard_analysis::reorder_sam::{reorder_header, reorder_record, Options};
use picard_analysis::sort_sam::SortingCollection;
use picard_analysis::stream_output::Output;

fn flag(args: &[String], key: &str) -> bool {
    arg(args, key)
        .map(|v| picard_analysis::args::is_true(&v))
        .unwrap_or(false)
}

/// The text of a SAM or BAM input's header, for a dictionary read from one. Its records are read
/// through and checked as the whole text was, and none is kept.
fn read_as_sam(path: &str) -> Result<String, Box<dyn std::error::Error>> {
    let (header, records) = Input::open_with(path, ValidationStringency::Lenient)?.into_parts();
    for record in records {
        picard_analysis::bam_stream::sam_line(&header, &record)
            .ok_or("records failed to re-encode as SAM")?;
    }
    Ok(header.encode())
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`: a FASTA's dictionary is the `.dict` file
/// beside it, and any other file is read for its own header.
fn read_dictionary(path: &str) -> Result<String, Box<dyn std::error::Error>> {
    for suffix in [".fasta", ".fa"] {
        if let Some(stem) = path.strip_suffix(suffix) {
            let beside = format!("{stem}.dict");
            if std::path::Path::new(&beside).exists() {
                return Ok(std::fs::read_to_string(beside)?);
            }
        }
    }
    read_as_sam(path)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("ReorderSam");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let dictionary = arg(&args, "SEQUENCE_DICTIONARY=")
        .or_else(|| arg(&args, "SD="))
        .ok_or("SEQUENCE_DICTIONARY= is required")?;

    let options = Options {
        allow_incomplete_dict_concordance: flag(&args, "ALLOW_INCOMPLETE_DICT_CONCORDANCE="),
        allow_contig_length_discordance: flag(&args, "ALLOW_CONTIG_LENGTH_DISCORDANCE="),
    };

    let opened = Input::open_with(&input, ValidationStringency::Lenient)?;
    let streamed = opened.streams();
    let (header, records) = opened.into_parts();
    let dict = read_dictionary(&dictionary)?;
    let (out_header, new_index, sort) =
        reorder_header(header, &dict, &options).map_err(|e| format!("{e:?}"))?;
    let reordered = |records: Records| {
        records.map(|mut record| {
            reorder_record(&mut record, &new_index);
            record
        })
    };

    // `SAMFileWriterFactory.makeWriter` decides by the output's name: `.sam` is text, a CRAM name
    // is CRAM, anything else is a BAM, including the `output.txt` the arrays pass.
    let mut out = Output::create(&output, &out_header, false)?;
    // The presorted=false writer re-sorts by the header's sort order, stably, so records already
    // in that order are written as they come. One that is not starts the output again, sorted.
    let mut previous: Option<BamRecord> = None;
    let mut in_order = true;
    for record in reordered(records) {
        if let (Some(compare), Some(previous)) = (sort, &previous) {
            if compare(previous, &record) == std::cmp::Ordering::Greater {
                in_order = false;
                break;
            }
        }
        out.write(&record)?;
        previous = Some(record);
    }
    if !in_order {
        drop(out);
        let compare = sort.expect("only a sorting writer finds records out of order");
        let mut sorting = SortingCollection::from_args(compare, &args)?;
        let again = if streamed {
            Input::reopen(&input, true)?
        } else {
            Input::open_with(&input, ValidationStringency::Lenient)?
        };
        for record in reordered(again.into_parts().1) {
            sorting.add(record)?;
        }
        let sorted = sorting.into_sorted()?;
        out = Output::create(&output, &out_header, false)?;
        for record in sorted {
            out.write(&record?)?;
        }
    }
    // A CRAM is indexed under CREATE_INDEX, a BAM never was here.
    out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
    Ok(())
}
