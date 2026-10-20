//! `RevertOriginalBaseQualitiesAndAddMateCigar` as a runnable binary: the covering array's port
//! side.
//!
//! Ports `picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar.doWork` at tag 3.4.0. The revert,
//! the mate-info pass and the sort orders live in
//! `picard_analysis::revert_original_quals_add_mate_cigar`.
//!
//! The tool's first act is to decide whether to do nothing. It reads the input twice: once to ask
//! whether any record carries an `OQ` and whether the first mapped pair already has its `MC`, and
//! again to do the work. When the answer is that there is nothing to do it returns zero having
//! written NO output file, which is a different thing from writing an unchanged one.
//!
//! `SORT_ORDER` decides two things: the `SO` the output header carries, and the order the records
//! come out in -- the writer is built unsorted, so it sorts. Two of the five orders name no
//! comparator, and the records then keep the query-name order the mate-info pass left them in.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::{coordinate, query_name};
use picard_analysis::args::arg;
use picard_analysis::revert_original_quals_add_mate_cigar::{
    add_mate_info, restore_original_qualities, skip_answer, sort_by_duplicate_order, CanSkip,
    Options, SortOrder,
};
use picard_analysis::sort_sam::SortingCollection;
use picard_analysis::stream_output::Output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("RevertOriginalBaseQualitiesAndAddMateCigar");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let sort_order = match arg(&args, "SORT_ORDER=") {
        None => None,
        Some(value) => {
            Some(SortOrder::parse(&value).ok_or_else(|| format!("unknown SORT_ORDER: {value}"))?)
        }
    };
    let options = Options {
        restore_original_qualities: arg(&args, "RESTORE_ORIGINAL_QUALITIES=")
            .map(|value| picard_analysis::args::is_true(&value))
            .unwrap_or(true),
        sort_order,
        max_records_to_examine: arg(&args, "MAX_RECORDS_TO_EXAMINE=")
            .and_then(|value| value.parse().ok())
            .unwrap_or(10_000),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let (input_header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();
    let mut records = records.peekable();

    // `canSkipSAMFile`, on its own pass over the first records of the input, which are then read
    // again. When it skips, nothing is written: the output file does not exist, which is the
    // answer the reference gives too.
    let mut head = Vec::new();
    let mut answer = CanSkip::FoundNoEvidence;
    while head.len() < options.max_records_to_examine.max(0) as usize {
        let Some(record) = records.next() else { break };
        let answered = skip_answer(&record, &options);
        head.push(record);
        if let Some(answered) = answered {
            answer = answered;
            break;
        }
    }
    if answer.skips() {
        return Ok(());
    }

    let output_order = options
        .sort_order
        .or_else(|| input_header.attributes.get("SO").and_then(SortOrder::parse));
    // `SO:unknown` names no comparator, and the writer -- built unsorted -- hands the records to a
    // `SortingCollection` that falls back to natural ordering. A SAM record is not `Comparable`,
    // so the run dies on a cast rather than on anything about the data. `unsorted` does NOT: the
    // writer has a case for it and never sorts.
    //
    // The class in the message is the record's, so it names `BAMRecord` for a BAM input and would
    // name `SAMRecord` for a text one.
    if output_order == Some(SortOrder::Unknown) {
        eprintln!(
            "Exception in thread \"main\" java.lang.ClassCastException: class \
             htsjdk.samtools.BAMRecord cannot be cast to class java.lang.Comparable \
             (htsjdk.samtools.BAMRecord is in unnamed module of loader 'app'; \
             java.lang.Comparable is in module java.base of loader 'bootstrap')"
        );
        std::process::exit(1);
    }
    let output_order = output_order.unwrap_or(SortOrder::Unsorted);
    let mut header = input_header.clone();
    header.set_sort_order(output_order.name());

    // The exception is thrown while the sorted records are walked, after the writer was built and
    // before it was handed a record (it sorts, so it holds them all until it is closed). A SAM
    // writer has flushed its header by then; a BAM writer has written only the blocks its header
    // filled. Neither is closed, so neither the index nor the MD5 file appears.
    let thrown = |thrown: String| -> Result<(), Box<dyn std::error::Error>> {
        let partial = if output.ends_with(".sam") {
            header.encode().into_bytes()
        } else {
            let writer = picard_analysis::output::bam_writer(Vec::new(), &header)
                .expect("in-memory BAM writer never fails");
            picard_analysis::java_writer::bam_on_disk(
                &writer.finish().expect("in-memory BAM writer never fails"),
            )
        };
        std::fs::write(&output, partial)?;
        eprintln!("Exception in thread \"main\" {thrown}");
        std::process::exit(1);
    };

    // Original qualities restored, then a query-name sort to group the templates.
    let mut by_name = SortingCollection::from_args(query_name::compare, &args)?;
    for mut record in head.into_iter().chain(records) {
        if options.restore_original_qualities {
            restore_original_qualities(&mut record);
        }
        by_name.add(record)?;
    }
    let mut by_name = by_name.into_sorted()?;
    // Each run of records sharing a name is a template, given its mate information as a whole.
    let mut next = by_name.next().transpose()?;
    let failed = std::cell::Cell::new(None);
    let mut templates = std::iter::from_fn(|| {
        let first = next.take()?;
        let mut template = vec![first];
        loop {
            match by_name.next() {
                Some(Ok(record)) if record.read_name == template[0].read_name => {
                    template.push(record)
                }
                Some(Ok(record)) => {
                    next = Some(record);
                    break;
                }
                Some(Err(e)) => return Some(vec![Err(e)]),
                None => break,
            }
        }
        if let Err(e) = add_mate_info(&mut template) {
            failed.set(Some(e));
            return None;
        }
        Some(template.into_iter().map(Ok).collect())
    })
    .flatten();

    // Then the output order: the writer sorts into the order the header names, and the two that
    // name no comparator leave the records in query-name order.
    let sorted: Box<dyn Iterator<Item = std::io::Result<BamRecord>>> = match output_order {
        SortOrder::Coordinate => {
            let mut sorting = SortingCollection::from_args(coordinate::compare, &args)?;
            for record in templates.by_ref() {
                sorting.add(record?)?;
            }
            Box::new(sorting.into_sorted()?)
        }
        // shortcut: SO:duplicate is sorted in memory, as its comparator reads the header, which
        // the SortingCollection's comparator cannot carry; stream it if a duplicate-sorted output
        // of a large file turns up.
        SortOrder::Duplicate => {
            let mut all = templates.by_ref().collect::<Result<Vec<_>, _>>()?;
            sort_by_duplicate_order(&header, &mut all);
            Box::new(all.into_iter().map(Ok))
        }
        SortOrder::Queryname | SortOrder::Unsorted | SortOrder::Unknown => Box::new(templates),
    };
    if let Some(e) = failed.take() {
        return thrown(e);
    }
    let mut out = Output::create(&output, &header, false)?;
    for record in sorted {
        out.write(&record?)?;
    }
    if let Some(e) = failed.take() {
        drop(out);
        return thrown(e);
    }
    // A CRAM is indexed under CREATE_INDEX, a BAM never was here.
    out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
    Ok(())
}
