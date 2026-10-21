//! `FixMateInformation` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.FixMateInformation.doWork` at tag 3.4.0. The mate fixing, the supplemental
//! records' mate info and the two sorts live in `picard_analysis::fix_mate_information`.
//!
//! `SORT_ORDER`'s default is the trap: it is `null`, not a sort order, and a null means "the input
//! header's order" rather than "coordinate". The tool sorts into queryname order to group
//! templates whatever the caller asked for, fixes each template, and re-sorts into the output
//! order.
//!
//! `ADD_MATE_CIGAR` decides whether `MC` is written on the primaries and copied onto the
//! supplementary records, or cleared on both. `IGNORE_MISSING_MATES` decides whether a template
//! with one end is passed through or refused.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::{coordinate, query_name};
use picard_analysis::args::arg;
use picard_analysis::bam_stream::Input;
use picard_analysis::fix_mate_information::{
    fix_templates, needs_grouping, output_order, FixMateError, Options, SortOrder,
};
use picard_analysis::sort_sam::SortingCollection;
use picard_analysis::stream_output::Output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("FixMateInformation");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    let sort_order = match arg(&args, "SORT_ORDER=").as_deref() {
        None => None,
        Some("coordinate") => Some(SortOrder::Coordinate),
        Some("queryname") => Some(SortOrder::Queryname),
        Some(other) => return Err(format!("unknown SORT_ORDER: {other}").into()),
    };
    let flag = |key: &str, default: bool| {
        arg(&args, key)
            .map(|v| picard_analysis::args::is_true(&v))
            .unwrap_or(default)
    };
    let options = Options {
        sort_order,
        assume_sorted: flag("ASSUME_SORTED=", false),
        add_mate_cigar: flag("ADD_MATE_CIGAR=", true),
        ignore_missing_mates: flag("IGNORE_MISSING_MATES=", true),
        create_index: flag("CREATE_INDEX=", false),
    };

    // The reference throws rather than exiting, so its handler prints the class before the
    // message: a row that refuses is a row this has to refuse the same way, text included.
    // Whether this run created the output, which a refusal then removes: the queryname order
    // streams into it while the templates are fixed.
    let created = std::cell::Cell::new(false);
    let refuse = |error: FixMateError| -> ! {
        if created.get() {
            let _ = std::fs::remove_file(&output);
        }
        eprintln!(
            "Exception in thread \"main\" {}: {}",
            error.java_class(),
            error.message()
        );
        std::process::exit(1);
    };

    let (mut header, records) =
        Input::open_with(&input, ValidationStringency::Lenient)?.into_parts();
    let output_order = output_order(&header, &options).unwrap_or_else(|e| refuse(e));
    let sort = |compare, records: &mut dyn Iterator<Item = std::io::Result<BamRecord>>| {
        let mut sorting = SortingCollection::from_args(compare, &args)?;
        for record in records {
            sorting.add(record?)?;
        }
        Ok::<_, Box<dyn std::error::Error>>(sorting.into_sorted()?)
    };
    let mut grouped: Box<dyn Iterator<Item = std::io::Result<BamRecord>>> =
        if needs_grouping(&header, &options) {
            Box::new(sort(query_name::compare, &mut records.map(Ok))?)
        } else {
            Box::new(records.map(Ok))
        };
    // Each run of records sharing a name is a template, fixed as a whole.
    let mut next = grouped.next().transpose()?;
    let mut fixed = std::iter::from_fn(|| {
        let first = next.take()?;
        let mut template = vec![first];
        loop {
            match grouped.next() {
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
        if let Err(error) = fix_templates(
            &mut template,
            options.add_mate_cigar,
            options.ignore_missing_mates,
        ) {
            refuse(error);
        }
        Some(template.into_iter().map(Ok).collect())
    })
    .flatten();
    header.set_sort_order(output_order.name());
    match output_order {
        // `createSamFileWriter` opens a queryname output presorted, so the fixed records are written
        // in the order they come, and the writer refuses one out of name order, as htsjdk's does.
        SortOrder::Queryname => {
            let mut out = Output::create(&output, &header, false)?;
            created.set(true);
            for record in fixed {
                out.write(&record?)?;
            }
            out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
            return Ok(());
        }
        // Any other order is sorted by the writer, after every record is fixed, so a refusal
        // leaves no output.
        SortOrder::Coordinate => {
            let sorted = sort(coordinate::compare, &mut fixed)?;
            let mut out = Output::create(&output, &header, false)?;
            for record in sorted {
                out.write(&record?)?;
            }
            // A CRAM is indexed under CREATE_INDEX, a BAM never was here.
            out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
        }
    }
    Ok(())
}
