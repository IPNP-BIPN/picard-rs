//! `CollectJumpingLibraryMetrics` as a runnable binary: the covering array's port side.
//!
//! Everything is read off the FIRST end of each pair, so a file's second ends never reach the
//! counters at all, and the chimera threshold is settled in a pass of its own before any pair is
//! bucketed. The order of the three chimera tests is the behaviour: an oversized insert counts as
//! oversized even when it is also tandem, and a tandem pair counts as tandem even when it also
//! spans two chromosomes.
//!
//! The refusal is the tool's, typo included: `SAM file must <name> must be sorted in coordintate
//! order`.

use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::args::arg;
use picard_analysis::jumping_library;
use picard_analysis::jumping_library::{unsorted_message, Arguments, Pair, Tally};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const DUPLICATE: u16 = 0x400;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CollectJumpingLibraryMetrics");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let inputs = picard_analysis::args::inputs(&args)?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let number = |key: &str, default: i64| -> Result<i64, Box<dyn std::error::Error>> {
        Ok(arg(&args, key)
            .map(|v| v.parse::<i64>())
            .transpose()?
            .unwrap_or(default))
    };
    let arguments = Arguments {
        minimum_mapping_quality: number("MINIMUM_MAPPING_QUALITY=", 0)? as i32,
        tail_limit: number("TAIL_LIMIT=", 10_000)? as i32,
        chimera_kb_min: number("CHIMERA_KB_MIN=", 100_000)?,
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    // Every file is read before anything is refused, as when the pairs were held: a file that does
    // not decode is answered first. Each is read once; what the two passes would have met is noted
    // on the way (`Tally`, `FileReads`).
    let mut tally = Tally::new(inputs.len());
    let read = inputs
        .iter()
        .map(|input| read_file(input, &mut tally, &arguments))
        .collect::<Result<Vec<_>, _>>()?;
    // `getFirstOfPairFlag()` throws for a read that is not paired, and both passes ask it of every
    // record they reach: `getOutieMode` over every file first, then the counting pass, file by
    // file, after each one's sort order is checked.
    let unpaired = || -> ! {
        eprintln!(
            "Exception in thread \"main\" java.lang.IllegalStateException: Inappropriate call if \
             not paired read"
        );
        std::process::exit(1);
    };
    if read.iter().any(|file| file.outward_reaches_unpaired) {
        unpaired();
    }
    // Every file is read in order into one set of counters; the first that is not coordinate
    // sorted is the refusal.
    for file in &read {
        if !file.sorted {
            eprintln!(
                "Exception in thread \"main\" picard.PicardException: {}",
                unsorted_message(&file.name)
            );
            std::process::exit(1);
        }
        if file.main_reaches_unpaired {
            unpaired();
        }
    }
    let metrics = tally.finish(&arguments);
    let mut file = MetricsFile::new();
    for header in picard_analysis::command_line::default_headers("CollectJumpingLibraryMetrics") {
        file.add_header(&header);
    }
    file.add_metric(&metrics);
    std::fs::write(&output, file.write())?;
    Ok(())
}

/// One file as the passes would have read it: its name, whether it is coordinate sorted, and
/// whether the sampling pass and the counting pass reach its first read that is not paired.
struct FileReads {
    name: String,
    sorted: bool,
    outward_reaches_unpaired: bool,
    main_reaches_unpaired: bool,
}

fn read_file(
    input: &str,
    tally: &mut Tally,
    arguments: &Arguments,
) -> Result<FileReads, Box<dyn std::error::Error>> {
    let (header, records) = picard_analysis::bam_stream::Input::open(input)?.into_parts();
    tally.start_file();
    let mut unpaired = None;
    // Whether a first end before the current record is unmapped at both ends with no reference,
    // which is where the counting pass stops.
    let mut main_stopped = false;
    for rec in records {
        if rec.flags & READ_PAIRED == 0 {
            unpaired.get_or_insert((tally.sampling_reaches(), !main_stopped));
            continue;
        }
        // "We're getting all our info from the first of each pair."
        if rec.flags & FIRST_OF_PAIR == 0 {
            continue;
        }
        let mate_quality = match rec.tags.get(Tag::new(b"MQ")) {
            Some(TagValue::Int(value)) => Some(*value as i32),
            _ => None,
        };
        let pair = Pair {
            reference_index: rec.reference_index,
            mate_reference_index: rec.mate_reference_index,
            reverse: rec.flags & READ_REVERSE != 0,
            mate_reverse: rec.flags & MATE_REVERSE != 0,
            insert_size: rec.inferred_insert_size as i64,
            duplicate: rec.flags & DUPLICATE != 0,
            mate_quality,
            mapping_quality: rec.mapping_quality as i32,
            unmapped: rec.flags & READ_UNMAPPED != 0,
            mate_unmapped: rec.flags & MATE_UNMAPPED != 0,
        };
        main_stopped |= pair.unmapped && pair.mate_unmapped && pair.reference_index < 0;
        tally.add(
            &pair,
            jumping_library::mate_is_to_the_right(&rec),
            arguments,
        );
    }
    let (outward, main) = unpaired.unwrap_or((false, false));
    Ok(FileReads {
        name: input.rsplit('/').next().unwrap_or(input).to_string(),
        sorted: header.attributes.get("SO") == Some("coordinate"),
        outward_reaches_unpaired: outward,
        main_reaches_unpaired: main,
    })
}
