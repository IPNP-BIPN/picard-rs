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
use picard_analysis::jumping_library::{
    collect_files, main_pass_reaches, outward_pass_reaches, unsorted_message, Arguments, Pair,
};

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

    let read = inputs
        .iter()
        .map(|input| pairs_of(input))
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
    let per_file = picard_analysis::jumping_library::SAMPLE_FOR_MODE / read.len().max(1);
    for file in &read {
        if let Some(before) = file.unpaired_at {
            if outward_pass_reaches(&file.pairs, before, per_file, &arguments) {
                unpaired();
            }
        }
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
        if file
            .unpaired_at
            .is_some_and(|before| main_pass_reaches(&file.pairs, before))
        {
            unpaired();
        }
    }
    let files: Vec<&[(Pair, bool)]> = read.iter().map(|file| file.pairs.as_slice()).collect();
    let metrics = collect_files(&files, &arguments);
    let mut file = MetricsFile::new();
    for header in picard_analysis::command_line::default_headers("CollectJumpingLibraryMetrics") {
        file.add_header(&header);
    }
    file.add_metric(&metrics);
    std::fs::write(&output, file.write())?;
    Ok(())
}

/// One file as the passes read it: its name and whether it is coordinate sorted, its first ends
/// in file order, each with whether its mate lies to the right, and how many of them come before
/// its first read that is not paired.
struct FileReads {
    name: String,
    sorted: bool,
    pairs: Vec<(Pair, bool)>,
    unpaired_at: Option<usize>,
}

fn pairs_of(input: &str) -> Result<FileReads, Box<dyn std::error::Error>> {
    let (header, records) = picard_analysis::bam_stream::Input::open(input)?.into_parts();
    let mut unpaired_at = None;
    let mut first_ends = 0usize;
    // "We're getting all our info from the first of each pair."
    let pairs = records
        .into_iter()
        .filter(|rec| {
            if rec.flags & READ_PAIRED == 0 {
                unpaired_at.get_or_insert(first_ends);
                return false;
            }
            let first = rec.flags & FIRST_OF_PAIR != 0;
            first_ends += usize::from(first);
            first
        })
        .map(|rec| {
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
            let mate_is_to_the_right = jumping_library::mate_is_to_the_right(&rec);
            (pair, mate_is_to_the_right)
        })
        .collect();
    Ok(FileReads {
        name: input.rsplit('/').next().unwrap_or(input).to_string(),
        sorted: header.attributes.get("SO") == Some("coordinate"),
        pairs,
        unpaired_at,
    })
}
