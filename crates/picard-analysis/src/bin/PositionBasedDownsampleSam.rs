//! `PositionBasedDownsampleSam` as a runnable binary: the covering array's port side.
//!
//! The library held the selection and said so: "reading the file and parsing the read names are
//! not ported". This is the rest, and the read names are where the surprise is.
//!
//! The fixtures' names are `read0314`, which the default `READ_NAME_REGEX` does not parse, so
//! every read gets `PhysicalLocationInt`'s defaults: tile `-1`, x `-1`, y `-1`. That is not a
//! degenerate case that skips the selection. `Coord` starts at zero on all four sides, so a tile
//! of unparsed reads spans `-1..0` rather than nothing at all, the widening is `1 / count` and
//! therefore zero, and every read's normalized position is exactly `0`. The mask then keeps all
//! of them or none of them depending on which side of a half the fraction falls, which is what
//! the array's two fractions measure.
//!
//! Ten rows, ten accepted, ten distinct outputs.

use picard_analysis::args::arg;
use picard_analysis::position_based_downsample_sam::{
    fraction_out_of_range_message, keep, PhysicalLocation, PG_PROGRAM_NAME,
};

use picard_analysis::stream_output::Output;

const DUPLICATE_READ: u16 = 0x400;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("PositionBasedDownsampleSam");
    picard_analysis::output::install();
    // `READ_NAME_REGEX`, which a name has to match whole for its tile and position to be read.
    let read_name_regex = picard_analysis::mark_duplicates::ReadNameRegex::from_arg(
        std::env::args()
            .find_map(|a| a.strip_prefix("READ_NAME_REGEX=").map(str::to_string))
            .as_deref(),
    )
    .unwrap_or_else(|thrown| {
        eprintln!("Exception in thread \"main\" {thrown}");
        std::process::exit(1);
    });
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let fraction: f64 = arg(&args, "FRACTION=")
        .or_else(|| arg(&args, "F="))
        .ok_or("FRACTION= is required")?
        .parse()?;
    let flag = |key: &str, default: bool| {
        arg(&args, key)
            .map(|v| picard_analysis::args::is_true(&v))
            .unwrap_or(default)
    };
    // `public boolean REMOVE_DUPLICATE_INFORMATION = true`.
    let remove_duplicate_information = flag("REMOVE_DUPLICATE_INFORMATION=", true);
    let allow_multiple = flag("ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS=", false);
    let stop_after = arg(&args, "STOP_AFTER=")
        .map(|v| v.parse::<usize>())
        .transpose()?;

    // `customCommandLineValidation`, which Barclay prints after the usage block.
    if !(0.0..=1.0).contains(&fraction) {
        eprintln!("{}", fraction_out_of_range_message(fraction));
        std::process::exit(1);
    }

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let opened = picard_analysis::bam_stream::Input::open(&input)?;
    let streamed = opened.streams();
    let (header, records) = opened.into_parts();
    // The first of the tool's two passes, which keeps only where each read was.
    let records: Vec<_> = records.map(|rec| rec.read_name).collect();

    // `checkProgramRecords`, before either pass: a file this tool has already downsampled is
    // refused outright unless the caller has said to go ahead anyway.
    if !allow_multiple {
        if let Some(previous) = header
            .programs
            .iter()
            .find(|pg| pg.attributes.get("PN") == Some(PG_PROGRAM_NAME))
        {
            eprintln!(
                "Exception in thread \"main\" picard.PicardException: Found previous Program \
                 Record that indicates that this file has been downsampled already with this \
                 program. Operation not supported! Previous PG: {}",
                previous.id
            );
            std::process::exit(1);
        }
    }

    // `getSamRecordLocation`: a name the regex does not parse leaves the location at its own
    // defaults, which are -1 and not zero.
    let locations: Vec<PhysicalLocation> = records
        .iter()
        .map(|read_name| {
            let parsed =
                picard_analysis::mark_duplicates::location_with(read_name, &read_name_regex);
            if parsed.known {
                PhysicalLocation {
                    tile: parsed.tile,
                    x: parsed.x,
                    y: parsed.y,
                }
            } else {
                PhysicalLocation {
                    tile: -1,
                    x: -1,
                    y: -1,
                }
            }
        })
        .collect();

    let mask = keep(&locations, fraction, stop_after);

    // `outputSamRecords`: the tool's own `@PG`, id made unique against the input's.
    let mut header = header;
    let pg = picard_analysis::command_line::pg_record(&header, PG_PROGRAM_NAME);
    header.programs.push(pg);
    // The second pass, written as it is read: a BAM whatever the name, and no file when a record
    // fails.
    let records = picard_analysis::bam_stream::Input::reopen(&input, streamed)?
        .into_parts()
        .1;
    let mut out = Output::bam(&output, &header, false)?;
    for (mut record, keep_it) in records.zip(mask) {
        if !keep_it {
            continue;
        }
        if remove_duplicate_information {
            record.flags &= !DUPLICATE_READ;
        }
        out.write(&record)?;
    }
    // `cram_output::write`: a CRAM takes `CREATE_INDEX` from the command line, a BAM no index.
    out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
    Ok(())
}
