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

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::mark_duplicates::location;
use picard_analysis::position_based_downsample_sam::{
    fraction_out_of_range_message, keep, PhysicalLocation, PG_PROGRAM_NAME,
};

const DUPLICATE_READ: u16 = 0x400;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);
    let remove_duplicate_information = flag("REMOVE_DUPLICATE_INFORMATION=", false);
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

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

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
        .map(|rec| {
            let parsed = location(&rec.read_name);
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

    let mut writer = BamWriter::new(Vec::new(), &header).map_err(|e| format!("{e:?}"))?;
    for (record, keep_it) in records.iter().zip(&mask) {
        if !keep_it {
            continue;
        }
        let mut record = record.clone();
        if remove_duplicate_information {
            record.flags &= !DUPLICATE_READ;
        }
        writer.write(&record).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;
    Ok(())
}
