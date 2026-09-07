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

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::jumping_library::{collect, unsorted_message, Arguments, Pair};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const DUPLICATE: u16 = 0x400;

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

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, records) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

    if header.attributes.get("SO") != Some("coordinate") {
        let name = input.rsplit('/').next().unwrap_or(&input);
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: {}",
            unsorted_message(name)
        );
        std::process::exit(1);
    }

    // "We're getting all our info from the first of each pair."
    let pairs: Vec<(Pair, bool)> = records
        .iter()
        .filter(|rec| rec.flags & READ_PAIRED != 0 && rec.flags & FIRST_OF_PAIR != 0)
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
            // `SamPairUtil.getPairOrientation` does not ask which START is further right: it
            // compares the two 5' ENDS, and which coordinate stands for which end depends on the
            // strand this read is on. Reading it as "the mate starts later" put three pairs in the
            // wrong bucket on this corpus.
            //
            //   read forward:  positive 5' is this read's start, negative 5' is the mate's start
            //                  plus the insert size;
            //   read reverse:  positive 5' is the mate's start, negative 5' is this read's END.
            //
            // `FR` is "positive before negative", and the library asks the question the other way
            // round -- whether the mate lies to the right -- so the two cases invert.
            let alignment_end = rec.alignment_start + rec.cigar.reference_length() as i32 - 1;
            let mate_is_to_the_right = if rec.flags & READ_REVERSE == 0 {
                rec.alignment_start < rec.mate_alignment_start + rec.inferred_insert_size
            } else {
                rec.mate_alignment_start >= alignment_end
            };
            (pair, mate_is_to_the_right)
        })
        .collect();

    let metrics = collect(&pairs, &arguments);
    let mut file = MetricsFile::new();
    file.add_header("CollectJumpingLibraryMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    file.add_metric(&metrics);
    std::fs::write(&output, file.write())?;
    Ok(())
}
