//! `MarkDuplicatesWithMateCigar` as a runnable binary: the covering array's port side.
//!
//! The sibling of `MarkDuplicates` that reads the mate's position off the `MC` tag instead of
//! waiting for the mate, which is why one pass over a coordinate-sorted file is enough for it and
//! why almost every row of its array is a refusal:
//!
//! * a file that is not coordinate-sorted is refused before a record is read;
//! * `DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES` is refused outright, whatever the file
//!   holds, because the two ends of a pair are scored separately here and summing base qualities
//!   would let them disagree;
//! * a pair with no `MC` is skipped, and `SKIP_PAIRS_WITH_NO_MATE_CIGAR=false` turns that skip
//!   into a refusal naming the read.
//!
//! The corpus gained `mate_cigar.bam` for this: without a file that carries the tag, every row
//! measured the refusal and nothing else.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::args::arg;
use picard_analysis::mark_duplicates::{Options, Record, ScoringStrategy, TaggingPolicy};
use picard_analysis::mate_cigar_duplicates::{
    mark_with_mate_cigar, MateCigarOptions, SortOrder as MateSortOrder,
};

const DUPLICATE_READ: u16 = 0x400;

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("MarkDuplicatesWithMateCigar");
    // `READ_NAME_REGEX`: `null` turns the locations, and with them the optical duplicates, off.
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
    let input = picard_analysis::args::single_input(&args)?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let metrics_file = arg(&args, "METRICS_FILE=").or_else(|| arg(&args, "M="));
    let flag = |key: &str, default: bool| {
        arg(&args, key)
            .map(|v| picard_analysis::args::is_true(&v))
            .unwrap_or(default)
    };

    let options = Options {
        scoring: match arg(&args, "DUPLICATE_SCORING_STRATEGY=").as_deref() {
            // `AbstractMarkDuplicatesCommandLineProgram.DUPLICATE_SCORING_STRATEGY =
            // TOTAL_MAPPED_REFERENCE_LENGTH`: an unset one is the one strategy this tool accepts.
            Some("SUM_OF_BASE_QUALITIES") => ScoringStrategy::SumOfBaseQualities,
            None | Some("TOTAL_MAPPED_REFERENCE_LENGTH") => {
                ScoringStrategy::TotalMappedReferenceLength
            }
            Some("RANDOM") => ScoringStrategy::Random,
            Some(other) => {
                return Err(format!("unknown DUPLICATE_SCORING_STRATEGY: {other}").into())
            }
        },
        remove_duplicates: flag("REMOVE_DUPLICATES=", false),
        remove_sequencing_duplicates: flag("REMOVE_SEQUENCING_DUPLICATES=", false),
        tagging_policy: match arg(&args, "TAGGING_POLICY=").as_deref() {
            None | Some("DontTag") => TaggingPolicy::DontTag,
            Some("OpticalOnly") => TaggingPolicy::OpticalOnly,
            Some("All") => TaggingPolicy::All,
            Some(other) => return Err(format!("unknown TAGGING_POLICY: {other}").into()),
        },
        clear_dt: flag("CLEAR_DT=", true),
        optical_duplicate_pixel_distance: arg(&args, "OPTICAL_DUPLICATE_PIXEL_DISTANCE=")
            .map(|v| v.parse::<i32>())
            .transpose()?
            .unwrap_or(100),
        parse_read_names: !read_name_regex.is_off(),
        read_name_regex,
        max_optical_duplicate_set_size: std::env::args()
            .find_map(|a| {
                a.strip_prefix("MAX_OPTICAL_DUPLICATE_SET_SIZE=")
                    .and_then(|v| v.parse().ok())
            })
            .unwrap_or(300_000),
        barcode_tag: arg(&args, "BARCODE_TAG="),
        assume_mate_cigar: true,
    };

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

    // ASSUME_SORT_ORDER is written into the header by `openInputs`, so it decides this check as
    // surely as the file does: a coordinate-sorted file assumed to be anything else is refused,
    // and the message says nothing about the assumption.
    let mut header = header;
    // `ASSUME_SORTED=true` stands for `ASSUME_SORT_ORDER=coordinate`.
    let assumed = arg(&args, "ASSUME_SORT_ORDER=")
        .or_else(|| arg(&args, "ASO="))
        .or_else(|| {
            arg(&args, "ASSUME_SORTED=")
                .or_else(|| arg(&args, "AS="))
                .is_some_and(|v| picard_analysis::args::is_true(&v))
                .then(|| "coordinate".to_string())
        });
    if let Some(order) = assumed {
        header.set_sort_order(&order);
    }
    let header_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    let order = match header_order.as_str() {
        "coordinate" => MateSortOrder::Coordinate,
        "queryname" => MateSortOrder::Queryname,
        _ => MateSortOrder::Unsorted,
    };

    // The read group's library, which is what a duplicate set is cut by, and the group's index in
    // the header, which is what `closeEnough` compares.
    let library_of = |record: &BamRecord| -> (String, i32) {
        let id = string_tag(record, b"RG");
        match id.and_then(|id| {
            header
                .read_groups
                .iter()
                .position(|group| group.id == id)
                .map(|position| (position, &header.read_groups[position]))
        }) {
            Some((position, group)) => (
                group
                    .attributes
                    .get("LB")
                    .unwrap_or("Unknown Library")
                    .to_string(),
                position as i32,
            ),
            None => ("Unknown Library".to_string(), -1),
        }
    };

    let marked_records: Vec<Record> = records
        .iter()
        .map(|record| {
            let (library, read_group) = library_of(record);
            Record {
                name: record.read_name.clone(),
                flags: record.flags,
                reference_index: record.reference_index,
                alignment_start: record.alignment_start,
                cigar: record.cigar.clone(),
                qualities: record.base_qualities.clone(),
                mate_reference_index: record.mate_reference_index,
                library,
                read_group,
                barcode: options.barcode_tag.as_ref().and_then(|tag| {
                    match record
                        .tags
                        .get(Tag::new(tag.as_bytes().try_into().unwrap_or(b"RX")))
                    {
                        Some(TagValue::Str(value)) => Some(value.clone()),
                        _ => None,
                    }
                }),
                existing_dt: string_tag(record, b"DT"),
                mate_cigar: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
                    _ => None,
                },
                mate_alignment_start: record.mate_alignment_start,
            }
        })
        .collect();

    let mate_options = MateCigarOptions {
        base: options.clone(),
        skip_pairs_with_no_mate_cigar: flag("SKIP_PAIRS_WITH_NO_MATE_CIGAR=", true),
        minimum_distance: arg(&args, "MINIMUM_DISTANCE=")
            .map(|v| v.parse::<i32>())
            .transpose()?
            .unwrap_or(-1),
    };
    let marking = match mark_with_mate_cigar(&marked_records, order, &mate_options) {
        Ok(marking) => marking,
        Err(refusal) => {
            eprintln!(
                "Exception in thread \"main\" {}: {}",
                refusal.exception(),
                refusal.message()
            );
            std::process::exit(1);
        }
    };

    let mut out_header = header.clone();
    // `COMMENT.forEach(outputHeader::addComment)`.
    for comment in args
        .iter()
        .filter_map(|a| a.strip_prefix("COMMENT=").or_else(|| a.strip_prefix("CO=")))
    {
        out_header.comments.push(format!("@CO\t{comment}"));
    }
    let add_pg_tag = flag("ADD_PG_TAG_TO_READS=", true);
    // `setPGIdsSeen`: one pass cannot know which program ids the reads carry, so the chain is
    // built from the header's, less every id another record names as its `PP`.
    let previous: Vec<&str> = out_header
        .programs
        .iter()
        .filter_map(|pg| pg.attributes.get("PP"))
        .collect();
    let mut seen: Vec<Option<String>> = Vec::new();
    for pg in &out_header.programs {
        let id = Some(pg.id.clone());
        if !previous.contains(&pg.id.as_str()) && !seen.contains(&id) {
            seen.push(id);
        }
    }
    let chained = picard_analysis::command_line::chained_pg_ids(
        &mut out_header,
        &picard_analysis::command_line::hash_set_order(seen),
        &picard_analysis::command_line::ProgramGroup::of("MarkDuplicatesWithMateCigar"),
        "MarkDuplicatesWithMateCigar",
    );

    let mut writer = BamWriter::new(Vec::new(), &out_header).map_err(|e| format!("{e:?}"))?;
    for (index, record) in records.iter().enumerate() {
        if !marking.written[index] {
            continue;
        }
        let mut written = record.clone();
        written.flags &= !DUPLICATE_READ;
        if marking.duplicate[index] {
            written.flags |= DUPLICATE_READ;
        }
        written.tags.remove(Tag::new(b"DT"));
        if let Some(code) = &marking.duplicate_type[index] {
            written
                .tags
                .insert(Tag::new(b"DT"), TagValue::Str(code.clone()));
        }
        // `updateProgramRecord` only CHAINS: a record with no `PG`, or with one the chain does
        // not hold (an id some other record follows), is left as it is, with a warning.
        if let (true, Some(chained)) = (add_pg_tag, &chained) {
            let old = string_tag(&written, b"PG");
            if let Some(id) = old
                .as_deref()
                .and_then(|old| picard_analysis::command_line::chained_id(chained, Some(old)))
            {
                written
                    .tags
                    .insert(Tag::new(b"PG"), TagValue::Str(id.to_string()));
            }
        }
        writer.write(&written).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;

    // `finalizeAndWriteMetrics`, shared with `MarkDuplicates`.
    if let Some(path) = metrics_file {
        let text = picard_analysis::mark_duplicates::duplication_metrics_text(
            "MarkDuplicatesWithMateCigar",
            &marking.metrics,
            &header
                .read_groups
                .iter()
                .map(|group| {
                    group
                        .attributes
                        .get("LB")
                        .unwrap_or("Unknown Library")
                        .to_string()
                })
                .collect::<Vec<_>>(),
            &marking.all_sets,
            &marking.optical_sets,
            &marking.non_optical_sets,
        );
        std::fs::write(path, text)?;
    }
    Ok(())
}
