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

use std::sync::Arc;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use picard_analysis::args::arg;
use picard_analysis::bam_stream::{Input, Records};
use picard_analysis::mark_duplicates::{
    base_quality_score, Options, Record, ScoringStrategy, Span, TaggingPolicy,
};
use picard_analysis::mate_cigar_duplicates::{
    mark_with_mate_cigar, MateCigarOptions, SortOrder as MateSortOrder,
};
use picard_analysis::stream_output::Output;

const DUPLICATE_READ: u16 = 0x400;

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("MarkDuplicatesWithMateCigar");
    picard_analysis::output::install();
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
    let inputs = picard_analysis::args::inputs(&args)?;
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

    // ASSUME_SORT_ORDER is written into the header by `openInputs`, so it decides this check as
    // surely as the file does: a coordinate-sorted file assumed to be anything else is refused,
    // and the message says nothing about the assumption.
    // `ASSUME_SORTED=true` stands for `ASSUME_SORT_ORDER=coordinate`.
    let assumed = arg(&args, "ASSUME_SORT_ORDER=")
        .or_else(|| arg(&args, "ASO="))
        .or_else(|| {
            arg(&args, "ASSUME_SORTED=")
                .or_else(|| arg(&args, "AS="))
                .is_some_and(|v| picard_analysis::args::is_true(&v))
                .then(|| "coordinate".to_string())
        });
    // Read twice, as MarkDuplicates reads it: once for the marking, which keeps only what it
    // needs of each record, and once for the writing. One input is streamed, and its first read
    // checks the whole file, so a damaged one is refused here as reading it whole refused it.
    // Several are merged in memory by `openInputs`, once per pass.
    let open =
        |streamed: Option<bool>| -> Result<(SamHeader, Records, bool), Box<dyn std::error::Error>> {
            if let [input] = inputs.as_slice() {
                let opened = match streamed {
                    None => Input::open(input)?,
                    Some(streamed) => Input::reopen(input, streamed)?,
                };
                let streams = opened.streams();
                let (mut header, records) = opened.into_parts();
                if let Some(order) = assumed.as_deref() {
                    header.set_sort_order(order);
                }
                return Ok((header, records, streams));
            }
            let read = |input: &String| picard_analysis::metrics_cli::read_input(input);
            let (header, records) = picard_analysis::merge_sam_files::open_inputs(
                inputs.iter().map(read).collect(),
                assumed.as_deref(),
            )
            .unwrap_or_else(|thrown| picard_analysis::metrics_cli::thrown(&thrown));
            // Only an order this tool refuses before reading a record makes reading one throw.
            let records = records.unwrap_or_default();
            Ok((header, Records::Decoded(records.into_iter(), None), false))
        };
    let (header, records, streamed) = open(None)?;
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
    // the header, which is what `closeEnough` compares. One shared string per read group, so the
    // records of a library hold a pointer each rather than a copy.
    let unknown_library: Arc<str> = Arc::from("Unknown Library");
    let libraries: Vec<Arc<str>> = header
        .read_groups
        .iter()
        .map(|group| match group.attributes.get("LB") {
            Some(library) => Arc::from(library),
            None => unknown_library.clone(),
        })
        .collect();
    let library_of = |record: &BamRecord| -> (Arc<str>, i32) {
        let id = match record.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(s)) => Some(s.as_str()),
            _ => None,
        };
        match id.and_then(|id| header.read_groups.iter().position(|group| group.id == id)) {
            Some(position) => (libraries[position].clone(), position as i32),
            None => (unknown_library.clone(), -1),
        }
    };

    let marked_records: Vec<Record> = records
        .map(|record| {
            let record = &record;
            let (library, read_group) = library_of(record);
            Record {
                name: record.read_name.clone(),
                flags: record.flags,
                reference_index: record.reference_index,
                alignment_start: record.alignment_start,
                span: Span::of(&record.cigar),
                read_length: record.read_bases.len() as u32,
                base_quality_score: base_quality_score(&record.base_qualities),
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
                mate_span: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text)
                        .ok()
                        .as_ref()
                        .map(Span::of),
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

    let mut out = Output::bam(&output, &out_header, false)?;
    for (index, mut written) in open(Some(streamed))?.1.enumerate() {
        if !marking.written[index] {
            continue;
        }
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
        out.write(&written)?;
    }
    // A CRAM is indexed under CREATE_INDEX, a BAM never was here.
    out.finish(picard_analysis::cram_output::create_index_argument(), false)?;

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
