//! `MarkDuplicates` as a runnable binary: the covering array's port side.
//!
//! The library decides which records are duplicates; this is the file around that decision, and
//! the file is where the arguments live.
//!
//! `ASSUME_SORT_ORDER` is the one that shapes everything. It is not a hint: the reference reads it
//! instead of the header, refuses anything that is neither coordinate nor queryname, and, for
//! queryname, rewrites the OUTPUT header to `SO:unknown GO:query` because a queryname-grouped file
//! is not a sorted one. The writer is always presorted, so assuming coordinate over a
//! queryname-sorted file does not re-sort anything -- it reaches htsjdk's own check and fails
//! there, with htsjdk's wording.
//!
//! What is not ported here is stated rather than implied: `TAG_DUPLICATE_SET_MEMBERS` needs the
//! representative-read index of every duplicate set, which is a second sorting collection in the
//! reference and a surface of its own, so the repository declares no value for it and the array
//! does not vary it.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use picard_analysis::args::arg;
use picard_analysis::bam_stream::{Input, Records};
use picard_analysis::mark_duplicates::{
    base_quality_score, Flags, Options, ReadEndsBuilder, Record, ScoringStrategy, Span,
    TaggingPolicy, Tally,
};
use picard_analysis::stream_output::Output;
use std::sync::Arc;

const DUPLICATE_READ: u16 = 0x400;

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    }
}

/// `UmiUtil.setMolecularIdentifier(rec, "", tag, duplex)`: the contig, the position of the end the
/// pair shares, a `/`, no UMI, and with `DUPLEX_UMI` the strand, `/A` on top and `/B` below.
///
/// The strand needs the mate's unclipped 5' end, which comes from its `MC` tag, so a pair on one
/// contig without one is htsjdk's refusal; and asking an unpaired read about its mate is Java's.
fn molecular_identifier(
    header: &SamHeader,
    record: &BamRecord,
    marked: &Record,
    duplex: bool,
) -> Result<String, String> {
    // `SAMRecord.getContig()` is null for an unmapped read, placed or not, and the builder
    // appends that as the word.
    let contig = if marked.unmapped() {
        "null".to_string()
    } else {
        usize::try_from(record.reference_index)
            .ok()
            .and_then(|i| header.sequences.get(i))
            .map_or("null".to_string(), |s| s.name.clone())
    };
    let mut out = picard_analysis::umi_duplicates::molecular_identifier(
        &contig,
        record.alignment_start,
        record.mate_alignment_start,
        marked.reverse_strand(),
        "",
    );
    if !duplex {
        return Ok(out);
    }
    // `getReadUnmappedFlag() || getMateUnmappedFlag()`: the second is only asked of a mapped read,
    // and asking it of an unpaired one throws.
    if marked.unmapped() {
        return Ok(out);
    }
    if !marked.paired() {
        return Err(
            "java.lang.IllegalStateException: Inappropriate call if not paired read".to_string(),
        );
    }
    if marked.mate_unmapped() {
        return Ok(out);
    }
    let top = if record.reference_index != record.mate_reference_index {
        marked.first_of_pair() == (record.reference_index < record.mate_reference_index)
    } else {
        let read_five = if marked.reverse_strand() {
            marked.unclipped_end()
        } else {
            marked.unclipped_start()
        };
        let mate_cigar = match record.tags.get(Tag::new(b"MC")) {
            Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
            _ => None,
        };
        let Some(mate_cigar) = mate_cigar else {
            return Err(format!(
                "htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found: {}",
                picard_analysis::view_sam::sam_string(header, record)
            ));
        };
        let mate = Span::of(&mate_cigar);
        let mate_five = if marked.mate_reverse_strand() {
            mate.unclipped_end(record.mate_alignment_start)
        } else {
            mate.unclipped_start(record.mate_alignment_start)
        };
        marked.first_of_pair() == (read_five <= mate_five)
    };
    out.push_str(if top { "/A" } else { "/B" });
    Ok(out)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("MarkDuplicates");
    picard_analysis::output::install();
    // `READ_NAME_REGEX`: `null` turns the locations, and with them the optical duplicates, off.
    let read_name_regex = picard_analysis::mark_duplicates::ReadNameRegex::from_arg(
        picard_analysis::args::process()
            .into_iter()
            .find_map(|a| a.strip_prefix("READ_NAME_REGEX=").map(str::to_string))
            .as_deref(),
    )
    .unwrap_or_else(|thrown| {
        eprintln!("Exception in thread \"main\" {thrown}");
        std::process::exit(1);
    });
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
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
            None | Some("SUM_OF_BASE_QUALITIES") => ScoringStrategy::SumOfBaseQualities,
            Some("TOTAL_MAPPED_REFERENCE_LENGTH") => ScoringStrategy::TotalMappedReferenceLength,
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
        max_optical_duplicate_set_size: picard_analysis::args::process()
            .into_iter()
            .find_map(|a| {
                a.strip_prefix("MAX_OPTICAL_DUPLICATE_SET_SIZE=")
                    .and_then(|v| v.parse().ok())
            })
            .unwrap_or(300_000),
        barcode_tag: arg(&args, "BARCODE_TAG="),
        assume_mate_cigar: false,
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    // `openInputs` sets the header's own sort order to ASSUME_SORT_ORDER, so by the time the check
    // below reads "the header", the header is already saying what it was told to say -- which is
    // why the refusal prints the assumed order twice.
    // `ASSUME_SORTED=true` is `ASSUME_SORT_ORDER=coordinate` by the time `openInputs` reads it.
    let assumed = arg(&args, "ASSUME_SORT_ORDER=")
        .or_else(|| arg(&args, "ASO="))
        .or_else(|| {
            arg(&args, "ASSUME_SORTED=")
                .or_else(|| arg(&args, "AS="))
                .is_some_and(|v| picard_analysis::args::is_true(&v))
                .then(|| "coordinate".to_string())
        });
    // Read twice, as Picard reads it: once for the marking, which keeps only what it needs of each
    // record, and once for the writing. One input is streamed, and its first read checks the whole
    // file, so a damaged one is refused here as reading it whole refused it. Several are merged in
    // memory by `openInputs`, once per pass.
    // shortcut: several inputs are still read whole and merged in memory; stream the merge once a
    // multi-input run on large files matters.
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
            // `buildSortedReadEndLists` reads every record before the order is checked.
            let records =
                records.unwrap_or_else(|thrown| picard_analysis::metrics_cli::thrown(&thrown));
            Ok((header, Records::Decoded(records.into_iter(), None), false))
        };
    let (header, records, streamed) = open(None)?;
    let header_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    if !matches!(header_order.as_str(), "coordinate" | "queryname") {
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: This program requires input that \
             are either coordinate or query sorted (according to the header, or at least \
             ASSUME_SORT_ORDER and the content.) Found ASSUME_SORT_ORDER={} and header \
             sortorder={header_order}",
            assumed.as_deref().unwrap_or("null")
        );
        std::process::exit(1);
    }

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

    // `getChainedPgIds`'s program ids are gathered on the same walk: one new `@PG` for every
    // program id the first pass saw on a record (no `PG` at all counting as one).
    let mut seen: Vec<Option<String>> = Vec::new();
    // `buildSortedReadEndLists` keys every primary mapped read whose mate is mapped on
    // `rec.getReadGroup().getReadGroupId()`, which is null with no RG tag or one the header does
    // not declare. The pass stops at the first unplaced unmapped read when coordinate-sorted.
    let mut past_unplaced = false;
    // What either pass reads of a record, built again on the second rather than kept.
    let summary = |record: &BamRecord| -> Record {
        let (library, read_group) = library_of(record);
        let barcode = options.barcode_tag.as_ref().and_then(|tag| {
            match record
                .tags
                .get(Tag::new(tag.as_bytes().try_into().unwrap_or(b"RX")))
            {
                Some(TagValue::Str(value)) => Some(value.clone()),
                _ => None,
            }
        });
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
            barcode,
            existing_dt: string_tag(record, b"DT"),
            mate_span: None,
            mate_alignment_start: record.mate_alignment_start,
        }
    };
    let mut ends = ReadEndsBuilder::new(&options, header_order == "queryname");
    for record in records {
        let id = string_tag(&record, b"PG");
        if !seen.contains(&id) {
            seen.push(id);
        }
        let read_group = library_of(&record).1;
        past_unplaced |=
            record.flags & 0x4 != 0 && record.reference_index == -1 && header_order == "coordinate";
        if !past_unplaced
            && record.flags & (0x1 | 0x4 | 0x8 | 0x100 | 0x800) == 0x1
            && read_group == -1
        {
            eprintln!(
                    "Exception in thread \"main\" java.lang.NullPointerException: Cannot invoke \
                     \"htsjdk.samtools.SAMReadGroupRecord.getReadGroupId()\" because the return value of \
                     \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null"
                );
            std::process::exit(1);
        }
        ends.push(&summary(&record));
    }
    let mut decisions = ends.decide();

    // `createOutHeader`: assuming queryname says the file is queryname GROUPED, which is not a
    // sort order, so the output says so.
    let mut out_header = header.clone();
    if assumed.as_deref() == Some("queryname") {
        out_header.set_sort_order("unknown");
        out_header.set_group_order("query");
    }
    // `COMMENT.forEach(outputHeader::addComment)`.
    for comment in args
        .iter()
        .filter_map(|a| a.strip_prefix("COMMENT=").or_else(|| a.strip_prefix("CO=")))
    {
        out_header.comments.push(format!("@CO\t{comment}"));
    }
    let tag_duplicate_set_members = flag("TAG_DUPLICATE_SET_MEMBERS=", false);
    let molecular_identifier_tag = arg(&args, "MOLECULAR_IDENTIFIER_TAG=")
        .filter(|tag| tag != "null" && options.barcode_tag.is_some());
    let duplex = flag("DUPLEX_UMI=", false);
    let add_pg_tag = flag("ADD_PG_TAG_TO_READS=", true);
    // `getChainedPgIds`: the ids seen above, in the `HashSet`'s order, each chained to the id it
    // replaces; none when `PROGRAM_RECORD_ID` is `null`.
    let group = picard_analysis::command_line::ProgramGroup::of("MarkDuplicates");
    let chained = picard_analysis::command_line::chained_pg_ids(
        &mut out_header,
        &picard_analysis::command_line::hash_set_order(seen),
        &group,
        "MarkDuplicates",
    );

    // htsjdk's own check in `SAMFileWriterImpl.addAlignment`, which is where a file assumed to be
    // coordinate-sorted and not is caught. The writer is presorted, so nothing re-sorts and the
    // record that goes backwards is the one that throws.
    let coordinate_out = out_header.attributes.get("SO") == Some("coordinate");
    let mut previous: Option<(i32, i32)> = None;

    // Picard's writer is opened before the second pass, and a failure in it leaves no file.
    let mut writer = Output::bam(&output, &out_header, false)?;
    let optical_pairs = std::mem::take(&mut decisions.optical_pairs);
    let all_sets = std::mem::take(&mut decisions.all_sets);
    let optical_sets = std::mem::take(&mut decisions.optical_sets);
    let non_optical_sets = std::mem::take(&mut decisions.non_optical_sets);
    let mut flags = Flags::new(decisions, header_order == "queryname");
    let mut tally = Tally::default();
    let (_, records, _) = open(Some(streamed))?;
    for record in records {
        let marked = summary(&record);
        let (duplicate, optical, duplicate_set) = flags.next(
            &marked.name,
            marked.secondary_or_supplementary() || marked.unmapped(),
        );
        let (duplicate_type, kept) = tally.add(&marked, &options, duplicate, optical);
        // `UmiUtil.setMolecularIdentifier`, reached for every record before the ones that are
        // removed are skipped, so a record that cannot be given one stops the run wherever it is.
        let molecular_identifier = match &molecular_identifier_tag {
            Some(_) => match molecular_identifier(&header, &record, &marked, duplex) {
                Ok(value) => Some(value),
                Err(thrown) => picard_analysis::metrics_cli::thrown(&thrown),
            },
            None => None,
        };
        if !kept {
            continue;
        }
        if coordinate_out {
            // `SAMSortOrderChecker.isSorted` under `coordinate`, whose key is the contig NAME and
            // the start, and whose message prints that key for both records.
            let here = (record.reference_index, record.alignment_start);
            if let Some(before) = previous {
                // `SAMRecordCoordinateComparator.fileOrderCompare`, in which an UNMAPPED record
                // (reference index -1) sorts LAST rather than first: comparing the pair as plain
                // tuples finds a violation one record early, and the message names the record it
                // found rather than the record htsjdk finds.
                let file_order = |left: (i32, i32), right: (i32, i32)| -> i32 {
                    if left.0 == -1 {
                        return if right.0 == -1 { 0 } else { 1 };
                    }
                    if right.0 == -1 {
                        return -1;
                    }
                    if left.0 != right.0 {
                        return left.0 - right.0;
                    }
                    left.1 - right.1
                };
                if file_order(before, here) > 0 {
                    let key = |reference_index: i32, start: i32| -> String {
                        let name = usize::try_from(reference_index)
                            .ok()
                            .and_then(|i| out_header.sequences.get(i))
                            .map(|s| s.name.as_str())
                            .unwrap_or("*");
                        format!("{name}:{start}")
                    };
                    // Both keys are the OFFENDING record, and that is htsjdk's, not a slip here:
                    // `isSorted` sets `prev = rec` before returning false, so the message's
                    // `getPreviousRecord()` returns the record that has just been rejected.
                    picard_analysis::metrics_cli::thrown(&format!(
                        "java.lang.IllegalArgumentException: Alignments added out of order in \
                         SAMFileWriterImpl.addAlignment for file://{output}. Sort order is \
                         coordinate. Offending records are at [{}] and [{}]",
                        key(here.0, here.1),
                        key(here.0, here.1)
                    ));
                }
            }
            previous = Some(here);
        }

        let flags = record.flags;
        let mut written = record;
        written.flags &= !DUPLICATE_READ;
        if duplicate {
            written.flags |= DUPLICATE_READ;
        }
        written.tags.remove(Tag::new(b"DT"));
        if let Some(code) = duplicate_type {
            written.tags.insert(Tag::new(b"DT"), TagValue::Str(code));
        }
        if tag_duplicate_set_members && flags & (0x4 | 0x100 | 0x800) == 0 {
            if let Some((representative, size)) = duplicate_set {
                written
                    .tags
                    .insert(Tag::new(b"DI"), TagValue::Int(representative as i64));
                written
                    .tags
                    .insert(Tag::new(b"DS"), TagValue::Int(size as i64));
            }
        }
        if let (Some(tag), Some(value)) = (&molecular_identifier_tag, molecular_identifier) {
            if let Ok(key) = <[u8; 2]>::try_from(tag.as_bytes()) {
                written.tags.insert(Tag::new(&key), TagValue::Str(value));
            }
        }
        // `rec.setAttribute(PG, chainedPgIds.get(rec.getStringAttribute(PG)))`.
        if let (true, Some(chained)) = (add_pg_tag, &chained) {
            let old = string_tag(&written, b"PG");
            match picard_analysis::command_line::chained_id(chained, old.as_deref()) {
                Some(id) => {
                    written
                        .tags
                        .insert(Tag::new(b"PG"), TagValue::Str(id.to_string()));
                }
                None => {
                    written.tags.remove(Tag::new(b"PG"));
                }
            }
        }
        writer.write(&written)?;
    }
    writer.finish(picard_analysis::cram_output::create_index_argument(), false)?;

    // The metrics file is written because the tool requires the argument, and its content is not
    // what this array compares: the OUTPUT is.
    if let Some(path) = metrics_file {
        let text = picard_analysis::mark_duplicates::duplication_metrics_text(
            "MarkDuplicates",
            &tally.finish(Some(&optical_pairs)),
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
            &all_sets,
            &optical_sets,
            &non_optical_sets,
        );
        std::fs::write(path, text)?;
    }
    Ok(())
}
