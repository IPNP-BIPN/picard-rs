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

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::mark_duplicates::{mark, Options, Record, ScoringStrategy, TaggingPolicy};

const DUPLICATE_READ: u16 = 0x400;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

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
        let mate = Record {
            alignment_start: record.mate_alignment_start,
            cigar: mate_cigar,
            ..marked.clone()
        };
        let mate_five = if marked.mate_reverse_strand() {
            mate.unclipped_end()
        } else {
            mate.unclipped_start()
        };
        marked.first_of_pair() == (read_five <= mate_five)
    };
    out.push_str(if top { "/A" } else { "/B" });
    Ok(out)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let metrics_file = arg(&args, "METRICS_FILE=").or_else(|| arg(&args, "M="));
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);

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
        max_optical_duplicate_set_size: std::env::args()
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

    // `openInputs` sets the header's own sort order to ASSUME_SORT_ORDER, so by the time the check
    // below reads "the header", the header is already saying what it was told to say -- which is
    // why the refusal prints the assumed order twice.
    // `ASSUME_SORTED=true` is `ASSUME_SORT_ORDER=coordinate` by the time `openInputs` reads it.
    let assumed = arg(&args, "ASSUME_SORT_ORDER=")
        .or_else(|| arg(&args, "ASO="))
        .or_else(|| {
            (arg(&args, "ASSUME_SORTED=")
                .or_else(|| arg(&args, "AS="))
                .as_deref()
                == Some("true"))
            .then(|| "coordinate".to_string())
        });
    let mut header = header;
    if let Some(order) = &assumed {
        header.set_sort_order(order);
    }
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
                mate_cigar: None,
                mate_alignment_start: record.mate_alignment_start,
            }
        })
        .collect();

    let marking = mark(&marked_records, &options);

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
    // The id every read is tagged with: `PROGRAM_RECORD_ID` made unique against the header's own
    // program records, or no tag at all when it is `null`. The corpus's reads carry no `PG`, so the
    // one chain `getChainedPgIds` builds is the one from `null`.
    let pg_id: Option<String> = match arg(&args, "PROGRAM_RECORD_ID=").as_deref() {
        Some("null") => None,
        Some(id) => Some(id.to_string()),
        None => Some("MarkDuplicates".to_string()),
    }
    .map(|id| picard_analysis::merge_sam_files::non_colliding_pg_id(&header, &id));

    // htsjdk's own check in `SAMFileWriterImpl.addAlignment`, which is where a file assumed to be
    // coordinate-sorted and not is caught. The writer is presorted, so nothing re-sorts and the
    // record that goes backwards is the one that throws.
    let coordinate_out = out_header.attributes.get("SO") == Some("coordinate");
    let mut previous: Option<(i32, i32)> = None;

    let mut writer = BamWriter::new(Vec::new(), &out_header).map_err(|e| format!("{e:?}"))?;
    for (index, record) in records.iter().enumerate() {
        // `UmiUtil.setMolecularIdentifier`, reached for every record before the ones that are
        // removed are skipped, so a record that cannot be given one stops the run wherever it is.
        let molecular_identifier = match &molecular_identifier_tag {
            Some(_) => {
                match molecular_identifier(&header, record, &marked_records[index], duplex) {
                    Ok(value) => Some(value),
                    Err(thrown) => {
                        eprintln!("Exception in thread \"main\" {thrown}");
                        std::process::exit(1);
                    }
                }
            }
            None => None,
        };
        if !marking.written[index] {
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
                    eprintln!(
                        "Exception in thread \"main\" java.lang.IllegalArgumentException: \
                         Alignments added out of order in SAMFileWriterImpl.addAlignment for \
                         file://{output}. Sort order is coordinate. Offending records are at [{}] \
                         and [{}]",
                        key(here.0, here.1),
                        key(here.0, here.1)
                    );
                    std::process::exit(1);
                }
            }
            previous = Some(here);
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
        if tag_duplicate_set_members && record.flags & (0x4 | 0x100 | 0x800) == 0 {
            if let Some((representative, size)) = marking.duplicate_set[index] {
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
        if let (true, Some(id)) = (add_pg_tag, &pg_id) {
            written
                .tags
                .insert(Tag::new(b"PG"), TagValue::Str(id.clone()));
        }
        writer.write(&written).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;

    // The metrics file is written because the tool requires the argument, and its content is not
    // what this array compares: the OUTPUT is.
    if let Some(path) = metrics_file {
        let text = picard_analysis::mark_duplicates::duplication_metrics_text(
            "MarkDuplicates",
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
