//! `SimpleMarkDuplicatesWithMateCigar` as a runnable binary: the covering array's port side.
//!
//! A `MarkDuplicates` subclass driven by htsjdk's duplicate-set iterator, which is what refuses a
//! record with no mate cigar: the refusal is htsjdk's `SAMException` and carries htsjdk's wording,
//! whatever `SKIP_PAIRS_WITH_NO_MATE_CIGAR` says. It also scores each end alone, as
//! `MarkDuplicates` does and unlike `MarkDuplicatesWithMateCigar`, because it inherits
//! `MarkDuplicates`'s scoring call.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::mark_duplicates::{Options, Record, ScoringStrategy, TaggingPolicy};
use picard_analysis::mate_cigar_duplicates::{
    simple_mark_with_mate_cigar, SortOrder as MateSortOrder,
};

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
        parse_read_names: true,
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

    // ASSUME_SORT_ORDER is written into the header by `openInputs`, so it decides this check as
    // surely as the file does: a coordinate-sorted file assumed to be anything else is refused,
    // and the message says nothing about the assumption.
    let mut header = header;
    if let Some(order) = arg(&args, "ASSUME_SORT_ORDER=").or_else(|| arg(&args, "ASO=")) {
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

    let marking = match simple_mark_with_mate_cigar(&marked_records, order, &options) {
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

    let out_header = header.clone();
    let add_pg_tag = flag("ADD_PG_TAG_TO_READS=", true);

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
        // `updateProgramRecord` only CHAINS: a record that already carries a `PG` gets the new id,
        // and a record with none gets a warning and nothing else. That is the opposite of
        // `MarkDuplicates`, which fills its chain from the records themselves and so stamps every
        // one of them; here the chain comes from the HEADER's program records, and this corpus has
        // none, so `ADD_PG_TAG_TO_READS` changes not a byte.
        if add_pg_tag {
            if let Some(TagValue::Str(existing)) = written.tags.get(Tag::new(b"PG")) {
                let chained = existing.clone();
                if header.programs.iter().any(|pg| pg.id == chained) {
                    written.tags.insert(Tag::new(b"PG"), TagValue::Str(chained));
                }
            }
        }
        writer.write(&written).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;

    // The metrics file is written because the tool requires the argument, and its content is not
    // what this array compares: the OUTPUT is.
    if let Some(path) = metrics_file {
        let mut text = String::from("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n");
        text.push_str(
            "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\t\
             UNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t\
             READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n",
        );
        for m in &marking.metrics {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                m.library,
                m.unpaired_reads_examined,
                m.read_pairs_examined,
                m.secondary_or_supplementary,
                m.unmapped_reads,
                m.unpaired_read_duplicates,
                m.read_pair_duplicates,
                m.read_pair_optical_duplicates,
                m.percent_duplication,
                m.estimated_library_size.unwrap_or(0),
            ));
        }
        std::fs::write(path, text)?;
    }
    Ok(())
}
