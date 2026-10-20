//! `UmiAwareMarkDuplicatesWithMateCigar` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.markduplicates.UmiAwareMarkDuplicatesWithMateCigar.doWork` at tag 3.4.0 over
//! `SimpleMarkDuplicatesWithMateCigar`: htsjdk's duplicate sets, each broken up again by UMI, the
//! records written with their duplicate flags and (where `MOLECULAR_IDENTIFIER_TAG` names one) the
//! molecular identifier of the UMI they were assigned, and a second metrics file about the UMIs.
//! The sets and the metrics are `picard_analysis::umi_aware`.
//!
//! The duplication metrics file the parent writes is accepted and written in the plain form the
//! other duplicate markers' binaries write it in: its content is not what this array compares.

use std::sync::Arc;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::bam_stream::{Input, Records};
use picard_analysis::mark_duplicates::{
    base_quality_score, Options, Record, ScoringStrategy, Span,
};
use picard_analysis::mate_cigar_duplicates::needs_mate_cigar;
use picard_analysis::metrics_cli::{fail, read_input, thrown, Args};
use picard_analysis::stream_output::Output;
use picard_analysis::umi_aware::{run, UmiAwareOptions};

const DUPLICATE_READ: u16 = 0x400;

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    }
}

fn tag_of(name: &str) -> Tag {
    let bytes = name.as_bytes();
    Tag::new(&[bytes[0], bytes[1]])
}

fn main() {
    picard_analysis::argument_check::check("UmiAwareMarkDuplicatesWithMateCigar");
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
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("M", "METRICS_FILE"),
        ("ASO", "ASSUME_SORT_ORDER"),
        ("UMI_METRICS", "UMI_METRICS_FILE"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let inputs = args.all("INPUT");
    if inputs.is_empty() {
        fail("Argument 'INPUT' is required");
    }
    let output = args.required("OUTPUT");
    let metrics_path = args.required("METRICS_FILE");
    let umi_metrics_path = args.required("UMI_METRICS_FILE");

    let scoring = match args
        .get("DUPLICATE_SCORING_STRATEGY")
        .unwrap_or("SUM_OF_BASE_QUALITIES")
    {
        "SUM_OF_BASE_QUALITIES" => ScoringStrategy::SumOfBaseQualities,
        "TOTAL_MAPPED_REFERENCE_LENGTH" => ScoringStrategy::TotalMappedReferenceLength,
        "RANDOM" => ScoringStrategy::Random,
        other => fail(&format!(
            "Argument 'DUPLICATE_SCORING_STRATEGY' cannot be set to '{other}': invalid value"
        )),
    };
    let base = Options {
        scoring,
        remove_duplicates: args.bool("REMOVE_DUPLICATES", false),
        parse_read_names: !read_name_regex.is_off(),
        read_name_regex,
        max_optical_duplicate_set_size: std::env::args()
            .find_map(|a| {
                a.strip_prefix("MAX_OPTICAL_DUPLICATE_SET_SIZE=")
                    .and_then(|v| v.parse().ok())
            })
            .unwrap_or(300_000),
        assume_mate_cigar: true,
        ..Options::default()
    };
    let umi_tag = args.get("UMI_TAG_NAME").unwrap_or("RX").to_string();
    let options = UmiAwareOptions {
        base: base.clone(),
        max_edit_distance_to_join: args.int("MAX_EDIT_DISTANCE_TO_JOIN", 1) as i32,
        molecular_identifier_tag: args.get("MOLECULAR_IDENTIFIER_TAG").is_some(),
        allow_missing_umis: args.bool("ALLOW_MISSING_UMIS", false),
        duplex_umi: args.bool("DUPLEX_UMI", false),
    };

    // ASSUME_SORT_ORDER is written into the header by `openInputs`, so it decides the check below
    // as surely as the file does.
    // `ASSUME_SORTED=true` stands for `ASSUME_SORT_ORDER=coordinate`.
    let assumed = args.get("ASSUME_SORT_ORDER").or_else(|| {
        args.get("ASSUME_SORTED")
            .is_some_and(picard_analysis::args::is_true)
            .then_some("coordinate")
    });
    // Read twice, as MarkDuplicates reads it: once for the marking, which keeps only what it
    // needs of each record, and once for the writing. One input is streamed, and its first read
    // checks the whole file, so a damaged one is refused here as reading it whole refused it.
    // Several are merged in memory by `openInputs`, once per pass.
    let open = |streamed: Option<bool>| -> (SamHeader, Records, bool) {
        if let [input] = inputs.as_slice() {
            let opened = match streamed {
                None => Input::open(input),
                Some(streamed) => Input::reopen(input, streamed),
            }
            .unwrap_or_else(|e| fail(&e.message()));
            let streams = opened.streams();
            let (mut header, records) = opened.into_parts();
            if let Some(order) = assumed {
                header.set_sort_order(order);
            }
            return (header, records, streams);
        }
        let (header, records) = picard_analysis::merge_sam_files::open_inputs(
            inputs.iter().map(|input| read_input(input)).collect(),
            assumed,
        )
        .unwrap_or_else(|message| thrown(&message));
        // Only an order this tool refuses before reading a record makes reading one throw.
        let records = records.unwrap_or_default();
        (header, Records::Decoded(records.into_iter(), None), false)
    };
    let (mut header, records, streamed) = open(None);
    if header.attributes.get("SO") != Some("coordinate") {
        thrown("picard.PicardException: This program requires inputs in coordinate SortOrder");
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
    let umi_tag_name: [u8; 2] = {
        let bytes = umi_tag.as_bytes();
        [bytes[0], bytes[1]]
    };
    // Each record's UMI is kept beside what the marking keeps of it.
    let (marked_records, umis): (Vec<Record>, Vec<Option<String>>) = records
        .map(|record| {
            let (library, read_group) = library_of(&record);
            let marked = Record {
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
                barcode: None,
                existing_dt: string_tag(&record, b"DT"),
                mate_span: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text)
                        .ok()
                        .as_ref()
                        .map(Span::of),
                    _ => None,
                },
                mate_alignment_start: record.mate_alignment_start,
            };
            (marked, string_tag(&record, &umi_tag_name))
        })
        .unzip();

    let contigs: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();

    // The refusal is htsjdk's, from `SAMUtils.getMateUnclippedStart`/`End`, and WHICH record it
    // names is decided by the sort `DuplicateSetIterator` runs before it cuts a set: Java's sort
    // begins by comparing `a[1]` against `a[0]`, so the record named is the second of the file.
    let examination = [1usize, 0]
        .into_iter()
        .chain(2..marked_records.len())
        .filter(|index| *index < marked_records.len());
    for index in examination {
        let record = &marked_records[index];
        if needs_mate_cigar(record) && record.mate_span.is_none() {
            // The record it names is read again, by its index.
            let named = open(Some(streamed))
                .1
                .nth(index)
                .expect("the record read the first time");
            thrown(&format!(
                "htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found: {}",
                picard_analysis::view_sam::sam_string(&header, &named)
            ));
        }
    }

    // `COMMENT.forEach(outputHeader::addComment)`, inherited from SimpleMarkDuplicatesWithMateCigar.
    for comment in std::env::args().filter_map(|a| a.strip_prefix("COMMENT=").map(str::to_string)) {
        header.comments.push(format!("@CO\t{comment}"));
    }
    let result = run(&marked_records, &umis, &contigs, &options).unwrap_or_else(|m| {
        // The output was opened before the sets were iterated and is never closed after the throw:
        // a SAM writer has flushed its header by then, a BAM writer's first block is still buffered.
        let partial = if output.ends_with(".sam") {
            write_sam(&header, &[]).unwrap_or_default()
        } else {
            String::new()
        };
        let _ = std::fs::write(&output, partial);
        thrown(&m)
    });

    let program_group =
        picard_analysis::command_line::ProgramGroup::of("UmiAwareMarkDuplicatesWithMateCigar");
    let molecular_identifier_tag = args.get("MOLECULAR_IDENTIFIER_TAG").map(tag_of);
    let mut out = Output::bam(&output, &header, false).unwrap_or_else(|e| fail(&format!("{e:?}")));
    for (index, mut written) in open(Some(streamed)).1.enumerate() {
        if base.remove_duplicates && result.duplicate[index] {
            continue;
        }
        written.flags &= !DUPLICATE_READ;
        if result.duplicate[index] {
            written.flags |= DUPLICATE_READ;
        }
        if result.umi_removed[index] {
            written.tags.remove(tag_of(&umi_tag));
        }
        if let (Some(tag), Some(identifier)) = (
            molecular_identifier_tag,
            &result.molecular_identifier[index],
        ) {
            written.tags.insert(tag, TagValue::Str(identifier.clone()));
        }
        // `record.setAttribute(PG, chainedPgIds.get(...))`, inherited from
        // SimpleMarkDuplicatesWithMateCigar: whenever PROGRAM_RECORD_ID is set, `ADD_PG_TAG_TO_READS`
        // unread. Nothing fills `pgIdsSeen` in this one-pass tool, so the chain is empty, no `@PG` is
        // added, and every record's `PG` is removed.
        if program_group.record_id.is_some() {
            written.tags.remove(Tag::new(b"PG"));
        }
        out.write(&written).unwrap_or_else(|e| fail(&e.to_string()));
    }
    // A CRAM is indexed under CREATE_INDEX, a BAM never was here.
    if let Err(e) = out.finish(picard_analysis::cram_output::create_index_argument(), false) {
        fail(&e.to_string());
    }

    // The UMI metrics: one row per library, in the order the reference's map iterates them.
    let mut file = MetricsFile::new();
    for header in
        picard_analysis::command_line::default_headers("UmiAwareMarkDuplicatesWithMateCigar")
    {
        file.add_header(&header);
    }
    for row in &result.metrics {
        file.add_metric(row);
    }
    if let Err(e) = std::fs::write(&umi_metrics_path, file.write()) {
        fail(&format!("{e}"));
    }

    // The duplication metrics, through `finalizeAndWriteMetrics` as every duplicate marker writes them.
    let tracked = picard_analysis::duplicate_set::track_sets(&result.sets, &marked_records, &base);
    let mut marking = picard_analysis::mark_duplicates::marking_from(
        &marked_records,
        &base,
        &result.duplicate,
        &tracked.optical,
    );
    marking.all_sets = tracked.all_sets;
    marking.optical_sets = tracked.optical_sets;
    marking.non_optical_sets = tracked.non_optical_sets;
    let text = picard_analysis::mark_duplicates::duplication_metrics_text(
        "UmiAwareMarkDuplicatesWithMateCigar",
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
    if let Err(e) = std::fs::write(&metrics_path, text) {
        fail(&format!("{e}"));
    }
}
