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

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::mark_duplicates::{Options, Record, ScoringStrategy};
use picard_analysis::mate_cigar_duplicates::{describe_with_contig, needs_mate_cigar};
use picard_analysis::metrics_cli::{fail, read_input, thrown, Args};
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
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("M", "METRICS_FILE"),
        ("ASO", "ASSUME_SORT_ORDER"),
        ("UMI_METRICS", "UMI_METRICS_FILE"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let input = args.required("INPUT");
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
        parse_read_names: true,
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

    let (mut header, records) = read_input(&input);
    // ASSUME_SORT_ORDER is written into the header by `openInputs`, so it decides the check below
    // as surely as the file does.
    if let Some(order) = args.get("ASSUME_SORT_ORDER") {
        header.set_sort_order(order);
    }
    if header.attributes.get("SO") != Some("coordinate") {
        thrown("picard.PicardException: This program requires inputs in coordinate SortOrder");
    }

    // The read group's library, which is what a duplicate set is cut by, and the group's index in
    // the header, which is what `closeEnough` compares.
    let library_of = |record: &BamRecord, header: &SamHeader| -> (String, i32) {
        match string_tag(record, b"RG").and_then(|id| {
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
            let (library, read_group) = library_of(record, &header);
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
                barcode: None,
                existing_dt: string_tag(record, b"DT"),
                mate_cigar: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
                    _ => None,
                },
                mate_alignment_start: record.mate_alignment_start,
            }
        })
        .collect();

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
        if needs_mate_cigar(record) && record.mate_cigar.is_none() {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found: {}",
                describe_with_contig(
                    record,
                    contigs
                        .get(record.reference_index.max(0) as usize)
                        .map(String::as_str),
                )
            ));
        }
    }

    let umi_tag_name: [u8; 2] = {
        let bytes = umi_tag.as_bytes();
        [bytes[0], bytes[1]]
    };
    let umis: Vec<Option<String>> = records
        .iter()
        .map(|record| string_tag(record, &umi_tag_name))
        .collect();
    let result = run(&marked_records, &umis, &contigs, &options).unwrap_or_else(|m| thrown(&m));

    let add_pg_tag = args.bool("ADD_PG_TAG_TO_READS", true);
    let molecular_identifier_tag = args.get("MOLECULAR_IDENTIFIER_TAG").map(tag_of);
    let mut writer =
        BamWriter::new(Vec::new(), &header).unwrap_or_else(|e| fail(&format!("{e:?}")));
    for (index, record) in records.iter().enumerate() {
        if base.remove_duplicates && result.duplicate[index] {
            continue;
        }
        let mut written = record.clone();
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
        // `updateProgramRecord` only CHAINS: a record that already carries a `PG` gets the new
        // id, and a record with none gets nothing.
        if add_pg_tag {
            if let Some(TagValue::Str(existing)) = written.tags.get(Tag::new(b"PG")) {
                let chained = existing.clone();
                if header.programs.iter().any(|pg| pg.id == chained) {
                    written.tags.insert(Tag::new(b"PG"), TagValue::Str(chained));
                }
            }
        }
        writer
            .write(&written)
            .unwrap_or_else(|e| fail(&format!("{e:?}")));
    }
    let bytes = writer.finish().unwrap_or_else(|e| fail(&format!("{e:?}")));
    if let Err(e) = std::fs::write(&output, bytes) {
        fail(&format!("{e}"));
    }

    // The UMI metrics: one row per library, in the order the reference's map iterates them.
    let mut file = MetricsFile::new();
    file.add_header("UmiAwareMarkDuplicatesWithMateCigar <command line>");
    file.add_header("Started on: <timestamp>");
    for row in &result.metrics {
        file.add_metric(row);
    }
    if let Err(e) = std::fs::write(&umi_metrics_path, file.write()) {
        fail(&format!("{e}"));
    }

    // The duplication metrics, in the plain form the other duplicate markers write them in.
    let marking = picard_analysis::mark_duplicates::marking_from(
        &marked_records,
        &base,
        &result.duplicate,
        &vec![false; marked_records.len()],
    );
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
    if let Err(e) = std::fs::write(&metrics_path, text) {
        fail(&format!("{e}"));
    }
}
