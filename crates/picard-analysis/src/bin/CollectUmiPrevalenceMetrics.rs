//! `CollectUmiPrevalenceMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectUmiPrevalenceMetrics.doWork` at tag 3.4.0. The six filters and
//! the counting live in `picard_analysis::collect_umi_prevalence_metrics`; the duplicate sets come
//! from `picard_analysis::duplicate_set`, which is the same iterator `MarkDuplicates` uses.
//!
//! The order is the whole tool: the filters run BEFORE the sets are cut, so a filtered read cannot
//! put two others in one set, and a set every one of whose reads was filtered is not a set at all.
//!
//! The barcode quality filter is the reverse of its name. It drops a read whose barcode has NO
//! base under the floor, so a file of well-formed barcodes reports nothing, and LOWERING the floor
//! drops more reads rather than fewer.

use std::collections::HashSet;

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::args::arg;
use picard_analysis::collect_umi_prevalence_metrics::{
    decode_barcode_qualities, filters_out, Arguments, Read, DEFAULT_BARCODE_QUALITY_TAG,
    DEFAULT_BARCODE_TAG, DEFAULT_MINIMUM_BARCODE_BASE_QUALITY, DEFAULT_MINIMUM_MAPPING_QUALITY,
};
use picard_analysis::duplicate_set::duplicate_sets;
use picard_analysis::mark_duplicates::{base_quality_score, Options, Record as SetRecord, Span};

fn string_tag(record: &BamRecord, name: &str) -> Option<String> {
    let bytes = name.as_bytes();
    if bytes.len() != 2 {
        return None;
    }
    match record.tags.get(Tag::new(&[bytes[0], bytes[1]])) {
        Some(TagValue::Str(value)) => Some(value.clone()),
        _ => None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("CollectUmiPrevalenceMetrics");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let number = |key: &str, default: i32| -> i32 {
        arg(&args, key)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let arguments = Arguments {
        minimum_mapping_quality: number("MINIMUM_MQ=", DEFAULT_MINIMUM_MAPPING_QUALITY),
        minimum_barcode_base_quality: number(
            "MINIMUM_BARCODE_BQ=",
            DEFAULT_MINIMUM_BARCODE_BASE_QUALITY,
        ),
        filter_unpaired_reads: arg(&args, "FILTER_UNPAIRED_READS=")
            .map(|value| picard_analysis::args::is_true(&value))
            .unwrap_or(true),
    };
    let barcode_tag = arg(&args, "BARCODE_TAG=").unwrap_or_else(|| DEFAULT_BARCODE_TAG.to_string());
    let barcode_quality_tag =
        arg(&args, "BARCODE_BQ=").unwrap_or_else(|| DEFAULT_BARCODE_QUALITY_TAG.to_string());

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();

    let library_of = |record: &BamRecord| -> (String, i32) {
        let id = string_tag(record, "RG");
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

    // `FilteringSamIterator` around the aggregate filter, which runs before the sets are cut. Each
    // kept record leaves what the sets are cut on and its barcode, and nothing else.
    let (set_records, barcodes): (Vec<SetRecord>, Vec<Option<String>>) = records
        .filter(|record| {
            let read = Read {
                unmapped: record.flags & 0x4 != 0,
                mapping_quality: i32::from(record.mapping_quality),
                secondary_or_supplementary: record.flags & 0x100 != 0 || record.flags & 0x800 != 0,
                paired: record.flags & 0x1 != 0,
                barcode: string_tag(record, &barcode_tag),
                barcode_qualities: string_tag(record, &barcode_quality_tag)
                    .map(|tag| decode_barcode_qualities(&tag)),
            };
            !filters_out(&read, &arguments)
        })
        .map(|record| {
            let (library, read_group) = library_of(&record);
            let set_record = SetRecord {
                name: record.read_name.clone(),
                flags: record.flags,
                reference_index: record.reference_index,
                alignment_start: record.alignment_start,
                span: Span::of(&record.cigar),
                read_length: record.read_bases.len() as u32,
                base_quality_score: base_quality_score(&record.base_qualities),
                mate_reference_index: record.mate_reference_index,
                library: library.into(),
                read_group,
                barcode: None,
                existing_dt: None,
                mate_span: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text)
                        .ok()
                        .as_ref()
                        .map(Span::of),
                    _ => None,
                },
                mate_alignment_start: record.mate_alignment_start,
            };
            (set_record, string_tag(&record, &barcode_tag))
        })
        .unzip();

    let options = Options::default();
    let mut bins: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
    for set in duplicate_sets(&set_records, &options) {
        // Every read in a set is already a survivor, so the count is the set's distinct barcodes.
        let barcodes: HashSet<Option<String>> =
            set.iter().map(|index| barcodes[*index].clone()).collect();
        *bins.entry(barcodes.len()).or_insert(0.0) += 1.0;
    }

    let mut file = MetricsFile::new();
    for header in picard_analysis::command_line::default_headers("CollectUmiPrevalenceMetrics") {
        file.add_header(&header);
    }
    file.histograms.push(Histogram {
        bin_label: "numUmis".to_string(),
        value_label: "duplicateSets".to_string(),
        key_class: "java.lang.Integer".to_string(),
        bins: bins
            .into_iter()
            .map(|(key, value)| (key.to_string(), value))
            .collect(),
    });
    std::fs::write(&output, file.write())?;
    Ok(())
}
