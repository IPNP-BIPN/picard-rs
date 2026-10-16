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
use std::io::Read as _;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::collect_umi_prevalence_metrics::{
    decode_barcode_qualities, filters_out, Arguments, Read, DEFAULT_BARCODE_QUALITY_TAG,
    DEFAULT_BARCODE_TAG, DEFAULT_MINIMUM_BARCODE_BASE_QUALITY, DEFAULT_MINIMUM_MAPPING_QUALITY,
};
use picard_analysis::duplicate_set::duplicate_sets;
use picard_analysis::mark_duplicates::{Options, Record as SetRecord};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

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
            .map(|value| value == "true")
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

    // `FilteringSamIterator` around the aggregate filter, which runs before the sets are cut.
    let kept: Vec<&BamRecord> = records
        .iter()
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
        .collect();

    let set_records: Vec<SetRecord> = kept
        .iter()
        .map(|record| {
            let (library, read_group) = library_of(record);
            SetRecord {
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
                existing_dt: None,
                mate_cigar: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
                    _ => None,
                },
                mate_alignment_start: record.mate_alignment_start,
            }
        })
        .collect();

    let options = Options::default();
    let mut bins: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
    for set in duplicate_sets(&set_records, &options) {
        // Every read in a set is already a survivor, so the count is the set's distinct barcodes.
        let barcodes: HashSet<Option<String>> = set
            .iter()
            .map(|index| string_tag(kept[*index], &barcode_tag))
            .collect();
        *bins.entry(barcodes.len()).or_insert(0.0) += 1.0;
    }

    let mut file = MetricsFile::new();
    file.add_header("CollectUmiPrevalenceMetrics <command line>");
    file.add_header("Started on: <timestamp>");
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
