//! `QualityScoreDistribution` as a runnable binary: the covering array's port side.
//!
//! Histograms only, like `MeanQualityByCycle`, and one argument more: `INCLUDE_NO_CALLS` decides
//! whether the quality at an `N` base is counted. It defaults to false, so a good quality on a
//! no-call is discarded, and the array varies it -- which is why this array has four distinct
//! outputs where the other two have two.
//!
//! `CHART_OUTPUT` is required by the reference and is accepted here, but nothing is written to it:
//! Picard renders that PDF by shelling out to R, which is a rendering of the metrics file rather
//! than the tool's answer. The array compares the metrics file, which is what both sides compute.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricsFile};
use picard_analysis::quality_score_distribution::{Options, QualityScoreDistribution};
use picard_analysis::single_pass_rejections::{
    check_sort_order, walk_reference, Rejection, SortOrder,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// The driver's own work, before and around the records: the sort-order refusal, and the
/// reference walker that refuses to rewind.
///
/// `SinglePassSamProgram.makeItSo` builds a `ReferenceSequenceFileWalker` whenever
/// `REFERENCE_SEQUENCE` is given -- whether or not the program reads a base from it -- and asks it
/// for every mapped record's contig. On a queryname-sorted input that walker is asked to go
/// backwards, which is why `ASSUME_SORTED=true` does not turn a refusal into a run: it moves the
/// refusal from Picard to htsjdk.
fn drive(
    input: &str,
    header: &SamHeader,
    records: &[BamRecord],
    assume_sorted: bool,
    with_reference: bool,
    mut accept: impl FnMut(&BamRecord),
) -> Result<(), Rejection> {
    let found = match header.attributes.get("SO") {
        Some("coordinate") => SortOrder::Coordinate,
        Some("queryname") => SortOrder::Queryname,
        Some("duplicate") => SortOrder::Duplicate,
        Some("unsorted") => SortOrder::Unsorted,
        _ => SortOrder::Unknown,
    };
    check_sort_order(input, found, assume_sorted)?;

    let mut current: Option<i32> = None;
    for record in records {
        if with_reference && record.reference_index != -1 {
            current = Some(walk_reference(current, record.reference_index)?);
        }
        accept(record);
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    // Required by the reference, so a row that omits it is refused before the input is opened.
    // Barclay refuses a missing required argument through the parser rather than by throwing, and
    // that is what the corpus records: `Argument 'CHART_OUTPUT' is required`.
    if arg(&args, "CHART_OUTPUT=").is_none() && arg(&args, "CHART=").is_none() {
        eprintln!("Argument CHART_OUTPUT was missing: Argument 'CHART_OUTPUT' is required");
        std::process::exit(1);
    }
    let flag = |key: &str| arg(&args, key).map(|v| v == "true").unwrap_or(false);

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

    let mut collector = QualityScoreDistribution::new(Options {
        pf_reads_only: flag("PF_READS_ONLY="),
        aligned_reads_only: flag("ALIGNED_READS_ONLY="),
        include_no_calls: flag("INCLUDE_NO_CALLS="),
    });
    let with_reference = arg(&args, "REFERENCE_SEQUENCE=").is_some() || arg(&args, "R=").is_some();
    // The reference throws rather than exiting, so its handler prints the class before the
    // message: a row that refuses is a row this has to refuse the same way, text included.
    if let Err(rejection) = drive(
        &input,
        &header,
        &records,
        flag("ASSUME_SORTED="),
        with_reference,
        |record| collector.accept(record),
    ) {
        eprintln!("Exception in thread \"main\" {}", rejection.thrown());
        std::process::exit(1);
    }

    let mut file = MetricsFile::new();
    file.add_header("QualityScoreDistribution <command line>");
    file.add_header("Started on: <timestamp>");
    let (q_bins, oq_bins) = collector.finish();
    let histogram = |label: &str, bins: &[(u8, f64)]| OutHistogram {
        bin_label: "QUALITY".to_string(),
        value_label: label.to_string(),
        key_class: "java.lang.Byte".to_string(),
        bins: bins.iter().map(|(k, v)| (k.to_string(), *v)).collect(),
    };
    file.histograms.push(histogram("COUNT_OF_Q", &q_bins));
    // `if (!oqHisto.isEmpty())`: the second histogram appears only when the input carried `OQ`,
    // and its absence is a different output shape rather than an empty column.
    if !oq_bins.is_empty() {
        file.histograms.push(histogram("COUNT_OF_OQ", &oq_bins));
    }
    std::fs::write(&output, file.write())?;
    Ok(())
}
