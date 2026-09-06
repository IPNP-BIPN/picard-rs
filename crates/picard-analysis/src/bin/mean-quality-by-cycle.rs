//! `MeanQualityByCycle` as a runnable binary: the covering array's port side.
//!
//! The tool writes no metric rows at all. Its whole body is one histogram, and a second one when
//! the input carries `OQ`, so what the arguments decide is which reads reach the histogram:
//! `PF_READS_ONLY` drops vendor-failed reads, `ALIGNED_READS_ONLY` drops unmapped ones, and the
//! array's fixtures make those two choices observable (its rows produce two distinct outputs).
//!
//! `CHART_OUTPUT` is required by the reference and is accepted here, but nothing is written to it:
//! Picard renders that PDF by shelling out to R, which is a rendering of the metrics file rather
//! than the tool's answer. The array compares the metrics file, which is what both sides compute.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::cycle::MeanQualityByCycle;
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

    let mut collector = MeanQualityByCycle::default();
    collector.pf_reads_only = flag("PF_READS_ONLY=");
    collector.aligned_reads_only = flag("ALIGNED_READS_ONLY=");
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
    file.add_header("MeanQualityByCycle <command line>");
    file.add_header("Started on: <timestamp>");
    file.histograms = collector.finish();
    std::fs::write(&output, file.write())?;
    Ok(())
}
