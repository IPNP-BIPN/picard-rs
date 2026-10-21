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

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use picard_analysis::args::arg;
use picard_analysis::quality_score_distribution::{Options, QualityScoreDistribution};
use picard_analysis::single_pass_rejections::{
    check_sort_order, walk_reference, Rejection, SortOrder,
};

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
    records: impl IntoIterator<Item = BamRecord>,
    assume_sorted: bool,
    with_reference: bool,
    stop_after: i64,
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
    let mut count = 0i64;
    for record in records {
        let record = &record;
        if with_reference && record.reference_index != -1 {
            current = Some(walk_reference(current, record.reference_index)?);
        }
        accept(record);
        // Counted after the record is handed over, as `makeItSo` counts it.
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("QualityScoreDistribution");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    // Optional since Picard 3.5.0, which refuses it only where there is no R.
    picard_analysis::metrics_cli::refuse_chart_without_r(
        "QualityScoreDistribution",
        arg(&args, "CHART_OUTPUT=").is_some() || arg(&args, "CHART=").is_some(),
    );
    let flag = |key: &str| {
        arg(&args, key)
            .map(|v| picard_analysis::args::is_true(&v))
            .unwrap_or(false)
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();

    let mut collector = QualityScoreDistribution::new(Options {
        pf_reads_only: flag("PF_READS_ONLY="),
        aligned_reads_only: flag("ALIGNED_READS_ONLY="),
        include_no_calls: flag("INCLUDE_NO_CALLS="),
    });
    let stop_after: i64 = arg(&args, "STOP_AFTER=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(0);
    let with_reference = arg(&args, "REFERENCE_SEQUENCE=").is_some() || arg(&args, "R=").is_some();
    // The reference throws rather than exiting, so its handler prints the class before the
    // message: a row that refuses is a row this has to refuse the same way, text included.
    if let Err(rejection) = drive(
        &input,
        &header,
        records,
        flag("ASSUME_SORTED="),
        with_reference,
        stop_after,
        |record| collector.accept(record),
    ) {
        eprintln!("Exception in thread \"main\" {}", rejection.thrown());
        std::process::exit(1);
    }

    let file = collector.metrics_file(&picard_analysis::command_line::default_headers(
        "QualityScoreDistribution",
    ));
    std::fs::write(&output, file.write())?;
    Ok(())
}
