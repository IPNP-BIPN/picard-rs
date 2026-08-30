//! `CollectQualityYieldMetrics` as a runnable binary: the covering array's port side.
//!
//! The array has existed since the harness landed and has only ever run the reference: the fuzzer
//! seeds from it, and the port half was empty because there was no binary to run. This is that
//! half.
//!
//! The tool counts bases and quality, so the arguments that matter are the ones deciding WHICH
//! bases are counted: whether the original qualities are used where a record carries them, and
//! whether secondary and supplementary records are counted at all.
//!
//! It reads a BAM or a SAM, telling the two apart by the file's own first bytes, because the
//! array's fixtures are one of each.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::quality_yield::{Options, QualityYieldMetricsCollector};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    let mut options = Options::default();
    if let Some(value) = arg(&args, "USE_ORIGINAL_QUALITIES=") {
        options.use_original_qualities = value == "true";
    }
    if let Some(value) = arg(&args, "INCLUDE_SECONDARY_ALIGNMENTS=") {
        options.include_secondary_alignments = value == "true";
    }
    if let Some(value) = arg(&args, "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=") {
        options.include_supplemental_alignments = value == "true";
    }
    // The flow-mode collector is a different tool in the reference, and asking this one for it is
    // refused rather than answered. The array holds the argument at `false` for that reason, and
    // a row that ever asks for `true` should be refused here too.
    if let Some(value) = arg(&args, "FLOW_MODE=") {
        if value == "true" {
            return Err("FLOW_MODE is obsolete. Flow support now provided by \
                 CollectQualityYieldMetricsFlow"
                .into());
        }
    }
    // Read so that a row naming them is not refused for naming them; neither changes what this
    // tool counts on the array's fixtures.
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (_header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
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

    let mut collector = QualityYieldMetricsCollector::new(options);
    for record in &records {
        collector.accept(record);
    }
    collector.finish();

    let mut file = MetricsFile::new();
    file.add_header("CollectQualityYieldMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    file.add_metric(collector.metrics());
    std::fs::write(&output, file.write())?;
    Ok(())
}
