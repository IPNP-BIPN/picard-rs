//! `DownsampleSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.DownsampleSam.doWork` at tag 3.4.0 for the default `ConstantMemory` strategy,
//! which is what `picard_analysis::downsample_sam` implements.
//!
//! The two arguments that decide the output are `PROBABILITY` and `RANDOM_SEED`, and the second is
//! not decoration: the keep decision is a `Murmur3` hash of the READ NAME seeded with it, so two
//! seeds give two different subsets of the same size distribution, and both mates of a template
//! share a decision because they share a name.
//!
//! `STRATEGY` picks how the keep decisions are made (see `picard_analysis::downsample_sam::keep_flags`):
//! `HighAccuracy` and `Chained` read in chunks of `ceil(1 / ACCURACY)` templates and shuffle each
//! chunk with `java.util.Random(RANDOM_SEED)`. `METRICS_FILE` is the quality yield of what was kept.
//!
//! The `@PG` record `doWork` adds carries the command line, temporary paths and all, and is
//! canonicalized away by the comparison exactly as the metrics tools' header is. This binary does
//! not write one, and the claim is over the surviving records and the rest of the header.

use std::io::{Read, Write};

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::downsample_sam::{keep_flags, Strategy, DEFAULT_SEED};
use picard_analysis::quality_yield::{Options, QualityYieldMetricsCollector};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let strategy = match arg(&args, "STRATEGY=").as_deref() {
        None | Some("ConstantMemory") => Strategy::ConstantMemory,
        Some("HighAccuracy") => Strategy::HighAccuracy,
        Some("Chained") => Strategy::Chained,
        Some(other) => return Err(format!("unknown STRATEGY: {other}").into()),
    };
    let probability: f64 = arg(&args, "PROBABILITY=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(1.0);
    let accuracy: f64 = arg(&args, "ACCURACY=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(0.0001);
    let seed: i32 = arg(&args, "RANDOM_SEED=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(DEFAULT_SEED);
    let metrics_file = arg(&args, "METRICS_FILE=").or_else(|| arg(&args, "M="));

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        (header, records)
    } else {
        htsjdk_bam::sam_file::read_sam_with(
            &String::from_utf8(raw)?,
            htsjdk_bam::text_parse::ValidationStringency::Lenient,
        )
        .map_err(|e| format!("{e:?}"))?
    };

    let names: Vec<String> = records.iter().map(|r| r.read_name.clone()).collect();
    let keep = match keep_flags(&names, strategy, probability, accuracy, seed) {
        Ok(keep) => keep,
        Err(thrown) => {
            eprintln!("Exception in thread \"main\" {thrown}");
            std::process::exit(1);
        }
    };
    let kept: Vec<BamRecord> = records
        .into_iter()
        .zip(keep)
        .filter_map(|(record, keep)| keep.then_some(record))
        .collect();

    if output.ends_with(".sam") {
        let sam = write_sam(&header, &kept).ok_or("records failed to re-encode as SAM")?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
    } else {
        // `SAMFileWriterFactory.makeWriter` again: anything that is not `.sam` is a BAM.
        let mut writer = picard_analysis::output::bam_writer(Vec::new(), &header)
            .map_err(|e| format!("{e:?}"))?;
        for record in &kept {
            writer.write(record).map_err(|e| format!("{e:?}"))?;
        }
        std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;
    }

    // `QualityYieldMetricsCollector(true, false, false)` over what was written.
    if let Some(path) = metrics_file {
        let mut collector = QualityYieldMetricsCollector::new(Options::default());
        for record in &kept {
            collector.accept(record);
        }
        collector.finish();
        let mut file = MetricsFile::new();
        file.add_header("DownsampleSam <command line>");
        file.add_header("Started on: <timestamp>");
        file.add_metric(collector.metrics());
        std::fs::write(path, file.write())?;
    }
    Ok(())
}
