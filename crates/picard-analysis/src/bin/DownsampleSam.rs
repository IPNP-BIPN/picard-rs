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
//! `doWork` adds its own `@PG` (`getPGRecord`) with the seed it used as `rs`, after moving the seed
//! off any an earlier run recorded in the header (`downsample_sam::reseed`).

use htsjdk_metrics::file::MetricsFile;
use picard_analysis::args::arg;
use picard_analysis::downsample_sam::{keep_flags, Strategy, DEFAULT_SEED};
use picard_analysis::quality_yield::{Options, QualityYieldMetricsCollector};
use picard_analysis::stream_output::Output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("DownsampleSam");
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

    // `SamReaderFactory` at LENIENT for a SAM; a BAM's records are validated as they are read.
    let opened = picard_analysis::bam_stream::Input::open_with(
        &input,
        htsjdk_bam::text_parse::ValidationStringency::Lenient,
    )?;
    let streamed = opened.streams();
    let (header, records) = opened.into_parts();
    // The decisions need every name first, so the records are read twice: once for the names, then
    // again to write the ones kept.
    let names: Vec<String> = records.map(|r| r.read_name).collect();

    // The seed is moved off any an earlier run recorded, and this run's `@PG` records the one used.
    let seed = picard_analysis::downsample_sam::reseed(&header, strategy, seed);
    let mut header = header;
    let mut pg = picard_analysis::command_line::pg_record(&header, "DownsampleSam");
    pg.attributes.set(
        picard_analysis::downsample_sam::RANDOM_SEED_TAG,
        &seed.to_string(),
    );
    header.programs.push(pg);

    let keep = match keep_flags(&names, strategy, probability, accuracy, seed) {
        Ok(keep) => keep,
        Err(thrown) => {
            eprintln!("Exception in thread \"main\" {thrown}");
            std::process::exit(1);
        }
    };
    let records = picard_analysis::bam_stream::Input::reopen(&input, streamed)?
        .into_parts()
        .1;
    // `SAMFileWriterFactory.makeWriter` again: `.sam` is SAM, and anything else a BAM (or a CRAM).
    let mut out = Output::create(&output, &header, false)?;
    // `QualityYieldMetricsCollector(true, false, false)` over what was written.
    let mut collector = metrics_file
        .as_ref()
        .map(|_| QualityYieldMetricsCollector::new(Options::default()));
    for (record, keep) in records.zip(keep) {
        if keep {
            out.write(&record)?;
            if let Some(collector) = &mut collector {
                collector.accept(&record);
            }
        }
    }
    // `cram_output::write`: a CRAM takes `CREATE_INDEX` from the command line, a BAM no index.
    out.finish(picard_analysis::cram_output::create_index_argument(), false)?;

    if let (Some(path), Some(mut collector)) = (metrics_file, collector) {
        collector.finish();
        let mut file = MetricsFile::new();
        for header in picard_analysis::command_line::default_headers("DownsampleSam") {
            file.add_header(&header);
        }
        file.add_metric(collector.metrics());
        std::fs::write(path, file.write())?;
    }
    Ok(())
}
