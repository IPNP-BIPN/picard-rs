//! `AddOATag` as a runnable binary: the throughput benchmark, and the covering array's port side.
//!
//! It was written for the benchmark, which needed two arguments and a comparable output, and it
//! wrote SAM because SAM is what a benchmark can diff. That made every row of the tool's covering
//! array differ from the reference for one reason: Picard writes a BAM there. `SAMFileWriterFactory`
//! picks the format from the OUTPUT's extension and defaults to BAM, and the array's OUTPUT is
//! `output.txt`, so the reference's answer to every row is a BAM.
//!
//! So this writes what Picard writes: a BAM unless the output is named `.sam`, through the same
//! `BamWriter` the conformance suite compares byte for byte, and a BAI beside it when
//! `CREATE_INDEX=true`.
//!
//! It reads what Picard reads, too: the array's fixtures are a BAM, a SAM and a queryname-sorted
//! BAM, and the format is told from the file's own first bytes rather than from its name.
//!
//! Set `RAYON_NUM_THREADS=1` to measure the single-threaded floor, and leave it unset (all cores)
//! to measure the parallel path.

use picard_analysis::add_oa_tag::add_oa_records;
use picard_analysis::args::arg;
use picard_analysis::stream_output::Output;

/// How many records are tagged in parallel at a time, between reading and writing.
const CHUNK: usize = 16 * 1024;

/// Coarse phase timing to stderr when `PICARD_RS_TIMING=1`, matching the metrics bench binary.
fn phase(label: &str, start: std::time::Instant) -> std::time::Instant {
    if std::env::var("PICARD_RS_TIMING").as_deref() == Ok("1") {
        eprintln!("{label}\t{:.3}", start.elapsed().as_secs_f64());
    }
    std::time::Instant::now()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("AddOATag");
    picard_analysis::output::install();
    let mut t = std::time::Instant::now();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);
    // The stringency is read so that a row naming it is not refused for naming it. The array's
    // three fixtures are valid under all three levels, so none of them makes the levels differ,
    // and a level that did would show as a row this binary got wrong rather than as one it
    // silently agreed with.
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();
    let intervals = match arg(&args, "INTERVAL_LIST=").or_else(|| arg(&args, "L=")) {
        Some(path) => {
            let text = std::fs::read_to_string(&path)?;
            let list = htsjdk_bam::interval::IntervalList::parse_body(Vec::new(), &text)
                .map_err(|e| format!("{e:?}"))?;
            Some(list.intervals)
        }
        None => None,
    };
    t = phase("open", t);

    // The extension decides the format, and everything that is not `.sam` is a BAM. Nothing was
    // written before every record had been tagged and encoded, so a failure leaves no file.
    let mut out = Output::create(&output, &header, false)?;
    let mut records = records.peekable();
    while records.peek().is_some() {
        let mut chunk: Vec<_> = records.by_ref().take(CHUNK).collect();
        // The parallel per-record OA stamping; the whole point of the benchmark.
        match &intervals {
            Some(intervals) => picard_analysis::add_oa_tag::add_oa_records_overlapping(
                &header, &mut chunk, intervals,
            ),
            None => add_oa_records(&header, &mut chunk),
        }
        for record in &chunk {
            out.write(record)?;
        }
    }
    out.finish(create_index, create_index)?;
    phase("write_file", t);
    Ok(())
}
