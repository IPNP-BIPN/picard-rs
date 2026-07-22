//! `AddOATag` as a runnable binary, for the multicore throughput benchmark.
//!
//! Deliberately minimal, like the `collect-alignment-summary-metrics` bench binary: it exists so the
//! per-record transform can be timed against Picard on the same BAM, and so the two SAM outputs can be
//! compared byte for byte in the same run. It is **not** the Barclay command line the program commits
//! to; the two arguments it accepts are the two the benchmark needs.
//!
//! It shares the work with the real thing: the same `BamReader`, the same `add_oa_records` transform
//! (which fans the per-record OA stamping across all cores via rayon), and the same SAM writer.
//! `AddOATag` adds no `@PG` and no timestamp, so the whole SAM is comparable raw, and the input is
//! coordinate-sorted with the output order preserved, so no sort enters the measurement.
//!
//! Set `RAYON_NUM_THREADS=1` to measure the single-threaded floor, and leave it unset (all cores) to
//! measure the parallel path; the two isolate the multicore win from the Rust-vs-JVM baseline.

use std::io::{Read, Write};

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use picard_analysis::add_oa_tag::add_oa_records;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// Coarse phase timing to stderr when `PICARD_RS_TIMING=1`, matching the metrics bench binary.
fn phase(label: &str, start: std::time::Instant) -> std::time::Instant {
    if std::env::var("PICARD_RS_TIMING").as_deref() == Ok("1") {
        eprintln!("{label}\t{:.3}", start.elapsed().as_secs_f64());
    }
    std::time::Instant::now()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut t = std::time::Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    t = phase("read_file", t);
    let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
    t = phase("bgzf_decompress", t);

    let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
    let header = reader.header.text.clone();
    let mut records: Vec<BamRecord> = reader
        .map(|r| r.map_err(|e| format!("{e:?}")))
        .collect::<Result<_, _>>()?;
    t = phase("decode", t);

    // The parallel per-record OA stamping; the whole point of the benchmark.
    add_oa_records(&header, &mut records);
    t = phase("add_oa", t);

    let sam = write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?;
    t = phase("encode_sam", t);

    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(sam.as_bytes())?;
    out.flush()?;
    phase("write_file", t);
    Ok(())
}
