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

use std::io::{Read, Write};

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::writer::BamWriter;
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

/// A BAM begins with a BGZF block, whose first two bytes are gzip's magic.
fn is_bgzf(raw: &[u8]) -> bool {
    raw.starts_with(&[0x1f, 0x8b])
}

/// The name htsjdk gives the index it writes beside an output.
///
/// `SAMFileWriterFactory` replaces a `.bam` suffix with `.bai` and otherwise appends one, so an
/// output called `output.txt` is indexed as `output.txt.bai`.
fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
    }
}

fn read_input(path: &str) -> Result<(SamHeader, Vec<BamRecord>), Box<dyn std::error::Error>> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    if !is_bgzf(&raw) {
        let text = String::from_utf8(raw)?;
        let (header, records) = read_sam(&text).map_err(|e| format!("{e:?}"))?;
        return Ok((header, records));
    }
    let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
    let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
    let header = reader.header.text.clone();
    let records: Vec<BamRecord> = reader
        .map(|r| r.map_err(|e| format!("{e:?}")))
        .collect::<Result<_, _>>()?;
    Ok((header, records))
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
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
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

    let (header, mut records) = read_input(&input)?;
    t = phase("read_input", t);

    // The parallel per-record OA stamping; the whole point of the benchmark.
    add_oa_records(&header, &mut records);
    t = phase("add_oa", t);

    // The extension decides the format, and everything that is not `.sam` is a BAM.
    if output.ends_with(".sam") {
        let sam = write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?;
        t = phase("encode_sam", t);
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        phase("write_file", t);
        return Ok(());
    }

    let mut writer = BamWriter::new(Vec::new(), &header)?;
    for record in &records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    let bam = writer.finish()?;
    t = phase("encode_bam", t);

    std::fs::write(&output, &bam)?;
    if create_index {
        let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
        std::fs::write(index_path(&output), index)?;
    }
    phase("write_file", t);
    Ok(())
}
