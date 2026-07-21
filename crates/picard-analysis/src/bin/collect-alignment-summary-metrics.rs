//! `CollectAlignmentSummaryMetrics` as a runnable binary, for the throughput benchmark.
//!
//! Deliberately minimal: it exists so the port can be timed against Picard on the same input,
//! and so the two outputs can be compared byte for byte in the same run. It is **not** the
//! Barclay-compatible command line the program commits to; that is generated from the inventory
//! and comes later. The three arguments it accepts are the three the benchmark needs.
//!
//! What it does share with the real thing is the work: the same reader, the same collector, the
//! same `MetricsFile` writer. A benchmark that measured a different code path would be worthless.

use std::io::Read;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::reader::BamReader;
use htsjdk_metrics::file::{Histogram as OutHistogram, MetricsFile};
use picard_analysis::alignment_summary::{GroupCollector, Options};

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
    let reference = arg(&args, "REFERENCE_SEQUENCE=").or_else(|| arg(&args, "R="));

    // Picard streams; this reads the whole file. That is a real difference in memory profile and
    // it is stated rather than hidden, because it is part of what any timing measures.
    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;

    let contigs = match &reference {
        Some(path) => read_fasta_file(path).map_err(|e| format!("{e:?}"))?,
        None => Vec::new(),
    };

    let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
    let mut collector = GroupCollector::new(Options::default());
    for record in reader {
        let rec = record.map_err(|e| format!("{e:?}"))?;
        // One contig in the benchmark corpus, as in the conformance corpus. A real reference
        // walker resolves per record; that is the walker's job and not this collector's.
        let bases = contigs.first().map(|c| c.bases.as_slice());
        collector.accept(&rec, bases);
    }
    collector.finish()?;

    let mut file = MetricsFile::new();
    file.add_header("CollectAlignmentSummaryMetrics <command line>");
    file.add_header("Started on: <timestamp>");
    for row in collector.rows() {
        file.add_metric(&row);
    }
    for (label, histogram) in collector.read_length_histograms() {
        file.histograms.push(OutHistogram {
            bin_label: "READ_LENGTH".to_string(),
            value_label: label.to_string(),
            key_class: "java.lang.Integer".to_string(),
            bins: histogram
                .bins()
                .map(|(id, count)| (format!("{}", id as i64), count))
                .collect(),
        });
    }
    std::fs::write(&output, file.write())?;
    Ok(())
}
