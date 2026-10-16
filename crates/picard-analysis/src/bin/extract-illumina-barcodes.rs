//! `ExtractIlluminaBarcodes` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.ExtractIlluminaBarcodes`, `ExtractBarcodesProgram`, `BarcodeExtractor`,
//! `DistanceMetric` and `SingleBarcodeDistanceMetric` at tag 3.4.0 over
//! `picard_analysis::illumina_reader` and `picard_analysis::barcode_extractor`.
//!
//! `BARCODE_FILE` replaces `INPUT_PARAMS_FILE` before it is read, and `MINIMUM_QUALITY` is read in
//! a field initialiser before the parser assigns it, so neither changes anything. `COMPRESS_OUTPUTS`
//! gzips each barcode file; the content is compared, not the deflater's bytes.

use htsjdk_metrics::file::MetricsFile;
use picard_analysis::barcode_extractor::{finalize, parse_input_file, Distance, Extractor, Metric};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::illumina_reader::{eamss, Run};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};

const TOOL: &str = "ExtractIlluminaBarcodes";

fn main() {
    picard_analysis::command_line::check("ExtractIlluminaBarcodes");
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("L", "LANE"),
        ("M", "METRICS_FILE"),
        ("RS", "READ_STRUCTURE"),
    ]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let metrics_file = args.required("METRICS_FILE");
    // BARCODE_FILE is required unless BARCODE is given, which Barclay checks first of all.
    if !args.given("BARCODE_FILE") && args.all("BARCODE").is_empty() {
        eprintln!("Argument BARCODE_FILE was missing: Argument 'BARCODE_FILE' is required unless one of {{[BARCODE]}} are provided");
        std::process::exit(1);
    }
    let output_dir = args
        .get("OUTPUT_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| basecalls.clone());
    let lanes: Vec<i32> = args
        .collection("LANE", &[])
        .iter()
        .filter_map(|v| v.parse().ok())
        .collect();
    let structure = parse_read_structure(&args.required("READ_STRUCTURE"))
        .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"));
    let barcode_file = args.get("BARCODE_FILE").map(str::to_string);
    let barcode_args = args.collection("BARCODE", &[]);
    // `COMPRESS_OUTPUTS` gzips each barcode file and names it `.gz`. The bytes are the deflater's;
    // the content is the tool's, and that is what is compared.
    let compress = args.bool("COMPRESS_OUTPUTS", false);
    let lengths: Vec<usize> = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Barcode)
        .map(|s| s.cycles)
        .collect();

    // `customCommandLineValidation`: this class's messages first, then its parent's.
    let mut messages = Vec::new();
    let mut declared: Vec<(String, Metric)> = Vec::new();
    if barcode_file.is_none() {
        let mut seen: Vec<&String> = Vec::new();
        for barcode in &barcode_args {
            if seen.contains(&barcode) {
                messages.push(format!("Barcode {barcode} specified more than once."));
            }
            seen.push(barcode);
            let mut seqs = Vec::new();
            let mut at = 0;
            for &l in &lengths {
                let piece = barcode.get(at..at + l).unwrap_or_else(|| {
                    thrown(&format!(
                        "java.lang.StringIndexOutOfBoundsException: begin {at}, end {}, length {}",
                        at + l,
                        barcode.len()
                    ))
                });
                seqs.push(piece.to_string());
                at += l;
            }
            let metric = Metric::new(None, None, &seqs);
            match declared.iter_mut().find(|(k, _)| k == barcode) {
                Some(slot) => slot.1 = metric,
                None => declared.push((barcode.clone(), metric)),
            }
        }
    }
    let mut parent = Vec::new();
    if let Some(path) = &barcode_file {
        declared = parse_input_file(path, &lengths, &mut parent);
        if declared.is_empty() {
            parent.push("No barcodes have been specified.".to_string());
        }
    }
    if (barcode_file.is_some() || !barcode_args.is_empty()) && declared.is_empty() {
        messages.push("No barcodes have been specified.".to_string());
    }
    messages.extend(parent);
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }

    // `doWork` opens with `assertDirectoryIsWritable(OUTPUT_DIR)`.
    if !output_dir.is_dir() {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Directory does not exist: {}",
            picard_analysis::fingerprint::file_uri(&output_dir.display().to_string())
        ));
    }
    let minimum = args.int("MINIMUM_BASE_QUALITY", 0) as i32;
    let mode = Distance::parse(args.get("DISTANCE_MODE"));
    let no_match_seqs: Vec<String> = lengths.iter().map(|l| "N".repeat(*l)).collect();
    let mut no_match = Metric::new(None, None, &no_match_seqs);
    let metrics_in_order: Vec<Metric> = declared.iter().map(|(_, m)| m.clone()).collect();
    let mut extractor = Extractor::new(
        &metrics_in_order,
        &lengths,
        (
            args.int("MAX_NO_CALLS", 2) as i32,
            args.int("MAX_MISMATCHES", 1) as i32,
            args.int("MIN_MISMATCH_DELTA", 1) as i32,
            minimum,
        ),
        mode,
    );

    // The provider reads every cycle the structure does not skip; the barcode reads are the `B`
    // segments among them.
    let mut reads: Vec<Vec<i32>> = Vec::new();
    let mut barcode_reads: Vec<usize> = Vec::new();
    let mut cycle = 1;
    for s in &structure {
        let cycles: Vec<i32> = (cycle..cycle + s.cycles as i32).collect();
        cycle += s.cycles as i32;
        if s.kind == SegmentKind::Skip {
            continue;
        }
        if s.kind == SegmentKind::Barcode {
            barcode_reads.push(reads.len());
        }
        reads.push(cycles);
    }
    let mut metrics = metrics_in_order;
    // `doWork` keeps one list of tile extractors across the lanes and submits the WHOLE list after
    // each lane, so a lane's tiles run once for every lane from it to the last. Each run counts
    // into the same extractor's metrics, and writes the same barcode file again.
    for (k, &lane) in lanes.iter().enumerate() {
        let times = (lanes.len() - k) as i64;
        let run = Run::new(&basecalls, lane);
        let tiles = run.available_tiles().unwrap_or_else(|e| thrown(&e));
        for tile in tiles {
            let clusters = run.clusters(tile, &reads, 2).unwrap_or_else(|e| thrown(&e));
            let mut text = String::new();
            for c in &clusters {
                let read: Vec<Vec<u8>> = barcode_reads
                    .iter()
                    .map(|i| c.reads[*i].0.clone())
                    .collect();
                // The provider runs EAMSS over every read it hands out, the barcode reads included.
                let quals: Vec<Vec<u8>> = barcode_reads
                    .iter()
                    .map(|i| {
                        let (bases, quals) = &c.reads[*i];
                        let mut quals = quals.clone();
                        eamss(bases, &mut quals);
                        quals
                    })
                    .collect();
                let m = extractor.find(&read, (minimum > 0).then_some(quals.as_slice()), false);
                // `updateMetrics`.
                let target = if m.matched {
                    let at = declared
                        .iter()
                        .position(|(k, _)| *k == m.barcode)
                        .unwrap_or_else(|| thrown("java.lang.NullPointerException"));
                    &mut metrics[at]
                } else {
                    &mut no_match
                };
                target.reads += times;
                if c.pf {
                    target.pf_reads += times;
                }
                if m.matched {
                    if m.mismatches == 0 {
                        target.perfect += times;
                        if c.pf {
                            target.pf_perfect += times;
                        }
                    } else if m.mismatches == 1 {
                        target.one += times;
                        if c.pf {
                            target.pf_one += times;
                        }
                    }
                }
                for r in &read {
                    text.push_str(&String::from_utf8_lossy(r));
                }
                text.push_str(&format!(
                    "\t{}\t{}\t{}\t{}\n",
                    if m.matched { "Y" } else { "N" },
                    m.barcode,
                    m.mismatches,
                    m.to_second
                ));
            }
            let suffix = if compress { ".gz" } else { "" };
            let path = output_dir.join(format!("s_{lane}_{tile:04}_barcode.txt{suffix}"));
            let written = if compress {
                std::fs::File::create(&path).and_then(|file| {
                    use std::io::Write;
                    let mut writer = htsjdk_bgzf::BgzfWriter::new(file);
                    writer.write_all(text.as_bytes())?;
                    writer.finish().map(|_| ())
                })
            } else {
                std::fs::write(&path, text)
            };
            if let Err(e) = written {
                thrown(&format!("htsjdk.samtools.SAMException: {e}"));
            }
        }
    }
    finalize(&mut metrics, &mut no_match);
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for m in &metrics {
        file.add_metric(m);
    }
    file.add_metric(&no_match);
    if let Err(e) = std::fs::write(&metrics_file, file.write()) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
