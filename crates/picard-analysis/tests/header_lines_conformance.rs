//! Conformance for the lines a tool writes about its own run, against Picard 3.5.0.
//!
//! Every other suite strips them (`strip_line_prefixes`, `strip_pg`): they hold paths and a clock
//! reading. This one compares them.
//!
//! # What this suite is for
//!
//!  * **the command line being the legacy parser's**, for 54 tools and four command lines each:
//!    set options in declaration order and by long name, three spaces, then every unset option
//!    with a non-null default; values as their type prints them (`+0042` is `42`, `0.250` is
//!    `0.25`, `TRUE` is `true`, `dir//x/` is `dir/x`); a collection as `[a, b]`, appended to its
//!    default elements, and emptied by `null`;
//!  * **the `@PG` records of a real run** of the seven tools that stamp one, on a header that
//!    already holds colliding ids and a chain: the ids `PgIdGenerator` hands out, one record per
//!    program id `MarkDuplicates` saw on the reads in `HashSet` order, one per unchained header id
//!    for `MarkDuplicatesWithMateCigar`, none for the two that never fill the set, the attributes
//!    in the order `doWork` sets them, DownsampleSam's `rs` after it reseeds, and IntervalListTools
//!    taking the first free integer;
//!  * **every read's `PG` tag** after the chain is applied;
//!  * **and the metrics file's two default headers**, the command line byte for byte and
//!    `Started on:` in `Date.toString()`'s format (its value is the clock).
//!
//! The runs read and wrote under `/work/hl`; the port runs in a temporary directory, which is
//! mapped back to `/work/hl` before anything is compared.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;
use std::process::Command;

use picard_analysis::command_line::{legacy_command_line, tool_options};

fn corpus() -> String {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/header_lines.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

/// The harness's escaping undone: `\\`, `\t` and `\n`.
fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// Every row of one kind, by case, in the dump's order.
fn rows(text: &str, kind: &str) -> Vec<(String, String)> {
    let prefix = format!("{kind}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| {
            let (case, payload) = rest.split_once('\t').expect("case and payload");
            (case.to_string(), unescape(payload))
        })
        .collect()
}

/// EstimateLibraryComplexity's unset `MAX_RECORDS_IN_RAM` is the reference's heap over a record
/// size: the dump ran with `-Xmx2g` and the port assumes the JVM's default heap, so the number is
/// compared by its shape here and by its formula in `command_line`'s unit test.
fn heap_masked(tool: &str, line: &str) -> String {
    if tool != "EstimateLibraryComplexity" {
        return line.to_string();
    }
    let Some((set, defaults)) = line.split_once("    ") else {
        return line.to_string();
    };
    let defaults: Vec<String> = defaults
        .split(' ')
        .map(|pair| match pair.strip_prefix("MAX_RECORDS_IN_RAM=") {
            Some(n) if !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => {
                "MAX_RECORDS_IN_RAM=<heap>".to_string()
            }
            _ => pair.to_string(),
        })
        .collect();
    format!("{set}    {}", defaults.join(" "))
}

#[test]
fn every_command_line_is_the_one_picard_records() {
    let text = corpus();
    let argvs: BTreeMap<String, String> = rows(&text, "argv").into_iter().collect();
    let lines = rows(&text, "cl");
    assert_eq!(lines.len(), 216, "four command lines for each of 54 tools");
    let mut failures = Vec::new();
    for (case, expected) in &lines {
        assert_ne!(
            expected, "FAILED",
            "{case}: the reference refused its own command line"
        );
        let tool = case.split('/').next().unwrap();
        let options = tool_options(tool).unwrap_or_else(|| panic!("{tool} has an option table"));
        let args: Vec<String> = argvs[case].split('\n').map(str::to_string).collect();
        let got = legacy_command_line(options, &args);
        let (got, expected) = (heap_masked(tool, &got), heap_masked(tool, expected));
        if got != expected {
            failures.push(format!("{case}\n  picard: {expected}\n  port:   {got}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} differ:\n{}",
        failures.len(),
        lines.len(),
        failures.join("\n")
    );
}

/// The binary that ports a tool, built alongside this test.
fn binary(tool: &str) -> &'static str {
    match tool {
        "MarkDuplicates" => env!("CARGO_BIN_EXE_MarkDuplicates"),
        "MarkDuplicatesWithMateCigar" => env!("CARGO_BIN_EXE_MarkDuplicatesWithMateCigar"),
        "SimpleMarkDuplicatesWithMateCigar" => {
            env!("CARGO_BIN_EXE_SimpleMarkDuplicatesWithMateCigar")
        }
        "UmiAwareMarkDuplicatesWithMateCigar" => {
            env!("CARGO_BIN_EXE_UmiAwareMarkDuplicatesWithMateCigar")
        }
        "DownsampleSam" => env!("CARGO_BIN_EXE_DownsampleSam"),
        "PositionBasedDownsampleSam" => env!("CARGO_BIN_EXE_PositionBasedDownsampleSam"),
        "IntervalListTools" => env!("CARGO_BIN_EXE_IntervalListTools"),
        "CollectQualityYieldMetrics" => env!("CARGO_BIN_EXE_CollectQualityYieldMetrics"),
        other => panic!("no binary for {other}"),
    }
}

const DIR: &str = "/work/hl";

/// A BAM's header `@PG` lines and each record's name and `PG` tag, as the dump prints them.
fn bam_program_lines(path: &Path) -> (String, String) {
    let raw = std::fs::read(path).expect("output");
    let plain = htsjdk_bgzf::decompress_all(&raw).expect("bgzf");
    let reader = htsjdk_bam::reader::BamReader::new(&plain).expect("bam");
    let header = reader.header.text.clone();
    let mut pg = String::new();
    for line in header.encode().lines() {
        if line.starts_with("@PG") {
            pg.push_str(line);
            pg.push('\n');
        }
    }
    let mut tags = String::new();
    for record in reader {
        let record = record.expect("record");
        let value = match record.tags.get(htsjdk_bam::tag::Tag::new(b"PG")) {
            Some(htsjdk_bam::tag::TagValue::Str(s)) => s.clone(),
            _ => "null".to_string(),
        };
        tags.push_str(&format!("{} {value}\n", record.read_name));
    }
    (pg, tags)
}

/// `Date.toString()`: `EEE MMM dd HH:mm:ss zzz yyyy`.
fn is_java_date(s: &str) -> bool {
    let parts: Vec<&str> = s.split(' ').collect();
    parts.len() == 6
        && ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"].contains(&parts[0])
        && [
            "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
        ]
        .contains(&parts[1])
        && parts[2].len() == 2
        && parts[2].bytes().all(|b| b.is_ascii_digit())
        && parts[3].len() == 8
        && parts[3]
            .split(':')
            .all(|f| f.len() == 2 && f.bytes().all(|b| b.is_ascii_digit()))
        && !parts[4].is_empty()
        && parts[5].len() == 4
        && parts[5].bytes().all(|b| b.is_ascii_digit())
}

#[test]
fn the_runs_write_the_records_and_headers_picard_writes() {
    let text = corpus();
    let dir = std::env::temp_dir().join(format!("header-lines-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let local = dir.to_str().unwrap().to_string();

    let inputs: BTreeMap<String, String> = rows(&text, "input").into_iter().collect();
    let (header, records) = htsjdk_bam::sam_file::read_sam(&inputs["in.sam"]).expect("the input");
    let mut writer = htsjdk_bam::writer::BamWriter::new(Vec::new(), &header).unwrap();
    for record in &records {
        writer.write(record).unwrap();
    }
    std::fs::write(dir.join("in.bam"), writer.finish().unwrap()).unwrap();
    std::fs::write(dir.join("in.interval_list"), &inputs["in.interval_list"]).unwrap();

    let exits: BTreeMap<String, String> = rows(&text, "exit").into_iter().collect();
    let pgs: BTreeMap<String, String> = rows(&text, "pg").into_iter().collect();
    let read_pgs: BTreeMap<String, String> = rows(&text, "readpg").into_iter().collect();
    let metrics = rows(&text, "metrics");
    let runs = rows(&text, "run");
    assert_eq!(runs.len(), 12);

    let mut failures = Vec::new();
    for (id, argv) in &runs {
        assert_eq!(exits[id], "0", "{id}: the reference run failed");
        let mut words = argv.split('\n');
        let tool = words.next().unwrap();
        let args: Vec<String> = words.map(|w| w.replace(DIR, &local)).collect();
        let out = Command::new(binary(tool))
            .args(&args)
            .output()
            .expect("run");
        if !out.status.success() {
            failures.push(format!(
                "{id}: exit {:?}: {}",
                out.status.code(),
                String::from_utf8_lossy(&out.stderr)
            ));
            continue;
        }
        let back = |s: String| s.replace(&local, DIR);
        if let Some(expected) = pgs.get(id) {
            let output = args
                .iter()
                .find_map(|a| a.strip_prefix("O="))
                .expect("an output");
            let (pg, tags) = if output.ends_with(".interval_list") {
                let written = std::fs::read_to_string(output).unwrap();
                let pg: String = written
                    .lines()
                    .filter(|l| l.starts_with("@PG"))
                    .map(|l| format!("{l}\n"))
                    .collect();
                (pg, String::new())
            } else {
                bam_program_lines(Path::new(output))
            };
            if &back(pg.clone()) != expected {
                failures.push(format!(
                    "{id} @PG\n  picard:\n{expected}  port:\n{}",
                    back(pg)
                ));
            }
            // By name: the tag is what is measured here, and SimpleMarkDuplicatesWithMateCigar
            // writes records of one position in an order of its own.
            let by_name = |text: &str| {
                let mut lines: Vec<&str> = text.lines().collect();
                lines.sort_unstable();
                lines.join("\n")
            };
            if let Some(expected) = read_pgs.get(id) {
                if by_name(&tags) != by_name(expected) {
                    failures.push(format!(
                        "{id} read PG\n  picard:\n{expected}  port:\n{tags}"
                    ));
                }
            }
        }
        for (case, expected) in metrics
            .iter()
            .filter(|(c, _)| c.starts_with(&format!("{id}/")))
        {
            let file = case.split('/').nth(1).unwrap();
            let written = std::fs::read_to_string(dir.join(file)).expect("the metrics file");
            let head: Vec<String> = written
                .lines()
                .take(4)
                .map(|l| back(l.to_string()))
                .collect();
            let want: Vec<&str> = expected.split('\n').collect();
            for (i, (got, want)) in head.iter().zip(&want).enumerate() {
                let same = match (
                    got.strip_prefix("# Started on: "),
                    want.strip_prefix("# Started on: "),
                ) {
                    (Some(g), Some(w)) => is_java_date(g) && is_java_date(w),
                    _ => got == want,
                };
                if !same {
                    failures.push(format!(
                        "{case} line {}\n  picard: {want}\n  port:   {got}",
                        i + 1
                    ));
                }
            }
            if head.len() != want.len() {
                failures.push(format!(
                    "{case}: {} header lines, picard {}",
                    head.len(),
                    want.len()
                ));
            }
        }
    }
    std::fs::remove_dir_all(&dir).ok();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
