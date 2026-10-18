//! Conformance for the input shapes the tools' own suites never carried, against Picard 3.5.0.
//!
//! Every case of `tools/inputshapes-conformance/InputShapesDump.java` is run through the port's own
//! binary, on the same SAM and with the same arguments, and must end the same way: the same exit
//! code or the same exception (class and message), and the same bytes in the output file, which
//! for a run that throws is what the writer's buffers had let through.

use std::collections::BTreeMap;
use std::io::Read;
use std::process::Command;

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/inputshapes.txt.gz");
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("corpus"))
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
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

/// The rows of each case, by kind, in the order the dump wrote the cases.
fn cases() -> Vec<(String, BTreeMap<String, String>)> {
    let mut cases: Vec<(String, BTreeMap<String, String>)> = Vec::new();
    for line in corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
    {
        let mut fields = line.splitn(3, '\t');
        let (kind, case) = (fields.next().unwrap(), fields.next().unwrap());
        let payload = unescape(fields.next().unwrap_or(""));
        if cases.last().map(|(name, _)| name.as_str()) != Some(case) {
            cases.push((case.to_string(), BTreeMap::new()));
        }
        cases
            .last_mut()
            .unwrap()
            .1
            .insert(kind.to_string(), payload);
    }
    cases
}

fn binary(tool: &str) -> &'static str {
    match tool {
        "CleanSam" => env!("CARGO_BIN_EXE_CleanSam"),
        "CollectSamErrorMetrics" => env!("CARGO_BIN_EXE_CollectSamErrorMetrics"),
        "RevertOriginalBaseQualitiesAndAddMateCigar" => {
            env!("CARGO_BIN_EXE_RevertOriginalBaseQualitiesAndAddMateCigar")
        }
        other => panic!("no binary for {other}"),
    }
}

#[test]
fn every_shape_ends_as_picard_ends_it() {
    let dir = std::env::temp_dir().join(format!("inputshapes-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let cases = cases();
    assert!(cases.len() >= 19, "the corpus has {} cases", cases.len());
    let mut failures = Vec::new();
    for (case, rows) in &cases {
        let ext = if rows.contains_key("output_bytes") {
            "bam"
        } else {
            "sam"
        };
        let input = dir.join(format!("{case}.in.sam"));
        // A tool that writes one file per suffix gets a prefix in a directory of its own.
        let per_suffix = rows.keys().any(|kind| kind.starts_with("file:"))
            || rows["tool"] == "CollectSamErrorMetrics";
        let case_dir = dir.join(case);
        std::fs::create_dir_all(&case_dir).unwrap();
        let output = if per_suffix {
            case_dir.join("out")
        } else {
            dir.join(format!("{case}.out.{ext}"))
        };
        std::fs::write(&input, &rows["input"]).unwrap();
        let _ = std::fs::remove_file(&output);
        let reference = dir.join(format!("{case}.fasta"));
        if let Some(fasta) = rows.get("reference") {
            std::fs::write(&reference, fasta).unwrap();
        }
        let mut command = Command::new(binary(&rows["tool"]));
        command.arg(format!("INPUT={}", input.display()));
        command.arg(format!("OUTPUT={}", output.display()));
        command.args(
            rows["args"]
                .split_whitespace()
                .map(|arg| arg.replace("<ref>", &reference.display().to_string())),
        );
        let run = command.output().unwrap();
        let stderr = String::from_utf8_lossy(&run.stderr);
        let thrown = stderr
            .lines()
            .find_map(|l| l.strip_prefix("Exception in thread \"main\" "));
        match (rows.get("exit"), rows.get("thrown")) {
            (Some(code), None) => {
                if run.status.code().map(|c| c.to_string()).as_deref() != Some(code) {
                    failures.push(format!(
                        "{case}: exit {:?}, Picard {code}; {stderr}",
                        run.status.code()
                    ));
                }
            }
            (None, Some(expected)) => {
                if thrown != Some(expected.as_str()) || run.status.code() != Some(1) {
                    failures.push(format!("{case}: threw {thrown:?}, Picard {expected:?}"));
                }
            }
            other => panic!("{case}: no outcome row: {other:?}"),
        }
        if per_suffix {
            // Every file the reference left, and no other, with the same table.
            let expected: BTreeMap<String, String> = rows
                .iter()
                .filter_map(|(kind, body)| {
                    Some((kind.strip_prefix("file:")?.to_string(), body.clone()))
                })
                .collect();
            let mut ours = BTreeMap::new();
            for entry in std::fs::read_dir(&case_dir).unwrap() {
                let path = entry.unwrap().path();
                let name = path.file_name().unwrap().to_string_lossy().into_owned();
                let body: String = std::fs::read_to_string(&path)
                    .unwrap()
                    .lines()
                    .filter(|line| !line.starts_with("# "))
                    .map(|line| format!("{line}\n"))
                    .collect();
                ours.insert(name.trim_start_matches("out.").to_string(), body);
            }
            if ours != expected {
                let names = |m: &BTreeMap<String, String>| m.keys().cloned().collect::<Vec<_>>();
                let differing: Vec<&String> = expected
                    .iter()
                    .filter(|(name, body)| ours.get(*name) != Some(body))
                    .map(|(name, _)| name)
                    .collect();
                failures.push(format!(
                    "{case}: files differ: ours {:?}, picard {:?}, differing {differing:?}",
                    names(&ours),
                    names(&expected)
                ));
            }
            continue;
        }
        let ours = std::fs::read(&output).ok();
        if let Some(expected) = rows.get("output") {
            if ours.as_deref() != Some(expected.as_bytes()) {
                failures.push(format!(
                    "{case}: output differs\n  ours:   {:?}\n  picard: {expected:?}",
                    ours.map(|o| String::from_utf8_lossy(&o).into_owned())
                ));
            }
        } else if let Some(bytes) = rows.get("output_bytes") {
            if ours.as_ref().map(|o| o.len().to_string()).as_deref() != Some(bytes.as_str()) {
                failures.push(format!(
                    "{case}: {:?} output bytes, Picard {bytes}",
                    ours.map(|o| o.len())
                ));
            }
        } else if rows.contains_key("output_absent") && ours.is_some() {
            failures.push(format!("{case}: wrote an output Picard did not"));
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
