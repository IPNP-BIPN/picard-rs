//! Conformance for `ReplaceSamHeader`'s block copy on damaged BAMs, against Picard 3.5.0.
//!
//! Every case of `tools/replacesamheader-damaged-conformance/ReplaceHeaderDamagedDump.java` is run
//! through the port's binary on the same HEADER and INPUT bytes, and must end the same way: the
//! same exit code or the same exception, and the same output bytes, or no output.

use std::collections::BTreeMap;
use std::io::Read;
use std::process::Command;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/replace_sam_header_damaged.txt.gz");
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("corpus"))
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
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
        if cases.last().map(|(name, _)| name.as_str()) != Some(case) {
            cases.push((case.to_string(), BTreeMap::new()));
        }
        let payload = fields.next().unwrap_or("").to_string();
        cases
            .last_mut()
            .unwrap()
            .1
            .insert(kind.to_string(), payload);
    }
    cases
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
}

#[test]
fn every_damaged_bam_ends_as_picard_ends_it() {
    let dir = std::env::temp_dir().join(format!("rshdamaged-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let cases = cases();
    assert_eq!(cases.len(), 10, "the corpus has {} cases", cases.len());
    let mut failures = Vec::new();
    for (case, rows) in &cases {
        let (header, ext) = match rows.get("header_bam") {
            Some(hex) => (hex, "bam"),
            None => (&rows["header_sam"], "sam"),
        };
        let header_path = dir.join(format!("{case}.header.{ext}"));
        let input = dir.join(format!("{case}.in.bam"));
        let output = dir.join(format!("{case}.out.bam"));
        std::fs::write(&header_path, unhex(header)).unwrap();
        std::fs::write(&input, unhex(&rows["input_bam"])).unwrap();
        let _ = std::fs::remove_file(&output);
        let run = Command::new(env!("CARGO_BIN_EXE_ReplaceSamHeader"))
            .arg(format!("INPUT={}", input.display()))
            .arg(format!("HEADER={}", header_path.display()))
            .arg(format!("OUTPUT={}", output.display()))
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&run.stderr).replace(&*dir.to_string_lossy(), "<dir>");
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
        let ours = std::fs::read(&output).ok();
        match rows.get("output") {
            Some(expected) if ours.as_deref() != Some(&unhex(expected)[..]) => {
                failures.push(format!("{case}: output differs from Picard's"))
            }
            None if ours.is_some() => {
                failures.push(format!("{case}: wrote an output Picard did not"))
            }
            _ => {}
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
