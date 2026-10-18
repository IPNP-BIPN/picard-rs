//! Conformance for what a screen of the non-BAM tools on large inputs found the port answering
//! differently from Picard 3.5.0, against `tools/nonbamshapes-conformance/NonBamShapesDump.java`.
//!
//! The dump carries every input it hands a tool (base64, under its path in the run's directory)
//! and every command line (with that directory written `<dir>`), so this test rebuilds the inputs
//! byte for byte, runs the port's binary on the same command line, and renders the same
//! observations: an exit code or an exception, metrics files without their `# ` lines, what an
//! index says (its kind, its properties, its sequences), whether a date has `Date.toString()`'s
//! shape, bait and FASTA names, and a dictionary without its `UR` field.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR};

fn golden() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/nonbam_shapes.txt.gz");
    let file = std::fs::File::open(&path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden is gzip");
    text
}

/// The dump's escaping undone: `\\`, `\t` and `\n`.
fn unescape(payload: &str) -> String {
    let mut out = String::with_capacity(payload.len());
    let mut chars = payload.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

fn base64(text: &str) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::with_capacity(text.len() * 3 / 4);
    let mut buffer = 0u32;
    let mut bits = 0;
    for byte in text.bytes().filter(|b| *b != b'=') {
        let value = ALPHABET.iter().position(|a| *a == byte).expect("base64") as u32;
        buffer = buffer << 6 | value;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((buffer >> bits) as u8);
        }
    }
    out
}

/// `GtcToVcf` is `gtc-to-vcf`: the class name's words, lower-cased and joined by hyphens.
fn binary(class: &str) -> PathBuf {
    let mut name = String::new();
    for (at, c) in class.chars().enumerate() {
        if c.is_ascii_uppercase() && at > 0 {
            name.push('-');
        }
        name.push(c.to_ascii_lowercase());
    }
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_GtcToVcf"));
    exe.with_file_name(name)
}

/// `exit` with the code, or `thrown` with the exception the port printed for the JVM's.
fn run(dir: &Path, class: &str, args: &str) -> (String, String) {
    let args: Vec<String> = args
        .split(' ')
        .map(|a| a.replace("<dir>", &dir.display().to_string()))
        .collect();
    let output = Command::new(binary(class))
        .args(&args)
        .current_dir(dir)
        .output()
        .expect("the binary runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    let thrown = stderr
        .lines()
        .find_map(|l| l.strip_prefix("Exception in thread \"main\" "));
    match thrown {
        Some(exception) if !output.status.success() => (
            "thrown".to_string(),
            exception.replace(&dir.display().to_string(), "<dir>"),
        ),
        _ => (
            "exit".to_string(),
            output.status.code().unwrap_or(-1).to_string(),
        ),
    }
}

fn metrics(kind: &str, path: &Path) -> (String, String) {
    match std::fs::read_to_string(path) {
        Err(_) => (format!("metrics_absent_{kind}"), String::new()),
        Ok(text) => (
            format!("metrics_{kind}"),
            text.lines()
                .filter(|l| !l.starts_with("# "))
                .collect::<Vec<_>>()
                .join("\n"),
        ),
    }
}

fn index(vcf: &Path) -> Vec<(String, String)> {
    let Ok(bytes) = std::fs::read(format!("{}.idx", vcf.display())) else {
        return vec![("idx_absent".to_string(), String::new())];
    };
    let index = TribbleIndex::read(&bytes).expect("an index the port wrote reads back");
    let kind = match index.index_type {
        LINEAR => "LinearIndex",
        INTERVAL_TREE => "IntervalTreeIndex",
        other => panic!("index type {other}"),
    };
    vec![
        ("idx_type".to_string(), kind.to_string()),
        (
            "idx_properties".to_string(),
            index
                .properties
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect::<Vec<_>>()
                .join("\n"),
        ),
        (
            "idx_sequences".to_string(),
            index.sequence_names().join(","),
        ),
    ]
}

/// `EEE MMM dd HH:mm:ss zzz yyyy`.
fn is_date(text: &str) -> bool {
    let parts: Vec<&str> = text.split(' ').collect();
    let word = |s: &str| {
        s.len() == 3
            && s.as_bytes()[0].is_ascii_uppercase()
            && s.as_bytes()[1..].iter().all(u8::is_ascii_lowercase)
    };
    let digits = |s: &str, n: usize| s.len() == n && s.bytes().all(|b| b.is_ascii_digit());
    parts.len() == 6
        && word(parts[0])
        && word(parts[1])
        && digits(parts[2], 2)
        && parts[3].split(':').count() == 3
        && parts[3].split(':').all(|p| digits(p, 2))
        && !parts[4].is_empty()
        && parts[4]
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"+-:".contains(&b))
        && digits(parts[5], 4)
}

fn date(path: &Path, prefix: &str) -> String {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    match text.lines().find_map(|l| l.strip_prefix(prefix)) {
        None => format!("no line starting {prefix}"),
        Some(rest) if is_date(rest) => "a date".to_string(),
        Some(_) => "not a date".to_string(),
    }
}

/// What the dump observed after each case's run, from the files the case wrote.
fn observe(dir: &Path, case: &str) -> Vec<(String, String)> {
    let at = |p: &str| dir.join(p);
    let mut rows = Vec::new();
    match case {
        "cvcm_types" => {
            rows.push(metrics(
                "detail",
                &at("cvcm/cvcm.variant_calling_detail_metrics"),
            ));
            rows.push(metrics(
                "summary",
                &at("cvcm/cvcm.variant_calling_summary_metrics"),
            ));
        }
        "fp_overlap" => {
            rows.push(metrics(
                "summary",
                &at("fp/check.fingerprinting_summary_metrics"),
            ));
            rows.push(metrics(
                "detail",
                &at("fp/check.fingerprinting_detail_metrics"),
            ));
        }
        "idx_sortvcf" => rows.extend(index(&at("idx/sorted.vcf"))),
        "idx_gtctovcf" => {
            rows.extend(index(&at("idx/gtc.vcf")));
            rows.push(("date".to_string(), date(&at("idx/gtc.vcf"), "##fileDate=")));
        }
        "idx_combine_create_index_true" => rows.extend(index(&at("idx/combined_true.vcf"))),
        "idx_combine_create_index_false" => rows.extend(index(&at("idx/combined_false.vcf"))),
        "idx_combine_default" => rows.extend(index(&at("idx/combined_default.vcf"))),
        "idx_haplotypes" => rows.extend(index(&at("idx/haplotypes.vcf"))),
        "idx_liftover" => rows.extend(index(&at("idx/lifted.vcf"))),
        "idx_liftover_disable_sort" => rows.extend(index(&at("idx/lifted_unsorted.vcf"))),
        "idx_extract" | "idx_contaminant" => {
            let vcf = if case == "idx_extract" {
                "idx/extracted.vcf"
            } else {
                "idx/contaminant.vcf"
            };
            rows.extend(index(&at(vcf)));
            rows.push(("date".to_string(), date(&at(vcf), "##fileDate=")));
        }
        "date_extended_manifest" => rows.push((
            "date".to_string(),
            date(&at("idx/report.txt"), "Generated on: "),
        )),
        "baits_thousands" => {
            let text = std::fs::read_to_string(at("baits/out/design.baits.interval_list"))
                .unwrap_or_default();
            rows.push((
                "baits".to_string(),
                text.lines()
                    .filter(|l| !l.starts_with('@'))
                    .collect::<Vec<_>>()
                    .join("\n"),
            ));
            let fasta =
                std::fs::read_to_string(at("baits/out/design.design.fasta")).unwrap_or_default();
            rows.push((
                "fasta_names".to_string(),
                fasta
                    .lines()
                    .filter(|l| l.starts_with('>'))
                    .collect::<Vec<_>>()
                    .join("\n"),
            ));
        }
        "dict_default_output" => match std::fs::read_to_string(at("dict/r.dict")) {
            Ok(text) => rows.push((
                "dict".to_string(),
                text.lines()
                    .map(|l| {
                        l.split('\t')
                            .filter(|f| !f.starts_with("UR:"))
                            .collect::<Vec<_>>()
                            .join("\t")
                    })
                    .collect::<Vec<_>>()
                    .join("\n"),
            )),
            Err(_) => rows.push(("dict_absent".to_string(), String::new())),
        },
        "barcodes_gzip" => rows.push((
            "gzip_bytes".to_string(),
            std::fs::read(at("barcodes/s_1_1101_barcode.txt.gz"))
                .unwrap_or_default()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect(),
        )),
        "dict_not_fasta_output" => rows.push((
            "dict_written".to_string(),
            at("dict/other.dict").exists().to_string(),
        )),
        _ => {}
    }
    rows
}

#[test]
fn every_case_answers_as_the_reference_did() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("nonbam-shapes-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let rows: Vec<(&str, &str, String)> = text
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let mut f = l.splitn(3, '\t');
            let kind = f.next().unwrap();
            let case = f.next().unwrap();
            (kind, case, unescape(f.next().unwrap_or("")))
        })
        .collect();
    // The directories the dump made for each group of cases, outputs and all.
    for group in ["cvcm", "idx", "baits/out", "fp", "dict", "cirm", "barcodes"] {
        std::fs::create_dir_all(dir.join(group)).unwrap();
    }
    for (kind, name, payload) in &rows {
        if *kind == "input" {
            let path = dir.join(name);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(&path, base64(payload)).unwrap();
        }
    }

    // The cases in the dump's order, each one's rows as the port produces them.
    let expected: Vec<(String, String, String)> = rows
        .iter()
        .filter(|(kind, _, _)| *kind != "input")
        .map(|(k, c, p)| (k.to_string(), c.to_string(), p.clone()))
        .collect();
    let mut produced: Vec<(String, String, String)> = Vec::new();
    let mut at = 0;
    while at < expected.len() {
        let (kind, case, class) = &expected[at];
        assert_eq!(kind, "tool", "each case opens with its tool");
        let args = &expected[at + 1];
        assert_eq!(args.0, "args");
        produced.push(expected[at].clone());
        produced.push(args.clone());
        let (outcome, value) = run(&dir, class, &args.2);
        produced.push((outcome, case.clone(), value));
        for (k, v) in observe(&dir, case) {
            produced.push((k, case.clone(), v));
        }
        // The liftover's reject is observed under a case name of its own.
        if case == "idx_liftover" {
            for (k, v) in index(&dir.join("idx/rejected.vcf")) {
                produced.push((k, "idx_liftover_reject".to_string(), v));
            }
        }
        at += 1;
        while at < expected.len() && expected[at].0 != "tool" {
            at += 1;
        }
    }
    for (e, p) in expected.iter().zip(&produced) {
        assert_eq!(e, p);
    }
    assert_eq!(expected.len(), produced.len(), "the same number of rows");
    let _ = std::fs::remove_dir_all(&dir);
}
