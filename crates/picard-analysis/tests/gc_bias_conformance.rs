//! Conformance for `CollectGcBiasMetrics` against Picard, through the port's own binary.
//!
//! Each case carries the input BAM, the reference bases, and **both** metrics files Picard
//! produced: the detail file of 101 GC bins and the summary. The binary runs on the same BAM and
//! the same FASTA (with the `.fai` and `.dict` the dump wrote beside it) and both of its files are
//! compared, two outputs from one run rather than the convenient one.
//!
//! The cases exist to hold the findings in `src/gc.rs` in place. `forward_and_reverse_same_bases`
//! pins that the two strands are charged to different windows; `past_last_window_start` pins that
//! a read whose window was never computed lands in GC bin 0; `reverse_near_contig_start` pins
//! that such a read is dropped from the bins while still counting as aligned.

use std::io::Read;

use std::process::Command;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gc_bias.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn rows(kind: &str) -> Vec<(String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let name = it.next()?.to_string();
            let payload = it.next()?.to_string();
            (k == kind).then_some((name, payload))
        })
        .collect()
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

fn from_hex(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
}

const CANONICALIZED: &[&str] = &["# CollectGcBiasMetrics ", "# Started on:"];

fn canonicalize(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        if CANONICALIZED.iter().any(|p| line.starts_with(*p)) {
            out.push_str("<canonicalized>\n");
        } else {
            out.push_str(line);
        }
    }
    out
}

#[test]
fn every_detail_and_summary_file_matches_picards() {
    let bams = rows("bam");
    let fastas = rows("fasta");
    let details = rows("detail");
    let summaries = rows("summary");
    assert_eq!(bams.len(), details.len());
    assert_eq!(bams.len(), summaries.len());
    assert!(
        bams.len() >= 9,
        "expected at least 9 cases, got {}",
        bams.len()
    );

    let dir = std::env::temp_dir().join(format!("gcbias-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut failures = Vec::new();
    for (i, (name, bam_hex)) in bams.iter().enumerate() {
        let reference = &fastas[i].1;
        let bam = dir.join(format!("{name}.bam"));
        std::fs::write(&bam, from_hex(bam_hex)).unwrap();
        // The FASTA, its `.fai` and its `.dict`, as `GcBiasDump.emit` writes them.
        let fasta = dir.join(format!("{name}.fasta"));
        let mut text = String::from(">chr1\n");
        for line in reference.as_bytes().chunks(60) {
            text.push_str(std::str::from_utf8(line).unwrap());
            text.push('\n');
        }
        std::fs::write(&fasta, text).unwrap();
        std::fs::write(
            dir.join(format!("{name}.fasta.fai")),
            format!("chr1\t{}\t6\t60\t61\n", reference.len()),
        )
        .unwrap();
        std::fs::write(
            dir.join(format!("{name}.dict")),
            format!(
                "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{}\n",
                reference.len()
            ),
        )
        .unwrap();
        let detail = dir.join(format!("{name}.detail.txt"));
        let summary = dir.join(format!("{name}.summary.txt"));
        let ran = Command::new(env!("CARGO_BIN_EXE_CollectGcBiasMetrics"))
            .arg(format!("INPUT={}", bam.display()))
            .arg(format!("OUTPUT={}", detail.display()))
            .arg(format!("SUMMARY_OUTPUT={}", summary.display()))
            .arg(format!(
                "CHART_OUTPUT={}",
                dir.join(format!("{name}.pdf")).display()
            ))
            .arg(format!("REFERENCE_SEQUENCE={}", fasta.display()))
            .arg("ASSUME_SORTED=true")
            .output()
            .unwrap();
        if !ran.status.success() {
            failures.push(format!(
                "{name}: exit {:?}, Picard 0: {}",
                ran.status.code(),
                String::from_utf8_lossy(&ran.stderr)
            ));
            continue;
        }

        for (label, path, theirs_raw) in [
            ("detail", &detail, &details[i].1),
            ("summary", &summary, &summaries[i].1),
        ] {
            let ours = canonicalize(&std::fs::read_to_string(path).unwrap_or_default());
            let theirs = canonicalize(&unescape(theirs_raw));
            if ours != theirs {
                let at = ours
                    .lines()
                    .zip(theirs.lines())
                    .position(|(a, b)| a != b)
                    .unwrap_or(0);
                failures.push(format!(
                    "{name} [{label}]: first difference at line {at}\n  picard: {}\n  ours  : {}",
                    theirs.lines().nth(at).unwrap_or("<end>"),
                    ours.lines().nth(at).unwrap_or("<end>")
                ));
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        failures.is_empty(),
        "{} divergences over {} cases:\n{}",
        failures.len(),
        bams.len(),
        failures.join("\n")
    );
}
