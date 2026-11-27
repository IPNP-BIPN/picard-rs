//! Conformance for `CollectIlluminaBasecallingMetrics` and `CollectHiSeqXPfFailMetrics` against
//! Picard, through the port's own binaries.
//!
//! Every run of `tools/illumina-conformance/BasecallingMetricsDump.java` and `PfFailMetricsDump.java`
//! is replayed over the same run directory with the same arguments, and the rows the dump prints
//! are rebuilt from the metrics files the binary wrote (comments dropped, as the dump drops them)
//! or from its error. The rebuilt text must equal the golden line for line.

mod common;

use common::{extract_barcodes, run, table, write_run, written, Dump};

/// `BasecallingMetricsDump.run`.
fn basecalling(dump: &mut Dump, case: &str, with_barcodes: bool, extra: &[&str]) {
    let dir = dump.dir();
    let basecalls = write_run(&dir.join("run"), 4);
    let out = dir.join("metrics.txt");
    // `LANE=2` is appended after the dump's own `LANE=1`, as the dump does.
    let mut argv = vec![
        format!("BASECALLS_DIR={}", basecalls.display()),
        "LANE=1".to_string(),
        format!("OUTPUT={}", out.display()),
    ];
    argv.extend(extra.iter().map(|argument| argument.to_string()));
    if !extra
        .iter()
        .any(|argument| argument.starts_with("READ_STRUCTURE="))
    {
        argv.push("READ_STRUCTURE=4T".to_string());
    }
    if with_barcodes {
        let barcodes = extract_barcodes(&dir, &basecalls);
        argv.push(format!("BARCODES_DIR={}", basecalls.display()));
        argv.push(format!("INPUT={}", barcodes.display()));
    }
    let ran = run(
        env!("CARGO_BIN_EXE_CollectIlluminaBasecallingMetrics"),
        &argv,
    );
    if dump.failed(case, &dir, &ran, false) {
        return;
    }
    dump.emit("metrics", case, &table(&out));
}

/// `PfFailMetricsDump.run`, over twenty-four cycles: the tool reads `N_CYCLES` before the parser
/// assigns it, so it always looks at the default's twenty-four.
fn pf_fail(dump: &mut Dump, case: &str, extra: &[&str]) {
    let dir = dump.dir();
    let basecalls = write_run(&dir.join("run"), 24);
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let mut argv = vec![
        format!("BASECALLS_DIR={}", basecalls.display()),
        format!("OUTPUT={}", out.join("pf").display()),
    ];
    argv.extend(extra.iter().map(|argument| argument.to_string()));
    if !extra.iter().any(|argument| argument.starts_with("LANE=")) {
        argv.push("LANE=1".to_string());
    }
    let ran = run(env!("CARGO_BIN_EXE_CollectHiSeqXPfFailMetrics"), &argv);
    if dump.failed(case, &dir, &ran, true) {
        return;
    }
    let files = written(&out);
    dump.emit("files", case, &files.join(" "));
    for file in files {
        dump.emit(
            "metrics",
            &format!("{case}.{file}"),
            &table(&out.join(&file)),
        );
    }
}

#[test]
fn basecalling_metrics_match_the_golden() {
    let mut dump = Dump::new("basecalling_metrics");
    basecalling(&mut dump, "the-whole-lane", false, &[]);
    basecalling(&mut dump, "by-barcode", true, &["READ_STRUCTURE=2T2B"]);
    basecalling(
        &mut dump,
        "a-barcode-segment-without-barcodes",
        false,
        &["READ_STRUCTURE=2T2B"],
    );
    basecalling(
        &mut dump,
        "a-skipped-segment",
        false,
        &["READ_STRUCTURE=2T2S"],
    );
    basecalling(&mut dump, "a-lane-that-is-not-there", false, &["LANE=2"]);
    dump.check("basecalling_metrics");
}

#[test]
fn pf_fail_metrics_match_the_golden() {
    let mut dump = Dump::new("pf_fail_metrics");
    pf_fail(&mut dump, "the-whole-lane", &[]);
    pf_fail(&mut dump, "two-cycles", &["N_CYCLES=2"]);
    pf_fail(&mut dump, "forty-cycles", &["N_CYCLES=40"]);
    pf_fail(&mut dump, "with-explicit-reads", &["PROB_EXPLICIT_READS=1"]);
    pf_fail(&mut dump, "a-lane-that-is-not-there", &["LANE=2"]);
    dump.check("pf_fail_metrics");
}
