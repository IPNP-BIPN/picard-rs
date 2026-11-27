//! Conformance for `IlluminaBasecallsToFastq` and `IlluminaBasecallsToSam` against Picard,
//! through the port's own binaries.
//!
//! Every run of `tools/illumina-conformance/BasecallsToFastqDump.java` and `BasecallsToSamDump.java`
//! is replayed over the same four-cluster run directory with the same arguments, and the rows the
//! dump prints are rebuilt from what the binary wrote: the files, their text (the BAM through
//! `ViewSam`), or the error. The rebuilt text must equal the golden line for line.

mod common;

use common::{extract_barcodes, run, write_run, written, Dump};

/// `BasecallsToFastqDump.run`.
fn fastq(dump: &mut Dump, case: &str, multiplex_params: Option<&str>, extra: &[&str]) {
    let dir = dump.dir();
    let basecalls = write_run(&dir.join("run"), 4);
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let mut argv = vec![
        format!("BASECALLS_DIR={}", basecalls.display()),
        "LANE=1".to_string(),
        "RUN_BARCODE=run17".to_string(),
        "MACHINE_NAME=machine".to_string(),
        "FLOWCELL_BARCODE=flowcell".to_string(),
        "COMPRESS_OUTPUTS=false".to_string(),
    ];
    argv.extend(extra.iter().map(|argument| argument.to_string()));
    if !extra
        .iter()
        .any(|argument| argument.starts_with("READ_STRUCTURE="))
    {
        argv.push("READ_STRUCTURE=4T".to_string());
    }
    match multiplex_params {
        None => argv.push(format!("OUTPUT_PREFIX={}", out.join("reads").display())),
        Some(params) => {
            extract_barcodes(&dir, &basecalls);
            argv.push(format!("BARCODES_DIR={}", basecalls.display()));
            let file = dir.join("params.tsv");
            std::fs::write(&file, params.replace("<out>", &out.display().to_string())).unwrap();
            argv.push(format!("MULTIPLEX_PARAMS={}", file.display()));
        }
    }
    let ran = run(env!("CARGO_BIN_EXE_IlluminaBasecallsToFastq"), &argv);
    if dump.failed(case, &dir, &ran, false) {
        return;
    }
    let files = written(&out);
    dump.emit("files", case, &files.join(" "));
    for file in files {
        let text = std::fs::read_to_string(out.join(&file)).unwrap();
        dump.emit("fastq", &format!("{case}.{file}"), &text);
    }
}

/// `BasecallsToSamDump.run`.
fn sam(dump: &mut Dump, case: &str, library_params: Option<&str>, extra: &[&str]) {
    let dir = dump.dir();
    let basecalls = write_run(&dir.join("run"), 4);
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let mut argv = vec![
        format!("BASECALLS_DIR={}", basecalls.display()),
        "LANE=1".to_string(),
        "RUN_BARCODE=run17".to_string(),
        "SEQUENCING_CENTER=centre".to_string(),
        "PLATFORM=ILLUMINA".to_string(),
        "SORT=false".to_string(),
    ];
    argv.extend(extra.iter().map(|argument| argument.to_string()));
    if !extra
        .iter()
        .any(|argument| argument.starts_with("READ_STRUCTURE="))
    {
        argv.push("READ_STRUCTURE=4T".to_string());
    }
    match library_params {
        None => {
            argv.push(format!("OUTPUT={}", out.join("reads.bam").display()));
            argv.push("SAMPLE_ALIAS=sample1".to_string());
            argv.push("LIBRARY_NAME=libraryA".to_string());
        }
        Some(params) => {
            extract_barcodes(&dir, &basecalls);
            let file = dir.join("library.params");
            std::fs::write(&file, params.replace("<out>", &out.display().to_string())).unwrap();
            argv.push(format!("LIBRARY_PARAMS={}", file.display()));
            argv.push(format!("BARCODES_DIR={}", basecalls.display()));
        }
    }
    let ran = run(env!("CARGO_BIN_EXE_IlluminaBasecallsToSam"), &argv);
    if dump.failed(case, &dir, &ran, true) {
        return;
    }
    let files = written(&out);
    dump.emit("files", case, &files.join(" "));
    for file in files {
        let viewed = run(
            env!("CARGO_BIN_EXE_ViewSam"),
            &[format!("INPUT={}", out.join(&file).display())],
        );
        assert_eq!(viewed.code, 0, "ViewSam {case}.{file}: {}", viewed.said);
        // `getSAMString()` of each read group, which carries no newline, then the records.
        let lines: Vec<&str> = viewed.stdout.lines().collect();
        let groups: String = lines
            .iter()
            .filter(|line| line.starts_with("@RG\t"))
            .copied()
            .collect();
        let records: String = lines
            .iter()
            .filter(|line| !line.starts_with('@') && !line.is_empty())
            .map(|line| format!("{line}\n"))
            .collect();
        dump.emit("header", &format!("{case}.{file}"), &groups);
        dump.emit("sam", &format!("{case}.{file}"), &records);
    }
}

#[test]
fn basecalls_to_fastq_matches_the_golden() {
    let mut dump = Dump::new("basecalls_to_fastq");
    fastq(&mut dump, "one-read-of-four", None, &[]);
    fastq(
        &mut dump,
        "two-reads-of-two",
        None,
        &["READ_STRUCTURE=2T2T"],
    );
    fastq(
        &mut dump,
        "a-barcode-segment",
        None,
        &["READ_STRUCTURE=2T2B"],
    );
    fastq(
        &mut dump,
        "a-skipped-segment",
        None,
        &["READ_STRUCTURE=2T2S"],
    );
    fastq(
        &mut dump,
        "without-the-non-pf-reads",
        None,
        &["INCLUDE_NON_PF_READS=false"],
    );
    fastq(
        &mut dump,
        "the-other-read-name-format",
        None,
        &["READ_NAME_FORMAT=ILLUMINA"],
    );
    fastq(
        &mut dump,
        "split-by-barcode",
        Some("OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n<out>/second\tCT\n"),
        &["READ_STRUCTURE=2T2B"],
    );
    fastq(
        &mut dump,
        "an-undeclared-barcode",
        Some("OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n"),
        &["READ_STRUCTURE=2T2B"],
    );
    fastq(
        &mut dump,
        "an-undeclared-barcode-ignored",
        Some("OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n"),
        &["READ_STRUCTURE=2T2B", "IGNORE_UNEXPECTED_BARCODES=true"],
    );
    fastq(
        &mut dump,
        "a-row-for-the-rest",
        Some("OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n<out>/rest\t\n"),
        &["READ_STRUCTURE=2T2B"],
    );
    dump.check("basecalls_to_fastq");
}

#[test]
fn basecalls_to_sam_matches_the_golden() {
    const TWO: &str = "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n\
                       AG\t<out>/first.bam\tsampleA\tlibraryA\n\
                       CT\t<out>/second.bam\tsampleB\tlibraryB\n";
    const ONE: &str = "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n\
                       AG\t<out>/first.bam\tsampleA\tlibraryA\n";
    let mut dump = Dump::new("basecalls_to_sam");
    sam(&mut dump, "one-read-of-four", None, &[]);
    sam(
        &mut dump,
        "two-reads-of-two",
        None,
        &["READ_STRUCTURE=2T2T"],
    );
    sam(
        &mut dump,
        "a-barcode-segment-without-params",
        None,
        &["READ_STRUCTURE=2T2B"],
    );
    sam(
        &mut dump,
        "a-skipped-segment",
        None,
        &["READ_STRUCTURE=2T2S"],
    );
    sam(
        &mut dump,
        "without-the-non-pf-reads",
        None,
        &["INCLUDE_NON_PF_READS=false"],
    );
    sam(
        &mut dump,
        "the-barcode-in-the-read-group",
        Some(TWO),
        &["READ_STRUCTURE=2T2B", "INCLUDE_BC_IN_RG_TAG=true"],
    );
    sam(
        &mut dump,
        "split-by-barcode",
        Some(TWO),
        &["READ_STRUCTURE=2T2B"],
    );
    sam(
        &mut dump,
        "an-undeclared-barcode",
        Some(ONE),
        &["READ_STRUCTURE=2T2B"],
    );
    sam(
        &mut dump,
        "an-undeclared-barcode-ignored",
        Some(ONE),
        &["READ_STRUCTURE=2T2B", "IGNORE_UNEXPECTED_BARCODES=true"],
    );
    dump.check("basecalls_to_sam");
}
