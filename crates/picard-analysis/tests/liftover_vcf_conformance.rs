//! Conformance for `LiftoverVcf` against Picard, through the port's own binary.
//!
//! Golden from `tools/liftovervcf-conformance/LiftoverVcfDump.java`: ten runs over a chain with
//! one forward block and one reversed one. Every run is replayed here with the same references,
//! chain, VCF and arguments, and the rows the dump prints (the exit status, the lifted and the
//! rejected variant lines, or the exception) are rebuilt from what the binary wrote and must equal
//! the golden line for line.
//!
//! # What this suite is for
//!
//!  * **a variant inside a block being renumbered by the block's offset**;
//!  * **a variant outside every block being rejected rather than dropped**;
//!  * **a reference allele the target does not carry being a different rejection**, carrying the
//!    locus and the alleles that were attempted;
//!  * **`--RECOVER_SWAPPED_REF_ALT` turning one of those into a lift**, genotypes and allele
//!    frequency and all;
//!  * **a reversed block complementing the alleles and counting from the other end**;
//!  * **the original position and alleles being recorded only where they are worth recording**;
//!  * **and the output being sorted by the target's coordinates rather than the input's.**

use std::io::Read;
use std::path::Path;
use std::process::Command;

const LENGTH: usize = 400;

fn golden() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/liftover_vcf.txt.gz");
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).expect("the golden"))
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn emit(out: &mut String, kind: &str, case: &str, payload: &str) {
    let payload = payload
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n");
    out.push_str(&format!("{kind}\t{case}\t{payload}\n"));
}

/// A contig of `ACGT` repeating, sixty bases a line, with its `.fai`.
fn write_fasta(path: &Path, contig: &str) {
    let bases: String = "ACGT".chars().cycle().take(LENGTH).collect();
    let mut text = format!(">{contig}\n");
    for line in bases.as_bytes().chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    std::fs::write(path, text).unwrap();
    std::fs::write(
        path.with_extension("fasta.fai"),
        format!("{contig}\t{LENGTH}\t{}\t60\t61\n", contig.len() + 2),
    )
    .unwrap();
    let made = Command::new(env!("CARGO_BIN_EXE_CreateSequenceDictionary"))
        .arg(format!("R={}", path.display()))
        .arg(format!("O={}", path.with_extension("dict").display()))
        .output()
        .unwrap();
    assert!(made.status.success(), "no dictionary for {contig}");
}

/// The source VCF, with the records given as `pos ref alt gt`.
fn vcf(records: &[[&str; 4]]) -> String {
    let mut text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length={LENGTH}>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
    );
    for [position, reference, alternate, genotype] in records {
        text.push_str(&format!(
            "chr1\t{position}\t.\t{reference}\t{alternate}\t100\tPASS\tAF=0.25\tGT\t{genotype}\n"
        ));
    }
    text
}

/// A VCF's variant lines, without its header.
fn records(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .split('\n')
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect::<Vec<_>>()
        .join("\n")
}

fn run(out: &mut String, base: &Path, case: &str, variants: &[[&str; 4]], extra: &[&str]) {
    let dir = base.join(case);
    std::fs::create_dir_all(&dir).unwrap();
    let dir_name = dir.display().to_string();
    write_fasta(&dir.join("source.fasta"), "chr1");
    let target = dir.join("target.fasta");
    write_fasta(&target, "chrT");
    let chain = dir.join("chain.txt");
    std::fs::write(
        &chain,
        format!(
            "chain 100 chr1 {LENGTH} + 0 100 chrT {LENGTH} + 12 112 1\n100\n\n\
             chain 100 chr1 {LENGTH} + 200 300 chrT {LENGTH} - 100 200 2\n100\n\n"
        ),
    )
    .unwrap();
    let input = dir.join("in.vcf");
    std::fs::write(&input, vcf(variants)).unwrap();
    let lifted = dir.join("out.vcf");
    let rejected = dir.join("reject.vcf");

    let ran = Command::new(env!("CARGO_BIN_EXE_LiftoverVcf"))
        .arg(format!("I={}", input.display()))
        .arg(format!("O={}", lifted.display()))
        .arg(format!("CHAIN={}", chain.display()))
        .arg(format!("REJECT={}", rejected.display()))
        .arg(format!("R={}", target.display()))
        .args(extra)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&ran.stderr);
    // The dump prints the root cause as `class:message`, the binary as Java's uncaught throwable.
    if let Some(thrown) = stderr
        .lines()
        .find_map(|line| line.strip_prefix("Exception in thread \"main\" "))
    {
        let thrown = thrown.replacen(": ", ":", 1).replace(&dir_name, "<dir>");
        emit(out, "error", case, &thrown);
        return;
    }
    emit(
        out,
        "code",
        case,
        &ran.status.code().unwrap_or(-1).to_string(),
    );
    if lifted.exists() {
        emit(out, "lifted", case, &records(&lifted));
    }
    if rejected.exists() {
        emit(out, "rejected", case, &records(&rejected));
    }
}

#[test]
fn every_run_ends_as_picards() {
    let base = std::env::temp_dir().join(format!("liftovervcf-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    let mut out = String::new();
    let inside = [["21", "A", "C", "0/1"]];
    let swapped = [["21", "C", "A", "0/1"]];

    run(&mut out, &base, "inside-a-block", &inside, &[]);
    run(
        &mut out,
        &base,
        "between-the-blocks",
        &[["150", "C", "A", "0/1"]],
        &[],
    );
    run(
        &mut out,
        &base,
        "inside-the-reversed-block",
        &[["250", "C", "A", "0/1"]],
        &[],
    );
    run(
        &mut out,
        &base,
        "a-reference-mismatch",
        &[["21", "C", "G", "0/1"]],
        &[],
    );
    run(&mut out, &base, "a-swapped-ref-and-alt", &swapped, &[]);
    run(
        &mut out,
        &base,
        "a-swapped-ref-and-alt-recovered",
        &swapped,
        &["RECOVER_SWAPPED_REF_ALT=true"],
    );
    run(
        &mut out,
        &base,
        "with-the-original-position",
        &inside,
        &["WRITE_ORIGINAL_POSITION=true"],
    );
    run(
        &mut out,
        &base,
        "with-the-original-alleles",
        &inside,
        &["WRITE_ORIGINAL_ALLELES=true"],
    );
    run(
        &mut out,
        &base,
        "with-the-original-alleles-after-a-swap",
        &swapped,
        &[
            "RECOVER_SWAPPED_REF_ALT=true",
            "WRITE_ORIGINAL_ALLELES=true",
        ],
    );
    run(
        &mut out,
        &base,
        "three-variants",
        &[
            ["250", "C", "A", "0/1"],
            ["21", "A", "C", "1/1"],
            ["150", "C", "A", "0/0"],
        ],
        &[],
    );
    let _ = std::fs::remove_dir_all(&base);

    let golden = golden();
    let ours: Vec<&str> = out.lines().collect();
    let picard: Vec<&str> = golden.lines().collect();
    let differing: Vec<String> = (0..ours.len().max(picard.len()))
        .filter(|&index| ours.get(index) != picard.get(index))
        .map(|index| {
            format!(
                "line {}:\n  ours:   {:?}\n  picard: {:?}",
                index + 1,
                ours.get(index),
                picard.get(index)
            )
        })
        .collect();
    assert!(differing.is_empty(), "{}", differing.join("\n"));
}
