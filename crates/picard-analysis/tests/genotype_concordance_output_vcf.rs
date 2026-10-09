//! `GenotypeConcordance OUTPUT_VCF=true`: the fourth file, run through the binary.
//!
//! The expected records are read off `writeVcfTuple` and `addToGenotypes` at Picard 3.4.0. The
//! covering array measures the same file against the reference; these cases pin the three rules
//! it is easiest to get wrong: a symbolic site is left out, a cell with no contingency writes `.`,
//! and `MISSING_SITES_HOM_REF` turns the truth's no-call into a hom-ref.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

const HEADER: &str = "##fileformat=VCFv4.2\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n\
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##ALT=<ID=DEL,Description=\"Deletion\">\n\
    ##contig=<ID=chr1,length=1000>\n";

/// A directory of its own per test, so the cases can run in parallel.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "genotype-concordance-output-vcf-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(dir: &Path, args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_genotype-concordance"))
        .current_dir(dir)
        .args(args)
        .status()
        .unwrap();
    assert!(status.success());
}

/// The VCF's text, decompressed across every BGZF member.
fn output_vcf(dir: &Path) -> String {
    let raw = std::fs::read(dir.join("out.genotype_concordance.vcf.gz")).unwrap();
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(&raw[..])
        .read_to_string(&mut text)
        .unwrap();
    text
}

fn records(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

#[test]
fn one_record_per_site_with_both_genotypes_and_the_contingency() {
    let dir = scratch("sites");
    std::fs::write(
        dir.join("v.vcf"),
        format!(
            "{HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n\
             chr1\t100\trs1\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:30:10\t0/0:20:5\n\
             chr1\t200\t.\tC\tT\t.\t.\t.\tGT:GQ:DP\t1/1:40:12\t./.:.:.\n\
             chr1\t300\t.\tG\t<DEL>\t.\t.\t.\tGT\t0/1\t0/0\n\
             chr1\t400\t.\tTCA\tT\t12.5\t.\t.\tGT:GQ:DP\t0/1:9:3\t0/1:50:20\n"
        ),
    )
    .unwrap();
    run(
        &dir,
        &[
            "TRUTH_VCF=v.vcf",
            "CALL_VCF=v.vcf",
            "TRUTH_SAMPLE=sample1",
            "CALL_SAMPLE=sample2",
            "OUTPUT=out",
            "OUTPUT_VCF=true",
        ],
    );
    let text = output_vcf(&dir);
    assert!(
        text.contains(
            "##INFO=<ID=CONC_ST,Number=.,Type=String,\
             Description=\"The genotype concordance contingency state(s)\">\n"
        ),
        "{text}"
    );
    // The call is the first sample column, whatever order the genotypes were added in.
    assert!(text.contains("\tFORMAT\tcall\ttruth\n"), "{text}");
    // The ID and the FILTER are not carried over, the quality is; the symbolic site is absent.
    assert_eq!(
        records(&text),
        [
            "chr1\t100\t.\tA\tG\t50\t.\tCONC_ST=TN,FN\tGT:DP:GQ\t0/0:5:20\t0/1:10:30",
            "chr1\t200\t.\tC\tT\t.\t.\tCONC_ST=.\tGT:DP:GQ\t./.\t1/1:12:40",
            "chr1\t400\t.\tTCA\tT\t12.50\t.\tCONC_ST=TP,TN\tGT:DP:GQ\t0/1:20:50\t0/1:3:9",
        ]
    );
    assert!(dir.join("out.genotype_concordance.vcf.gz.tbi").exists());
}

#[test]
fn missing_sites_hom_ref_writes_a_missing_truth_as_hom_ref() {
    let dir = scratch("hom-ref");
    let one = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n";
    std::fs::write(
        dir.join("truth.vcf"),
        format!("{HEADER}{one}chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"),
    )
    .unwrap();
    std::fs::write(
        dir.join("call.vcf"),
        format!(
            "{HEADER}{one}chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n\
             chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/1\n"
        ),
    )
    .unwrap();
    std::fs::write(
        dir.join("all.interval_list"),
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nchr1\t1\t1000\t+\tall\n",
    )
    .unwrap();
    let common = [
        "TRUTH_VCF=truth.vcf",
        "CALL_VCF=call.vcf",
        "OUTPUT=out",
        "OUTPUT_VCF=true",
        "INTERVALS=all.interval_list",
    ];

    run(&dir, &common);
    let without = output_vcf(&dir);
    assert!(
        records(&without)[1].ends_with("\tGT\t0/1\t./."),
        "{without}"
    );

    let mut with = common.to_vec();
    with.push("MISSING_SITES_HOM_REF=true");
    run(&dir, &with);
    let with = output_vcf(&dir);
    assert!(records(&with)[1].ends_with("\tGT\t0/1\t0/0"), "{with}");
}

#[test]
fn without_the_flag_there_is_no_fourth_file() {
    let dir = scratch("off");
    std::fs::write(
        dir.join("v.vcf"),
        format!(
            "{HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n\
             chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"
        ),
    )
    .unwrap();
    run(&dir, &["TRUTH_VCF=v.vcf", "CALL_VCF=v.vcf", "OUTPUT=out"]);
    assert!(!dir.join("out.genotype_concordance.vcf.gz").exists());
    assert!(dir
        .join("out.genotype_concordance_summary_metrics")
        .exists());
}
