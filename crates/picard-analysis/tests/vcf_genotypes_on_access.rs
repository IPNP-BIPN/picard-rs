//! A malformed genotype is refused only by what looks at it, as `LazyGenotypesContext` does.
//! The expected answers were measured on the Picard 3.5.0 oracle.

use std::process::{Command, Output};

const VCF: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t10\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:3\t0/0:4\n\
chr1\t20\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:3:9:9\t0/0:4\n\
chr1\t30\t.\tA\tC\t.\t.\t.\tGT:DP\t0/7:x\t0/0:4\n\
chr1\t40\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:3\n\
chr1\t50\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:3\t0/0:4\n";

fn records(path: &std::path::Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("chr1"))
        .map(|line| format!("{line}\n"))
        .collect()
}

fn scratch(name: &str) -> std::path::PathBuf {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn input(dir: &std::path::Path) -> String {
    let path = dir.join("in.vcf");
    std::fs::write(&path, VCF).unwrap();
    format!("I={}", path.display())
}

fn run(tool: &str, args: &[String]) -> Output {
    Command::new(tool).args(args).output().unwrap()
}

#[test]
fn a_copy_passes_a_malformed_block_through() {
    let body: String = VCF.lines().skip(5).map(|l| format!("{l}\n")).collect();
    for tool in [
        env!("CARGO_BIN_EXE_SortVcf"),
        env!("CARGO_BIN_EXE_MergeVcfs"),
        env!("CARGO_BIN_EXE_VcfFormatConverter"),
    ] {
        let dir = scratch(&format!("on_access_{}", tool.rsplit('/').next().unwrap()));
        let out = dir.join("out.vcf");
        let mut args = vec![
            input(&dir),
            format!("O={}", out.display()),
            "CREATE_INDEX=false".to_string(),
        ];
        if tool.ends_with("VcfFormatConverter") {
            args.push("REQUIRE_INDEX=false".to_string());
        }
        let output = run(tool, &args);
        assert!(output.status.success(), "{tool}: {output:?}");
        assert_eq!(records(&out), body, "{tool}");
    }
}

#[test]
fn encoding_the_genotypes_refuses_at_that_record() {
    let dir = scratch("on_access_bcf");
    let out = dir.join("out.bcf");
    let output = run(
        env!("CARGO_BIN_EXE_VcfFormatConverter"),
        &[
            input(&dir),
            format!("O={}", out.display()),
            "CREATE_INDEX=false".to_string(),
            "REQUIRE_INDEX=false".to_string(),
        ],
    );
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("There are too many keys"), "{stderr}");
}

#[test]
fn naming_a_sample_looks_and_refuses_while_none_does_not() {
    let dir = scratch("on_access_sites");
    let out = dir.join("out.vcf");
    let args = vec![
        input(&dir),
        format!("O={}", out.display()),
        "CREATE_INDEX=false".to_string(),
    ];
    let tool = env!("CARGO_BIN_EXE_MakeSitesOnlyVcf");
    assert!(run(tool, &args).status.success());
    let mut named = args;
    named.push("SAMPLE=A".to_string());
    assert_eq!(run(tool, &named).status.code(), Some(1));
}
