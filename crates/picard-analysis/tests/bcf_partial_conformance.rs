//! Conformance for what VcfFormatConverter leaves on disk when its writer refuses a record
//! partway through, against Picard 3.5.0.
//!
//! Golden from `tools/bcf-conformance/PartialOutputDump.java`. Its inputs are generated from a
//! recipe that [`vcf`] repeats, and the golden carries each input's length and MD5, so a drift in
//! the recipe fails as such rather than as a difference in the output.

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::process::Command;

use md5::{Digest, Md5};

fn golden() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bcf_partial.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn digest(bytes: &[u8]) -> String {
    let hash: String = Md5::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    format!("{} {hash}", bytes.len())
}

/// `PartialOutputDump.vcf`.
fn vcf(n: usize, bad: usize, kind: &str, pad: usize) -> String {
    let mut b = String::from(
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
         ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n\
         ##contig=<ID=chr1,length=100000000>\n",
    );
    for c in 0..pad {
        b.push_str(&format!("##contig=<ID=pad{c},length=1000>\n"));
    }
    b.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n");
    for i in 0..n {
        let (mut chrom, mut alt, mut info) = ("chr1", "G", format!("DP={};AC=1", i % 97));
        if i == bad {
            match kind {
                "undeclared" => info = "DP=3;XX=1".into(),
                "count" => info = "DP=3;AC=1,2".into(),
                "contig" => chrom = "chrZ",
                "missing" => {
                    alt = "G,T";
                    info = "DP=3;AC=1,.".into();
                }
                _ => panic!("{kind}"),
            }
        }
        b.push_str(&format!(
            "{chrom}\t{}\trs{i}\tA\t{alt}\t{}.5\tPASS\t{info}\tGT:DP\t0/1:{}\t1/1:7\t0/0:12\n",
            i * 10 + 1,
            i % 60,
            i % 30
        ));
    }
    b
}

#[test]
fn partial_outputs_match_picard() {
    let rows: HashMap<(String, String), String> = golden()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let parts: Vec<&str> = l.splitn(3, '\t').collect();
            ((parts[0].into(), parts[1].into()), parts[2].into())
        })
        .collect();
    assert_eq!(rows.len(), 24);

    let dir = std::env::temp_dir().join(format!("bcf-partial-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("a directory");
    let converter = env!("CARGO_BIN_EXE_VcfFormatConverter");
    let mut failures = Vec::new();
    let cases: [(&str, usize, usize, &str, usize); 6] = [
        ("late-undeclared", 20000, 18000, "undeclared", 0),
        ("late-count", 20000, 18000, "count", 0),
        ("late-contig", 20000, 18000, "contig", 0),
        ("late-missing", 20000, 18000, "missing", 0),
        ("big-header", 10, 5, "undeclared", 4000),
        ("small", 10, 5, "undeclared", 0),
    ];
    for (name, n, bad, kind, pad) in cases {
        let text = vcf(n, bad, kind, pad);
        let input = dir.join(format!("{name}.vcf"));
        std::fs::write(&input, &text).expect("written");
        let key = |step: &str| (step.to_string(), name.to_string());
        assert_eq!(rows[&key("input")], digest(text.as_bytes()), "{name} input");
        let mut steps = vec![("bcf", ".bcf", false), ("vcf", ".out.vcf", false)];
        if kind == "undeclared" {
            steps.push(("bcf-indexed", ".indexed.bcf", true));
            steps.push(("vcf-indexed", ".indexed.vcf", true));
        }
        for (step, suffix, index) in steps {
            let out = dir.join(format!("{name}{suffix}"));
            let output = Command::new(converter)
                .arg(format!("I={}", input.display()))
                .arg(format!("O={}", out.display()))
                .arg("REQUIRE_INDEX=false")
                .arg(format!("CREATE_INDEX={index}"))
                .output()
                .expect("runs");
            let outcome = if output.status.success() {
                "ok".to_string()
            } else {
                String::from_utf8_lossy(&output.stderr)
                    .lines()
                    .rev()
                    .find_map(|l| l.strip_prefix("Exception in thread \"main\" "))
                    .map(|e| e.replace(&dir.display().to_string(), "<dir>"))
                    .unwrap_or_else(|| format!("exit {}", output.status.code().unwrap_or(-1)))
            };
            let idx = Path::new(&format!("{}.idx", out.display())).exists();
            let port = format!(
                "{} idx={} {outcome}",
                digest(&std::fs::read(&out).unwrap_or_default()),
                if idx { "yes" } else { "no" }
            );
            if rows[&key(step)] != port {
                failures.push(format!(
                    "{name}/{step}:\n  oracle {}\n  port   {port}",
                    rows[&key(step)]
                ));
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
