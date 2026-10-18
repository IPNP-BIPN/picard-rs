//! Conformance for `CollectVariantCallingMetrics` on an edge-case VCF against Picard 3.5.0.
//!
//! The golden prints the input VCF, the dbSNP one, and both tables, the detail rows in the order
//! the reference wrote them. The binary is run on the printed inputs and its tables compared whole:
//! which counter each site reaches (`VariantContext.getType()`, a half-called genotype being
//! MIXED, the dbSNP indel bitset covering a deletion's span) and the `HashMap` order of the rows.

use std::io::Read;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_variant_calling_metrics_edge.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn field(text: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("{kind}/{name}"));
    line[prefix.len()..]
        .replace("\\n", "\n")
        .replace("\\t", "\t")
        .replace("\\\\", "\\")
}

fn table(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn the_edge_cases_match_the_reference() {
    let golden = corpus();
    let dir = std::env::temp_dir().join(format!("cvcm-edge-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("in.vcf");
    let dbsnp = dir.join("dbsnp.vcf");
    std::fs::write(&input, field(&golden, "in", "edge")).unwrap();
    std::fs::write(&dbsnp, field(&golden, "in", "dbsnp")).unwrap();
    // The binary only asks that an index exist, as the reference's processor does.
    std::fs::write(dir.join("in.vcf.idx"), b"").unwrap();
    let out = dir.join("out");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_CollectVariantCallingMetrics"))
        .arg(format!("I={}", input.display()))
        .arg(format!("DBSNP={}", dbsnp.display()))
        .arg(format!("O={}", out.display()))
        .status()
        .unwrap();
    assert!(status.success());
    for (kind, extension) in [
        ("detail", "variant_calling_detail_metrics"),
        ("summary", "variant_calling_summary_metrics"),
    ] {
        let written = std::fs::read_to_string(format!("{}.{extension}", out.display())).unwrap();
        assert_eq!(table(&written), field(&golden, kind, "edge"), "{kind}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
