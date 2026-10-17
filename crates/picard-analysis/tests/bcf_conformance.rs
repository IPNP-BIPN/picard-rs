//! Conformance for the VCF tools with a `.bcf` on either side, against Picard 3.5.0.
//!
//! Golden from `tools/bcf-conformance/BcfDump.java`, which prints each case's VCF beside what
//! the tools made of it, so this runs the same binaries on the same bytes:
//!
//!  * **VcfFormatConverter** VCF to BCF, that BCF back to VCF, and to BCF again, which copies
//!    every genotype block it never decoded;
//!  * **SortVcf** and **MergeVcfs** reading that BCF and writing another;
//!  * **MakeSitesOnlyVcf** writing a BCF from the VCF;
//!  * **and the writer's refusals**, which name a Java exception and its message.

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::process::Command;

fn golden() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bcf_tools.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).expect("hex"))
        .collect()
}

/// One tool run: `(kind, payload)` as the dump would print it.
fn run(binary: &str, out: &Path, args: &[String]) -> (String, String) {
    let output = Command::new(binary).args(args).output().expect("runs");
    if output.status.success() {
        return ("bytes".into(), hex(&std::fs::read(out).expect("an output")));
    }
    let stderr = String::from_utf8_lossy(&output.stderr);
    let exception = stderr
        .lines()
        .rev()
        .find_map(|l| l.strip_prefix("Exception in thread \"main\" "))
        .map(str::to_string)
        .unwrap_or_else(|| format!("exit {}", output.status.code().unwrap_or(-1)));
    ("error".into(), exception)
}

#[test]
fn bcf_tools_match_picard() {
    let text = golden();
    let mut rows: HashMap<(String, String), String> = HashMap::new();
    let mut cases: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let parts: Vec<&str> = line.splitn(3, '\t').collect();
        if parts[0] == "vcf" {
            cases.push((parts[1].to_string(), unhex(parts[2])));
        } else {
            rows.insert(
                (parts[0].to_string(), parts[1].to_string()),
                parts[2].to_string(),
            );
        }
    }
    assert!(cases.len() >= 10, "{} cases", cases.len());

    let dir = std::env::temp_dir().join(format!("bcf-conformance-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("a directory");
    let converter = env!("CARGO_BIN_EXE_VcfFormatConverter");
    let mut failures = Vec::new();
    for (name, vcf_bytes) in &cases {
        let path = |suffix: &str| dir.join(format!("{name}{suffix}"));
        let vcf = path(".vcf");
        std::fs::write(&vcf, vcf_bytes).expect("written");
        let s = |p: &Path| p.display().to_string();
        let mut check = |step: &str, binary: &str, out: &Path, args: Vec<String>| -> bool {
            let (kind, payload) = run(binary, out, &args);
            let (want_key, ok) = if kind == "bytes" {
                ((step.to_string(), name.clone()), true)
            } else {
                (("error".to_string(), format!("{name}/{step}")), false)
            };
            if rows.get(&want_key) != Some(&payload) {
                let oracle = rows
                    .get(&(step.to_string(), name.clone()))
                    .map(|p| format!("{} bytes", p.len() / 2))
                    .or_else(|| {
                        rows.get(&("error".into(), format!("{name}/{step}")))
                            .cloned()
                    })
                    .unwrap_or_else(|| "nothing".into());
                let port = if ok {
                    format!("{} bytes", payload.len() / 2)
                } else {
                    payload
                };
                failures.push(format!("{name}/{step}: oracle {oracle}, port {port}"));
            }
            ok
        };
        let no_index = "CREATE_INDEX=false".to_string();
        let bcf = path(".bcf");
        if !check(
            "converter-bcf",
            converter,
            &bcf,
            vec![
                format!("I={}", s(&vcf)),
                format!("O={}", s(&bcf)),
                "REQUIRE_INDEX=false".into(),
                no_index.clone(),
            ],
        ) {
            continue;
        }
        for (step, suffix) in [
            ("converter-vcf", ".back.vcf"),
            ("converter-bcf2", ".again.bcf"),
        ] {
            let out = path(suffix);
            check(
                step,
                converter,
                &out,
                vec![
                    format!("I={}", s(&bcf)),
                    format!("O={}", s(&out)),
                    "REQUIRE_INDEX=false".into(),
                    no_index.clone(),
                ],
            );
        }
        let sorted = path(".sorted.bcf");
        check(
            "sortvcf-bcf",
            env!("CARGO_BIN_EXE_SortVcf"),
            &sorted,
            vec![
                format!("I={}", s(&bcf)),
                format!("O={}", s(&sorted)),
                no_index.clone(),
            ],
        );
        let merged = path(".merged.bcf");
        check(
            "mergevcfs-bcf",
            env!("CARGO_BIN_EXE_MergeVcfs"),
            &merged,
            vec![
                format!("I={}", s(&bcf)),
                format!("O={}", s(&merged)),
                no_index.clone(),
            ],
        );
        let sites = path(".sites.bcf");
        check(
            "sitesonly-bcf",
            env!("CARGO_BIN_EXE_MakeSitesOnlyVcf"),
            &sites,
            vec![
                format!("I={}", s(&vcf)),
                format!("O={}", s(&sites)),
                no_index.clone(),
            ],
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
