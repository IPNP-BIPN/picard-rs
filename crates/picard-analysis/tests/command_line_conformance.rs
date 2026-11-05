//! Conformance for what every binary does before its tool runs, against Picard 3.5.0: refuse an
//! argument the tool does not declare, and answer `-h`, `--help` and `--version`.
//!
//! Goldens from `tools/commandline-conformance/`: `CommandLineDump.java` records what each probe
//! command line printed on each stream and exited with, and `CommandLineDeclarationsDump.java` the
//! names each of the two parsers accepts.
//!
//! The binaries print the probe golden's texts with the user's words put in, so the first test is
//! the one that matters: for every tool and every probe, the port's answer to the probe's own
//! command line is the reference's, byte for byte. The other two hold the ends that test cannot
//! see: that every binary asks for its own tool, and that no row of a covering array is refused.

use std::io::Read;
use std::path::Path;

use picard_analysis::argument_check::{answer, declared, probe_args, unescape, Answer};

fn gunzip(name: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

#[test]
fn every_probe_is_answered_as_the_reference_answered_it() {
    let golden = gunzip("command_line.txt.gz");
    let mut compared = 0;
    for line in golden.lines().filter(|l| l.starts_with("run\t")) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (tool, probe) = (fields[1], fields[2]);
        let expected = Answer {
            exit: fields[3].parse().unwrap(),
            stdout: unescape(fields[4]),
            stderr: unescape(fields[5]),
        };
        let words: Vec<String> = probe_args(probe).iter().map(|w| w.to_string()).collect();
        assert!(
            !words.is_empty(),
            "a probe this test does not know: {probe}"
        );
        let port = answer(tool, &words).unwrap_or_else(|| panic!("{tool} {probe}: not refused"));
        assert_eq!(port.exit, expected.exit, "{tool} {probe}: exit");
        assert_eq!(port.stdout, expected.stdout, "{tool} {probe}: stdout");
        if port.stderr != expected.stderr {
            let at = port
                .stderr
                .chars()
                .zip(expected.stderr.chars())
                .take_while(|(a, b)| a == b)
                .count();
            let near = |t: &str| {
                t.chars()
                    .skip(at.saturating_sub(40))
                    .take(120)
                    .collect::<String>()
            };
            panic!(
                "{tool} {probe}: stderr differs at character {at}\n port:      {:?}\n reference: {:?}",
                near(&port.stderr),
                near(&expected.stderr)
            );
        }
        compared += 1;
    }
    // 121 tools, ten probes each, less the two bare-word probes of CompareSAMs, which takes
    // positional words.
    assert_eq!(compared, 121 * 10 - 2);
}

/// The binary each tool runs as, from the covering arrays' manifest.
fn ports() -> Vec<(String, String)> {
    let manifest = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/conformance/manifest.json"),
    )
    .unwrap();
    let pair = regex::Regex::new(r#""tool": "([^"]+)",\s*"port": "([^"]+)""#).unwrap();
    pair.captures_iter(&manifest)
        .map(|c| (c[1].to_string(), c[2].to_string()))
        .collect()
}

#[test]
fn every_binary_checks_its_command_line_as_its_own_tool() {
    let ports = ports();
    assert_eq!(ports.len(), 121);
    for (tool, port) in &ports {
        let source = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("src/bin/{port}.rs")),
        )
        .unwrap();
        // `SetNmAndUqTags` is the deprecated name of `SetNmMdAndUqTags`, run by the same binary.
        let own = if tool == "SetNmAndUqTags" {
            "SetNmMdAndUqTags"
        } else {
            tool
        };
        let call = format!("picard_analysis::argument_check::check(\"{own}\");");
        let body = source.split_once("fn main()").expect("a main").1;
        let first = body.lines().nth(1).unwrap_or_default().trim();
        assert_eq!(
            first, call,
            "{port} must check its command line before anything else"
        );
        assert!(
            declared(own).is_some(),
            "{own} is in the declarations golden"
        );
    }
}

#[test]
fn no_covering_array_row_is_refused() {
    let label = regex::Regex::new(r#""(--[^"]+)": "#).unwrap();
    let held =
        regex::Regex::new(r#"\{[^{}]*?"argument": "(--[^"]+)"[^{}]*?"source": "fixture"[^{}]*?\}"#)
            .unwrap();
    for (tool, _) in ports() {
        let own = if tool == "SetNmAndUqTags" {
            "SetNmMdAndUqTags".to_string()
        } else {
            tool.clone()
        };
        let array = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(format!("../../tools/coverage/arrays/{tool}.t2.json")),
        )
        .unwrap();
        let mut names: Vec<String> = Vec::new();
        for labels in array.split("\"labels\": {").skip(1) {
            let block = labels.split('}').next().unwrap();
            names.extend(label.captures_iter(block).map(|c| c[1].to_string()));
        }
        names.extend(held.captures_iter(&array).map(|c| c[1].to_string()));
        names.sort();
        names.dedup();
        // The runner hands a port `NAME=value`, the dashes stripped; the positional pseudo-argument
        // goes last as bare words, which only CompareSAMs has.
        let args: Vec<String> = names
            .iter()
            .filter(|n| !n.starts_with("--[NA"))
            .map(|n| format!("{}=x", n.trim_start_matches('-')))
            .collect();
        // `x` stands for any value, which the option's type may refuse: only the name is checked.
        for arg in &args {
            let refused = answer(&own, std::slice::from_ref(arg));
            assert!(
                refused.is_none_or(|a| !a.stderr.contains("ERROR: Unrecognized option")),
                "{tool}: {arg}"
            );
        }
    }
}

/// What Picard 3.5.0 printed for each of these command lines, checked against the reference.
#[test]
fn a_value_the_legacy_parser_refuses_is_refused_with_its_message() {
    let cases = [
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt O=m2.txt",
            "ERROR: Option 'O' cannot be specified more than once.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam OUTPUT=m1.txt o=m2.txt",
            "ERROR: Option 'O' cannot be specified more than once.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt I=x.sam",
            "ERROR: Option 'I' cannot be specified more than once.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt O=m2.txt BAD=1",
            "ERROR: Option 'O' cannot be specified more than once.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt BAD=1 O=m2.txt",
            "ERROR: Unrecognized option: BAD",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=null",
            "ERROR: non-null value must be provided for 'O'.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt VALIDATION_STRINGENCY=strict",
            "ERROR: 'strict' is not a valid value for ValidationStringency.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt VALIDATION_STRINGENCY=",
            "ERROR: '' is not a valid value for ValidationStringency.",
        ),
        (
            "CollectQualityYieldMetrics",
            "I=in.sam O=m1.txt MAX_RECORDS_IN_RAM=99999999999",
            "ERROR: Problem constructing Integer from the string '99999999999'.",
        ),
        (
            "MarkDuplicates",
            "I=in.sam O=d.sam M=dm.txt MAX_OPTICAL_DUPLICATE_SET_SIZE=abc",
            "ERROR: Problem constructing Long from the string 'abc'.",
        ),
        (
            "DownsampleSam",
            "I=in.sam O=d.sam PROBABILITY=abc",
            "ERROR: Problem constructing Double from the string 'abc'.",
        ),
    ];
    for (tool, line, message) in cases {
        let words: Vec<String> = line.split(' ').map(str::to_string).collect();
        let refused = answer(tool, &words).unwrap_or_else(|| panic!("{tool} {line}: not refused"));
        assert_eq!(refused.exit, 1, "{tool} {line}");
        assert!(
            refused
                .stderr
                .contains(&format!("\n{message}\n\nUSAGE: {tool} ")),
            "{tool} {line}:\n{}",
            refused.stderr
        );
        assert!(
            refused.stderr.contains(&format!(
                "{tool} -{}",
                line.split(' ').next().unwrap().replace('=', " ")
            )),
            "{tool} {line}: the banner prints the words as typed"
        );
    }
    for (tool, line) in [
        ("SortSam", "I= in.sam O= s.sam SORT_ORDER=coordinate"),
        ("DownsampleSam", "I=in.sam O=d.sam PROBABILITY=1e-1d"),
        ("CollectAlignmentSummaryMetrics", "I=in.sam O=a.txt R=null"),
        (
            "CollectAlignmentSummaryMetrics",
            "I=in.sam O=a.txt ADAPTER_SEQUENCE=A ADAPTER_SEQUENCE=C",
        ),
    ] {
        let words: Vec<String> = line.split(' ').map(str::to_string).collect();
        assert_eq!(answer(tool, &words), None, "{tool} {line}");
    }
}
