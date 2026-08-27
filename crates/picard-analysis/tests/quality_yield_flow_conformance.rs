//! Conformance for `CollectQualityYieldMetricsFlow` against Picard 3.4.0.
//!
//! Each case carries the file the tool read, as SAM without its header, and the metrics table it
//! wrote. The flows of these fixtures are their homopolymer runs and every flow of a read carries
//! that read's own base quality, so the port is handed the same numbers the tool derived from the
//! flow matrices without those matrices being ported.
//!
//! # What this suite is for
//!
//!  * **the unit being the flow and not the base**;
//!  * **a skipped secondary or supplementary read being counted nowhere, TOTAL_READS included**;
//!  * **a vendor-failed read being counted in TOTAL_READS and left out of PF_READS**;
//!  * **the two include arguments being independent**;
//!  * **PF_Q20_FLOWS counting the 30s as well, so it is never smaller than PF_Q30_FLOWS**;
//!  * **PF_Q20_EQUIVALENT_YIELD being the sum of the qualities divided by twenty**;
//!  * **the mean flow count being an integer division**;
//!  * **and an empty file producing a table of zeros rather than none.**

use std::io::Read as _;

use picard_analysis::quality_yield_flow::{
    collect, cycle_count, flow_quality, outcome, Metrics, Outcome, Record, MAX_QUAL,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("quality_yield_flow.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("{kind} {name} is in the corpus"));
    unescape(&line[prefix.len()..])
}

/// One case's records, with the flow qualities the fixture gives them: one per homopolymer run,
/// each carrying the read's own base quality, which is uniform across every fixture read.
fn reads(text: &str, case: &str) -> Vec<(bool, bool, bool, Vec<u8>)> {
    field(text, "sam", case)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("a flag word");
            let bases = columns[9].as_bytes();
            let quality = columns[10].as_bytes()[0] - 33;
            let mut runs = 0;
            for (index, base) in bases.iter().enumerate() {
                if index == 0 || bases[index - 1] != *base {
                    runs += 1;
                }
            }
            (
                flags & 0x100 != 0,
                flags & 0x800 != 0,
                flags & 0x200 != 0,
                vec![quality; runs],
            )
        })
        .collect()
}

/// The metrics table of one case, as its nine values.
fn metrics(text: &str, case: &str) -> Vec<String> {
    let table = field(text, "metrics", case);
    let mut lines = table.lines().filter(|line| !line.is_empty());
    let header = lines.next().expect("a header line");
    assert_eq!(
        header.split('\t').next(),
        Some("TOTAL_READS"),
        "{case} header"
    );
    lines
        .next()
        .expect("a value line")
        .split('\t')
        .map(str::to_string)
        .collect()
}

/// The reference writes a whole number without a decimal point and a fraction with one.
fn format(value: f64) -> String {
    if value == value.trunc() {
        format!("{}", value as i64)
    } else {
        format!("{value}")
    }
}

fn as_written(m: &Metrics) -> Vec<String> {
    vec![
        m.total_reads.to_string(),
        m.pf_reads.to_string(),
        m.mean_pf_read_number_of_flows().to_string(),
        m.pf_flows.to_string(),
        m.pf_q20_flows.to_string(),
        format(m.pct_pf_q20_flows()),
        m.pf_q30_flows.to_string(),
        format(m.pct_pf_q30_flows()),
        m.pf_q20_equivalent_yield.to_string(),
    ]
}

/// The include flags each case ran under, which the dump names on the command line.
const CASES: &[(&str, bool, bool)] = &[
    ("one-read", false, false),
    ("homopolymers-and-singles", false, false),
    ("low-quality", false, false),
    ("between-thresholds", false, false),
    ("vendor-failed", false, false),
    ("secondary-out", false, false),
    ("secondary-in", true, false),
    ("supplementary-out", false, false),
    ("supplementary-in", false, true),
    ("secondary-in-supplementary-out", true, false),
    ("with-histogram", false, false),
    ("empty", false, false),
];

fn collected(text: &str, case: &str, secondary: bool, supplementary: bool) -> Metrics {
    let reads = reads(text, case);
    let records: Vec<Record> = reads
        .iter()
        .map(
            |(is_secondary, is_supplementary, fails, qualities)| Record {
                read_length: qualities.len(),
                secondary: *is_secondary,
                supplementary: *is_supplementary,
                fails_vendor_quality: *fails,
                flow_qualities: qualities,
            },
        )
        .collect();
    collect(&records, secondary, supplementary)
}

/// Every case's nine metrics are what the port reaches.
#[test]
fn every_case_reaches_the_same_metrics() {
    let text = corpus();
    for (case, secondary, supplementary) in CASES {
        let ours = collected(&text, case, *secondary, *supplementary);
        assert_eq!(as_written(&ours), metrics(&text, case), "{case}");
    }
}

/// The unit is the flow: eight bases in four homopolymers are four flows and eight singles are
/// eight, so two reads of eight bases each give twelve flows and not sixteen.
#[test]
fn the_unit_is_the_flow_and_not_the_base() {
    let text = corpus();
    assert_eq!(collected(&text, "one-read", false, false).pf_flows, 4);
    let both = collected(&text, "homopolymers-and-singles", false, false);
    assert_eq!(both.pf_flows, 12);
    assert_eq!(both.pf_reads, 2);
    // And the mean is an integer division, so twelve over two is six exactly.
    assert_eq!(both.mean_pf_read_number_of_flows(), 6);
}

/// A skipped read is counted nowhere, TOTAL_READS included, which a vendor-failed read is not.
#[test]
fn a_skipped_read_is_counted_nowhere() {
    let text = corpus();
    let skipped = collected(&text, "secondary-out", false, false);
    assert_eq!(skipped.total_reads, 1);
    assert_eq!(skipped.pf_reads, 1);
    let failed = collected(&text, "vendor-failed", false, false);
    assert_eq!(failed.total_reads, 2);
    assert_eq!(failed.pf_reads, 1);
    // Its flows are counted nowhere either.
    assert_eq!(failed.pf_flows, 8);
}

/// The two include arguments are independent: naming the secondary one leaves the supplementary
/// records out all the same.
#[test]
fn the_two_include_arguments_are_independent() {
    let text = corpus();
    let one = collected(&text, "secondary-in-supplementary-out", true, false);
    assert_eq!(one.total_reads, 2);
    assert_eq!(one.pf_flows, 16);
    // With both named the third read would be counted as well.
    let both = collected(&text, "secondary-in-supplementary-out", true, true);
    assert_eq!(both.total_reads, 3);
    assert_eq!(both.pf_flows, 24);
}

/// PF_Q20_FLOWS counts the 30s too, so it is never smaller than PF_Q30_FLOWS, and a quality
/// between the two thresholds counts for one and not the other.
#[test]
fn the_twenties_include_the_thirties() {
    let text = corpus();
    for (case, _, _) in CASES {
        let ours = collected(&text, case, false, false);
        assert!(ours.pf_q20_flows >= ours.pf_q30_flows, "{case}");
    }
    let low = collected(&text, "low-quality", false, false);
    assert_eq!(low.pf_q20_flows, 16);
    assert_eq!(low.pf_q30_flows, 8);
    let middle = collected(&text, "between-thresholds", false, false);
    assert_eq!(middle.pf_q20_flows, 8);
    assert_eq!(middle.pf_q30_flows, 0);
}

/// The yield is the sum of the qualities divided by twenty, so it moves with the qualities and
/// not only with the flows: four flows at 40 give 8 and eight flows at 25 give 10.
#[test]
fn the_yield_is_the_sum_over_twenty() {
    let text = corpus();
    let one = collected(&text, "one-read", false, false);
    assert_eq!(one.pf_flows, 4);
    assert_eq!(one.pf_q20_equivalent_yield, 8);
    let middle = collected(&text, "between-thresholds", false, false);
    assert_eq!(middle.pf_flows, 8);
    assert_eq!(middle.pf_q20_equivalent_yield, 10);
    assert!(middle.pf_flows > one.pf_flows);
    assert!(middle.pf_q20_equivalent_yield > one.pf_q20_equivalent_yield);
}

/// An empty file writes a table of zeros, the derived fields included.
#[test]
fn an_empty_file_writes_zeros() {
    let text = corpus();
    let empty = collected(&text, "empty", false, false);
    assert_eq!(empty, Metrics::default());
    assert_eq!(empty.mean_pf_read_number_of_flows(), 0);
    assert_eq!(empty.pct_pf_q20_flows(), 0.0);
    assert_eq!(metrics(&text, "empty"), vec!["0"; 9]);
}

/// A read of no bases is skipped before anything is counted, as a secondary one is.
#[test]
fn a_read_of_no_bases_is_skipped() {
    let empty = Record {
        read_length: 0,
        secondary: false,
        supplementary: false,
        fails_vendor_quality: false,
        flow_qualities: &[],
    };
    assert_eq!(outcome(&empty, true, true), Outcome::Skipped);
    assert_eq!(collect(&[empty], true, true), Metrics::default());
}

/// The per-flow quality is the phred of the error probability, and a probability of exactly zero
/// answers the ceiling rather than an infinity.
#[test]
fn the_flow_quality_is_a_clamped_phred() {
    assert_eq!(flow_quality(0.0), MAX_QUAL as u8);
    assert_eq!(flow_quality(0.01), 20);
    assert_eq!(flow_quality(0.001), 30);
    // Below the ceiling's probability the answer is clamped rather than growing.
    assert_eq!(flow_quality(1e-30), MAX_QUAL as u8);
    assert_eq!(flow_quality(1.0), 0);
}

/// The histogram's cycles are four flows each, rounded up.
#[test]
fn the_histogram_cycles_are_four_flows() {
    assert_eq!(cycle_count(4), 1);
    assert_eq!(cycle_count(5), 2);
    assert_eq!(cycle_count(8), 2);
    assert_eq!(cycle_count(0), 0);
}
