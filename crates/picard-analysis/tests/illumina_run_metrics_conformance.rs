//! Conformance for the two run summaries against Picard 3.4.0.
//!
//! Goldens from `tools/illumina-conformance/`: `CollectIlluminaBasecallingMetrics` and
//! `CollectHiSeqXPfFailMetrics`, both over the same four-cluster fixture.
//!
//! # What this suite is for
//!
//!  * **the read structure deciding how many cycles count as bases**;
//!  * **the PF columns coming off the filter file rather than the basecalls**;
//!  * **a lane row counting every cluster whatever its barcode**;
//!  * **the failure classes partitioning the filter's own answer**;
//!  * **and `--N_CYCLES` doing nothing at all.**

use std::io::Read;

use picard_analysis::illumina_basecalls::{position_in_name, Cluster};
use picard_analysis::illumina_files::{decode_basecall, parse_read_structure, BaseCall};
use picard_analysis::illumina_run_metrics::{
    base_cycles, basecalling_metrics, classify, pf_fail_metrics, PfFailure, CYCLES_JUDGED,
};

fn corpus(name: &str) -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.txt.gz"));
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

fn table(text: &str, kind: &str, case: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let body = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = body.lines();
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .map(str::to_string)
        .collect();
    let rows = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    (header, rows)
}

fn column(header: &[String], rows: &[Vec<String>], row: usize, name: &str) -> String {
    let index = header
        .iter()
        .position(|column| column == name)
        .unwrap_or_else(|| panic!("{name}"));
    rows[row][index].clone()
}

/// The fixture's clusters, with the barcode `2T2B` would read off each.
fn clusters() -> Vec<(Cluster, Option<String>)> {
    let cycles = ["ACGT", "ACGT", "AACC", "GGTT"];
    (0..4)
        .map(|index| {
            let calls: Vec<BaseCall> = cycles
                .iter()
                .map(|cycle| {
                    decode_basecall(
                        (30 << 2)
                            | match cycle.as_bytes()[index] {
                                b'A' => 0,
                                b'C' => 1,
                                b'G' => 2,
                                _ => 3,
                            },
                    )
                })
                .collect();
            let barcode: String = calls[2..4].iter().map(|call| call.base as char).collect();
            (
                Cluster {
                    calls,
                    passed_filter: index != 3,
                    x: position_in_name(100.0 * (index as f32 + 1.0), 0.0).0,
                    y: position_in_name(0.0, 200.0 * (index as f32 + 1.0)).1,
                },
                Some(barcode),
            )
        })
        .collect()
}

/// The whole lane as one row, and the read structure deciding what a base is.
#[test]
fn the_basecalling_metrics_are_the_goldens() {
    let text = corpus("basecalling_metrics");
    for (case, structure, bases) in [
        ("the-whole-lane", "4T", 16),
        ("a-barcode-segment-without-barcodes", "2T2B", 8),
        ("a-skipped-segment", "2T2S", 8),
    ] {
        let parsed = parse_read_structure(structure).expect("a structure");
        let rows = basecalling_metrics(1, &clusters(), &parsed, &[]);
        let (header, recorded) = table(&text, "metrics", case);
        assert_eq!(recorded.len(), 1, "{case}");
        assert_eq!(
            column(&header, &recorded, 0, "TOTAL_BASES"),
            bases.to_string(),
            "{case}"
        );
        assert_eq!(rows[0].total_bases, bases, "{case}");
        // The PF columns are the filter file's answer: three of the four clusters passed.
        assert_eq!(rows[0].pf_clusters, 3, "{case}");
        assert_eq!(column(&header, &recorded, 0, "PF_CLUSTERS"), "3", "{case}");
        assert_eq!(rows[0].total_clusters, 4, "{case}");
    }
    assert_eq!(
        base_cycles(&parse_read_structure("2T2B").expect("a structure")),
        2
    );
}

/// Split by barcode: a row apiece and a row for the lane.
#[test]
fn the_barcodes_are_rows_of_their_own() {
    let text = corpus("basecalling_metrics");
    let structure = parse_read_structure("2T2B").expect("a structure");
    let rows = basecalling_metrics(
        1,
        &clusters(),
        &structure,
        &["AG".to_string(), "CT".to_string()],
    );
    let (header, recorded) = table(&text, "metrics", "by-barcode");
    assert_eq!(rows.len(), recorded.len());
    for (row, index) in rows.iter().zip(0..recorded.len()) {
        assert_eq!(
            column(&header, &recorded, index, "MOLECULAR_BARCODE_SEQUENCE_1"),
            row.barcode
        );
        assert_eq!(
            column(&header, &recorded, index, "TOTAL_CLUSTERS"),
            row.total_clusters.to_string()
        );
        assert_eq!(
            column(&header, &recorded, index, "PF_CLUSTERS"),
            row.pf_clusters.to_string()
        );
    }
    // The lane's own row counts every cluster, which is why it is four where each barcode is two.
    assert_eq!(rows.last().expect("the lane row").total_clusters, 4);
}

/// The failure classes partition the filter's answer, and `--N_CYCLES` changes nothing.
#[test]
fn the_pf_failures_are_classified() {
    let text = corpus("pf_fail_metrics");
    // The fixture for this tool has twenty-four cycles, which is what the tool always looks at.
    let clusters: Vec<Cluster> = clusters()
        .into_iter()
        .map(|(cluster, _)| {
            let mut calls = Vec::new();
            for index in 0..CYCLES_JUDGED {
                calls.push(cluster.calls[index % cluster.calls.len()]);
            }
            Cluster { calls, ..cluster }
        })
        .collect();
    let row = pf_fail_metrics("1101", &clusters, CYCLES_JUDGED);
    let (header, recorded) = table(&text, "metrics", "the-whole-lane.pf.pffail_summary_metrics");
    // The golden writes an `All` row and a row per tile; the port's is the tile's.
    let tile = recorded
        .iter()
        .position(|values| values[0] == "1101")
        .expect("the tile row");
    assert_eq!(
        column(&header, &recorded, tile, "READS"),
        row.reads.to_string()
    );
    assert_eq!(
        column(&header, &recorded, tile, "PF_FAIL_READS"),
        row.pf_fail_reads.to_string()
    );
    assert_eq!(
        column(&header, &recorded, tile, "PF_FAIL_POLYCLONAL"),
        row.pf_fail_polyclonal.to_string()
    );
    assert_eq!(row.pf_fail_reads, 1);
    assert_eq!(row.pf_fail_polyclonal, 1);
    // The one failing cluster is polyclonal: its bases are called and its qualities are good.
    assert_eq!(classify(&clusters[3], CYCLES_JUDGED), PfFailure::Polyclonal);

    // And the three `--N_CYCLES` cases of the golden are the same output three times, because the
    // read structure was built from the argument's default before the parser could assign it.
    let two = field(&text, "metrics", "two-cycles.pf.pffail_summary_metrics").expect("the golden");
    let forty =
        field(&text, "metrics", "forty-cycles.pf.pffail_summary_metrics").expect("the golden");
    let default =
        field(&text, "metrics", "the-whole-lane.pf.pffail_summary_metrics").expect("the golden");
    assert_eq!(two, default);
    assert_eq!(forty, default);
}
