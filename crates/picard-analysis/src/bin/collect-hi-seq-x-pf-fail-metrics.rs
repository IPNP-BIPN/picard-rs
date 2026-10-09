//! `CollectHiSeqXPfFailMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.quality.CollectHiSeqXPfFailMetrics` at tag 3.4.0 over
//! `picard_analysis::illumina_reader`.
//!
//! The read structure is `N_CYCLES + "T"` built in a field initialiser, which runs before the
//! command line is parsed, so it is always twenty-four cycles. A failing cluster is classified by
//! its no-calls (counted as `.`, which the reader never produces: a no-call is `N`) and by how many
//! of its qualities exceed two; each one is written to the detailed file when an unseeded
//! `Random` draws below `PROB_EXPLICIT_READS`, which is only reproducible at 0 and 1.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::illumina_reader::Run;
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};

const TOOL: &str = "CollectHiSeqXPfFailMetrics";

const SUMMARY: &[&str] = &[
    "TILE",
    "READS",
    "PF_FAIL_READS",
    "PCT_PF_FAIL_READS",
    "PF_FAIL_EMPTY",
    "PCT_PF_FAIL_EMPTY",
    "PF_FAIL_POLYCLONAL",
    "PCT_PF_FAIL_POLYCLONAL",
    "PF_FAIL_MISALIGNED",
    "PCT_PF_FAIL_MISALIGNED",
    "PF_FAIL_UNKNOWN",
    "PCT_PF_FAIL_UNKNOWN",
];

struct Row {
    class: &'static str,
    columns: &'static [&'static str],
    values: Vec<Value>,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        self.class
    }
    fn columns(&self) -> &[&'static str] {
        self.columns
    }
    fn values(&self) -> Vec<Value> {
        self.values.clone()
    }
}

#[derive(Default, Clone)]
struct Summary {
    reads: i64,
    fail: i64,
    empty: i64,
    polyclonal: i64,
    misaligned: i64,
    unknown: i64,
}

impl Summary {
    fn add(&mut self, o: &Summary) {
        self.reads += o.reads;
        self.fail += o.fail;
        self.empty += o.empty;
        self.polyclonal += o.polyclonal;
        self.misaligned += o.misaligned;
        self.unknown += o.unknown;
    }

    fn row(&self, tile: &str) -> Row {
        let pct = |n: i64| {
            if self.reads != 0 {
                n as f64 / self.reads as f64
            } else {
                0.0
            }
        };
        Row {
            class: "picard.illumina.quality.CollectHiSeqXPfFailMetrics$PFFailSummaryMetric",
            columns: SUMMARY,
            values: vec![
                Value::Str(tile.to_string()),
                Value::Long(self.reads),
                Value::Long(self.fail),
                Value::Double(pct(self.fail)),
                Value::Long(self.empty),
                Value::Double(pct(self.empty)),
                Value::Long(self.polyclonal),
                Value::Double(pct(self.polyclonal)),
                Value::Long(self.misaligned),
                Value::Double(pct(self.misaligned)),
                Value::Long(self.unknown),
                Value::Double(pct(self.unknown)),
            ],
        }
    }
}

fn main() {
    let args = Args::from_env(&[("B", "BASECALLS_DIR"), ("O", "OUTPUT"), ("L", "LANE")]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let output = args.required("OUTPUT");
    let lane: i32 = args.required("LANE").parse().unwrap_or(0);
    let n_cycles = args.int("N_CYCLES", 24);
    let probability = args.double("PROB_EXPLICIT_READS", 0.0);
    let mut errors = Vec::new();
    if n_cycles < 0 {
        errors.push("Number of Cycles to look at must be greater than 0".to_string());
    }
    if !(0.0..=1.0).contains(&probability) {
        errors.push(
            "PROB_EXPLICIT_READS must be a probability, i.e., 0 <= PROB_EXPLICIT_READS <= 1"
                .to_string(),
        );
    }
    if !errors.is_empty() {
        refuse_validation(TOOL, &errors);
    }

    let run = Run::new(&basecalls, lane);
    let tiles = run.available_tiles().unwrap_or_else(|e| thrown(&e));
    let read: Vec<i32> = (1..=24).collect();
    let mut per_tile: Vec<(i32, Summary, Vec<Row>)> = Vec::new();
    for &tile in &tiles {
        let clusters = run
            .clusters(tile, std::slice::from_ref(&read), 2)
            .unwrap_or_else(|e| thrown(&e));
        let mut s = Summary::default();
        let mut detailed = Vec::new();
        for c in &clusters {
            s.reads += 1;
            if c.pf {
                continue;
            }
            s.fail += 1;
            let (bases, quals) = &c.reads[0];
            let length = bases.len() as i64;
            let ns = bases.iter().filter(|b| **b == b'.').count() as i64;
            let q = quals.iter().filter(|q| **q > 2).count() as i64;
            let class = if ns >= length - 1 {
                "MISALIGNED"
            } else if ns <= 1 && q <= length / 3 {
                "EMPTY"
            } else if ns <= 1 && q >= length / 2 {
                "POLYCLONAL"
            } else {
                "UNKNOWN"
            };
            match class {
                "MISALIGNED" => s.misaligned += 1,
                "EMPTY" => s.empty += 1,
                "POLYCLONAL" => s.polyclonal += 1,
                _ => s.unknown += 1,
            }
            if probability >= 1.0 {
                detailed.push(Row {
                    class:
                        "picard.illumina.quality.CollectHiSeqXPfFailMetrics$PFFailDetailedMetric",
                    columns: &["TILE", "X", "Y", "NUM_N", "NUM_Q_GT_TWO", "CLASSIFICATION"],
                    values: vec![
                        Value::Long(i64::from(tile)),
                        Value::Long(i64::from(c.x)),
                        Value::Long(i64::from(c.y)),
                        Value::Long(ns),
                        Value::Long(q),
                        Value::Str(class.to_string()),
                    ],
                });
            }
        }
        per_tile.push((tile, s, detailed));
    }
    let header = |file: &mut MetricsFile| {
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
    };
    if probability > 0.0 {
        let mut detail = MetricsFile::new();
        header(&mut detail);
        for (_, _, rows) in &per_tile {
            for r in rows {
                detail.add_metric(r);
            }
        }
        if let Err(e) = std::fs::write(format!("{output}.pffail_detailed_metrics"), detail.write())
        {
            thrown(&format!("htsjdk.samtools.SAMException: {e}"));
        }
    }
    let mut summary = MetricsFile::new();
    header(&mut summary);
    let mut total = Summary::default();
    for (_, s, _) in &per_tile {
        total.add(s);
    }
    summary.add_metric(&total.row("All"));
    for (tile, s, _) in &per_tile {
        summary.add_metric(&s.row(&tile.to_string()));
    }
    if let Err(e) = std::fs::write(format!("{output}.pffail_summary_metrics"), summary.write()) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
