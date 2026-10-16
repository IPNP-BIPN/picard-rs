//! `CreateBafRegressMetricsFile` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.CreateBafRegressMetricsFile.doWork` at tag 3.4.0. The parse is
//! `picard_analysis::create_baf_regress_metrics`; this writes `OUTPUT.bafregress_metrics`.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::create_baf_regress_metrics::{output_name, parse, ParseError};
use picard_analysis::fingerprint::reference_path;
use picard_analysis::metrics_cli::{thrown, Args};

const COLUMNS: &[&str] = &[
    "SAMPLE",
    "ESTIMATE",
    "STDERR",
    "TVAL",
    "PVAL",
    "LOG10_PVAL",
    "CALL_RATE",
    "NHOM",
];

struct Row(Vec<Value>);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.arrays.BafRegressMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        self.0.clone()
    }
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT")]);
    let input = args.required("INPUT");
    let output = output_name(&args.required("OUTPUT"));
    let text = std::fs::read_to_string(&input)
        .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMException: {e}")));
    let rows = match parse(&text) {
        Ok(rows) => rows,
        Err(ParseError::Header(line)) => thrown(&format!(
            "picard.PicardException: Unrecognized header line: '{line}' in {}",
            reference_path(&input)
        )),
        Err(ParseError::EntryCount { .. }) => {
            thrown("picard.PicardException: Error parsing bafRegress Output")
        }
        Err(ParseError::Number(text)) => thrown(&format!(
            "java.lang.NumberFormatException: For input string: \"{text}\""
        )),
        Err(ParseError::EndedEarly) => thrown("java.lang.NullPointerException"),
    };
    let mut file = MetricsFile::new();
    file.add_header("CreateBafRegressMetricsFile <command line>");
    file.add_header("Started on: <timestamp>");
    for r in rows {
        file.add_metric(&Row(vec![
            Value::Str(r.sample),
            Value::Double(r.estimate),
            Value::Double(r.standard_error),
            Value::Double(r.t_value),
            Value::Double(r.p_value),
            Value::Double(r.log10_p_value),
            Value::Double(r.call_rate),
            Value::Long(i64::from(r.number_homozygous)),
        ]));
    }
    if let Err(e) = std::fs::write(&output, file.write()) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
