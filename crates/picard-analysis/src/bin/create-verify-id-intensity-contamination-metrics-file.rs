//! `CreateVerifyIDIntensityContaminationMetricsFile` as a runnable binary: the covering array's
//! port side.
//!
//! Ports `picard.arrays.CreateVerifyIDIntensityContaminationMetricsFile.doWork` at tag 3.4.0. The
//! parse is `picard_analysis::create_verify_id_intensity_metrics`; this writes
//! `OUTPUT.verifyidintensity_metrics`.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::create_verify_id_intensity_metrics::{output_name, parse, ParseError};
use picard_analysis::fingerprint::reference_path;
use picard_analysis::metrics_cli::{thrown, Args};

const COLUMNS: &[&str] = &["ID", "PCT_MIX", "LLK", "LLK0"];

struct Row(Vec<Value>);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.arrays.VerifyIDIntensityContaminationMetrics"
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
        Err(ParseError::Unrecognised(line)) => thrown(&format!(
            "picard.PicardException: Unrecognized line: {line} in {}",
            reference_path(&input)
        )),
        Err(ParseError::EndedEarly) => thrown("java.lang.NullPointerException"),
    };
    let mut file = MetricsFile::new();
    file.add_header("CreateVerifyIDIntensityContaminationMetricsFile <command line>");
    file.add_header("Started on: <timestamp>");
    for r in rows {
        file.add_metric(&Row(vec![
            Value::Long(i64::from(r.id)),
            Value::Double(r.percent_mix),
            Value::Double(r.log_likelihood),
            Value::Double(r.log_likelihood_zero),
        ]));
    }
    if let Err(e) = std::fs::write(&output, file.write()) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
