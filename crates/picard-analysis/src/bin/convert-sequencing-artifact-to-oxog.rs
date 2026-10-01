//! `ConvertSequencingArtifactToOxoG` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.artifacts.ConvertSequencingArtifactToOxoG.doWork` at tag 3.4.0. The
//! arithmetic is the library's (`convert_artifact_to_oxog`); what this adds is what decides the
//! file around it:
//!
//! * `customCommandLineValidation`, which derives the three file names and refuses with every
//!   missing pair at once;
//! * the rows in the order Picard writes them, libraries and contexts each in the iteration order
//!   of the `HashSet` they were collected into;
//! * the two lookups that are not checked, which a pair of tables that do not line up turns into
//!   the `NullPointerException` the JVM describes, message and all.

use std::collections::HashMap;

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::convert_artifact_to_oxog::{
    is_oxo_g, reverse_complement, CpcgMetrics, BAIT_BIAS_DETAILS_EXT, NO_BAIT_BIAS_MESSAGE,
    NO_OXOG_OUT_MESSAGE, NO_PRE_ADAPTER_MESSAGE, OXOG_METRICS_EXT, PRE_ADAPTER_DETAILS_EXT,
};
use picard_analysis::java_hash_map::JavaHashMap;
use picard_analysis::metrics_cli::{fail, refuse_validation, thrown, Args};

const TOOL: &str = "ConvertSequencingArtifactToOxoG";

/// `CollectOxoGMetrics.CpcgMetrics`, in its declaration order.
pub const CPCG_COLUMNS: [&str; 20] = [
    "SAMPLE_ALIAS",
    "LIBRARY",
    "CONTEXT",
    "TOTAL_SITES",
    "TOTAL_BASES",
    "REF_NONOXO_BASES",
    "REF_OXO_BASES",
    "REF_TOTAL_BASES",
    "ALT_NONOXO_BASES",
    "ALT_OXO_BASES",
    "OXIDATION_ERROR_RATE",
    "OXIDATION_Q",
    "C_REF_REF_BASES",
    "G_REF_REF_BASES",
    "C_REF_ALT_BASES",
    "G_REF_ALT_BASES",
    "C_REF_OXO_ERROR_RATE",
    "C_REF_OXO_Q",
    "G_REF_OXO_ERROR_RATE",
    "G_REF_OXO_Q",
];

struct Row(CpcgMetrics);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.analysis.CollectOxoGMetrics$CpcgMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &CPCG_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let m = &self.0;
        vec![
            Value::Str(m.sample_alias.clone()),
            Value::Str(m.library.clone()),
            Value::Str(m.context.clone()),
            Value::Long(m.total_sites),
            Value::Long(m.total_bases),
            Value::Long(m.ref_nonoxo_bases),
            Value::Long(m.ref_oxo_bases),
            Value::Long(m.ref_total_bases),
            Value::Long(m.alt_nonoxo_bases),
            Value::Long(m.alt_oxo_bases),
            Value::Double(m.oxidation_error_rate),
            Value::Double(m.oxidation_q),
            Value::Long(m.c_ref_ref_bases),
            Value::Long(m.g_ref_ref_bases),
            Value::Long(m.c_ref_alt_bases),
            Value::Long(m.g_ref_alt_bases),
            Value::Double(m.c_ref_oxo_error_rate),
            Value::Double(m.c_ref_oxo_q),
            Value::Double(m.g_ref_oxo_error_rate),
            Value::Double(m.g_ref_oxo_q),
        ]
    }
}

/// `MetricsFile.readBeans`, reduced to what this tool reads: each row as column name to text.
fn read_beans(path: &str) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{path}"
        ))
    });
    let mut lines = text.lines();
    for line in lines.by_ref() {
        if line.starts_with("## METRICS CLASS") {
            break;
        }
    }
    let columns: Vec<String> = match lines.next() {
        Some(header) => header.split('\t').map(str::to_string).collect(),
        None => return Vec::new(),
    };
    let mut rows = Vec::new();
    for line in lines {
        if line.is_empty() {
            break;
        }
        rows.push(
            columns
                .iter()
                .cloned()
                .zip(line.split('\t').map(str::to_string))
                .collect(),
        );
    }
    rows
}

fn long(row: &HashMap<String, String>, column: &str) -> i64 {
    row.get(column)
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| fail(&format!("could not read {column}")))
}

fn base(row: &HashMap<String, String>, column: &str) -> char {
    row.get(column)
        .and_then(|v| v.chars().next())
        .unwrap_or('\0')
}

fn text(row: &HashMap<String, String>, column: &str) -> String {
    row.get(column).cloned().unwrap_or_default()
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT_BASE"), ("O", "OUTPUT_BASE")]);
    let input_base = args.get("INPUT_BASE").map(str::to_string);
    let mut output_base = args.get("OUTPUT_BASE").map(str::to_string);
    let mut pre_adapter_in = args.get("PRE_ADAPTER_IN").map(str::to_string);
    let mut bait_bias_in = args.get("BAIT_BIAS_IN").map(str::to_string);
    let mut oxog_out = args.get("OXOG_OUT").map(str::to_string);

    // customCommandLineValidation. A null base concatenates as the text "null".
    let name = |base: &Option<String>| base.clone().unwrap_or_else(|| "null".to_string());
    let mut errors = Vec::new();
    if output_base.is_none() {
        output_base = input_base.clone();
    }
    if pre_adapter_in.is_none() {
        if input_base.is_none() {
            errors.push(NO_PRE_ADAPTER_MESSAGE.to_string());
        }
        pre_adapter_in = Some(format!("{}{PRE_ADAPTER_DETAILS_EXT}", name(&input_base)));
    }
    if bait_bias_in.is_none() {
        if input_base.is_none() {
            errors.push(NO_BAIT_BIAS_MESSAGE.to_string());
        }
        bait_bias_in = Some(format!("{}{BAIT_BIAS_DETAILS_EXT}", name(&input_base)));
    }
    if oxog_out.is_none() {
        if output_base.is_none() {
            errors.push(NO_OXOG_OUT_MESSAGE.to_string());
        }
        oxog_out = Some(format!("{}{OXOG_METRICS_EXT}", name(&output_base)));
    }
    if !errors.is_empty() {
        refuse_validation(TOOL, &errors);
    }
    let (pre_adapter_in, bait_bias_in, oxog_out) = (
        pre_adapter_in.expect("set"),
        bait_bias_in.expect("set"),
        oxog_out.expect("set"),
    );

    let pre_rows = read_beans(&pre_adapter_in);
    let bait_rows = read_beans(&bait_bias_in);

    let sample_alias = match pre_rows.first() {
        Some(row) => text(row, "SAMPLE_ALIAS"),
        None => thrown("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"),
    };
    let mut libraries: JavaHashMap<()> = JavaHashMap::new();
    let mut contexts: JavaHashMap<()> = JavaHashMap::new();
    for row in &pre_rows {
        libraries.put(&text(row, "LIBRARY"), ());
        if base(row, "REF_BASE") == 'C' {
            contexts.put(&text(row, "CONTEXT"), ());
        }
    }

    let mut pre_map: HashMap<String, HashMap<String, &HashMap<String, String>>> = HashMap::new();
    let mut bait_map: HashMap<String, HashMap<String, &HashMap<String, String>>> = HashMap::new();
    for (library, _) in libraries.iter() {
        pre_map.insert(library.to_string(), HashMap::new());
        bait_map.insert(library.to_string(), HashMap::new());
    }
    let unput =
        "java.lang.NullPointerException: Cannot invoke \"java.util.Map.put(Object, Object)\" \
                 because the return value of \"java.util.Map.get(Object)\" is null";
    for row in &pre_rows {
        if is_oxo_g(base(row, "REF_BASE"), base(row, "ALT_BASE")) {
            match pre_map.get_mut(&text(row, "LIBRARY")) {
                Some(by_context) => {
                    by_context.insert(text(row, "CONTEXT"), row);
                }
                None => thrown(unput),
            }
        }
    }
    for row in &bait_rows {
        if is_oxo_g(base(row, "REF_BASE"), base(row, "ALT_BASE")) {
            match bait_map.get_mut(&text(row, "LIBRARY")) {
                Some(by_context) => {
                    by_context.insert(text(row, "CONTEXT"), row);
                }
                None => thrown(unput),
            }
        }
    }

    let mut rows = Vec::new();
    for (library, _) in libraries.iter() {
        for (context, _) in contexts.iter() {
            let pre = pre_map[library]
                .get(&reverse_complement(context))
                .copied()
                .unwrap_or_else(|| {
                    thrown(
                        "java.lang.NullPointerException: Cannot read field \"PRO_REF_BASES\" \
                         because \"preAdapter\" is null",
                    )
                });
            let (pro_ref, pro_alt, con_ref, con_alt) = (
                long(pre, "PRO_REF_BASES"),
                long(pre, "PRO_ALT_BASES"),
                long(pre, "CON_REF_BASES"),
                long(pre, "CON_ALT_BASES"),
            );
            let total_bases = pro_ref + pro_alt + con_ref + con_alt;
            let oxidation_error_rate = (pro_alt - con_alt).max(1) as f64 / total_bases as f64;
            let bait = bait_map[library].get(context).copied().unwrap_or_else(|| {
                thrown(
                    "java.lang.NullPointerException: Cannot read field \"FWD_CXT_REF_BASES\" \
                     because \"baitBiasFwd\" is null",
                )
            });
            let c_ref_ref = long(bait, "FWD_CXT_REF_BASES");
            let g_ref_ref = long(bait, "REV_CXT_REF_BASES");
            let c_ref_alt = long(bait, "FWD_CXT_ALT_BASES");
            let g_ref_alt = long(bait, "REV_CXT_ALT_BASES");
            let c_rate = c_ref_alt as f64 / (c_ref_alt + c_ref_ref) as f64;
            let g_rate = g_ref_alt as f64 / (g_ref_alt + g_ref_ref) as f64;
            // `Math.max` answers NaN when either side is NaN; `f64::max` would not.
            let java_max = |a: f64, b: f64| if a.is_nan() { a } else { a.max(b) };
            let c_ref_oxo_error_rate = java_max(c_rate - g_rate, 1e-10);
            let g_ref_oxo_error_rate = java_max(g_rate - c_rate, 1e-10);
            rows.push(Row(CpcgMetrics {
                sample_alias: sample_alias.clone(),
                library: library.to_string(),
                context: context.to_string(),
                total_sites: 0,
                total_bases,
                ref_total_bases: pro_ref + con_ref,
                ref_nonoxo_bases: con_ref,
                ref_oxo_bases: pro_ref,
                alt_nonoxo_bases: con_alt,
                alt_oxo_bases: pro_alt,
                oxidation_error_rate,
                oxidation_q: -10.0 * oxidation_error_rate.log10(),
                c_ref_ref_bases: c_ref_ref,
                g_ref_ref_bases: g_ref_ref,
                c_ref_alt_bases: c_ref_alt,
                g_ref_alt_bases: g_ref_alt,
                c_ref_oxo_error_rate,
                g_ref_oxo_error_rate,
                c_ref_oxo_q: -10.0 * c_ref_oxo_error_rate.log10(),
                g_ref_oxo_q: -10.0 * g_ref_oxo_error_rate.log10(),
            }));
        }
    }

    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for row in &rows {
        file.add_metric(row);
    }
    if let Err(e) = std::fs::write(&oxog_out, file.write()) {
        fail(&format!("{e}"));
    }
}
