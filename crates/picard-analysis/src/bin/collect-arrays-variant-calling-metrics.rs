//! `CollectArraysVariantCallingMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.CollectArraysVariantCallingMetrics.doWork` and
//! `picard.arrays.ArraysCallingMetricAccumulator` at tag 3.4.0. The control codes and the counting
//! rules are `picard_analysis::collect_arrays_variant_calling_metrics`; this reads the header,
//! walks the records and writes the three files.
//!
//! # Two dates, two zones
//!
//! The AutoCall date is parsed as `MM/dd/yyyy HH:mm` in the JVM's zone, which the oracle leaves at
//! UTC, and the imaging date as `MM/dd/yyyy hh:mm:ss a` in America/New_York; both are written as
//! `Iso8601Date` in the JVM's zone. So the imaging date moves by four or five hours, by the US
//! daylight-saving rule of its year.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::variant::VariantContext;
use picard_analysis::fingerprint::file_uri;
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};
use picard_analysis::vcf_io::read_path;

const TOOL: &str = "CollectArraysVariantCallingMetrics";

const CONTROLS: [&str; 23] = [
    "DNP(High)",
    "DNP(Bgnd)",
    "Biotin(High)",
    "Biotin(Bgnd)",
    "Extension(A)",
    "Extension(T)",
    "Extension(C)",
    "Extension(G)",
    "TargetRemoval",
    "Hyb(High)",
    "Hyb(Medium)",
    "Hyb(Low)",
    "String(PM)",
    "String(MM)",
    "NSB(Bgnd)Red",
    "NSB(Bgnd)Purple",
    "NSB(Bgnd)Blue",
    "NSB(Bgnd)Green",
    "NP(A)",
    "NP(T)",
    "NP(C)",
    "NP(G)",
    "Restore",
];

const SUMMARY: &[&str] = &[
    "NUM_ASSAYS",
    "NUM_NON_FILTERED_ASSAYS",
    "NUM_FILTERED_ASSAYS",
    "NUM_ZEROED_OUT_ASSAYS",
    "NUM_SNPS",
    "NUM_INDELS",
    "NUM_CALLS",
    "NUM_AUTOCALL_CALLS",
    "NUM_NO_CALLS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "PCT_DBSNP",
    "CALL_RATE",
    "AUTOCALL_CALL_RATE",
    "NUM_SINGLETONS",
];

const DETAIL_OWN: &[&str] = &[
    "CHIP_WELL_BARCODE",
    "SAMPLE_ALIAS",
    "ANALYSIS_VERSION",
    "CHIP_TYPE",
    "AUTOCALL_PF",
    "AUTOCALL_DATE",
    "IMAGING_DATE",
    "IS_ZCALLED",
    "GTC_CALL_RATE",
    "AUTOCALL_GENDER",
    "FP_GENDER",
    "REPORTED_GENDER",
    "GENDER_CONCORDANCE_PF",
    "HET_PCT",
    "CLUSTER_FILE_NAME",
    "P95_GREEN",
    "P95_RED",
    "AUTOCALL_VERSION",
    "ZCALL_VERSION",
    "EXTENDED_MANIFEST_VERSION",
    "HET_HOMVAR_RATIO",
    "SCANNER_NAME",
    "PIPELINE_VERSION",
];

struct Row {
    class: &'static str,
    columns: Vec<&'static str>,
    values: Vec<Value>,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        self.class
    }
    fn columns(&self) -> &[&'static str] {
        &self.columns
    }
    fn values(&self) -> Vec<Value> {
        self.values.clone()
    }
}

/// The summary counters, which the detail metric extends.
#[derive(Default, Clone)]
struct Counts {
    assays: i64,
    non_filtered: i64,
    filtered: i64,
    zeroed: i64,
    snps: i64,
    indels: i64,
    calls: i64,
    autocall_calls: i64,
    no_calls: i64,
    in_dbsnp: i64,
    singletons: i64,
    hets: i64,
    hom_vars: i64,
}

impl Counts {
    fn add(&mut self, o: &Counts) {
        self.assays += o.assays;
        self.non_filtered += o.non_filtered;
        self.filtered += o.filtered;
        self.zeroed += o.zeroed;
        self.snps += o.snps;
        self.indels += o.indels;
        self.calls += o.calls;
        self.autocall_calls += o.autocall_calls;
        self.no_calls += o.no_calls;
        self.in_dbsnp += o.in_dbsnp;
        self.singletons += o.singletons;
        self.hets += o.hets;
        self.hom_vars += o.hom_vars;
    }

    /// The summary columns, with `calculateDerivedFields`' float divisions.
    fn values(&self) -> Vec<Value> {
        let rate = |n: i64, d: i64| Value::Double(f64::from(n as f32 / d as f32));
        vec![
            Value::Long(self.assays),
            Value::Long(self.non_filtered),
            Value::Long(self.filtered),
            Value::Long(self.zeroed),
            Value::Long(self.snps),
            Value::Long(self.indels),
            Value::Long(self.calls),
            Value::Long(self.autocall_calls),
            Value::Long(self.no_calls),
            Value::Long(self.in_dbsnp),
            Value::Long(self.snps - self.in_dbsnp),
            rate(self.in_dbsnp, self.snps),
            rate(self.calls, self.non_filtered),
            rate(self.autocall_calls, self.non_filtered),
            Value::Long(self.singletons),
        ]
    }
}

fn header_value(lines: &[HeaderLine], key: &str) -> Option<String> {
    lines.iter().find_map(|l| match l {
        HeaderLine::Unstructured { key: k, value } if k == key => Some(value.clone()),
        _ => None,
    })
}

fn required(lines: &[HeaderLine], key: &str) -> String {
    header_value(lines, key).unwrap_or_else(|| {
        thrown(&format!(
            "java.lang.IllegalArgumentException: Input VCF file is missing header line of type '{key}'"
        ))
    })
}

/// `Sex.fromString(s).toSymbol()`.
fn sex_symbol(s: &str) -> String {
    let symbol = match s.to_ascii_lowercase().as_str() {
        "m" | "male" => "M",
        "f" | "female" => "F",
        "u" | "unknown" => "U",
        "n" | "notreported" => "N",
        _ => thrown(&format!(
            "picard.PicardException: Unrecognized Sex string: {s}"
        )),
    };
    symbol.to_string()
}

/// Days from 1970-01-01 to a civil date.
fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = (m + 9) % 12;
    let doy = (153 * mp + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// The day of the month of the `n`th Sunday (`n` = 0 for the last) of a month.
fn sunday(y: i64, m: i64, n: i64) -> i64 {
    // 1970-01-01 was a Thursday: day 0 has weekday 4 (Sunday = 0).
    let weekday = |d: i64| (days_from_civil(y, m, d) + 4).rem_euclid(7);
    if n > 0 {
        let first = (7 - weekday(1)) % 7 + 1;
        first + 7 * (n - 1)
    } else {
        let days = days_from_civil(
            if m == 12 { y + 1 } else { y },
            if m == 12 { 1 } else { m + 1 },
            1,
        ) - days_from_civil(y, m, 1);
        days - weekday(days)
    }
}

/// America/New_York's offset in hours at a local wall-clock time, by the US rule of the year.
fn new_york_offset(y: i64, mo: i64, d: i64, h: i64) -> i64 {
    let (start, end) = if y >= 2007 {
        ((3, sunday(y, 3, 2)), (11, sunday(y, 11, 1)))
    } else {
        ((4, sunday(y, 4, 1)), (10, sunday(y, 10, 0)))
    };
    let at = (mo, d, h);
    let dst = at >= (start.0, start.1, 2) && at < (end.0, end.1, 2);
    if dst {
        -4
    } else {
        -5
    }
}

/// Seconds since the epoch to `yyyy-MM-dd'T'HH:mm:ss+0000`.
fn iso(seconds: i64) -> String {
    let days = seconds.div_euclid(86_400);
    let secs = seconds.rem_euclid(86_400);
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}+0000",
        secs / 3600,
        secs % 3600 / 60,
        secs % 60
    )
}

fn number(s: &str) -> Option<i64> {
    s.trim().parse().ok()
}

/// `MM/dd/yyyy HH:mm` (AutoCall) or `MM/dd/yyyy hh:mm:ss a` (imaging), to the ISO text.
fn parse_date(value: &str, key: &str, imaging: bool) -> String {
    let fail = || -> ! {
        thrown(&format!(
            "java.lang.IllegalArgumentException: Unrecognized date for '{key}' in VCF header ({value})"
        ))
    };
    let (date, rest) = value.split_once(' ').unwrap_or_else(|| fail());
    let parts: Vec<i64> = date
        .split('/')
        .map(|p| number(p).unwrap_or_else(|| fail()))
        .collect();
    if parts.len() != 3 {
        fail();
    }
    let (mo, d, y) = (parts[0], parts[1], parts[2]);
    let (time, ampm) = match rest.split_once(' ') {
        Some((t, a)) => (t, Some(a)),
        None => (rest, None),
    };
    let t: Vec<i64> = time
        .split(':')
        .map(|p| number(p).unwrap_or_else(|| fail()))
        .collect();
    let (mut h, mi, s) = (t[0], *t.get(1).unwrap_or(&0), *t.get(2).unwrap_or(&0));
    if imaging {
        match ampm {
            Some(a) if a.eq_ignore_ascii_case("PM") => h = h % 12 + 12,
            Some(a) if a.eq_ignore_ascii_case("AM") => h %= 12,
            _ => fail(),
        }
    }
    let mut seconds = days_from_civil(y, mo, d) * 86_400 + h * 3600 + mi * 60 + s;
    if imaging {
        seconds -= new_york_offset(y, mo, d, h) * 3600;
    }
    iso(seconds)
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("SD", "SEQUENCE_DICTIONARY"),
    ]);
    let input = args.required("INPUT");
    let dbsnp_path = args.required("DBSNP");
    let output = args.required("OUTPUT");
    let threshold = args.double("CALL_RATE_PF_THRESHOLD", 0.98);
    if threshold <= 0.0 || threshold > 1.0 {
        refuse_validation(
            TOOL,
            &["The parameter CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0".to_string()],
        );
    }
    if !std::path::Path::new(&format!("{input}.idx")).exists()
        && !std::path::Path::new(&format!("{input}.tbi")).exists()
    {
        thrown(&format!(
            "htsjdk.tribble.TribbleException: An index is required, but none found with file ending .idx, for input source: {}",
            file_uri(&input)
        ));
    }
    let vcf = read_path(&input).unwrap_or_else(|e| thrown(&e));
    let lines = &vcf.file.header.lines;
    let dbsnp = read_path(&dbsnp_path).unwrap_or_else(|e| thrown(&e));
    let known: Vec<(String, i64)> = dbsnp
        .records
        .iter()
        .map(|r| &r.variant)
        .filter(|v| v.alleles.iter().all(|a| a.len() == 1))
        .map(|v| (v.contig.clone(), v.start))
        .collect();

    // The control codes, each `control|category|red|green`.
    let mut controls = MetricsFile::new();
    controls.add_header(&format!("{TOOL} <command line>"));
    controls.add_header("Started on: <timestamp>");
    for control in CONTROLS {
        let value = required(lines, control);
        let t: Vec<&str> = value.split('|').collect();
        controls.add_metric(&Row {
            class:
                "picard.arrays.CollectArraysVariantCallingMetrics$ArraysControlCodesSummaryMetrics",
            columns: vec!["CONTROL", "CATEGORY", "RED", "GREEN"],
            values: vec![
                Value::Str(t[0].to_string()),
                Value::Str(t.get(1).unwrap_or(&"").to_string()),
                Value::Long(t.get(2).and_then(|v| number(v)).unwrap_or(0)),
                Value::Long(t.get(3).and_then(|v| number(v)).unwrap_or(0)),
            ],
        });
    }

    // `ArraysCallingMetricAccumulator.setup`.
    let sample_alias = required(lines, "sampleAlias");
    let pipeline_version = header_value(lines, "pipelineVersion");
    let analysis_version = header_value(lines, "analysisVersionNumber").and_then(|v| number(&v));
    let chip_type = required(lines, "arrayType");
    let gender_or_not = |key: &str| {
        header_value(lines, key)
            .map(|g| sex_symbol(&g))
            .unwrap_or_else(|| "N".to_string())
    };
    let reported = gender_or_not("expectedGender");
    let fingerprint = gender_or_not("fingerprintGender");
    let gtc_call_rate =
        header_value(lines, "gtcCallRate").and_then(|v| v.trim().parse::<f64>().ok());
    let autocall_gender = sex_symbol(&required(lines, "autocallGender"));
    let autocall_version = required(lines, "autocallVersion");
    let autocall_date = parse_date(&required(lines, "autocallDate"), "autocallDate", false);
    let imaging_date = parse_date(&required(lines, "imagingDate"), "imagingDate", true);
    let manifest_version = required(lines, "extendedIlluminaManifestVersion");
    let zcall_version = header_value(lines, "zcallVersion");
    let zcall_thresholds = header_value(lines, "zcallThresholds");
    let cluster = required(lines, "clusterFile");
    let p95_green = number(&required(lines, "p95Green")).unwrap_or(0);
    let p95_red = number(&required(lines, "p95Red")).unwrap_or(0);
    let scanner = required(lines, "scannerName");

    // `accumulate`, per sample.
    let samples = vcf.file.header.samples.clone();
    let mut per_sample: Vec<Counts> = vec![Counts::default(); samples.len()];
    for record in &vcf.records {
        let vc: &VariantContext = &record.variant;
        let variants: Vec<&htsjdk_vcf::variant::Genotype> = vc
            .genotypes
            .iter()
            .filter(|g| g.is_het() || g.is_hom_var())
            .take(2)
            .collect();
        let chromosomes: i32 = variants
            .iter()
            .map(|g| if g.is_het() { 1 } else { 2 })
            .sum();
        let singleton = (chromosomes == 1)
            .then(|| variants.last().map(|g| g.sample_name.clone()))
            .flatten();
        let filters = vc.filters.clone().unwrap_or_default();
        let is_filtered = !filters.is_empty();
        let is_snp =
            vc.alleles.len() > 1 && vc.alleles.iter().all(|a| a.len() == 1 && !a.is_symbolic());
        let is_indel = !is_snp
            && vc.alleles.len() > 1
            && vc.alleles.iter().all(|a| !a.is_symbolic())
            && vc.alleles[1..]
                .iter()
                .any(|a| a.len() != vc.alleles[0].len());
        for (n, sample) in samples.iter().enumerate() {
            let Some(g) = vc.genotypes.iter().find(|g| &g.sample_name == sample) else {
                continue;
            };
            let m = &mut per_sample[n];
            m.assays += 1;
            if !is_filtered || filters.iter().any(|f| f == "DUPE") {
                m.non_filtered += 1;
                if g.is_called() {
                    m.calls += 1;
                    // `getExtendedAttribute("GTA", genotypeString)`: absent, the call itself.
                    let gta = match g.get("GTA") {
                        Some(htsjdk_vcf::variant::Value::Str(s)) => Some(s.clone()),
                        _ => None,
                    };
                    if gta.as_deref() != Some("./.") {
                        m.autocall_calls += 1;
                    }
                } else {
                    m.no_calls += 1;
                }
                if is_snp {
                    m.snps += 1;
                    if known.iter().any(|(c, p)| *c == vc.contig && *p == vc.start) {
                        m.in_dbsnp += 1;
                    }
                } else if is_indel {
                    m.indels += 1;
                }
                if singleton.as_deref() == Some(sample.as_str()) {
                    m.singletons += 1;
                }
                if g.is_het() {
                    m.hets += 1;
                } else if g.is_hom_var() {
                    m.hom_vars += 1;
                }
            } else {
                m.filtered += 1;
                if filters.iter().any(|f| f == "ZEROED_OUT_ASSAY") {
                    m.zeroed += 1;
                }
            }
        }
    }

    let header = |file: &mut MetricsFile| {
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
    };
    let mut detail = MetricsFile::new();
    header(&mut detail);
    let mut summary_counts = Counts::default();
    let text = |v: &Option<String>| v.clone().map_or(Value::Null, Value::Str);
    for (sample, c) in samples.iter().zip(&per_sample) {
        summary_counts.add(c);
        let autocall_rate = c.autocall_calls as f32 / c.non_filtered as f32;
        let concordance = {
            let all = [&reported, &fingerprint, &autocall_gender];
            let count = |s: &str| all.iter().filter(|x| x.as_str() == s).count();
            if count("U") + count("N") == 3 {
                false
            } else {
                (count("F") > 1 && count("M") == 0) || (count("M") > 1 && count("F") == 0)
            }
        };
        let mut values = vec![
            Value::Str(sample.clone()),
            Value::Str(sample_alias.clone()),
            analysis_version.map_or(Value::Null, Value::Long),
            Value::Str(chip_type.clone()),
            Value::Bool(f64::from(autocall_rate) > threshold),
            Value::Str(autocall_date.clone()),
            Value::Str(imaging_date.clone()),
            Value::Bool(zcall_thresholds.as_deref().is_some_and(|z| !z.is_empty())),
            gtc_call_rate.map_or(Value::Null, Value::Double),
            Value::Str(autocall_gender.clone()),
            Value::Str(fingerprint.clone()),
            Value::Str(reported.clone()),
            Value::Bool(concordance),
            Value::Double(c.hets as f64 / c.calls as f64),
            Value::Str(cluster.clone()),
            Value::Long(p95_green),
            Value::Long(p95_red),
            Value::Str(autocall_version.clone()),
            text(&zcall_version),
            Value::Str(manifest_version.clone()),
            Value::Double(c.hets as f64 / c.hom_vars as f64),
            Value::Str(scanner.clone()),
            text(&pipeline_version),
        ];
        values.extend(c.values());
        let mut columns: Vec<&'static str> = DETAIL_OWN.to_vec();
        columns.extend_from_slice(SUMMARY);
        detail.add_metric(&Row {
            class: "picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingDetailMetrics",
            columns,
            values,
        });
    }
    let mut summary = MetricsFile::new();
    header(&mut summary);
    summary.add_metric(&Row {
        class:
            "picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingSummaryMetrics",
        columns: SUMMARY.to_vec(),
        values: summary_counts.values(),
    });
    let prefix = picard_analysis::fingerprint::reference_path(&output);
    let prefix = std::path::absolute(&output)
        .map(|p| p.display().to_string())
        .unwrap_or(prefix);
    for (ext, file) in [
        ("arrays_variant_calling_detail_metrics", &detail),
        ("arrays_variant_calling_summary_metrics", &summary),
        ("arrays_control_code_summary_metrics", &controls),
    ] {
        if let Err(e) = std::fs::write(format!("{prefix}.{ext}"), file.write()) {
            thrown(&format!("htsjdk.samtools.SAMException: {e}"));
        }
    }
}
