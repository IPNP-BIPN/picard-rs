//! `CrosscheckFingerprints` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.fingerprint.CrosscheckFingerprints.doWork` at tag 3.4.0:
//!
//! * the mutex pairs (`SECOND_INPUT`/`MATRIX_OUTPUT`, `INPUT_SAMPLE_MAP`/`INPUT_SAMPLE_FILE_MAP`)
//!   and `customCommandLineValidation`;
//! * every file fingerprinted in order (`fingerprintFiles` with one thread) into a
//!   `ConcurrentHashMap`, each fingerprint then capped at `-MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK`
//!   (`CappedHaplotypeProbabilities`) into a `HashMap`, and the sample maps applied;
//! * `crossCheckGrouped`: both sides merged by `CROSSCHECK_BY` and EVERY pair compared, a group
//!   with itself included, in the `HashMap` order of the merged identities; or, with
//!   `SECOND_INPUT` and `CHECK_SAME_SAMPLE`, `checkFingerprintsBySample`, which forces `SAMPLE`;
//! * the verdict from `LOD_THRESHOLD` used with both signs, `OUTPUT_ERRORS_ONLY` keeping only the
//!   unexpected and inconclusive rows, and the exit code: `EXIT_CODE_WHEN_NO_VALID_CHECKS` when no
//!   row kept has a non-zero LOD (before anything is written), else `EXIT_CODE_WHEN_MISMATCH` when
//!   a comparison was unexpected;
//! * `MATRIX_OUTPUT` through `NumberFormat.getInstance()` with four fraction digits.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprinting::{
    absolute_path, calculate_match_results, fingerprint_files, log, merge_entries_by, DataType,
    Fingerprint, FingerprintMap, HaplotypeMap, IdDetails, JavaMap, Probabilities, SamOptions,
    Stringency,
};
use picard_analysis::metrics_cli::{fail, thrown, Args};

const TOOL: &str = "CrosscheckFingerprints";

/// The row keys, the column keys and the LODs of `MATRIX_OUTPUT`.
type Matrix = (Vec<String>, Vec<String>, Vec<Vec<f64>>);

#[derive(Clone, Copy, PartialEq, Eq)]
enum Verdict {
    ExpectedMatch,
    ExpectedMismatch,
    UnexpectedMatch,
    UnexpectedMismatch,
    Inconclusive,
}

impl Verdict {
    fn name(self) -> &'static str {
        match self {
            Verdict::ExpectedMatch => "EXPECTED_MATCH",
            Verdict::ExpectedMismatch => "EXPECTED_MISMATCH",
            Verdict::UnexpectedMatch => "UNEXPECTED_MATCH",
            Verdict::UnexpectedMismatch => "UNEXPECTED_MISMATCH",
            Verdict::Inconclusive => "INCONCLUSIVE",
        }
    }

    fn is_expected(self) -> Option<bool> {
        match self {
            Verdict::ExpectedMatch | Verdict::ExpectedMismatch => Some(true),
            Verdict::UnexpectedMatch | Verdict::UnexpectedMismatch => Some(false),
            Verdict::Inconclusive => None,
        }
    }
}

/// `getMatchResults`.
fn verdict(expected_to_match: bool, lod: f64, threshold: f64) -> Verdict {
    if expected_to_match {
        if lod < threshold {
            Verdict::UnexpectedMismatch
        } else if lod > -threshold {
            Verdict::ExpectedMatch
        } else {
            Verdict::Inconclusive
        }
    } else if lod > -threshold {
        Verdict::UnexpectedMatch
    } else if lod < threshold {
        Verdict::ExpectedMismatch
    } else {
        Verdict::Inconclusive
    }
}

struct Metric {
    left: IdDetails,
    right: IdDetails,
    verdict: Verdict,
    data_type: DataType,
    lod: f64,
    lod_tn: f64,
    lod_nt: f64,
}

fn opt(value: &Option<String>) -> Value {
    match value {
        Some(s) => Value::Str(s.clone()),
        None => Value::Null,
    }
}

fn lane(value: Option<i32>) -> Value {
    match value {
        Some(v) => Value::Long(i64::from(v)),
        None => Value::Null,
    }
}

impl MetricBean for Metric {
    fn class_name(&self) -> &str {
        "picard.fingerprint.CrosscheckMetric"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "LEFT_GROUP_VALUE",
            "RIGHT_GROUP_VALUE",
            "RESULT",
            "DATA_TYPE",
            "LOD_SCORE",
            "LOD_SCORE_TUMOR_NORMAL",
            "LOD_SCORE_NORMAL_TUMOR",
            "LEFT_RUN_BARCODE",
            "LEFT_LANE",
            "LEFT_MOLECULAR_BARCODE_SEQUENCE",
            "LEFT_LIBRARY",
            "LEFT_SAMPLE",
            "LEFT_FILE",
            "RIGHT_RUN_BARCODE",
            "RIGHT_LANE",
            "RIGHT_MOLECULAR_BARCODE_SEQUENCE",
            "RIGHT_LIBRARY",
            "RIGHT_SAMPLE",
            "RIGHT_FILE",
        ]
    }
    fn values(&self) -> Vec<Value> {
        let l = &self.left;
        let r = &self.right;
        vec![
            opt(&l.group),
            opt(&r.group),
            Value::Str(self.verdict.name().into()),
            Value::Str(self.data_type.name().into()),
            Value::Double(self.lod),
            Value::Double(self.lod_tn),
            Value::Double(self.lod_nt),
            opt(&l.run_barcode),
            lane(l.run_lane),
            opt(&l.molecular_barcode),
            opt(&l.library),
            opt(&l.sample),
            opt(&l.file),
            opt(&r.run_barcode),
            lane(r.run_lane),
            opt(&r.molecular_barcode),
            opt(&r.library),
            opt(&r.sample),
            opt(&r.file),
        ]
    }
}

fn throw(t: (String, String)) -> ! {
    thrown(&format!("{}: {}", t.0, t.1))
}

/// The usage dump Barclay prints before a command-line refusal, and the refusal.
fn refuse(message: &str) -> ! {
    eprintln!("USAGE: {TOOL} [arguments]\n");
    eprintln!("{message}");
    std::process::exit(1);
}

/// `NumberFormat.getInstance()` (en_US) with at most four fraction digits: grouping, half-even.
fn number_format(value: f64) -> String {
    if value.is_nan() {
        return "NaN".into();
    }
    if value.is_infinite() {
        return if value > 0.0 {
            "∞".into()
        } else {
            "-∞".into()
        };
    }
    let text = format!("{:.4}", value);
    let (int_part, frac) = text.split_once('.').unwrap_or((&text, ""));
    let negative = int_part.starts_with('-');
    let digits = int_part.trim_start_matches('-');
    let mut grouped = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(c);
    }
    let frac = frac.trim_end_matches('0');
    let mut out = String::new();
    if negative {
        out.push('-');
    }
    out.push_str(&grouped);
    if !frac.is_empty() {
        out.push('.');
        out.push_str(frac);
    }
    out
}

/// A two-column tab-separated map, as `getStringStringMap` reads it.
fn read_string_map(path: &str, field: &str) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&e.to_string()));
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            thrown(&format!(
                "java.lang.IllegalArgumentException: Each line of the {field} must have exactly two strings separated by a tab."
            ));
        }
        if out.iter().any(|(k, _)| k == parts[0]) {
            thrown(&format!(
                "java.lang.IllegalArgumentException: Strings in first column of the {field} must be unique. found [{}] twice.",
                parts[0]
            ));
        }
        out.push((parts[0].to_string(), parts[1].to_string()));
    }
    out
}

/// `remapFingerprints`: every identity whose sample the map names is re-keyed under the new
/// sample, removed and put back in the `LinkedHashSet` order of the keys as they were.
fn remap(fp_map: &mut FingerprintMap, sample_map: &[(String, String)], field: &str) {
    let samples: Vec<String> = fp_map
        .iter()
        .filter_map(|(id, _)| id.sample.clone())
        .collect();
    let mut resulting: Vec<String> = Vec::new();
    for s in &samples {
        if !resulting.contains(s) {
            resulting.push(s.clone());
        }
    }
    for (from, to) in sample_map {
        if let Some(pos) = resulting.iter().position(|s| s == from) {
            resulting.remove(pos);
            resulting.push(to.clone());
        }
    }
    let mut unique: Vec<&String> = Vec::new();
    for s in &resulting {
        if !unique.contains(&s) {
            unique.push(s);
        }
    }
    if unique.len() != resulting.len() {
        thrown(&format!(
            "java.lang.IllegalArgumentException: After applying the mapping found in the {field} the resulting sample names must be unique when taken together with the remaining unmapped samples."
        ));
    }
    let ids: Vec<IdDetails> = fp_map.iter().map(|(id, _)| id.clone()).collect();
    for id in ids {
        let Some(sample) = &id.sample else { continue };
        let Some((_, to)) = sample_map.iter().find(|(k, _)| k == sample) else {
            continue;
        };
        let fp = fp_map.remove(id.hash(), &id).expect("listed above");
        let mut renamed = id.clone();
        renamed.sample = Some(to.clone());
        fp_map.put(renamed.hash(), renamed, fp);
    }
}

#[allow(clippy::too_many_arguments)]
fn cross_check_grouped(
    lhs: &FingerprintMap,
    rhs: &FingerprintMap,
    by: DataType,
    map: &HaplotypeMap,
    opts: &Options,
    individuals: &Option<Vec<(String, String)>>,
    metrics: &mut Vec<Metric>,
    matrix: &mut Option<Matrix>,
) -> i64 {
    let to_entries = |m: &FingerprintMap| -> Vec<(IdDetails, Fingerprint)> {
        m.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
    };
    let lhs_by = merge_entries_by(&to_entries(lhs), by, map).unwrap_or_else(|t| throw(t));
    let rhs_by = merge_entries_by(&to_entries(rhs), by, map).unwrap_or_else(|t| throw(t));
    let lhs_ids: Vec<(&IdDetails, &Fingerprint)> = lhs_by.iter().collect();
    let rhs_ids: Vec<(&IdDetails, &Fingerprint)> = rhs_by.iter().collect();
    if opts.matrix_output {
        *matrix = Some((
            lhs_ids
                .iter()
                .map(|(k, _)| k.group.clone().unwrap_or_default())
                .collect(),
            rhs_ids
                .iter()
                .map(|(k, _)| k.group.clone().unwrap_or_default())
                .collect(),
            vec![vec![0.0; rhs_ids.len()]; lhs_ids.len()],
        ));
    }
    let resolve = |sample: &Option<String>| -> Option<String> {
        if let Some(ind) = individuals {
            if let Some(s) = sample {
                if let Some((_, v)) = ind.iter().find(|(k, _)| k == s) {
                    return Some(v.clone());
                }
            }
        }
        sample.clone()
    };
    let mut unexpected = 0;
    for (row, (lhs_id, lhs_fp)) in lhs_ids.iter().enumerate() {
        for (col, (rhs_id, rhs_fp)) in rhs_ids.iter().enumerate() {
            let l = resolve(&lhs_id.sample);
            let r = resolve(&rhs_id.sample);
            let Some(l) = l else {
                thrown("java.lang.NullPointerException");
            };
            let expected = opts.expect_all || Some(&l) == r.as_ref();
            let results = calculate_match_results(
                lhs_fp,
                rhs_fp,
                map,
                opts.loss_of_het,
                false,
                opts.tumor_aware,
            )
            .unwrap_or_else(|t| throw(t));
            let v = verdict(expected, results.lod, opts.lod_threshold);
            if !opts.errors_only || v == Verdict::Inconclusive || v.is_expected() == Some(false) {
                metrics.push(Metric {
                    left: (*lhs_id).clone(),
                    right: (*rhs_id).clone(),
                    verdict: v,
                    data_type: by,
                    lod: results.lod,
                    lod_tn: results.lod_tn,
                    lod_nt: results.lod_nt,
                });
            }
            if v != Verdict::Inconclusive && v.is_expected() == Some(false) {
                unexpected += 1;
            }
            if let Some((_, _, cells)) = matrix.as_mut() {
                cells[row][col] = results.lod;
            }
        }
    }
    unexpected
}

fn check_by_sample(
    first: &FingerprintMap,
    second: &FingerprintMap,
    map: &HaplotypeMap,
    opts: &Options,
    metrics: &mut Vec<Metric>,
) -> i64 {
    let to_entries = |m: &FingerprintMap| -> Vec<(IdDetails, Fingerprint)> {
        m.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
    };
    let one =
        merge_entries_by(&to_entries(first), DataType::Sample, map).unwrap_or_else(|t| throw(t));
    let two =
        merge_entries_by(&to_entries(second), DataType::Sample, map).unwrap_or_else(|t| throw(t));
    // `Collectors.toMap(id -> id.group, id -> id)`: a HashMap keyed by the group string.
    let by_group = |m: &FingerprintMap| {
        let mut out: picard_analysis::java_hash_map::JavaHashMap<IdDetails> =
            picard_analysis::java_hash_map::JavaHashMap::new();
        for (id, _) in m.iter() {
            out.put(id.group.as_deref().unwrap_or("null"), id.clone());
        }
        out
    };
    let s1 = by_group(&one);
    let s2 = by_group(&two);
    let mut samples: Vec<String> = Vec::new();
    for (k, _) in s1.iter().chain(s2.iter()) {
        if !samples.iter().any(|s| s == k) {
            samples.push(k.to_string());
        }
    }
    let mut unexpected = 0;
    for sample in samples {
        let (Some(lhs_id), Some(rhs_id)) = (s1.get(&sample), s2.get(&sample)) else {
            let side = if s1.get(&sample).is_none() {
                "LEFT"
            } else {
                "RIGHT"
            };
            log(
                "ERROR",
                TOOL,
                &format!("sample {sample} is missing from {side} group"),
            );
            unexpected += 1;
            continue;
        };
        let lhs_fp = one.get(lhs_id.hash(), lhs_id).expect("keyed above");
        let rhs_fp = two.get(rhs_id.hash(), rhs_id).expect("keyed above");
        if lhs_fp.map.is_empty() || rhs_fp.map.is_empty() {
            unexpected += 1;
        }
        let results = calculate_match_results(
            lhs_fp,
            rhs_fp,
            map,
            opts.loss_of_het,
            false,
            opts.tumor_aware,
        )
        .unwrap_or_else(|t| throw(t));
        let v = verdict(true, results.lod, opts.lod_threshold);
        let keep = if !opts.errors_only {
            true
        } else {
            match v.is_expected() {
                Some(e) => !e,
                // `!result.isExpected()` unboxes the null an inconclusive verdict carries.
                None => thrown(
                    "java.lang.NullPointerException: Cannot invoke \"java.lang.Boolean.booleanValue()\" because the return value of \"picard.fingerprint.CrosscheckMetric$FingerprintResult.isExpected()\" is null",
                ),
            }
        };
        if keep {
            metrics.push(Metric {
                left: lhs_id.clone(),
                right: rhs_id.clone(),
                verdict: v,
                data_type: DataType::Sample,
                lod: results.lod,
                lod_tn: results.lod_tn,
                lod_nt: results.lod_nt,
            });
        }
        if v != Verdict::Inconclusive && v.is_expected() == Some(false) {
            unexpected += 1;
        }
        if results.lod == 0.0 {
            unexpected += 1;
        }
    }
    unexpected
}

struct Options {
    lod_threshold: f64,
    expect_all: bool,
    errors_only: bool,
    loss_of_het: f64,
    tumor_aware: bool,
    matrix_output: bool,
}

/// `capFingerprints`: every block's evidence capped, into a `HashMap` built in the source's order.
fn cap(
    entries: &[(IdDetails, Fingerprint)],
    map: &HaplotypeMap,
    max_effect: f64,
) -> FingerprintMap {
    let mut out: FingerprintMap = JavaMap::new();
    for (id, fp) in entries {
        let mut capped = Fingerprint::new(fp.sample.clone(), fp.source.clone(), fp.info.clone());
        for p in fp.map.values() {
            capped.add(map, Probabilities::capped(p, map, -max_effect));
        }
        out.put(id.hash(), id.clone(), capped);
    }
    out
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("SI", "SECOND_INPUT"),
        ("O", "OUTPUT"),
        ("MO", "MATRIX_OUTPUT"),
        ("H", "HAPLOTYPE_MAP"),
        ("LOD", "LOD_THRESHOLD"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let inputs = args.all("INPUT");
    let second = args.all("SECOND_INPUT");
    let output = args.get("OUTPUT").map(str::to_string);
    let matrix_output = args.get("MATRIX_OUTPUT").map(str::to_string);
    let input_sample_map = args.get("INPUT_SAMPLE_MAP").map(str::to_string);
    let input_sample_file_map = args.get("INPUT_SAMPLE_FILE_MAP").map(str::to_string);
    let second_sample_map = args.get("SECOND_INPUT_SAMPLE_MAP").map(str::to_string);
    let individual_map = args.get("SAMPLE_INDIVIDUAL_MAP").map(str::to_string);

    // Barclay's mutex check, then the required arguments.
    if !second.is_empty() && matrix_output.is_some() {
        refuse(
            "Argument 'SECOND_INPUT' cannot be used in conjunction with argument(s) MATRIX_OUTPUT",
        );
    }
    if input_sample_map.is_some() && input_sample_file_map.is_some() {
        refuse("Argument 'INPUT_SAMPLE_MAP' cannot be used in conjunction with argument(s) INPUT_SAMPLE_FILE_MAP");
    }
    if inputs.is_empty() {
        refuse("Argument 'INPUT' is required");
    }
    let haplotype_map = args.required("HAPLOTYPE_MAP");
    let mode_same_sample = match args.get("CROSSCHECK_MODE").unwrap_or("CHECK_SAME_SAMPLE") {
        "CHECK_SAME_SAMPLE" => true,
        "CHECK_ALL_OTHERS" => false,
        other => fail(&format!(
            "Argument 'CROSSCHECK_MODE' cannot be set to '{other}'"
        )),
    };
    let mut crosscheck_by = DataType::parse(args.get("CROSSCHECK_BY").unwrap_or("READGROUP"))
        .unwrap_or_else(|| fail("Argument 'CROSSCHECK_BY' has an invalid value"));
    let stringency = match args.get("VALIDATION_STRINGENCY").unwrap_or("STRICT") {
        "STRICT" => Stringency::Strict,
        "LENIENT" => Stringency::Lenient,
        _ => Stringency::Silent,
    };
    let opts = Options {
        lod_threshold: args.double("LOD_THRESHOLD", 0.0),
        expect_all: args.bool("EXPECT_ALL_GROUPS_TO_MATCH", false),
        errors_only: args.bool("OUTPUT_ERRORS_ONLY", false),
        loss_of_het: args.double("LOSS_OF_HET_RATE", 0.5),
        tumor_aware: args.bool("CALCULATE_TUMOR_AWARE_RESULTS", true),
        matrix_output: matrix_output.is_some(),
    };
    let require_index = args.bool("REQUIRE_INDEX_FILES", false);
    let allow_duplicates = args.bool("ALLOW_DUPLICATE_READS", false);
    let exit_mismatch = args.int("EXIT_CODE_WHEN_MISMATCH", 1) as i32;
    let exit_no_valid = args.int("EXIT_CODE_WHEN_NO_VALID_CHECKS", 1) as i32;
    let max_effect = args.double("MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK", 3.0);
    if max_effect < 0.0 {
        refuse(&format!(
            "Argument 'MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK' has value {max_effect}, which is less than the minimum allowed value of 0.0"
        ));
    }

    // customCommandLineValidation. Its two sample-map checks test `SECOND_INPUT == null`, which a
    // collection argument never is (Barclay hands it an empty list), so they never refuse: a sample
    // map without SECOND_INPUT is applied (INPUT's) or ignored (SECOND_INPUT's).
    if args.get("REFERENCE_SEQUENCE").is_none()
        && inputs
            .iter()
            .chain(second.iter())
            .any(|i| i.ends_with(".cram"))
    {
        refuse("REFERENCE must be provided when using CRAM as input.");
    }

    // doWork.
    if !second.is_empty() && mode_same_sample && crosscheck_by != DataType::Sample {
        crosscheck_by = DataType::Sample;
    }
    let map_text = std::fs::read_to_string(&haplotype_map).unwrap_or_else(|e| fail(&e.to_string()));
    let map = HaplotypeMap::from_database(&map_text, &absolute_path(&haplotype_map))
        .unwrap_or_else(|t| throw(t));
    let options = SamOptions {
        allow_duplicate_reads: allow_duplicates,
        stringency,
        ..SamOptions::default()
    };

    let uncapped =
        fingerprint_files(&inputs, &map, &options, require_index).unwrap_or_else(|t| throw(t));
    let mut fp_map = cap(&uncapped, &map, max_effect);
    if let Some(path) = &input_sample_map {
        let m = read_string_map(path, "INPUT_SAMPLE_MAP");
        remap(&mut fp_map, &m, "INPUT_SAMPLE_MAP");
    }
    if let Some(path) = &input_sample_file_map {
        // Column one is the new sample, column two the file it applies to, matched by URI.
        let m: Vec<(String, String)> = read_string_map(path, "INPUT_SAMPLE_FILE_MAP")
            .into_iter()
            .map(|(sample, file)| (picard_analysis::fingerprinting::uri_of(&file), sample))
            .collect();
        // Every file fingerprinted, mapped or not, must hold one sample.
        let mut files: Vec<(Option<String>, Vec<Option<String>>)> = Vec::new();
        for (id, _) in fp_map.iter() {
            match files.iter_mut().find(|(f, _)| *f == id.file) {
                Some((_, samples)) => {
                    if !samples.contains(&id.sample) {
                        samples.push(id.sample.clone());
                    }
                }
                None => files.push((id.file.clone(), vec![id.sample.clone()])),
            }
        }
        if let Some((file, samples)) = files.iter().find(|(_, s)| s.len() > 1) {
            thrown(&format!(
                "java.lang.IllegalArgumentException: fingerprinting file ({}in INPUT_SAMPLE_FILE_MAP contains multiple samples: {}",
                file.as_deref().unwrap_or("null"),
                samples
                    .iter()
                    .map(|s| s.as_deref().unwrap_or("null"))
                    .collect::<String>()
            ));
        }
        let ids: Vec<IdDetails> = fp_map.iter().map(|(id, _)| id.clone()).collect();
        for id in ids {
            let Some(file) = &id.file else { continue };
            let Some((_, sample)) = m.iter().find(|(f, _)| f == file) else {
                continue;
            };
            let fp = fp_map.remove(id.hash(), &id).expect("listed above");
            let mut renamed = id.clone();
            renamed.sample = Some(sample.clone());
            fp_map.put(renamed.hash(), renamed, fp);
        }
    }
    let individuals = individual_map
        .as_ref()
        .map(|p| read_string_map(p, "SAMPLE_INDIVIDUAL_MAP"));

    let mut metrics: Vec<Metric> = Vec::new();
    let mut matrix = None;
    let unexpected = if second.is_empty() {
        cross_check_grouped(
            &fp_map,
            &fp_map,
            crosscheck_by,
            &map,
            &opts,
            &individuals,
            &mut metrics,
            &mut matrix,
        )
    } else {
        let uncapped2 =
            fingerprint_files(&second, &map, &options, require_index).unwrap_or_else(|t| throw(t));
        let mut fp_map2 = cap(&uncapped2, &map, max_effect);
        if let Some(path) = &second_sample_map {
            let m = read_string_map(path, "SECOND_INPUT_SAMPLE_MAP");
            remap(&mut fp_map2, &m, "SECOND_INPUT_SAMPLE_MAP");
        }
        if mode_same_sample {
            check_by_sample(&fp_map, &fp_map2, &map, &opts, &mut metrics)
        } else {
            cross_check_grouped(
                &fp_map,
                &fp_map2,
                crosscheck_by,
                &map,
                &opts,
                &individuals,
                &mut metrics,
                &mut matrix,
            )
        }
    };

    if metrics.iter().all(|m| m.lod == 0.0) {
        log(
            "ERROR",
            TOOL,
            "No non-zero results found. This is likely an error. Probable cause: there are no reads or variants at fingerprinting sites ",
        );
        std::process::exit(exit_no_valid);
    }
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for m in &metrics {
        file.add_metric(m);
    }
    match &output {
        Some(path) => std::fs::write(path, file.write()).unwrap_or_else(|e| fail(&e.to_string())),
        None => print!("{}", file.write()),
    }
    if let (Some(path), Some((lhs, rhs, cells))) = (&matrix_output, &matrix) {
        let mut text = String::from(crosscheck_by.name());
        for key in rhs {
            text.push('\t');
            text.push_str(key);
        }
        text.push('\n');
        for (row, key) in lhs.iter().enumerate() {
            text.push_str(key);
            for lod in &cells[row] {
                text.push('\t');
                text.push_str(&number_format(*lod));
            }
            text.push('\n');
        }
        std::fs::write(path, text).unwrap_or_else(|e| fail(&e.to_string()));
    }
    if unexpected > 0 {
        std::process::exit(exit_mismatch);
    }
}
