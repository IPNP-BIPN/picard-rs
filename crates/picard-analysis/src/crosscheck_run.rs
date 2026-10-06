//! `CrosscheckFingerprints.doWork`, shared with its deprecated subclass
//! `CrosscheckReadGroupFingerprints`, over [`crate::fingerprint`].
//!
//! Ported from `picard.fingerprint.CrosscheckFingerprints`, `FingerprintChecker.fingerprintFiles`,
//! `Fingerprint.mergeFingerprintsBy` and `CappedHaplotypeProbabilities` at tag 3.4.0.
//!
//! # Every order in the output is a hash table's
//!
//! The rows come out in the order of the merged fingerprints' `HashMap`, keyed by
//! `FingerprintIdDetails`; before that, `fingerprintFiles` collects into a `ConcurrentHashMap`
//! (whose table is sized from the number of files and whose resize reverses part of a bin), which
//! `capFingerprints` copies into a `HashMap`, which `mergeFingerprintsBy` groups through a
//! `HashMap` keyed by the group's string. Each of those is reproduced here, and the file name is
//! part of every key, so the port hashes the path the reference was given
//! ([`crate::fingerprint::reference_path`]).

use crate::fingerprint::{
    calculate_match_results, file_uri, fingerprint_sam_file, java_hash_order,
    load_vcf_fingerprints, Fingerprint, HaplotypeMap, IdDetails, MatchResults, Probs, SamOptions,
    SharedRandom,
};
use crate::java_hash_map::string_hash_code;

/// `CrosscheckMetric.DataType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    File,
    Sample,
    Library,
    ReadGroup,
}

impl DataType {
    pub fn parse(name: &str) -> Option<DataType> {
        match name {
            "FILE" => Some(DataType::File),
            "SAMPLE" => Some(DataType::Sample),
            "LIBRARY" => Some(DataType::Library),
            "READGROUP" => Some(DataType::ReadGroup),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            DataType::File => "FILE",
            DataType::Sample => "SAMPLE",
            DataType::Library => "LIBRARY",
            DataType::ReadGroup => "READGROUP",
        }
    }
}

/// The arguments `doWork` reads.
#[derive(Debug, Clone)]
pub struct Options {
    pub tool: String,
    pub inputs: Vec<String>,
    pub second_inputs: Vec<String>,
    pub haplotype_map: String,
    pub output: Option<String>,
    pub matrix_output: Option<String>,
    pub sample_individual_map: Option<String>,
    pub check_all_others: bool,
    pub lod_threshold: f64,
    pub crosscheck_by: DataType,
    pub tumor_aware: bool,
    pub allow_duplicate_reads: bool,
    pub output_errors_only: bool,
    pub loss_of_het_rate: f64,
    pub expect_all_groups_to_match: bool,
    pub exit_code_when_mismatch: i32,
    pub exit_code_when_no_valid_checks: i32,
    pub max_effect_of_each_haplotype_block: f64,
    pub require_index_files: bool,
    pub strict: bool,
}

/// `ConcurrentHashMap(initialCapacity)` filled in this order, then iterated: the table starts at
/// `tableSizeFor(c + c/2 + 1)`, doubles when the count reaches three quarters of it, and a resize
/// moves the nodes before a bin's last run to the FRONT of their new bin, in reverse.
pub fn concurrent_hash_order<T>(items: Vec<(i32, T)>, initial_capacity: usize) -> Vec<T> {
    let spread = |h: i32| ((h ^ ((h as u32) >> 16) as i32) & 0x7fff_ffff) as usize;
    let wanted = initial_capacity + (initial_capacity >> 1) + 1;
    let mut n = wanted.next_power_of_two().max(1);
    let mut size_ctl = n - (n >> 2);
    let mut table: Vec<Vec<(usize, T)>> = (0..n).map(|_| Vec::new()).collect();
    let mut count = 0usize;
    for (hash, item) in items {
        let h = spread(hash);
        table[h & (n - 1)].push((h, item));
        count += 1;
        if count >= size_ctl {
            let mut next: Vec<Vec<(usize, T)>> = (0..2 * n).map(|_| Vec::new()).collect();
            for (i, bin) in table.into_iter().enumerate() {
                if bin.is_empty() {
                    continue;
                }
                let bits: Vec<bool> = bin.iter().map(|(h, _)| h & n != 0).collect();
                let mut last_run = 0;
                for k in 1..bits.len() {
                    if bits[k] != bits[k - 1] {
                        last_run = k;
                    }
                }
                let mut low: Vec<(usize, T)> = Vec::new();
                let mut high: Vec<(usize, T)> = Vec::new();
                let mut nodes: Vec<(usize, T)> = bin;
                let tail: Vec<(usize, T)> = nodes.split_off(last_run);
                if bits[last_run] {
                    high.extend(tail);
                } else {
                    low.extend(tail);
                }
                for node in nodes {
                    if node.0 & n == 0 {
                        low.insert(0, node);
                    } else {
                        high.insert(0, node);
                    }
                }
                next[i] = low;
                next[i + n] = high;
            }
            table = next;
            n *= 2;
            size_ctl = n - (n >> 2);
        }
    }
    table.into_iter().flatten().map(|(_, item)| item).collect()
}

/// `Fingerprint.getFingerprintIdDetailsStringFunction`.
pub fn group_of(details: &IdDetails, by: DataType) -> String {
    let null = |s: &Option<String>| s.clone().unwrap_or_else(|| "null".to_string());
    let value = match by {
        DataType::ReadGroup => details.platform_unit.clone(),
        DataType::Library => Some(format!(
            "{}::{}",
            null(&details.sample),
            null(&details.library)
        )),
        DataType::File => Some(format!(
            "{}::{}",
            null(&details.file),
            null(&details.sample)
        )),
        DataType::Sample => details.sample.clone(),
    };
    value.unwrap_or_else(|| details.hash_code().to_string())
}

/// `Fingerprint.mergeFingerprintsBy`: one entry per group, its details merged when several share
/// it, with `group` set; in the `HashMap` order of the merged details.
pub fn merge_by(
    entries: &[(IdDetails, Fingerprint)],
    by: DataType,
) -> Vec<(IdDetails, Fingerprint)> {
    let mut groups: Vec<(String, Vec<&(IdDetails, Fingerprint)>)> = Vec::new();
    for entry in entries {
        let key = group_of(&entry.0, by);
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, list)) => list.push(entry),
            None => groups.push((key, vec![entry])),
        }
    }
    let groups = java_hash_order(
        groups
            .into_iter()
            .map(|(k, l)| (string_hash_code(&k), (k, l)))
            .collect(),
    );
    let mut merged: Vec<(IdDetails, Fingerprint)> = Vec::new();
    for (key, list) in groups {
        let mut id = if list.len() == 1 {
            list[0].0.clone()
        } else {
            let mut id = IdDetails::default();
            for (details, _) in &list {
                id.merge(details);
            }
            id
        };
        id.group = Some(key);
        let fp = if list.len() == 1 {
            list[0].1.clone()
        } else {
            let first = &list[0].0;
            let mut fp = Fingerprint::new(first.sample.clone(), None, Some(group_of(first, by)));
            for (_, f) in &list {
                fp.merge(f);
            }
            fp
        };
        merged.push((id, fp));
    }
    java_hash_order(
        merged
            .into_iter()
            .map(|(d, f)| (d.hash_code(), (d, f)))
            .collect(),
    )
}

/// `FingerprintChecker.fingerprintFiles` followed by `capFingerprints`.
fn fingerprint_files(
    files: &[String],
    map: &HaplotypeMap,
    options: &Options,
    random: &mut SharedRandom,
) -> Result<Vec<(IdDetails, Fingerprint)>, String> {
    let mut collected: Vec<(i32, (IdDetails, Fingerprint))> = Vec::new();
    for file in files {
        let reads = file.ends_with(".bam") || file.ends_with(".sam") || file.ends_with(".cram");
        let found: Vec<(IdDetails, Fingerprint)> = if reads {
            let sam_options = SamOptions {
                allow_duplicates: options.allow_duplicate_reads,
                strict: options.strict,
                ..SamOptions::default()
            };
            fingerprint_sam_file(
                file,
                map,
                &sam_options,
                random,
                &|m, b| Probs::sequence(m, b),
                &|p, snp, base, qual| p.add_base(snp, base, qual),
            )?
        } else {
            if options.require_index_files
                && !std::path::Path::new(&format!("{file}.idx")).exists()
                && !std::path::Path::new(&format!("{file}.tbi")).exists()
            {
                return Err(format!(
                    "picard.PicardException: Input VCF file {file} has no index while user required index to proceed."
                ));
            }
            let by_sample = load_vcf_fingerprints(file, map, None, 0.01)?;
            let uri = file_uri(file);
            java_hash_order(
                by_sample
                    .into_iter()
                    .map(|(sample, fp)| {
                        let details = IdDetails {
                            sample: Some(sample),
                            file: Some(uri.clone()),
                            ..IdDetails::default()
                        };
                        (details.hash_code(), (details, fp))
                    })
                    .collect(),
            )
        };
        for (details, fp) in found {
            match collected.iter().position(|(_, (d, _))| d.same(&details)) {
                Some(at) => collected[at].1 .1 = fp,
                None => collected.push((details.hash_code(), (details, fp))),
            }
        }
    }
    let concurrent = concurrent_hash_order(collected, files.len());
    let cap = -options.max_effect_of_each_haplotype_block;
    Ok(java_hash_order(
        concurrent
            .into_iter()
            .map(|(details, fp)| {
                let mut capped =
                    Fingerprint::new(fp.sample.clone(), fp.source.clone(), fp.info.clone());
                for (key, probs) in &fp.blocks {
                    capped.blocks.insert(key.clone(), probs.capped(cap));
                }
                (details.hash_code(), (details, capped))
            })
            .collect(),
    ))
}

/// `CrosscheckMetric.FingerprintResult`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// One row of the table, as `getMatchDetails` fills a `CrosscheckMetric`.
#[derive(Debug, Clone)]
pub struct Metric {
    pub left: IdDetails,
    pub right: IdDetails,
    pub result: &'static str,
    pub data_type: DataType,
    pub results: MatchResults,
}

/// A matrix's row keys, column keys and LODs.
pub type Matrix = (Vec<String>, Vec<String>, Vec<Vec<f64>>);

/// What a run produced.
pub struct Outcome {
    pub metrics: Vec<Metric>,
    /// The matrix's row keys, column keys and LODs, when one was asked for.
    pub matrix: Option<Matrix>,
    pub unexpected: usize,
    /// `log.error` lines the run printed, in order.
    pub errors: Vec<String>,
}

/// A `TabbedInputParser` file of two columns, as `getStringStringMap` reads it.
fn read_two_columns(path: &str, field: &str) -> Result<Vec<(String, String)>, String> {
    let text =
        std::fs::read_to_string(path).map_err(|e| format!("htsjdk.samtools.SAMException: {e}"))?;
    let mut pairs: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            return Err(format!("java.lang.IllegalArgumentException: Each line of the {field} must have exactly two strings separated by a tab. Found: [{}] right before [{}]", fields.join(", "), line));
        }
        if pairs.iter().any(|(k, _)| k == fields[0]) {
            return Err(format!("java.lang.IllegalArgumentException: Strings in first column of the {field} must be unique. found [{}] twice", fields[0]));
        }
        pairs.push((fields[0].to_string(), fields[1].to_string()));
    }
    Ok(pairs)
}

/// `doWork`, up to the writing.
pub fn run(options: &Options) -> Result<Outcome, String> {
    let map = HaplotypeMap::load(&options.haplotype_map)?;
    let mut by = options.crosscheck_by;
    let same_sample = !options.second_inputs.is_empty() && !options.check_all_others;
    if same_sample {
        by = DataType::Sample;
    }
    let mut random = SharedRandom(crate::theoretical_sensitivity::JavaRandom::new(42));
    let first = fingerprint_files(&options.inputs, &map, options, &mut random)?;
    let individuals = match &options.sample_individual_map {
        Some(path) => Some(read_two_columns(path, "SAMPLE_INDIVIDUAL_MAP")?),
        None => None,
    };
    let resolve = |sample: &Option<String>| -> Option<String> {
        let s = sample.clone()?;
        Some(
            individuals
                .as_ref()
                .and_then(|m| m.iter().find(|(k, _)| *k == s).map(|(_, v)| v.clone()))
                .unwrap_or(s),
        )
    };
    let mut metrics: Vec<Metric> = Vec::new();
    let mut errors: Vec<String> = Vec::new();
    let mut unexpected = 0usize;
    let mut matrix = None;
    let mut cross = |lhs: Vec<(IdDetails, Fingerprint)>,
                     rhs: Vec<(IdDetails, Fingerprint)>,
                     metrics: &mut Vec<Metric>|
     -> usize {
        let mut unexpected = 0;
        let mut lods = vec![vec![0.0; rhs.len()]; lhs.len()];
        for (row, (lid, lfp)) in lhs.iter().enumerate() {
            for (col, (rid, rfp)) in rhs.iter().enumerate() {
                let expected = options.expect_all_groups_to_match
                    || resolve(&lid.sample) == resolve(&rid.sample);
                let results = calculate_match_results(
                    &map,
                    lfp,
                    rfp,
                    options.loss_of_het_rate,
                    false,
                    options.tumor_aware,
                );
                let v = verdict(expected, results.lod, options.lod_threshold);
                if !options.output_errors_only
                    || v == Verdict::Inconclusive
                    || v.is_expected() == Some(false)
                {
                    metrics.push(Metric {
                        left: lid.clone(),
                        right: rid.clone(),
                        result: v.name(),
                        data_type: by,
                        results: results.clone(),
                    });
                }
                if v.is_expected() == Some(false) {
                    unexpected += 1;
                }
                lods[row][col] = results.lod;
            }
        }
        if options.matrix_output.is_some() {
            matrix = Some((
                lhs.iter()
                    .map(|(d, _)| d.group.clone().unwrap_or_default())
                    .collect(),
                rhs.iter()
                    .map(|(d, _)| d.group.clone().unwrap_or_default())
                    .collect(),
                lods,
            ));
        }
        unexpected
    };
    if options.second_inputs.is_empty() {
        let lhs = merge_by(&first, by);
        let rhs = merge_by(&first, by);
        unexpected += cross(lhs, rhs, &mut metrics);
    } else {
        let second = fingerprint_files(&options.second_inputs, &map, options, &mut random)?;
        if options.check_all_others {
            let lhs = merge_by(&first, by);
            let rhs = merge_by(&second, by);
            unexpected += cross(lhs, rhs, &mut metrics);
        } else {
            // `checkFingerprintsBySample`.
            let lhs = merge_by(&first, DataType::Sample);
            let rhs = merge_by(&second, DataType::Sample);
            let keyed = |list: &Vec<(IdDetails, Fingerprint)>| -> Vec<String> {
                java_hash_order(
                    list.iter()
                        .map(|(d, _)| {
                            let g = d.group.clone().unwrap_or_default();
                            (string_hash_code(&g), g)
                        })
                        .collect(),
                )
            };
            let mut samples = keyed(&lhs);
            for s in keyed(&rhs) {
                if !samples.contains(&s) {
                    samples.push(s);
                }
            }
            for sample in samples {
                let l = lhs
                    .iter()
                    .find(|(d, _)| d.group.as_deref() == Some(&sample));
                let r = rhs
                    .iter()
                    .find(|(d, _)| d.group.as_deref() == Some(&sample));
                let (Some((lid, lfp)), Some((rid, rfp))) = (l, r) else {
                    errors.push(format!(
                        "sample {sample} is missing from {} group",
                        if l.is_none() { "LEFT" } else { "RIGHT" }
                    ));
                    unexpected += 1;
                    continue;
                };
                if lfp.blocks.is_empty() || rfp.blocks.is_empty() {
                    errors.push(format!("sample {sample} from {} group was not fingerprinted.  Probably there are no reads/variants at fingerprinting sites.", if lfp.blocks.is_empty() { "LEFT" } else { "RIGHT" }));
                    unexpected += 1;
                }
                let results = calculate_match_results(
                    &map,
                    lfp,
                    rfp,
                    options.loss_of_het_rate,
                    false,
                    options.tumor_aware,
                );
                let v = verdict(true, results.lod, options.lod_threshold);
                if !options.output_errors_only || v.is_expected() != Some(true) {
                    metrics.push(Metric {
                        left: lid.clone(),
                        right: rid.clone(),
                        result: v.name(),
                        data_type: DataType::Sample,
                        results: results.clone(),
                    });
                }
                if v.is_expected() == Some(false) {
                    unexpected += 1;
                }
                if results.lod == 0.0 {
                    errors.push("LOD score of zero found when checking sample fingerprints.  Probably there are no reads/variants at fingerprinting sites for one of the samples".to_string());
                    unexpected += 1;
                }
            }
        }
    }
    Ok(Outcome {
        metrics,
        matrix,
        unexpected,
        errors,
    })
}

/// `NumberFormat.getInstance()` with at most four fraction digits, as `writeMatrix` uses it.
pub fn matrix_number(value: f64) -> String {
    let text = format!("{value:.4}");
    let (sign, digits) = match text.strip_prefix('-') {
        Some(rest) => ("-", rest.to_string()),
        None => ("", text),
    };
    let (int, frac) = digits.split_once('.').unwrap_or((&digits, ""));
    let frac = frac.trim_end_matches('0');
    let mut grouped = String::new();
    for (i, c) in int.chars().enumerate() {
        if i > 0 && (int.len() - i) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(c);
    }
    if frac.is_empty() {
        format!("{sign}{grouped}")
    } else {
        format!("{sign}{grouped}.{frac}")
    }
}

/// `writeMatrix`.
pub fn matrix_text(by: DataType, rows: &[String], cols: &[String], lods: &[Vec<f64>]) -> String {
    let mut out = String::from(by.name());
    for c in cols {
        out.push('\t');
        out.push_str(c);
    }
    out.push('\n');
    for (r, key) in rows.iter().enumerate() {
        out.push_str(key);
        for lod in &lods[r] {
            out.push('\t');
            out.push_str(&matrix_number(*lod));
        }
        out.push('\n');
    }
    out
}

/// The `CrosscheckMetric` columns, in `getFields()` order.
const COLUMNS: &[&str] = &[
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
];

impl htsjdk_metrics::file::MetricBean for Metric {
    fn class_name(&self) -> &str {
        "picard.fingerprint.CrosscheckMetric"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<htsjdk_metrics::file::Value> {
        use htsjdk_metrics::file::Value;
        let s = |v: &Option<String>| v.clone().map_or(Value::Null, Value::Str);
        let i = |v: &Option<i32>| v.map_or(Value::Null, |n| Value::Long(i64::from(n)));
        vec![
            s(&self.left.group),
            s(&self.right.group),
            Value::Str(self.result.to_string()),
            Value::Str(self.data_type.name().to_string()),
            Value::Double(self.results.lod),
            Value::Double(self.results.lod_tn),
            Value::Double(self.results.lod_nt),
            s(&self.left.run_barcode),
            i(&self.left.run_lane),
            s(&self.left.molecular_barcode),
            s(&self.left.library),
            s(&self.left.sample),
            s(&self.left.file),
            s(&self.right.run_barcode),
            i(&self.right.run_lane),
            s(&self.right.molecular_barcode),
            s(&self.right.library),
            s(&self.right.sample),
            s(&self.right.file),
        ]
    }
}

/// The end of `doWork`: refuse a run with nothing but zeroes, write the table and the matrix,
/// and the exit code.
pub fn write_outcome(options: &Options, outcome: &Outcome) -> i32 {
    for e in &outcome.errors {
        crate::vcf_io::log_error("CrosscheckFingerprints", e);
    }
    if outcome.metrics.iter().all(|m| m.results.lod == 0.0) {
        crate::vcf_io::log_error(
            "CrosscheckFingerprints",
            "No non-zero results found. This is likely an error. Probable cause: there are no reads or variants at fingerprinting sites ",
        );
        return options.exit_code_when_no_valid_checks;
    }
    let mut file = htsjdk_metrics::file::MetricsFile::new();
    file.add_header(&format!("{} <command line>", options.tool));
    file.add_header("Started on: <timestamp>");
    for m in &outcome.metrics {
        file.add_metric(m);
    }
    let text = file.write();
    match &options.output {
        Some(path) => {
            if let Err(e) = std::fs::write(path, &text) {
                crate::metrics_cli::thrown(&format!("htsjdk.samtools.SAMException: {e}"));
            }
        }
        None => print!("{text}"),
    }
    if let (Some(path), Some((rows, cols, lods))) = (&options.matrix_output, &outcome.matrix) {
        let by = if !options.second_inputs.is_empty() && !options.check_all_others {
            DataType::Sample
        } else {
            options.crosscheck_by
        };
        if let Err(e) = std::fs::write(path, matrix_text(by, rows, cols, lods)) {
            crate::metrics_cli::thrown(&format!("java.lang.RuntimeException: {e}"));
        }
    }
    if outcome.unexpected > 0 {
        options.exit_code_when_mismatch
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_numbers_are_grouped_and_trimmed() {
        assert_eq!(matrix_number(1.33714), "1.3371");
        assert_eq!(matrix_number(-6.09551), "-6.0955");
        assert_eq!(matrix_number(2.5), "2.5");
        assert_eq!(matrix_number(1234.0), "1,234");
    }

    #[test]
    fn a_concurrent_map_resizes_at_three_quarters() {
        // Two files size the table at four, so the third entry resizes it to eight: 5 and 1 shared
        // bucket 1, and the split sends 5 to bucket 5 while 1 stays.
        let order = concurrent_hash_order(vec![(2, "two"), (5, "five"), (1, "one")], 2);
        assert_eq!(order, vec!["one", "two", "five"]);
    }
}
