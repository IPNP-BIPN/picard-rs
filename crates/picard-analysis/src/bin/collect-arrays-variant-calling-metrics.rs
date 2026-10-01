//! `CollectArraysVariantCallingMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.CollectArraysVariantCallingMetrics.doWork`,
//! `picard.arrays.ArraysCallingMetricAccumulator` and the `MergeableMetricBase` folds the result
//! goes through, at tag 3.4.0.
//!
//! In the reference's order: the input is opened with an index REQUIRED (`VCFFileReader(INPUT,
//! true)`), the dictionary comes from `SEQUENCE_DICTIONARY` or the input, dbSNP is loaded into a
//! SNP bitset and an indel bitset after its own dictionary is asserted equal to that one, the
//! twenty-three control codes are read from the header (each `<control>|<category>|<red>|<green>`),
//! and only then does each accumulator's `setup` read the sample-level header lines, required and
//! optional, the genders through `Sex.fromString` and the two dates through their
//! `SimpleDateFormat`s (`picard_analysis::java_date`).
//!
//! The variants are walked per contig of the header's dictionary, through the index, so a record
//! on a contig the header does not declare is never seen. Each accumulator holds a detail row for
//! every sample; the rows are grouped by sample in a `HashMap` (which is the order they are
//! written in), folded (`@MergeByAdding` summed, `@MergeByAssertEquals` taken), and the summary is
//! the sum of every detail. The derived fields divide by a `float` where the reference does: the
//! call rates and `PCT_DBSNP` are floats widened for printing, the two het ratios doubles.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::variant::{Genotype, Value as VcfValue, VariantContext};
use picard_analysis::java_date::{iso8601, parse as parse_date, Zone};
use picard_analysis::java_hash_map::JavaHashMap;
use picard_analysis::java_number::{parse_double, parse_int};
use picard_analysis::metrics_cli::{refuse_validation, Args};
use picard_analysis::vcf_io::{
    die, extract_dictionary, header_dictionary, is_same_sequence, read_path, Sequence,
};

const TOOL: &str = "CollectArraysVariantCallingMetrics";

/// `ArraysControlInfo.CONTROL_INFO`, the header keys in the order the rows are written.
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

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

fn number_format(message: String) -> String {
    format!("java.lang.NumberFormatException: {message}")
}

/// `VCFHeader.getOtherHeaderLine(key).getValue()`.
fn other_line<'a>(lines: &'a [HeaderLine], key: &str) -> Option<&'a str> {
    lines.iter().rev().find_map(|line| match line {
        HeaderLine::Unstructured { key: k, value } if k == key => Some(value.as_str()),
        _ => None,
    })
}

fn required<'a>(lines: &'a [HeaderLine], key: &str) -> Result<&'a str, String> {
    other_line(lines, key).ok_or_else(|| {
        format!(
            "java.lang.IllegalArgumentException: Input VCF file is missing header line of type \
             '{key}'"
        )
    })
}

/// `Sex.fromString(value).toSymbol()`.
fn sex_symbol(value: &str) -> Result<&'static str, String> {
    let sexes = [
        ("M", "Male"),
        ("F", "Female"),
        ("U", "Unknown"),
        ("N", "NotReported"),
    ];
    let matches: Vec<&'static str> = sexes
        .iter()
        .filter(|(symbol, name)| {
            value.eq_ignore_ascii_case(symbol) || value.eq_ignore_ascii_case(name)
        })
        .map(|(symbol, _)| *symbol)
        .collect();
    if matches.len() == 1 {
        Ok(matches[0])
    } else {
        Err(format!(
            "picard.PicardException: Unrecognized Sex string: {value}"
        ))
    }
}

/// `getSexConcordance`, over the three symbols.
fn sex_concordance(reported: &str, fingerprint: &str, autocall: &str) -> bool {
    let votes = [reported, fingerprint, autocall];
    let count = |symbol: &str| votes.iter().filter(|v| **v == symbol).count();
    if count("U") + count("N") == 3 {
        return false;
    }
    (count("F") > 1 && count("M") == 0) || (count("M") > 1 && count("F") == 0)
}

/// What `setup` reads from the header, which every detail row starts from.
#[derive(Clone)]
struct SampleHeader {
    sample_alias: String,
    pipeline_version: Option<String>,
    analysis_version: Option<i32>,
    chip_type: String,
    reported_gender: &'static str,
    fingerprint_gender: &'static str,
    gtc_call_rate: Option<f64>,
    autocall_gender: &'static str,
    autocall_version: String,
    autocall_date: String,
    imaging_date: String,
    extended_manifest_version: String,
    zcall_version: Option<String>,
    zcall_thresholds: Option<String>,
    cluster_file: String,
    p95_green: i32,
    p95_red: i32,
    scanner_name: String,
}

/// `ArraysCallingMetricAccumulator.setup`, in its order.
fn setup(lines: &[HeaderLine]) -> Result<SampleHeader, String> {
    let optional = |key: &str| other_line(lines, key).map(str::to_string);
    let optional_gender = |key: &str| -> Result<&'static str, String> {
        match other_line(lines, key) {
            Some(value) => sex_symbol(value),
            None => Ok("N"),
        }
    };
    let date = |key: &str, pattern: &str, zone: Zone| -> Result<String, String> {
        let value = required(lines, key)?;
        parse_date(pattern, zone, value)
            .map(iso8601)
            .ok_or_else(|| {
                format!(
                "java.lang.IllegalArgumentException: Unrecognized date for '{key}' in VCF header \
                 ({value})"
            )
            })
    };
    let integer = |value: &str| parse_int(value).map_err(number_format);

    let sample_alias = required(lines, "sampleAlias")?.to_string();
    let pipeline_version = optional("pipelineVersion");
    let analysis_version = match other_line(lines, "analysisVersionNumber") {
        Some(value) => Some(integer(value)?),
        None => None,
    };
    let chip_type = required(lines, "arrayType")?.to_string();
    let reported_gender = optional_gender("expectedGender")?;
    let fingerprint_gender = optional_gender("fingerprintGender")?;
    let gtc_call_rate = match other_line(lines, "gtcCallRate") {
        Some(value) => Some(parse_double(value).map_err(number_format)?),
        None => None,
    };
    let autocall_gender = sex_symbol(required(lines, "autocallGender")?)?;
    let autocall_version = required(lines, "autocallVersion")?.to_string();
    let autocall_date = date("autocallDate", "MM/dd/yyyy HH:mm", Zone::Utc)?;
    let imaging_date = date("imagingDate", "MM/dd/yyyy hh:mm:ss a", Zone::NewYork)?;
    let extended_manifest_version = required(lines, "extendedIlluminaManifestVersion")?.to_string();
    let zcall_version = optional("zcallVersion");
    let zcall_thresholds = optional("zcallThresholds");
    let cluster_file = required(lines, "clusterFile")?.to_string();
    let p95_green = integer(required(lines, "p95Green")?)?;
    let p95_red = integer(required(lines, "p95Red")?)?;
    let scanner_name = required(lines, "scannerName")?.to_string();
    Ok(SampleHeader {
        sample_alias,
        pipeline_version,
        analysis_version,
        chip_type,
        reported_gender,
        fingerprint_gender,
        gtc_call_rate,
        autocall_gender,
        autocall_version,
        autocall_date,
        imaging_date,
        extended_manifest_version,
        zcall_version,
        zcall_thresholds,
        cluster_file,
        p95_green,
        p95_red,
        scanner_name,
    })
}

/// The `@MergeByAdding` counters of one row.
#[derive(Default, Clone, Copy)]
struct Counts {
    assays: i64,
    non_filtered_assays: i64,
    filtered_assays: i64,
    zeroed_out_assays: i64,
    snps: i64,
    indels: i64,
    calls: i64,
    autocall_calls: i64,
    no_calls: i64,
    in_db_snp: i64,
    singletons: i64,
    hets: i64,
    hom_vars: i64,
}

impl Counts {
    fn add(&mut self, other: &Counts) {
        self.assays += other.assays;
        self.non_filtered_assays += other.non_filtered_assays;
        self.filtered_assays += other.filtered_assays;
        self.zeroed_out_assays += other.zeroed_out_assays;
        self.snps += other.snps;
        self.indels += other.indels;
        self.calls += other.calls;
        self.autocall_calls += other.autocall_calls;
        self.no_calls += other.no_calls;
        self.in_db_snp += other.in_db_snp;
        self.singletons += other.singletons;
        self.hets += other.hets;
        self.hom_vars += other.hom_vars;
    }

    /// The summary's columns, derived fields included, in declaration order.
    fn summary_values(&self) -> Vec<Value> {
        let pct_dbsnp = self.in_db_snp as f32 / self.snps as f32;
        let call_rate = self.calls as f32 / self.non_filtered_assays as f32;
        let autocall_call_rate = self.autocall_calls as f32 / self.non_filtered_assays as f32;
        vec![
            Value::Long(self.assays),
            Value::Long(self.non_filtered_assays),
            Value::Long(self.filtered_assays),
            Value::Long(self.zeroed_out_assays),
            Value::Long(self.snps),
            Value::Long(self.indels),
            Value::Long(self.calls),
            Value::Long(self.autocall_calls),
            Value::Long(self.no_calls),
            Value::Long(self.in_db_snp),
            Value::Long(self.snps - self.in_db_snp),
            Value::Double(f64::from(pct_dbsnp)),
            Value::Double(f64::from(call_rate)),
            Value::Double(f64::from(autocall_call_rate)),
            Value::Long(self.singletons),
        ]
    }
}

const SUMMARY_COLUMNS: [&str; 15] = [
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

/// `Class.getFields()` on the detail class: its own public fields, then the summary's.
const DETAIL_COLUMNS: [&str; 38] = [
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

struct Summary(Counts);

impl MetricBean for Summary {
    fn class_name(&self) -> &str {
        "picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &SUMMARY_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        self.0.summary_values()
    }
}

struct Detail<'a> {
    sample: String,
    header: &'a SampleHeader,
    counts: Counts,
    threshold: f64,
}

impl MetricBean for Detail<'_> {
    fn class_name(&self) -> &str {
        "picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &DETAIL_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let h = self.header;
        let c = &self.counts;
        let optional = |value: &Option<String>| match value {
            Some(v) => Value::Str(v.clone()),
            None => Value::Null,
        };
        let autocall_call_rate = c.autocall_calls as f32 / c.non_filtered_assays as f32;
        let mut values = vec![
            Value::Str(self.sample.clone()),
            Value::Str(h.sample_alias.clone()),
            h.analysis_version
                .map_or(Value::Null, |v| Value::Long(i64::from(v))),
            Value::Str(h.chip_type.clone()),
            // `AUTOCALL_CALL_RATE > CALL_RATE_PF_THRESHOLD`: the float widened to the Double.
            Value::Bool(f64::from(autocall_call_rate) > self.threshold),
            Value::Str(h.autocall_date.clone()),
            Value::Str(h.imaging_date.clone()),
            Value::Bool(h.zcall_thresholds.as_deref().is_some_and(|t| !t.is_empty())),
            h.gtc_call_rate.map_or(Value::Null, Value::Double),
            Value::Str(h.autocall_gender.to_string()),
            Value::Str(h.fingerprint_gender.to_string()),
            Value::Str(h.reported_gender.to_string()),
            Value::Bool(sex_concordance(
                h.reported_gender,
                h.fingerprint_gender,
                h.autocall_gender,
            )),
            Value::Double(c.hets as f64 / c.calls as f64),
            Value::Str(h.cluster_file.clone()),
            Value::Long(i64::from(h.p95_green)),
            Value::Long(i64::from(h.p95_red)),
            Value::Str(h.autocall_version.clone()),
            optional(&h.zcall_version),
            Value::Str(h.extended_manifest_version.clone()),
            Value::Double(c.hets as f64 / c.hom_vars as f64),
            Value::Str(h.scanner_name.clone()),
            optional(&h.pipeline_version),
        ];
        values.extend(c.summary_values());
        values
    }
}

/// `VariantContext.getType()`, as far as `isSNP` and `isIndel` ask.
#[derive(PartialEq)]
enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

fn variant_type(vc: &VariantContext) -> VariantType {
    let reference = vc.reference();
    let mut kind: Option<VariantType> = None;
    for allele in vc.alternate_alleles() {
        let this = if allele.is_symbolic() {
            VariantType::Symbolic
        } else if reference.len() == allele.len() {
            if allele.len() == 1 {
                VariantType::Snp
            } else {
                VariantType::Mnp
            }
        } else {
            VariantType::Indel
        };
        match &kind {
            None => kind = Some(this),
            Some(existing) if *existing != this => return VariantType::Mixed,
            Some(_) => {}
        }
    }
    kind.unwrap_or(VariantType::NoVariation)
}

/// `vc.isFiltered()`: filters applied, and at least one.
fn filters(vc: &VariantContext) -> &[String] {
    vc.filters.as_deref().unwrap_or(&[])
}

/// `CallingMetricAccumulator.getSingletonSample`: the first two het or hom-var genotypes, their
/// variant chromosomes added up; exactly one is a singleton.
fn singleton_sample(vc: &VariantContext) -> Option<&str> {
    let carriers: Vec<&Genotype> = vc
        .genotypes
        .iter()
        .filter(|g| g.is_het() || g.is_hom_var())
        .take(2)
        .collect();
    let chromosomes: i32 = carriers
        .iter()
        .map(|g| if g.is_het() { 1 } else { 2 })
        .sum();
    if chromosomes == 1 {
        carriers.last().map(|g| g.sample_name.as_str())
    } else {
        None
    }
}

/// The dbSNP sites of one variant kind, by contig.
#[derive(Default)]
struct Sites(std::collections::HashMap<String, std::collections::BTreeSet<i64>>);

impl Sites {
    fn contains(&self, contig: &str, position: i64) -> bool {
        self.0
            .get(contig)
            .is_some_and(|set| set.contains(&position))
    }
}

/// `DbSnpBitSetUtil.createSnpAndIndelBitSets`: only the SNP bitset is ever asked.
fn load_db_snp(path: &str, dictionary: Option<&[Sequence]>) -> Result<Sites, String> {
    let vcf = read_path(path)?;
    if let (Some(own), Some(input)) = (header_dictionary(&vcf.file.header), dictionary) {
        if let Err(cause) = assert_dictionary_lists(&own, input) {
            return Err(format!(
                "picard.PicardException: Sequence dictionary for (DBSNP: file:{}) does not match \
                 sequence dictionary for (INPUT)\nCaused by: \
                 htsjdk.samtools.util.SequenceUtil$SequenceListsDifferException: {cause}",
                absolute(path)
            ));
        }
    }
    let mut snps = Sites::default();
    for record in &vcf.records {
        let vc = &record.variant;
        if variant_type(vc) == VariantType::Snp {
            let set = snps.0.entry(vc.contig.clone()).or_default();
            for position in vc.start..=vc.stop {
                set.insert(position);
            }
        }
    }
    Ok(snps)
}

/// `SequenceUtil.assertSequenceListsEqual`, as the message of its exception.
fn assert_dictionary_lists(first: &[Sequence], second: &[Sequence]) -> Result<(), String> {
    if first.len() != second.len() {
        return Err(format!(
            "Sequence dictionaries are not the same size ({}, {})",
            first.len(),
            second.len()
        ));
    }
    for (index, (a, b)) in first.iter().zip(second).enumerate() {
        if !is_same_sequence(a, b) {
            return Err(format!("Sequences at index {index} don't match"));
        }
    }
    Ok(())
}

/// One sample's counters over one variant: `updateDetailMetric`.
fn update(
    counts: &mut Counts,
    genotype: &Genotype,
    vc: &VariantContext,
    snps: &Sites,
    singleton: bool,
) {
    counts.assays += 1;
    let filters = filters(vc);
    if filters.is_empty() || filters.iter().any(|f| f == "DUPE") {
        counts.non_filtered_assays += 1;
        if genotype.is_called() {
            counts.calls += 1;
            // `getExtendedAttribute(GTA, getGenotypeString())`: a called genotype's own string is
            // never `./.`, so only a GTA that says so is not an autocall.
            let autocall = match genotype.get("GTA") {
                Some(VcfValue::Str(gta)) => gta != "./.",
                _ => true,
            };
            if autocall {
                counts.autocall_calls += 1;
            }
        } else {
            counts.no_calls += 1;
        }
        match variant_type(vc) {
            VariantType::Snp => {
                counts.snps += 1;
                if snps.contains(&vc.contig, vc.start) {
                    counts.in_db_snp += 1;
                }
            }
            VariantType::Indel => counts.indels += 1,
            _ => {}
        }
        if singleton {
            counts.singletons += 1;
        }
        if genotype.is_het() {
            counts.hets += 1;
        } else if genotype.is_hom_var() {
            counts.hom_vars += 1;
        }
    } else {
        counts.filtered_assays += 1;
        if filters.iter().any(|f| f == "ZEROED_OUT_ASSAY") {
            counts.zeroed_out_assays += 1;
        }
    }
}

struct Control {
    control: String,
    category: String,
    red: i32,
    green: i32,
}

impl MetricBean for Control {
    fn class_name(&self) -> &str {
        "picard.arrays.CollectArraysVariantCallingMetrics$ArraysControlCodesSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &["CONTROL", "CATEGORY", "RED", "GREEN"]
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.control.clone()),
            Value::Str(self.category.clone()),
            Value::Long(i64::from(self.red)),
            Value::Long(i64::from(self.green)),
        ]
    }
}

/// `parseControlHeaderString`: `split("\\|")`, which drops trailing empty tokens, and the four
/// tokens read in argument order, so a short value fails on the first index it lacks.
fn parse_control(value: &str) -> Result<Control, String> {
    let mut tokens: Vec<&str> = value.split('|').collect();
    while tokens.len() > 1 && tokens.last() == Some(&"") {
        tokens.pop();
    }
    if tokens.len() == 1 && tokens[0].is_empty() && value.contains('|') {
        tokens.clear();
    }
    let token = |index: usize| {
        tokens.get(index).copied().ok_or_else(|| {
            format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {index} out of bounds for length {}",
                tokens.len()
            )
        })
    };
    let control = token(0)?.to_string();
    let category = token(1)?.to_string();
    let red = parse_int(token(2)?).map_err(number_format)?;
    let green = parse_int(token(3)?).map_err(number_format)?;
    Ok(Control {
        control,
        category,
        red,
        green,
    })
}

fn new_metrics_file() -> MetricsFile {
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    file
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("SD", "SEQUENCE_DICTIONARY"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let db_snp = args.required("DBSNP");
    let sequence_dictionary = args.get("SEQUENCE_DICTIONARY").map(str::to_string);
    let threshold = args.double("CALL_RATE_PF_THRESHOLD", 0.98);
    let _processors = args.int("NUM_PROCESSORS", 0);

    // `customCommandLineValidation`.
    if threshold <= 0.0 || threshold > 1.0 {
        refuse_validation(
            TOOL,
            &["The parameter CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0".to_string()],
        );
    }

    for path in std::iter::once(&input)
        .chain(std::iter::once(&db_snp))
        .chain(sequence_dictionary.iter())
    {
        if !std::path::Path::new(path).is_file() {
            die(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(path)
            ));
        }
    }

    // `new VCFFileReader(INPUT, true)`: the index is required.
    if !std::path::Path::new(&format!("{input}.idx")).is_file() {
        die(&format!(
            "htsjdk.tribble.TribbleException: An index is required, but none found with file \
             ending .idx, for input source: file://{}",
            absolute(&input)
        ));
    }
    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    let lines = &vcf.file.header.lines;

    let dictionary = extract_dictionary(sequence_dictionary.as_deref().unwrap_or(&input))
        .unwrap_or_else(|e| die(&format!("htsjdk.samtools.SAMException: {e}")));
    let snps = load_db_snp(&db_snp, dictionary.as_deref()).unwrap_or_else(|e| die(&e));

    let mut control_file = new_metrics_file();
    for control in CONTROLS {
        let value = required(lines, control).unwrap_or_else(|e| die(&e));
        let parsed = parse_control(value).unwrap_or_else(|e| die(&e));
        control_file.add_metric(&parsed);
    }

    // The segments come from the header's own dictionary; without one the generator has nothing
    // to call `getSequences` on.
    let segments = header_dictionary(&vcf.file.header).unwrap_or_else(|| {
        die("java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \"dict\" is null")
    });
    let header = setup(lines).unwrap_or_else(|e| die(&e));

    // One detail row per sample, in the order the HashMap of rows hands them out.
    let samples = &vcf.file.header.samples;
    let mut counts: Vec<Counts> = vec![Counts::default(); samples.len()];
    for segment in &segments {
        for record in &vcf.records {
            let vc = &record.variant;
            let overlaps =
                vc.contig == segment.name && vc.start <= segment.length.max(1) && vc.stop >= 1;
            if !overlaps {
                continue;
            }
            let singleton = singleton_sample(vc);
            for genotype in vc.genotypes.iter() {
                let Some(index) = samples.iter().position(|s| *s == genotype.sample_name) else {
                    continue;
                };
                let is_singleton = singleton == Some(genotype.sample_name.as_str());
                update(&mut counts[index], genotype, vc, &snps, is_singleton);
            }
        }
    }

    let mut order: JavaHashMap<usize> = JavaHashMap::new();
    for (index, sample) in samples.iter().enumerate() {
        order.put(sample, index);
    }
    let mut detail_file = new_metrics_file();
    let mut total = Counts::default();
    for (sample, index) in order.iter() {
        total.add(&counts[*index]);
        detail_file.add_metric(&Detail {
            sample: sample.to_string(),
            header: &header,
            counts: counts[*index],
            threshold,
        });
    }
    let mut summary_file = new_metrics_file();
    summary_file.add_metric(&Summary(total));

    let prefix = format!("{}.", absolute(&output));
    let write = |extension: &str, file: &MetricsFile| {
        if let Err(e) = std::fs::write(format!("{prefix}{extension}"), file.write()) {
            die(&format!(
                "htsjdk.samtools.SAMException: Could not write metrics file: {e}"
            ));
        }
    };
    write("arrays_variant_calling_detail_metrics", &detail_file);
    write("arrays_variant_calling_summary_metrics", &summary_file);
    write("arrays_control_code_summary_metrics", &control_file);
}
