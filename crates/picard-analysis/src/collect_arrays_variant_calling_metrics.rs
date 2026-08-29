//! `CollectArraysVariantCallingMetrics`: what an array's metrics take from the header.
//!
//! Almost everything the tool reports comes out of the VCF's HEADER rather than out of its
//! variants: the chip, the sample, the two genders, the call rate the genotyping software already
//! computed, the dates, the cluster file and the twenty-three control codes. What the variants
//! decide is a handful of counters and one threshold.
//!
//! Ported from `picard.arrays.CollectArraysVariantCallingMetrics`,
//! `picard.arrays.ArraysCallingMetricAccumulator` and `picard.arrays.illumina.ArraysControlInfo`.

/// The header lines the accumulator refuses to run without, in the order it asks for them.
pub const REQUIRED_HEADER_LINES: [&str; 9] = [
    "sampleAlias",
    "arrayType",
    "autocallGender",
    "autocallDate",
    "imagingDate",
    "extendedIlluminaManifestVersion",
    "clusterFile",
    "p95Green",
    "p95Red",
];

/// The header lines it reads when they are there and does without when they are not.
///
/// Which of the two lists a line is on is not guessable from its name: `clusterFile` is required
/// and `zcallThresholds` is not, `autocallVersion` is required and `pipelineVersion` is not.
pub const OPTIONAL_HEADER_LINES: [&str; 7] = [
    "pipelineVersion",
    "analysisVersionNumber",
    "expectedGender",
    "fingerprintGender",
    "gtcCallRate",
    "zcallVersion",
    "zcallThresholds",
];

/// `scannerName` and `autocallVersion` are required too, and are named apart because they are read
/// through the plain accessor rather than through the optional one.
pub const REQUIRED_HEADER_LINES_ALSO: [&str; 2] = ["autocallVersion", "scannerName"];

/// The message a missing required line produces, which names the line.
pub fn missing_header_line_message(name: &str) -> String {
    format!("Input VCF file is missing header line of type '{name}'")
}

/// A gender that is not in the header, which is `NotReported` and not an error.
pub const NOT_REPORTED_GENDER: &str = "U";

/// The twenty-three control codes, as `ArraysControlInfo.CONTROL_INFO` declares them.
///
/// The pair beside each name is the category and the two intensities, and the intensities in the
/// class itself are zero: what a run reports comes from the VCF's own header line, whose value is
/// the four fields separated by a pipe.
pub const CONTROL_CODES: [(&str, &str); 23] = [
    ("DNP(High)", "Staining"),
    ("DNP(Bgnd)", "Staining"),
    ("Biotin(High)", "Staining"),
    ("Biotin(Bgnd)", "Staining"),
    ("Extension(A)", "Extension"),
    ("Extension(T)", "Extension"),
    ("Extension(C)", "Extension"),
    ("Extension(G)", "Extension"),
    ("TargetRemoval", "TargetRemoval"),
    ("Hyb(High)", "Hybridization"),
    ("Hyb(Medium)", "Hybridization"),
    ("Hyb(Low)", "Hybridization"),
    ("String(PM)", "Stringency"),
    ("String(MM)", "Stringency"),
    ("NSB(Bgnd)Red", "Non-SpecificBinding"),
    ("NSB(Bgnd)Purple", "Non-SpecificBinding"),
    ("NSB(Bgnd)Blue", "Non-SpecificBinding"),
    ("NSB(Bgnd)Green", "Non-SpecificBinding"),
    ("NP(A)", "Non-Polymorphic"),
    ("NP(T)", "Non-Polymorphic"),
    ("NP(C)", "Non-Polymorphic"),
    ("NP(G)", "Non-Polymorphic"),
    ("Restore", "Restoration"),
];

/// One control code as the header carries it: `<control>|<category>|<red>|<green>`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ControlCode {
    pub control: String,
    pub category: String,
    pub red: i32,
    pub green: i32,
}

/// `parseControlHeaderString`, which splits on a pipe and takes four fields.
pub fn parse_control_header(value: &str) -> Option<ControlCode> {
    let fields: Vec<&str> = value.split('|').collect();
    if fields.len() < 4 {
        return None;
    }
    Some(ControlCode {
        control: fields[0].to_string(),
        category: fields[1].to_string(),
        red: fields[2].parse().ok()?,
        green: fields[3].parse().ok()?,
    })
}

/// The three files, named from one prefix and a dot.
pub const DETAIL_FILE_EXTENSION: &str = "arrays_variant_calling_detail_metrics";
pub const SUMMARY_FILE_EXTENSION: &str = "arrays_variant_calling_summary_metrics";
pub const CONTROL_FILE_EXTENSION: &str = "arrays_control_code_summary_metrics";

/// `OUTPUT.getAbsolutePath() + "."` and then each extension.
pub fn file_names(prefix: &str) -> [String; 3] {
    [
        format!("{prefix}.{DETAIL_FILE_EXTENSION}"),
        format!("{prefix}.{SUMMARY_FILE_EXTENSION}"),
        format!("{prefix}.{CONTROL_FILE_EXTENSION}"),
    ]
}

/// What one assay contributes, before the counters are added up.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct AssayCounts {
    pub assays: i64,
    pub non_filtered_assays: i64,
    pub filtered_assays: i64,
    pub zeroed_out_assays: i64,
    pub calls: i64,
    pub autocall_calls: i64,
    pub no_calls: i64,
}

/// `DUPE`, the filter that is counted as if it had passed.
pub const DUPE_FILTER: &str = "DUPE";
/// `ZEROED_OUT_ASSAY`, the filter that means the chip could not read the assay.
pub const ZEROED_OUT_FILTER: &str = "ZEROED_OUT_ASSAY";
/// `VCFConstants.EMPTY_GENOTYPE`, which the autocall test compares against.
pub const EMPTY_GENOTYPE: &str = "./.";

/// Whether an assay is counted as non-filtered, which a duplicate one is.
pub fn is_counted_as_non_filtered(filters: &[&str]) -> bool {
    filters.is_empty() || filters.contains(&DUPE_FILTER)
}

/// Whether a called genotype is an autocall.
///
/// The comparison is against `./.` and NOT against a dot, and the default the reference asks for
/// is the genotype's own string, so the two spellings do opposite things: a `GTA` of `./.` is not
/// an autocall, and a `GTA` of `.` is one, a single dot being a missing attribute.
pub fn is_autocall(gta: Option<&str>, genotype_string: &str) -> bool {
    gta.unwrap_or(genotype_string) != EMPTY_GENOTYPE
}

/// `updateDetailMetric` as far as the assay counters go.
pub fn assay_counts(
    filters: &[&str],
    called: bool,
    gta: Option<&str>,
    genotype: &str,
) -> AssayCounts {
    let mut counts = AssayCounts {
        assays: 1,
        ..AssayCounts::default()
    };
    if is_counted_as_non_filtered(filters) {
        counts.non_filtered_assays = 1;
        if called {
            counts.calls = 1;
            if is_autocall(gta, genotype) {
                counts.autocall_calls = 1;
            }
        } else {
            counts.no_calls = 1;
        }
    }
    if !filters.is_empty() && !filters.contains(&DUPE_FILTER) {
        counts.filtered_assays = 1;
        if filters.contains(&ZEROED_OUT_FILTER) {
            counts.zeroed_out_assays = 1;
        }
    }
    counts
}

/// `CALL_RATE_PF_THRESHOLD`, whose default is the reference's own, and whose range is checked.
pub const DEFAULT_CALL_RATE_PF_THRESHOLD: f64 = 0.98;

/// The message an out-of-range threshold produces, which is a validation error and not an
/// exception.
pub const CALL_RATE_PF_THRESHOLD_MESSAGE: &str =
    "The parameter CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0";

/// Whether the threshold itself is accepted.
pub fn call_rate_threshold_is_valid(threshold: f64) -> bool {
    threshold > 0.0 && threshold <= 1.0
}

/// `CALL_RATE`, which is over the NON-FILTERED assays and not over every assay.
///
/// The division is by a `float` in the reference, so the value is a float widened to a double
/// rather than a double computed from two longs.
pub fn call_rate(calls: i64, non_filtered_assays: i64) -> f64 {
    f64::from(calls as f32 / non_filtered_assays as f32)
}

/// `HET_PCT`, which is over the CALLS and not over the assays, and is a double throughout.
pub fn het_pct(hets: i64, calls: i64) -> f64 {
    hets as f64 / calls as f64
}

/// `IS_ZCALLED`, which is the presence of the thresholds file and not of the version.
pub fn is_zcalled(zcall_thresholds: Option<&str>) -> bool {
    zcall_thresholds
        .map(|value| !value.is_empty())
        .unwrap_or(false)
}

/// Whether the sample passes, which the threshold decides against the AUTOCALL call rate and not
/// against the call rate the header reported.
pub fn passes_call_rate(autocall_call_rate: f64, threshold: f64) -> bool {
    autocall_call_rate > threshold
}

/// `getSexConcordance`, which is a vote of three and not a comparison of two.
///
/// Three unknowns are a failure. Otherwise it passes when one sex has more than one vote and the
/// other has none, so a sample reported female, fingerprinted female and autocalled unknown
/// passes, and one reported female and autocalled male does not.
pub fn sex_concordance(reported: &str, fingerprint: &str, autocall: &str) -> bool {
    let count = |symbol: &str| {
        [reported, fingerprint, autocall]
            .iter()
            .filter(|value| **value == symbol)
            .count()
    };
    let females = count("F");
    let males = count("M");
    let unknown = count("U") + count("N");
    if unknown == 3 {
        return false;
    }
    (females > 1 && males == 0) || (males > 1 && females == 0)
}
