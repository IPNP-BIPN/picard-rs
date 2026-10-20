//! `CompareGtcFiles` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.illumina.CompareGtcFiles.doWork` at tag 3.4.0. The files are read by
//! `picard_analysis::infinium`; this is the comparison, which the reference makes by reflection
//! over every public zero-argument getter of `InfiniumGTCFile`.
//!
//! # What reflection decides
//!
//! The getters are compared in the order `Class.getMethods` lists them, which is the JVM's and
//! not the declaration's; they are taken here in declaration order, which only shows when two
//! fields differ at once. A getter whose value is null on either side is skipped with a warning,
//! and so is one where exactly one side is the INTEGER 0 or the FLOAT 0.0 (a `long` or a `double`
//! zero never is, being boxed to another class). Float arrays are compared on their decimal
//! digits cut to three places, through `BigDecimal.valueOf(float)`.

use picard_analysis::extract_run::{java_double_to_string, java_float_to_string};
use picard_analysis::infinium::{Bpm, Gtc, ReadError, Reader};
use picard_analysis::metrics_cli::{thrown, Args};
use picard_analysis::vcf_io::log_error;

/// A getter's boxed value.
enum Got<'a> {
    Text(Option<&'a str>),
    Int(i32),
    Long(i64),
    Double(f64),
    Float(f32),
    Ints(Option<&'a [i32]>),
    Floats(Option<&'a [f32]>),
    Bytes(Option<&'a [i8]>),
    Pairs(Option<&'a [[i8; 2]]>),
}

fn getters(g: &Gtc) -> Result<Vec<(&'static str, Got<'_>)>, String> {
    let percentile = |p: &Option<[i32; 3]>, i: usize| -> Result<Got, String> {
        p.map(|v| Got::Int(v[i]))
            .ok_or_else(|| "java.lang.NullPointerException".to_string())
    };
    Ok(vec![
        (
            "getHetPercent",
            Got::Double(f64::from(g.ab_calls) / f64::from(g.num_calls)),
        ),
        ("getSampleName", Got::Text(g.sample_name.as_deref())),
        ("getSamplePlate", Got::Text(g.sample_plate.as_deref())),
        ("getSampleWell", Got::Text(g.sample_well.as_deref())),
        ("getClusterFile", Got::Text(g.cluster_file.as_deref())),
        ("getSnpManifest", Got::Text(g.snp_manifest.as_deref())),
        ("getImagingDate", Got::Text(g.imaging_date.as_deref())),
        ("getAutoCallDate", Got::Text(g.auto_call_date.as_deref())),
        (
            "getAutoCallVersion",
            Got::Text(g.auto_call_version.as_deref()),
        ),
        (
            "getRawControlXIntensities",
            Got::Ints(g.raw_control_x.as_deref()),
        ),
        (
            "getRawControlYIntensities",
            Got::Ints(g.raw_control_y.as_deref()),
        ),
        ("getScannerName", Got::Text(g.scanner_name.as_deref())),
        ("getPmtGreen", Got::Int(g.pmt_green)),
        ("getPmtRed", Got::Int(g.pmt_red)),
        ("getScannerVersion", Got::Text(g.scanner_version.as_deref())),
        ("getImagingUser", Got::Text(g.imaging_user.as_deref())),
        ("getCallRate", Got::Double(g.call_rate)),
        ("getGender", Got::Text(g.gender.as_deref())),
        ("getNumberOfSnps", Got::Int(g.number_of_snps)),
        ("getNumCalls", Got::Int(g.num_calls)),
        ("getNumNoCalls", Got::Int(g.num_no_calls)),
        ("getPloidy", Got::Int(g.ploidy)),
        ("getPloidyType", Got::Int(g.ploidy_type)),
        ("getP05Red", percentile(&g.red_percentiles, 0)?),
        ("getP50Red", percentile(&g.red_percentiles, 1)?),
        ("getP95Red", percentile(&g.red_percentiles, 2)?),
        ("getP05Green", percentile(&g.green_percentiles, 0)?),
        ("getP50Green", percentile(&g.green_percentiles, 1)?),
        ("getP95Green", percentile(&g.green_percentiles, 2)?),
        ("getLogRDev", Got::Float(g.log_r_dev)),
        ("getP10GC", Got::Float(g.p10_gc)),
        ("getP50GC", Got::Float(g.p50_gc)),
        ("getNumIntensityOnly", Got::Int(g.num_intensity_only)),
        ("getAaCalls", Got::Long(g.aa_calls)),
        ("getBbCalls", Got::Long(g.bb_calls)),
        ("getSentrixBarcode", Got::Text(g.sentrix_barcode.as_deref())),
        ("getDx", Got::Int(g.dx)),
        ("getBaseCalls", Got::Pairs(g.base_calls.as_deref())),
        ("getAbCalls", Got::Int(g.ab_calls)),
        ("getRawXIntensities", Got::Ints(g.raw_x.as_deref())),
        ("getRawYIntensities", Got::Ints(g.raw_y.as_deref())),
        (
            "getNormalizedXIntensities",
            Got::Floats(Some(&g.normalized_x)),
        ),
        (
            "getNormalizedYIntensities",
            Got::Floats(Some(&g.normalized_y)),
        ),
        ("getbAlleleFreqs", Got::Floats(g.b_allele_freqs.as_deref())),
        ("getLogRRatios", Got::Floats(g.log_r_ratios.as_deref())),
        ("getRIlmn", Got::Floats(Some(&g.r_ilmn))),
        ("getThetaIlmn", Got::Floats(Some(&g.theta_ilmn))),
        ("getGenotypeBytes", Got::Bytes(g.genotypes.as_deref())),
        (
            "getGenotypeScores",
            Got::Floats(g.genotype_scores.as_deref()),
        ),
        ("getIdentifier", Got::Text(Some(&g.identifier))),
    ])
}

const IGNORED: &[&str] = &[
    "getClass",
    "getAutoCallDate",
    "getImagingDate",
    "getNumberOfEntries",
    "getSampleName",
    "getSamplePlate",
    "getSampleWell",
];

/// `BigDecimal.valueOf(float).setScale(3, ROUND_DOWN)`, as its digits.
fn three_places(value: f32) -> String {
    let text = java_double_to_string(f64::from(value));
    if text.contains('E') {
        // A tiny or huge value: the rounding is on the plain expansion.
        let plain = format!("{:.20}", f64::from(value));
        let (whole, frac) = plain.split_once('.').unwrap_or((&plain, ""));
        return format!("{whole}.{:0<3}", &frac[..frac.len().min(3)]);
    }
    let (whole, frac) = text.split_once('.').unwrap_or((&text, ""));
    format!("{whole}.{:0<3}", &frac[..frac.len().min(3)])
}

/// `arrayDifferences` over two lists of comparable keys.
fn differences<T: PartialEq>(name: &str, a: &[T], b: &[T], errors: &mut Vec<String>) -> usize {
    if a.len() != b.len() {
        errors.push(format!(
            "{name} do not match. Arrays of different lengths. ( {} vs {} )",
            a.len(),
            b.len()
        ));
        return 0;
    }
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

fn float_keys(values: &[f32]) -> Vec<String> {
    values
        .iter()
        .map(|v| {
            if v.is_nan() {
                "NaN".to_string()
            } else {
                three_places(*v)
            }
        })
        .collect()
}

fn compare(one: &Gtc, two: &Gtc) -> Result<Vec<String>, String> {
    let mut errors = Vec::new();
    for ((name, a), (_, b)) in getters(one)?.into_iter().zip(getters(two)?) {
        if IGNORED.contains(&name) {
            continue;
        }
        let mut found: Vec<String> = Vec::new();
        match (a, b) {
            (Got::Text(a), Got::Text(b)) => {
                let (Some(a), Some(b)) = (a, b) else { continue };
                if a != b {
                    found.push(format!("{name} does not match ( {a} vs {b} )"));
                }
            }
            (Got::Int(a), Got::Int(b)) => {
                if (a == 0) != (b == 0) {
                    continue;
                }
                if a != b {
                    found.push(format!("{name} does not match ( {a} vs {b} )"));
                }
            }
            (Got::Long(a), Got::Long(b)) => {
                if a != b {
                    found.push(format!("{name} does not match ( {a} vs {b} )"));
                }
            }
            (Got::Double(a), Got::Double(b)) => {
                if a.to_bits() != b.to_bits() && !(a.is_nan() && b.is_nan()) {
                    found.push(format!(
                        "{name} does not match ( {} vs {} )",
                        java_double_to_string(a),
                        java_double_to_string(b)
                    ));
                }
            }
            (Got::Float(a), Got::Float(b)) => {
                if (a == 0.0 && a.is_sign_positive()) != (b == 0.0 && b.is_sign_positive()) {
                    continue;
                }
                if a.to_bits() != b.to_bits() && !(a.is_nan() && b.is_nan()) {
                    found.push(format!(
                        "{name} does not match ( {} vs {} )",
                        java_float_to_string(a),
                        java_float_to_string(b)
                    ));
                }
            }
            (Got::Ints(Some(a)), Got::Ints(Some(b))) => {
                let n = differences(name, a, b, &mut found);
                if n > 0 {
                    found.push(format!(
                        "{name} do not match. {n} elements of the array differ."
                    ));
                }
            }
            (Got::Bytes(Some(a)), Got::Bytes(Some(b))) => {
                let n = differences(name, a, b, &mut found);
                if n > 0 {
                    found.push(format!(
                        "{name} do not match. {n} elements of the array differ."
                    ));
                }
            }
            (Got::Floats(Some(a)), Got::Floats(Some(b))) => {
                let n = differences(name, &float_keys(a), &float_keys(b), &mut found);
                if n > 0 {
                    found.push(format!(
                        "{name} do not match. {n} elements of the array differ."
                    ));
                }
            }
            (Got::Pairs(Some(a)), Got::Pairs(Some(b))) => {
                if a.len() != b.len() {
                    found.push(format!(
                        "{name} do not match. Arrays of different lengths. ( {} vs {} )",
                        a.len(),
                        b.len()
                    ));
                } else {
                    let mut n = 0;
                    for (x, y) in a.iter().zip(b) {
                        n += differences(name, x, y, &mut found);
                    }
                    if n > 0 {
                        found.push(format!(
                            "{name} do not match. {n} elements of the array differ."
                        ));
                    }
                }
            }
            _ => continue,
        }
        errors.extend(found);
    }
    Ok(errors)
}

fn file_error() -> ! {
    thrown("picard.PicardException: File error: ")
}

fn main() {
    picard_analysis::argument_check::check("CompareGtcFiles");
    picard_analysis::output::install();
    let args = Args::from_env(&[("I", "INPUT"), ("BPM", "ILLUMINA_BEAD_POOL_MANIFEST_FILE")]);
    let inputs = args.collection("INPUT", &[]);
    let bpm_path = args.required("ILLUMINA_BEAD_POOL_MANIFEST_FILE");
    let open = |p: &str| Reader::open(p).unwrap_or_else(|_| file_error());
    let bpm = Bpm::parse(&mut open(&bpm_path)).unwrap_or_else(|_| file_error());
    let parse = |p: &str| -> Gtc {
        match Gtc::parse(&mut open(p), &bpm) {
            Ok(g) => g,
            Err(ReadError::Io(_)) | Err(ReadError::Picard(_)) => file_error(),
        }
    };
    let one = parse(&inputs[0]);
    let two = parse(&inputs[1]);
    let errors = compare(&one, &two).unwrap_or_else(|_| file_error());
    if !errors.is_empty() {
        for e in &errors {
            log_error("CompareGtcFiles", e);
        }
        std::process::exit(1);
    }
}
