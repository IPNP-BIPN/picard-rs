//! `GtcToVcf`: a chip's genotype calls as a VCF.
//!
//! Four files decide one record: the calls themselves, the bead pool manifest that says what each
//! locus is, the cluster file that says how the call was made, and the extended manifest that says
//! where the locus sits on build 37 and which alleles it has there.
//!
//! A call is a letter pair on the chip, `AA`, `AB` or `BB`, and it becomes a genotype against the
//! BUILD's alleles rather than the chip's. Where the reference base is neither of the chip's two
//! alleles, the record is written with both of them as alternates and filtered `TRIALLELIC`.
//!
//! Ported from `picard.arrays.GtcToVcf`, `picard.arrays.illumina.InfiniumGTCFile`,
//! `picard.arrays.illumina.InfiniumTransformation`,
//! `picard.arrays.illumina.Build37ExtendedIlluminaManifestRecord` and
//! `picard.arrays.illumina.InfiniumVcfFields` in Picard 3.4.0.

use crate::liftover_vcf::format_vcf_double;

/// A call file's genotype codes.
pub const NO_CALL: u8 = 0;
pub const AA_CALL: u8 = 1;
pub const AB_CALL: u8 = 2;
pub const BB_CALL: u8 = 3;

/// The filters the tool can add.
pub const DUPE: &str = "DUPE";
pub const TRIALLELIC: &str = "TRIALLELIC";
pub const ZEROED_OUT_ASSAY: &str = "ZEROED_OUT_ASSAY";

/// What the extended manifest says happened to a locus when it was built.
///
/// The two that survive are `PASS` and `DUPE`: a duplicate is written out with a filter, and every
/// other flag drops the locus from the file entirely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flag {
    Pass,
    Dupe,
    IlluminaFlagged,
    LiftoverFailed,
    UnsupportedGenomeBuild,
    ProbeSequenceMismatch,
    MissingAlleleBProbeseq,
    SourceSequenceInvalid,
    IndelNotFound,
    IndelConflict,
}

impl Flag {
    pub fn parse(name: &str) -> Option<Flag> {
        Some(match name {
            "PASS" => Flag::Pass,
            "DUPE" => Flag::Dupe,
            "ILLUMINA_FLAGGED" => Flag::IlluminaFlagged,
            "LIFTOVER_FAILED" => Flag::LiftoverFailed,
            "UNSUPPORTED_GENOME_BUILD" => Flag::UnsupportedGenomeBuild,
            "PROBE_SEQUENCE_MISMATCH" => Flag::ProbeSequenceMismatch,
            "MISSING_ALLELE_B_PROBESEQ" => Flag::MissingAlleleBProbeseq,
            "SOURCE_SEQUENCE_INVALID" => Flag::SourceSequenceInvalid,
            "INDEL_NOT_FOUND" => Flag::IndelNotFound,
            "INDEL_CONFLICT" => Flag::IndelConflict,
            _ => return None,
        })
    }

    /// A flagged locus is one that failed, and a duplicate did not.
    pub fn is_fail(&self) -> bool {
        !matches!(self, Flag::Pass | Flag::Dupe)
    }

    pub fn is_dupe(&self) -> bool {
        matches!(self, Flag::Dupe)
    }
}

/// One row of the extended manifest, as far as the VCF needs it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ManifestRecord {
    pub name: String,
    /// The chip's own contig and position, which are written to the record as they were.
    pub chr: String,
    pub position: i32,
    pub genome_build: String,
    pub b37_chr: String,
    pub b37_pos: i32,
    pub ref_allele: String,
    pub allele_a: String,
    pub allele_b: String,
    pub rs_id: String,
    pub ilmn_strand: String,
    pub probe_a: String,
    pub probe_b: String,
    pub bead_set_id: i32,
    pub source: String,
    pub flag: Flag,
}

/// One locus's cluster, from the `.egt`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cluster {
    pub total_score: f32,
    /// The counts of the three genotypes in the training set.
    pub n: [i32; 3],
    pub dev_r: [f32; 3],
    pub mean_r: [f32; 3],
    pub dev_theta: [f32; 3],
    pub mean_theta: [f32; 3],
}

/// One locus's call, from the `.gtc`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Call {
    pub genotype: u8,
    pub score: f32,
    pub raw_x: i32,
    pub raw_y: i32,
    pub b_allele_freq: f32,
    pub log_r_ratio: f32,
}

/// The normalization one locus's intensities go through.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transformation {
    pub offset_x: f32,
    pub offset_y: f32,
    pub scale_x: f32,
    pub scale_y: f32,
    pub shear: f32,
    pub theta: f32,
}

/// The intensities the record reports, which are not the ones the chip measured.
///
/// The offsets come off first, then the rotation, then the shear, and a negative result is taken
/// to nought before either is scaled.
pub fn normalize(raw_x: i32, raw_y: i32, transformation: &Transformation) -> (f32, f32) {
    let temp_x = raw_x as f32 - transformation.offset_x;
    let temp_y = raw_y as f32 - transformation.offset_y;
    let theta = f64::from(transformation.theta);
    let temp_x2 = theta.cos() * f64::from(temp_x) + theta.sin() * f64::from(temp_y);
    let temp_y2 = -theta.sin() * f64::from(temp_x) + theta.cos() * f64::from(temp_y);
    let mut temp_x3 = temp_x2 - f64::from(transformation.shear) * temp_y2;
    let temp_y2 = temp_y2.max(0.0);
    if temp_x3 < 0.0 {
        temp_x3 = 0.0;
    }
    (
        (temp_x3 / f64::from(transformation.scale_x)) as f32,
        (temp_y2 / f64::from(transformation.scale_y)) as f32,
    )
}

/// The chip's two summary numbers: the total signal, and how it is split between the two alleles.
///
/// The distance is Manhattan rather than Euclidean, which is why the total is the sum of the two
/// intensities and not the hypotenuse.
pub fn r_and_theta(normalized_x: f32, normalized_y: f32) -> (f32, f32) {
    let x = f64::from(normalized_x);
    let y = f64::from(normalized_y);
    let theta = 2.0 * ((y / x).atan() / std::f64::consts::PI);
    ((x + y) as f32, theta as f32)
}

/// A cluster's position in the intensity plane, converted from the polar form the file holds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EuclideanValues {
    pub mean_x: f32,
    pub mean_y: f32,
    pub dev_x: f32,
    pub dev_y: f32,
}

/// `polarToEuclidean`, the same Manhattan conversion with the deviations propagated through it.
pub fn polar_to_euclidean(r: f32, dev_r: f32, theta: f32, dev_theta: f32) -> EuclideanValues {
    let theta_variance = f64::from(dev_theta).powi(2);
    let r_variance = f64::from(dev_r).powi(2);
    let half_pi = std::f64::consts::PI / 2.0;
    let normalized_theta = half_pi * f64::from(theta);
    let r_over_x = 1.0 + normalized_theta.tan();

    let theta_variance_factor_x =
        -(half_pi * f64::from(r)) * (r_over_x * normalized_theta.cos()).powi(-2);
    let r_variance_factor_x = 1.0 / r_over_x;
    let variance_x = (theta_variance_factor_x.powi(2) * theta_variance)
        + (r_variance_factor_x.powi(2) * r_variance);
    let theta_variance_factor_y = -theta_variance_factor_x;
    let r_variance_factor_y = 1.0 - r_variance_factor_x;
    let variance_y = (theta_variance_factor_y.powi(2) * theta_variance)
        + (r_variance_factor_y.powi(2) * r_variance);

    let mean_x = f64::from(r) / r_over_x;
    let mean_y = f64::from(r) - mean_x;
    EuclideanValues {
        mean_x: mean_x as f32,
        mean_y: mean_y as f32,
        dev_x: variance_x.powf(0.5) as f32,
        dev_y: variance_y.powf(0.5) as f32,
    }
}

/// `formatFloatForVcf`: three decimal places at most, no grouping, and a dot for anything that is
/// not a number.
pub fn format_float_for_vcf(value: f32) -> String {
    if value.is_nan() || value.is_infinite() {
        return ".".to_string();
    }
    let scaled = f64::from(value) * 1000.0;
    // Java's DecimalFormat rounds half to even.
    let rounded = round_half_even(scaled) / 1000.0;
    let text = format!("{rounded:.3}");
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

fn round_half_even(value: f64) -> f64 {
    let rounded = value.round();
    if (value - value.trunc()).abs() == 0.5 && rounded % 2.0 != 0.0 {
        rounded - value.signum()
    } else {
        rounded
    }
}

/// The alleles of one record: the reference first, then whichever of the chip's two differ from
/// it.
pub fn assay_alleles(record: &ManifestRecord) -> Vec<String> {
    let mut alleles = vec![record.ref_allele.clone()];
    if record.allele_a != record.ref_allele {
        alleles.push(record.allele_a.clone());
    }
    if record.allele_b != record.ref_allele {
        alleles.push(record.allele_b.clone());
    }
    alleles
}

/// The sample's two alleles, from the chip's call.
pub fn called_alleles(record: &ManifestRecord, call: &Call) -> Option<[String; 2]> {
    match call.genotype {
        NO_CALL => None,
        AA_CALL => Some([record.allele_a.clone(), record.allele_a.clone()]),
        AB_CALL => Some([record.allele_a.clone(), record.allele_b.clone()]),
        BB_CALL => Some([record.allele_b.clone(), record.allele_b.clone()]),
        _ => None,
    }
}

/// The genotype as the file writes it: the indices of the called alleles among the record's own.
pub fn genotype_field(record: &ManifestRecord, call: &Call) -> String {
    let Some(called) = called_alleles(record, call) else {
        return "./.".to_string();
    };
    let alleles = assay_alleles(record);
    let index = |allele: &String| {
        alleles
            .iter()
            .position(|candidate| candidate == allele)
            .expect("a called allele is one of the record's")
            .to_string()
    };
    format!("{}/{}", index(&called[0]), index(&called[1]))
}

/// The allele counts a record carries: how many of each alternate were called, and how many
/// alleles were called at all.
pub fn chromosome_counts(record: &ManifestRecord, call: &Call) -> (Vec<i32>, i32) {
    let alleles = assay_alleles(record);
    let mut counts = vec![0; alleles.len().saturating_sub(1)];
    let mut total = 0;
    if let Some(called) = called_alleles(record, call) {
        for allele in called {
            total += 1;
            if let Some(position) = alleles.iter().position(|candidate| *candidate == allele) {
                if position > 0 {
                    counts[position - 1] += 1;
                }
            }
        }
    }
    (counts, total)
}

/// One INFO entry.
fn attribute(key: &str, value: String) -> (String, String) {
    (key.to_string(), value)
}

/// The three genotypes' suffixes, in the order the cluster file holds them.
pub const GENOTYPE_VALUES: [&str; 3] = ["AA", "AB", "BB"];

/// The whole INFO field of one record, sorted the way the encoder sorts it.
pub fn info_fields(
    record: &ManifestRecord,
    call: &Call,
    cluster: &Cluster,
) -> Vec<(String, String)> {
    let (counts, total) = chromosome_counts(record, call);
    let mut fields = Vec::new();
    if !counts.is_empty() {
        fields.push(attribute(
            "AC",
            counts
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join(","),
        ));
        fields.push(attribute(
            "AF",
            counts
                .iter()
                .map(|count| {
                    format_vcf_double(if total == 0 {
                        0.0
                    } else {
                        f64::from(*count) / f64::from(total)
                    })
                })
                .collect::<Vec<_>>()
                .join(","),
        ));
    }
    fields.push(attribute("AN", total.to_string()));
    // The allele that is the reference is written with a star, which is how the file says which of
    // the chip's two the build agrees with.
    let star = |allele: &String| {
        if *allele == record.ref_allele {
            format!("{allele}*")
        } else {
            allele.clone()
        }
    };
    fields.push(attribute("ALLELE_A", star(&record.allele_a)));
    fields.push(attribute("ALLELE_B", star(&record.allele_b)));
    fields.push(attribute("BEADSET_ID", record.bead_set_id.to_string()));
    fields.push(attribute(
        "GC_SCORE",
        format_float_for_vcf(cluster.total_score),
    ));
    fields.push(attribute("ILLUMINA_BUILD", record.genome_build.clone()));
    fields.push(attribute("ILLUMINA_CHR", record.chr.clone()));
    fields.push(attribute("ILLUMINA_POS", record.position.to_string()));
    fields.push(attribute("ILLUMINA_STRAND", record.ilmn_strand.clone()));
    for (ordinal, suffix) in GENOTYPE_VALUES.iter().enumerate() {
        fields.push(attribute(
            &format!("N_{suffix}"),
            cluster.n[ordinal].to_string(),
        ));
    }
    fields.push(attribute("PROBE_A", record.probe_a.clone()));
    fields.push(attribute("PROBE_B", record.probe_b.clone()));
    // A source with semicolons in it would end the INFO field early, so they become commas and
    // spaces become underscores.
    fields.push(attribute(
        "SOURCE",
        record.source.replace(';', ",").replace(' ', "_"),
    ));
    for (ordinal, suffix) in GENOTYPE_VALUES.iter().enumerate() {
        let euclidean = polar_to_euclidean(
            cluster.mean_r[ordinal],
            cluster.dev_r[ordinal],
            cluster.mean_theta[ordinal],
            cluster.dev_theta[ordinal],
        );
        for (key, value) in [
            (format!("devR_{suffix}"), cluster.dev_r[ordinal]),
            (format!("devTHETA_{suffix}"), cluster.dev_theta[ordinal]),
            (format!("devX_{suffix}"), euclidean.dev_x),
            (format!("devY_{suffix}"), euclidean.dev_y),
            (format!("meanR_{suffix}"), cluster.mean_r[ordinal]),
            (format!("meanTHETA_{suffix}"), cluster.mean_theta[ordinal]),
            (format!("meanX_{suffix}"), euclidean.mean_x),
            (format!("meanY_{suffix}"), euclidean.mean_y),
        ] {
            fields.push(attribute(&key, format_float_for_vcf(value)));
        }
    }
    if !record.rs_id.is_empty() {
        fields.push(attribute("refSNP", record.rs_id.clone()));
    }
    fields.sort_by(|left, right| left.0.cmp(&right.0));
    fields
}

/// The FORMAT keys a record carries, in the order they are written.
pub const FORMAT: [&str; 10] = [
    "GT", "BAF", "IGC", "LRR", "NORMX", "NORMY", "R", "THETA", "X", "Y",
];

/// The sample's column.
pub fn sample_field(
    record: &ManifestRecord,
    call: &Call,
    transformation: &Transformation,
) -> String {
    let (normalized_x, normalized_y) = normalize(call.raw_x, call.raw_y, transformation);
    let (r, theta) = r_and_theta(normalized_x, normalized_y);
    [
        genotype_field(record, call),
        format_float_for_vcf(call.b_allele_freq),
        format_float_for_vcf(call.score),
        format_float_for_vcf(call.log_r_ratio),
        format_float_for_vcf(normalized_x),
        format_float_for_vcf(normalized_y),
        format_float_for_vcf(r),
        format_float_for_vcf(theta),
        call.raw_x.to_string(),
        call.raw_y.to_string(),
    ]
    .join(":")
}

/// The filters one record carries.
///
/// The triallelic filter is added when the record is written rather than when it is built, which
/// is why it is decided by the alternates and not by the manifest.
pub fn filters(record: &ManifestRecord, cluster: &Cluster) -> Vec<String> {
    let mut filters = Vec::new();
    if cluster.total_score == 0.0 {
        filters.push(ZEROED_OUT_ASSAY.to_string());
    }
    if record.flag.is_dupe() {
        filters.push(DUPE.to_string());
    }
    if assay_alleles(record).len() > 2 {
        filters.push(TRIALLELIC.to_string());
    }
    filters
}

/// One variant line.
pub fn variant_line(
    record: &ManifestRecord,
    call: &Call,
    cluster: &Cluster,
    transformation: &Transformation,
) -> String {
    let alleles = assay_alleles(record);
    let info = info_fields(record, call, cluster)
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<_>>()
        .join(";");
    let filters = filters(record, cluster);
    format!(
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t{}",
        record.b37_chr,
        record.b37_pos,
        record.name,
        alleles[0],
        alleles[1..].join(","),
        if filters.is_empty() {
            ".".to_string()
        } else {
            filters.join(";")
        },
        info,
        FORMAT.join(":"),
        sample_field(record, call, transformation)
    )
}

/// Every record of one run, in the order the file writes them.
///
/// A locus the manifest flagged is not written at all; a duplicate is, with a filter. The output
/// is sorted by the TARGET build's coordinates, which is not the order the chip lists its loci in.
pub fn records(
    manifest: &[ManifestRecord],
    calls: &[Call],
    clusters: &[Cluster],
    transformations: &[Transformation],
) -> Vec<String> {
    let mut rows: Vec<(String, i32, String)> = manifest
        .iter()
        .enumerate()
        .filter(|(_, record)| !record.flag.is_fail())
        .map(|(index, record)| {
            (
                record.b37_chr.clone(),
                record.b37_pos,
                variant_line(
                    record,
                    &calls[index],
                    &clusters[index],
                    &transformations[index],
                ),
            )
        })
        .collect();
    rows.sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    rows.into_iter().map(|(_, _, line)| line).collect()
}

/// The header lines the command line decides, in the order the file sorts them.
pub fn header_lines(
    analysis_version_number: Option<i32>,
    pipeline_version: Option<&str>,
    sample_alias: &str,
    sample_column: &str,
) -> Vec<String> {
    let mut lines = vec!["##fileformat=VCFv4.2".to_string()];
    let mut declared = Vec::new();
    if let Some(version) = analysis_version_number {
        declared.push(format!("##analysisVersionNumber={version}"));
    }
    if let Some(version) = pipeline_version {
        declared.push(format!("##pipelineVersion={version}"));
    }
    declared.push(format!("##sampleAlias={sample_alias}"));
    declared.sort();
    lines.extend(declared);
    lines.push(format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_column}"
    ));
    lines
}
