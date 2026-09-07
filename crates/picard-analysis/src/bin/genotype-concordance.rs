//! `GenotypeConcordance` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.GenotypeConcordance.doWork` at tag 3.4.0. The states, the GA4GH schemes, the
//! contingency arithmetic and the allele normalisation live in
//! `picard_analysis::genotype_concordance`.
//!
//! The tool walks two VCFs in tandem, one sample from each, files every site under a pair of
//! states, and writes three metrics files from one basename: a summary of ratios, the whole
//! contingency table row by row, and the four counters summed.
//!
//! Three things about it are easy to miss.
//!
//! A site the other file does not have is not skipped: it is `MISSING`, which is a state like any
//! other and contributes to the table. `MISSING_SITES_HOM_REF` goes further and adds the WHOLE
//! uncovered region -- the interval list's base count minus the sites seen -- to the
//! missing/missing cell, so the counts stop being a number of variants.
//!
//! Each side is subset to its own sample before anything is read, and the subset REDERIVES the
//! alleles: a site with two samples where this one is hom-ref keeps only the reference allele, so
//! its type becomes `NO_VARIATION` and it is no longer a SNP at all.
//!
//! `IGNORE_FILTER_STATUS` does not skip a filtered call. It clears the state the filter would have
//! given, so the call is read from its alleles like any other -- and a filtered TRUTH site has no
//! such escape.
//!
//! Not ported: `OUTPUT_VCF`, which writes a fourth file holding both samples' genotypes and the
//! contingency state as an INFO attribute. The array compares the summary metrics, so a row that
//! asks for it is measured on everything but that file.

use std::collections::HashMap;

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use picard_analysis::genotype_concordance::{
    contingency, contingency_string, determine_state, file_names, is_var, CallState, Cell, Counts,
    GenotypeView, SiteView, TruthState, CALL_DECLARATION_ORDER, GA4GH, GA4GH_MISSING_AS_HOM_REF,
    HET_CALL_STATES, HET_TRUTH_STATES, HOM_VAR_CALL_STATES, HOM_VAR_TRUTH_STATES,
    TRUTH_DECLARATION_ORDER, VAR_CALL_STATES, VAR_TRUTH_STATES,
};

/// `VariantContext.Type`, as far as the counter cares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantType {
    NoVariation,
    Snp,
    Indel,
    Mixed,
    Other,
}

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn throw(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
    std::process::exit(1);
}

/// One interval of an interval list, which is all this tool reads out of one.
#[derive(Debug, Clone)]
struct Interval {
    contig: String,
    start: i64,
    end: i64,
}

/// `IntervalList.fromPath` followed by `uniqued()`: the intervals sorted by the sequence
/// dictionary's order and merged where they touch.
///
/// The merge is what the base count is taken over, and the base count is what
/// `MISSING_SITES_HOM_REF` fills the missing/missing cell from, so two intervals that abut are one
/// interval and not two.
fn read_intervals(text: &str) -> (Vec<String>, Vec<Interval>) {
    let mut contigs = Vec::new();
    let mut intervals = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("@SQ\t") {
            if let Some(name) = rest.split('\t').find_map(|f| f.strip_prefix("SN:")) {
                contigs.push(name.to_string());
            }
            continue;
        }
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        intervals.push(Interval {
            contig: fields[0].to_string(),
            start: fields[1].parse().unwrap_or(0),
            end: fields[2].parse().unwrap_or(0),
        });
    }
    let order = |contig: &str| {
        contigs
            .iter()
            .position(|name| name == contig)
            .unwrap_or(usize::MAX)
    };
    intervals.sort_by(|a, b| {
        order(&a.contig)
            .cmp(&order(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    let mut merged: Vec<Interval> = Vec::new();
    for interval in intervals {
        match merged.last_mut() {
            // `IntervalList.getUniqueIntervals`: overlapping OR abutting intervals become one.
            Some(last) if last.contig == interval.contig && interval.start <= last.end + 1 => {
                last.end = last.end.max(interval.end);
            }
            _ => merged.push(interval),
        }
    }
    (contigs, merged)
}

/// `ByIntervalListVariantContextIterator`: the records of each interval in turn, minus the ones
/// the PREVIOUS interval already returned.
fn within_intervals(records: &[VariantContext], intervals: &[Interval]) -> Vec<VariantContext> {
    let overlaps = |record: &VariantContext, interval: &Interval| {
        record.contig == interval.contig
            && record.start <= interval.end
            && record.stop >= interval.start
    };
    let mut kept = Vec::new();
    let mut previous: Option<&Interval> = None;
    for interval in intervals {
        for record in records {
            if overlaps(record, interval) && !previous.is_some_and(|last| overlaps(record, last)) {
                kept.push(record.clone());
            }
        }
        previous = Some(interval);
    }
    kept
}

/// `subContextFromSample`, which is a subset AND a rederivation.
///
/// The alleles are rebuilt from the one genotype's called alleles, with the reference added back
/// if the genotype did not call it, and then filtered through the record's own allele order. A
/// hom-ref genotype therefore leaves a one-allele site, whose type is `NO_VARIATION`.
fn sub_context(record: &VariantContext, sample: &str) -> Option<(VariantType, SiteView)> {
    let genotype = record
        .genotypes
        .iter()
        .find(|genotype| genotype.sample_name == sample)?;

    let mut called: Vec<&htsjdk_vcf::allele::Allele> = Vec::new();
    let mut added_reference = false;
    for allele in &genotype.alleles {
        added_reference = added_reference || allele.is_reference();
        if !allele.is_no_call() {
            called.push(allele);
        }
    }
    let reference = record
        .alleles
        .iter()
        .find(|allele| allele.is_reference())
        .expect("a reference allele");
    let mut kept: Vec<&htsjdk_vcf::allele::Allele> = record
        .alleles
        .iter()
        .filter(|allele| called.iter().any(|called| *called == *allele))
        .collect();
    if !added_reference && !kept.iter().any(|allele| allele.is_reference()) {
        // `allelesOfGenotypes` adds the reference when no genotype allele was the reference, and
        // the rederived list keeps the record's own order, so it lands first.
        kept.insert(0, reference);
    }

    let variant_type = variant_type(&kept, reference);
    Some((
        variant_type,
        SiteView {
            is_mixed: variant_type == VariantType::Mixed,
            // `VariantContext.isFiltered()`: filters applied and not empty. A record that never
            // had filters applied prints `.` and is not filtered.
            is_filtered: record.filters.as_ref().is_some_and(|f| !f.is_empty()),
            reference: reference.base_string(),
            genotype: GenotypeView {
                alleles: genotype
                    .alleles
                    .iter()
                    .map(|allele| allele.base_string())
                    .collect(),
                is_filtered: genotype.filters.as_ref().is_some_and(|f| !f.is_empty()),
                gq: genotype.gq.unwrap_or(-1),
                dp: genotype.dp.unwrap_or(-1),
            },
        },
    ))
}

/// `determineType`: one allele is no variation, and alternates that disagree are a mixed site.
fn variant_type(
    alleles: &[&htsjdk_vcf::allele::Allele],
    reference: &htsjdk_vcf::allele::Allele,
) -> VariantType {
    if alleles.len() <= 1 {
        return VariantType::NoVariation;
    }
    let mut found: Option<VariantType> = None;
    for allele in alleles {
        if allele.is_reference() {
            continue;
        }
        let biallelic = if allele.is_symbolic() {
            VariantType::Other
        } else if reference.len() == allele.len() {
            if allele.len() == 1 {
                VariantType::Snp
            } else {
                VariantType::Other // MNP
            }
        } else {
            VariantType::Indel
        };
        match found {
            None => found = Some(biallelic),
            Some(previous) if previous != biallelic => return VariantType::Mixed,
            Some(_) => {}
        }
    }
    found.unwrap_or(VariantType::NoVariation)
}

struct SummaryMetrics {
    variant_type: &'static str,
    truth_sample: String,
    call_sample: String,
    values: [f64; 11],
}

const SUMMARY_COLUMNS: [&str; 14] = [
    "VARIANT_TYPE",
    "TRUTH_SAMPLE",
    "CALL_SAMPLE",
    "HET_SENSITIVITY",
    "HET_PPV",
    "HET_SPECIFICITY",
    "HOMVAR_SENSITIVITY",
    "HOMVAR_PPV",
    "HOMVAR_SPECIFICITY",
    "VAR_SENSITIVITY",
    "VAR_PPV",
    "VAR_SPECIFICITY",
    "GENOTYPE_CONCORDANCE",
    "NON_REF_GENOTYPE_CONCORDANCE",
];

impl MetricBean for SummaryMetrics {
    fn class_name(&self) -> &str {
        "picard.vcf.GenotypeConcordanceSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &SUMMARY_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let mut out = vec![
            Value::Str(self.variant_type.to_string()),
            Value::Str(self.truth_sample.clone()),
            Value::Str(self.call_sample.clone()),
        ];
        out.extend(self.values.iter().map(|value| Value::Double(*value)));
        out
    }
}

struct DetailMetrics {
    variant_type: &'static str,
    truth_sample: String,
    call_sample: String,
    truth_state: TruthState,
    call_state: CallState,
    count: i64,
    contingency_values: String,
}

const DETAIL_COLUMNS: [&str; 7] = [
    "VARIANT_TYPE",
    "TRUTH_SAMPLE",
    "CALL_SAMPLE",
    "TRUTH_STATE",
    "CALL_STATE",
    "COUNT",
    "CONTINGENCY_VALUES",
];

impl MetricBean for DetailMetrics {
    fn class_name(&self) -> &str {
        "picard.vcf.GenotypeConcordanceDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &DETAIL_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.variant_type.to_string()),
            Value::Str(self.truth_sample.clone()),
            Value::Str(self.call_sample.clone()),
            Value::Str(self.truth_state.name().to_string()),
            Value::Str(self.call_state.name().to_string()),
            Value::Long(self.count),
            Value::Str(self.contingency_values.clone()),
        ]
    }
}

struct ContingencyMetrics {
    variant_type: &'static str,
    truth_sample: String,
    call_sample: String,
    counts: [i64; 5],
}

const CONTINGENCY_COLUMNS: [&str; 8] = [
    "VARIANT_TYPE",
    "TRUTH_SAMPLE",
    "CALL_SAMPLE",
    "TP_COUNT",
    "TN_COUNT",
    "FP_COUNT",
    "FN_COUNT",
    "EMPTY_COUNT",
];

impl MetricBean for ContingencyMetrics {
    fn class_name(&self) -> &str {
        "picard.vcf.GenotypeConcordanceContingencyMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &CONTINGENCY_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let mut out = vec![
            Value::Str(self.variant_type.to_string()),
            Value::Str(self.truth_sample.clone()),
            Value::Str(self.call_sample.clone()),
        ];
        out.extend(self.counts.iter().map(|count| Value::Long(*count)));
        out
    }
}

/// The sequence dictionary a VCF header carries, as `(name, length)` in `##contig` order.
///
/// It is the whole genome's length that `MISSING_SITES_HOM_REF` counts against when no interval
/// list was given, so a missing `length` is a zero rather than a guess.
fn dictionary(header: &htsjdk_vcf::header::VcfHeader) -> Vec<(String, i64)> {
    header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Contig { fields, .. } => {
                let id = fields.iter().find(|(key, _)| key == "ID")?.1.clone();
                let length = fields
                    .iter()
                    .find(|(key, _)| key == "length")
                    .and_then(|(_, value)| value.parse().ok())
                    .unwrap_or(0);
                Some((id, length))
            }
            _ => None,
        })
        .collect()
}

fn metrics_file<B: MetricBean>(beans: Vec<B>) -> String {
    let mut file = MetricsFile::new();
    file.add_header("GenotypeConcordance <command line>");
    file.add_header("Started on: <timestamp>");
    for bean in beans {
        file.add_metric(&bean);
    }
    file.write()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let truth_vcf = arg(&args, "TRUTH_VCF=").ok_or("TRUTH_VCF= is required")?;
    let call_vcf = arg(&args, "CALL_VCF=").ok_or("CALL_VCF= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let flag = |key: &str, default: bool| arg(&args, key).map(|v| v == "true").unwrap_or(default);
    let number = |key: &str, default: i32| -> i32 {
        arg(&args, key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let min_gq = number("MIN_GQ=", 0);
    let min_dp = number("MIN_DP=", 0);
    let output_all_rows = flag("OUTPUT_ALL_ROWS=", false);
    let missing_sites_hom_ref = flag("MISSING_SITES_HOM_REF=", false);
    let ignore_filter_status = flag("IGNORE_FILTER_STATUS=", false);

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let truth = read_vcf(&std::fs::read_to_string(&truth_vcf)?).map_err(|e| format!("{e:?}"))?;
    let call = read_vcf(&std::fs::read_to_string(&call_vcf)?).map_err(|e| format!("{e:?}"))?;

    // `TRUTH_SAMPLE is required when the TRUTH_VCF has more than one sample`, and the same for the
    // call side: a one-sample file names its own sample, and a two-sample file refuses.
    let sample_of = |given: Option<String>,
                     file: &htsjdk_vcf::reader::VcfFile,
                     which: &str,
                     path: &str|
     -> String {
        match given {
            Some(name) => {
                if !file.header.samples.contains(&name) {
                    throw(&format!(
                        "File {path} does not contain genotypes for sample {name}"
                    ));
                }
                name
            }
            None => {
                if file.header.samples.len() > 1 {
                    throw(&format!(
                        "{which}_SAMPLE is required when the {which}_VCF has more than one sample"
                    ));
                }
                file.header.samples.first().cloned().unwrap_or_default()
            }
        }
    };
    let truth_sample = sample_of(arg(&args, "TRUTH_SAMPLE="), &truth, "TRUTH", &truth_vcf);
    let call_sample = sample_of(arg(&args, "CALL_SAMPLE="), &call, "CALL", &call_vcf);

    // The intervals decide two things: which sites are compared, and -- through the base count --
    // how many missing/missing sites `MISSING_SITES_HOM_REF` adds.
    let intervals: Vec<String> = args
        .iter()
        .filter_map(|a| a.strip_prefix("INTERVALS=").map(str::to_string))
        .collect();
    let (contigs, merged) = if intervals.is_empty() {
        (Vec::new(), Vec::new())
    } else {
        // One list is taken as it lies; INTERSECT_INTERVALS only decides how a SECOND list would
        // be combined with the first.
        read_intervals(&std::fs::read_to_string(&intervals[0])?)
    };
    let truth_records = if merged.is_empty() {
        truth.records.clone()
    } else {
        within_intervals(&truth.records, &merged)
    };
    let call_records = if merged.is_empty() {
        call.records.clone()
    } else {
        within_intervals(&call.records, &merged)
    };

    // `VariantContextComparator`: the sequence dictionary's order, then the start.
    let contig_order: HashMap<String, usize> = if contigs.is_empty() {
        dictionary(&truth.header)
            .into_iter()
            .enumerate()
            .map(|(index, (name, _))| (name, index))
            .collect()
    } else {
        contigs
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), index))
            .collect()
    };
    let key = |record: &VariantContext| {
        (
            *contig_order.get(&record.contig).unwrap_or(&usize::MAX),
            record.start,
        )
    };

    let mut snp_counter = Counts::new();
    let mut indel_counter = Counts::new();
    let (mut left, mut right) = (0usize, 0usize);
    while left < truth_records.len() || right < call_records.len() {
        // `PairedVariantSubContextIterator`: whichever side is behind is taken alone, and a tie
        // takes both.
        let (truth_site, call_site) = match (truth_records.get(left), call_records.get(right)) {
            (Some(l), Some(r)) => match key(l).cmp(&key(r)) {
                std::cmp::Ordering::Equal => {
                    left += 1;
                    right += 1;
                    (Some(l), Some(r))
                }
                std::cmp::Ordering::Less => {
                    left += 1;
                    (Some(l), None)
                }
                std::cmp::Ordering::Greater => {
                    right += 1;
                    (None, Some(r))
                }
            },
            (Some(l), None) => {
                left += 1;
                (Some(l), None)
            }
            (None, Some(r)) => {
                right += 1;
                (None, Some(r))
            }
            (None, None) => break,
        };

        let truth_view = truth_site.and_then(|record| sub_context(record, &truth_sample));
        let call_view = call_site.and_then(|record| sub_context(record, &call_sample));
        let truth_type = truth_view
            .as_ref()
            .map_or(VariantType::NoVariation, |(kind, _)| *kind);
        let call_type = call_view
            .as_ref()
            .map_or(VariantType::NoVariation, |(kind, _)| *kind);

        let states = determine_state(
            truth_view.as_ref().map(|(_, view)| view),
            call_view.as_ref().map(|(_, view)| view),
            min_gq,
            min_dp,
            ignore_filter_status,
        );
        let (truth_state, call_state) = match states {
            Ok(states) => states,
            Err(message) => throw(&message),
        };

        // `classifyVariants`: which counter a pair lands in, and the pairs that land in neither.
        // A MIXED site is filed under the OTHER side's type, which is why a truth SNP against a
        // mixed call is a SNP row.
        let counter = match (truth_type, call_type) {
            (
                VariantType::Snp,
                VariantType::Snp | VariantType::Mixed | VariantType::NoVariation,
            ) => Some(&mut snp_counter),
            (
                VariantType::Indel,
                VariantType::Indel | VariantType::Mixed | VariantType::NoVariation,
            ) => Some(&mut indel_counter),
            (VariantType::Mixed, VariantType::Snp) => Some(&mut snp_counter),
            (VariantType::Mixed, VariantType::Indel) => Some(&mut indel_counter),
            (VariantType::NoVariation, VariantType::Snp) => Some(&mut snp_counter),
            (VariantType::NoVariation, VariantType::Indel) => Some(&mut indel_counter),
            _ => None,
        };
        if let Some(counter) = counter {
            counter.increment(truth_state, call_state);
        }
    }

    if missing_sites_hom_ref {
        // The whole region minus what was seen, as a count of sites that are missing on both
        // sides. It is a difference of doubles and can be negative if there were more variants
        // than bases, which the reference does not guard against either.
        let base_count = if merged.is_empty() {
            dictionary(&truth.header)
                .iter()
                .map(|(_, length)| *length)
                .sum::<i64>()
        } else {
            merged.iter().map(|i| i.end - i.start + 1).sum::<i64>()
        } as f64;
        let snp_seen = snp_counter.size();
        let indel_seen = indel_counter.size();
        snp_counter.increment_by(
            TruthState::Missing,
            CallState::Missing,
            base_count - snp_seen,
        );
        indel_counter.increment_by(
            TruthState::Missing,
            CallState::Missing,
            base_count - indel_seen,
        );
    }

    let scheme: &[(CallState, [Cell; 11])] = if missing_sites_hom_ref {
        &GA4GH_MISSING_AS_HOM_REF
    } else {
        &GA4GH
    };

    for (kind, counter) in [("SNP", &snp_counter), ("INDEL", &indel_counter)] {
        if let Err((truth_state, call_state)) = counter.validate_against(scheme) {
            throw(&format!(
                "Found counts for an illegal set of states: [{}, {}]",
                truth_state.name(),
                call_state.name()
            ));
        }
        let _ = kind;
    }

    let summary: Vec<SummaryMetrics> = [("SNP", &snp_counter), ("INDEL", &indel_counter)]
        .into_iter()
        .map(|(kind, counter)| SummaryMetrics {
            variant_type: kind,
            truth_sample: truth_sample.clone(),
            call_sample: call_sample.clone(),
            values: [
                counter.sensitivity(scheme, &HET_TRUTH_STATES),
                counter.ppv(scheme, &HET_CALL_STATES),
                // "The specificity for all heterozygous variants cannot be calculated."
                f64::NAN,
                counter.sensitivity(scheme, &HOM_VAR_TRUTH_STATES),
                counter.ppv(scheme, &HOM_VAR_CALL_STATES),
                f64::NAN,
                counter.sensitivity(scheme, &VAR_TRUTH_STATES),
                counter.ppv(scheme, &VAR_CALL_STATES),
                counter.specificity(scheme, &VAR_TRUTH_STATES),
                counter.genotype_concordance(missing_sites_hom_ref, true),
                counter.genotype_concordance(missing_sites_hom_ref, false),
            ],
        })
        .collect();

    let mut detail: Vec<DetailMetrics> = Vec::new();
    for (kind, counter) in [("SNP", &snp_counter), ("INDEL", &indel_counter)] {
        for truth_state in TRUTH_DECLARATION_ORDER {
            for call_state in CALL_DECLARATION_ORDER {
                let count = counter.count(truth_state, call_state);
                if count > 0 || output_all_rows {
                    detail.push(DetailMetrics {
                        variant_type: kind,
                        truth_sample: truth_sample.clone(),
                        call_sample: call_sample.clone(),
                        truth_state,
                        call_state,
                        count,
                        contingency_values: contingency_string(
                            contingency(scheme, call_state, truth_state).expect("a cell"),
                        ),
                    });
                }
            }
        }
    }

    let contingency_metrics: Vec<ContingencyMetrics> =
        [("SNP", &snp_counter), ("INDEL", &indel_counter)]
            .into_iter()
            .map(|(kind, counter)| {
                let counts = counter.contingency_counts(scheme);
                ContingencyMetrics {
                    variant_type: kind,
                    truth_sample: truth_sample.clone(),
                    call_sample: call_sample.clone(),
                    counts: [counts.tp, counts.tn, counts.fp, counts.fn_, counts.empty],
                }
            })
            .collect();

    let [summary_path, detail_path, contingency_path] = file_names(&output);
    std::fs::write(summary_path, metrics_file(summary))?;
    std::fs::write(detail_path, metrics_file(detail))?;
    std::fs::write(contingency_path, metrics_file(contingency_metrics))?;
    let _ = is_var;
    Ok(())
}
