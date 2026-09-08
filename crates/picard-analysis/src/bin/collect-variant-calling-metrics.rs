//! `CollectVariantCallingMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.CollectVariantCallingMetrics.doWork` and
//! `picard.vcf.CallingMetricAccumulator` at tag 3.4.0. Which variant reaches which counter, and
//! the columns the counters produce, live in `picard_analysis::collect_variant_calling_metrics`.
//!
//! The tool writes two files from one prefix: a summary row for the file, and a detail row per
//! sample. They are not each other's sums. The summary counts a site ONCE however many samples
//! carry it, and a detail row counts it only for a sample whose call is not homozygous reference.
//!
//! Three things decide what is counted at all. `TARGET_INTERVALS` restricts the sites, and the
//! dbSNP membership is built over those same intervals. A site where every call is homozygous
//! reference is excluded before any counter is touched. And a FILTERED site increments its
//! filtered column and nothing else.
//!
//! The dbSNP answer is two answers: a bitset of SNP sites and another of indel sites, each asked
//! for the site's own type.

use std::collections::HashSet;

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use picard_analysis::collect_variant_calling_metrics::{
    accumulate_detail, accumulate_site, file_names, Call, Counts, HiddenCounts, Site,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// One interval of an interval list, uniqued the way `IntervalList.uniqued` does.
struct Interval {
    contig: String,
    start: i64,
    end: i64,
}

fn read_intervals(text: &str) -> Vec<Interval> {
    let mut contigs: Vec<String> = Vec::new();
    let mut intervals: Vec<Interval> = Vec::new();
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
            Some(last) if last.contig == interval.contig && interval.start <= last.end + 1 => {
                last.end = last.end.max(interval.end);
            }
            _ => merged.push(interval),
        }
    }
    merged
}

fn within(record: &VariantContext, intervals: &[Interval]) -> bool {
    intervals.iter().any(|interval| {
        record.contig == interval.contig
            && record.start <= interval.end
            && record.stop >= interval.start
    })
}

/// `VariantContext.isSNP()` over a decoded record.
fn record_is_snp(record: &VariantContext) -> bool {
    let reference = record
        .alleles
        .iter()
        .find(|allele| allele.is_reference())
        .map(|allele| allele.base_string())
        .unwrap_or_default();
    reference.len() == 1
        && record
            .alleles
            .iter()
            .filter(|allele| !allele.is_reference())
            .all(|allele| allele.base_string().len() == 1)
}

fn site_of(
    record: &VariantContext,
    samples: &[String],
    snps: &HashSet<(String, i64)>,
    indels: &HashSet<(String, i64)>,
) -> Site {
    let reference = record
        .alleles
        .iter()
        .find(|allele| allele.is_reference())
        .map(|allele| allele.base_string())
        .unwrap_or_default();
    let alternates: Vec<String> = record
        .alleles
        .iter()
        .filter(|allele| !allele.is_reference() && !allele.is_no_call())
        .map(|allele| allele.base_string())
        .collect();
    let index_of = |allele: &htsjdk_vcf::allele::Allele| -> usize {
        record
            .alleles
            .iter()
            .position(|other| other == allele)
            .unwrap_or(0)
    };
    let calls = samples
        .iter()
        .map(|sample| {
            let genotype = record
                .genotypes
                .iter()
                .find(|genotype| genotype.sample_name == *sample);
            match genotype {
                None => Call {
                    alleles: None,
                    gq: -1,
                    allele_depths: None,
                },
                Some(genotype) => Call {
                    alleles: if genotype.alleles.iter().any(|allele| allele.is_no_call()) {
                        None
                    } else {
                        Some(genotype.alleles.iter().map(index_of).collect())
                    },
                    gq: genotype.gq.unwrap_or(-1),
                    allele_depths: genotype.ad.clone(),
                },
            }
        })
        .collect();
    let key = (record.contig.clone(), record.start);
    Site {
        reference,
        alternates,
        filtered: record.filters.as_ref().is_some_and(|f| !f.is_empty()),
        in_db_snp_snps: snps.contains(&key),
        in_db_snp_indels: indels.contains(&key),
        calls,
    }
}

/// The twenty columns every row carries, and the five a detail row carries first.
const SUMMARY_COLUMNS: [&str; 20] = [
    "TOTAL_SNPS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "FILTERED_SNPS",
    "PCT_DBSNP",
    "DBSNP_TITV",
    "NOVEL_TITV",
    "TOTAL_INDELS",
    "NOVEL_INDELS",
    "FILTERED_INDELS",
    "PCT_DBSNP_INDELS",
    "NUM_IN_DB_SNP_INDELS",
    "DBSNP_INS_DEL_RATIO",
    "NOVEL_INS_DEL_RATIO",
    "TOTAL_MULTIALLELIC_SNPS",
    "NUM_IN_DB_SNP_MULTIALLELIC",
    "TOTAL_COMPLEX_INDELS",
    "NUM_IN_DB_SNP_COMPLEX_INDELS",
    "SNP_REFERENCE_BIAS",
    "NUM_SINGLETONS",
];

const DETAIL_COLUMNS: [&str; 25] = [
    "SAMPLE_ALIAS",
    "HET_HOMVAR_RATIO",
    "PCT_GQ0_VARIANTS",
    "TOTAL_GQ0_VARIANTS",
    "TOTAL_HET_DEPTH",
    "TOTAL_SNPS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "FILTERED_SNPS",
    "PCT_DBSNP",
    "DBSNP_TITV",
    "NOVEL_TITV",
    "TOTAL_INDELS",
    "NOVEL_INDELS",
    "FILTERED_INDELS",
    "PCT_DBSNP_INDELS",
    "NUM_IN_DB_SNP_INDELS",
    "DBSNP_INS_DEL_RATIO",
    "NOVEL_INS_DEL_RATIO",
    "TOTAL_MULTIALLELIC_SNPS",
    "NUM_IN_DB_SNP_MULTIALLELIC",
    "TOTAL_COMPLEX_INDELS",
    "NUM_IN_DB_SNP_COMPLEX_INDELS",
    "SNP_REFERENCE_BIAS",
    "NUM_SINGLETONS",
];

/// `calculateDerivedFields`, which is where the ratios come from -- and which does NOT guard its
/// divisions: a file with no known SNPs writes `?` for `PCT_DBSNP`, and the TI/TV ratios stay at
/// zero rather than becoming NaN because they are only assigned when the denominator is positive.
fn summary_values(counts: &Counts, hidden: &HiddenCounts) -> Vec<Value> {
    let novel_snps = counts.total_snps - counts.num_in_db_snp;
    let novel_indels = counts.total_indels - counts.num_in_db_snp_indels;
    let db_snp_titv = if hidden_transversions(counts, true) > 0 {
        counts.db_snp_transitions as f64 / counts.db_snp_transversions as f64
    } else {
        0.0
    };
    let novel_titv = if hidden_transversions(counts, false) > 0 {
        counts.novel_transitions as f64 / counts.novel_transversions as f64
    } else {
        0.0
    };
    let db_snp_ins_del = if hidden.db_snp_deletions > 0 {
        hidden.db_snp_insertions as f64 / hidden.db_snp_deletions as f64
    } else {
        0.0
    };
    let novel_ins_del = if hidden.novel_deletions > 0 {
        hidden.novel_insertions as f64 / hidden.novel_deletions as f64
    } else {
        0.0
    };
    vec![
        Value::Long(counts.total_snps),
        Value::Long(counts.num_in_db_snp),
        Value::Long(novel_snps),
        Value::Long(counts.filtered_snps),
        // PCT_DBSNP is a FLOAT in the reference, which is why it prints six digits of a float's
        // precision rather than a double's.
        Value::Double(f64::from(
            counts.num_in_db_snp as f32 / counts.total_snps as f32,
        )),
        Value::Double(db_snp_titv),
        Value::Double(novel_titv),
        Value::Long(counts.total_indels),
        Value::Long(novel_indels),
        Value::Long(counts.filtered_indels),
        Value::Double(f64::from(
            counts.num_in_db_snp_indels as f32 / counts.total_indels as f32,
        )),
        Value::Long(counts.num_in_db_snp_indels),
        Value::Double(db_snp_ins_del),
        Value::Double(novel_ins_del),
        Value::Long(counts.total_multiallelic_snps),
        Value::Long(counts.num_in_db_snp_multiallelic),
        Value::Long(hidden.total_complex_indels),
        Value::Long(hidden.num_in_db_snp_complex_indels),
        Value::Double(
            hidden.reference_allele_observations as f64
                / (hidden.reference_allele_observations + hidden.alternate_allele_observations)
                    as f64,
        ),
        Value::Long(counts.num_singletons),
    ]
}

fn hidden_transversions(counts: &Counts, known: bool) -> i64 {
    if known {
        counts.db_snp_transversions
    } else {
        counts.novel_transversions
    }
}

struct SummaryRow(Counts, HiddenCounts);

impl MetricBean for SummaryRow {
    fn class_name(&self) -> &str {
        "picard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &SUMMARY_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        summary_values(&self.0, &self.1)
    }
}

struct DetailRow {
    sample: String,
    counts: Counts,
    hidden: HiddenCounts,
}

impl MetricBean for DetailRow {
    fn class_name(&self) -> &str {
        "picard.vcf.CollectVariantCallingMetrics$VariantCallingDetailMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &DETAIL_COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let hets = self.hidden.number_of_hets;
        let hom_var = self.hidden.number_of_hom_var;
        let mut out = vec![
            Value::Str(self.sample.clone()),
            Value::Double(hets as f64 / hom_var as f64),
            Value::Double(self.hidden.total_gq0_variants as f64 / (hets + hom_var) as f64),
            Value::Long(self.hidden.total_gq0_variants),
            Value::Long(
                self.hidden.reference_allele_observations
                    + self.hidden.alternate_allele_observations,
            ),
        ];
        out.extend(summary_values(&self.counts, &self.hidden));
        out
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let dbsnp = arg(&args, "DBSNP=").ok_or("DBSNP= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if arg(&args, "GVCF_INPUT=").as_deref() == Some("true") {
        // `GvcfMetricAccumulator.setup` refuses a file with more than one sample, and the message
        // is the tool's answer for those rows.
        let vcf = read_vcf(&std::fs::read_to_string(&input)?).map_err(|e| format!("{e:?}"))?;
        if vcf.header.samples.len() != 1 {
            eprintln!(
                "Exception in thread \"main\" java.lang.IllegalArgumentException: Expected to have exactly 1 sample in a GVCF, found {}",
                vcf.header.samples.len()
            );
            std::process::exit(1);
        }
    }

    let vcf = read_vcf(&std::fs::read_to_string(&input)?).map_err(|e| format!("{e:?}"))?;
    let known = read_vcf(&std::fs::read_to_string(&dbsnp)?).map_err(|e| format!("{e:?}"))?;
    let intervals = match arg(&args, "TARGET_INTERVALS=") {
        None => Vec::new(),
        Some(path) => read_intervals(&std::fs::read_to_string(path)?),
    };

    // `createSnpAndIndelBitSets`: one set of SNP sites and one of indel sites, over the intervals.
    let mut snp_sites: HashSet<(String, i64)> = HashSet::new();
    let mut indel_sites: HashSet<(String, i64)> = HashSet::new();
    for record in &known.records {
        if !intervals.is_empty() && !within(record, &intervals) {
            continue;
        }
        let key = (record.contig.clone(), record.start);
        if record_is_snp(record) {
            snp_sites.insert(key);
        } else {
            indel_sites.insert(key);
        }
    }

    let samples = vcf.header.samples.clone();
    let mut summary = Counts::default();
    let mut summary_hidden = HiddenCounts::default();
    let mut details: Vec<(Counts, HiddenCounts)> =
        samples.iter().map(|_| Default::default()).collect();

    for record in &vcf.records {
        if !intervals.is_empty() && !within(record, &intervals) {
            continue;
        }
        let site = site_of(record, &samples, &snp_sites, &indel_sites);
        if site.is_excluded() {
            continue;
        }
        let singleton = site.singleton_sample();
        accumulate_site(
            &mut summary,
            &mut summary_hidden,
            None,
            &site,
            None,
            singleton.is_some(),
        );
        for (index, call) in site.calls.iter().enumerate() {
            if call.is_hom_ref() {
                continue;
            }
            let (counts, hidden) = &mut details[index];
            accumulate_detail(
                counts,
                hidden,
                &mut summary_hidden,
                &site,
                call,
                singleton == Some(index),
            );
        }
    }

    let (detail_path, summary_path) = file_names(&output);
    let mut summary_file = MetricsFile::new();
    summary_file.add_header("CollectVariantCallingMetrics <command line>");
    summary_file.add_header("Started on: <timestamp>");
    summary_file.add_metric(&SummaryRow(summary, summary_hidden));
    std::fs::write(summary_path, summary_file.write())?;

    let mut detail_file = MetricsFile::new();
    detail_file.add_header("CollectVariantCallingMetrics <command line>");
    detail_file.add_header("Started on: <timestamp>");
    for (index, sample) in samples.iter().enumerate() {
        let (counts, hidden) = details[index].clone();
        detail_file.add_metric(&DetailRow {
            sample: sample.clone(),
            counts,
            hidden,
        });
    }
    std::fs::write(detail_path, detail_file.write())?;
    Ok(())
}
