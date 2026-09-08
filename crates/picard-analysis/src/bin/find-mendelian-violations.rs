//! `FindMendelianViolations` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.MendelianViolations.FindMendelianViolations.doWork` at tag 3.4.0. The
//! violation classes, the filters and the sex-chromosome rules live in
//! `picard_analysis::find_mendelian_violations`.
//!
//! A trio comes out of the pedigree file, and every row of one names a trio: a row whose parents
//! are the placeholder `0` builds a trio for two samples the VCF does not have, and that trio is
//! dropped rather than refused. Only the rows whose three samples are all in the file survive.
//!
//! The genotype quality the filters read is the one the likelihoods give, not the `GQ` column on
//! its own, so a call whose `PL` disagrees with its `GQ` is judged on the `PL`.

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use htsjdk_vcf::reader::read_vcf;
use picard_analysis::find_mendelian_violations::{
    collect, Genotype, Interval, Metrics, Options, Sex, Site, Trio,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn args_of(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

const COLUMNS: [&str; 14] = [
    "FAMILY_ID",
    "MOTHER",
    "FATHER",
    "OFFSPRING",
    "OFFSPRING_SEX",
    "NUM_VARIANT_SITES",
    "NUM_DIPLOID_DENOVO",
    "NUM_HOMVAR_HOMVAR_HET",
    "NUM_HOMREF_HOMVAR_HOM",
    "NUM_HOM_HET_HOM",
    "NUM_HAPLOID_DENOVO",
    "NUM_HAPLOID_OTHER",
    "NUM_OTHER",
    "TOTAL_MENDELIAN_VIOLATIONS",
];

struct Row(Metrics);

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.vcf.MendelianViolations.MendelianViolationMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let m = &self.0;
        vec![
            Value::Str(m.family_id.clone()),
            Value::Str(m.mother.clone()),
            Value::Str(m.father.clone()),
            Value::Str(m.offspring.clone()),
            Value::Str(m.offspring_sex.clone()),
            Value::Long(m.num_variant_sites as i64),
            Value::Long(m.num_diploid_denovo as i64),
            Value::Long(m.num_homvar_homvar_het as i64),
            Value::Long(m.num_homref_homvar_hom as i64),
            Value::Long(m.num_hom_het_hom as i64),
            Value::Long(m.num_haploid_denovo as i64),
            Value::Long(m.num_haploid_other as i64),
            Value::Long(m.num_other as i64),
            Value::Long(m.total_mendelian_violations as i64),
        ]
    }
}

/// `PedFile.fromFile`: six whitespace-separated columns, the fifth being the sex.
fn read_pedigree(text: &str) -> Vec<Trio> {
    text.lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .filter_map(|line| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 5 {
                return None;
            }
            Some(Trio {
                family_id: fields[0].to_string(),
                father: fields[2].to_string(),
                mother: fields[3].to_string(),
                offspring: fields[1].to_string(),
                offspring_sex: match fields[4] {
                    "1" => Sex::Male,
                    "2" => Sex::Female,
                    _ => Sex::Unknown,
                },
            })
        })
        .collect()
}

/// `X:60001-2699520` and the three others the default names.
fn parse_interval(text: &str) -> Option<Interval> {
    let (contig, range) = text.split_once(':')?;
    let (start, end) = range.split_once('-')?;
    Some(Interval {
        contig: contig.to_string(),
        start: start.parse().ok()?,
        end: end.parse().ok()?,
    })
}

fn genotype_of(record: &htsjdk_vcf::variant::VariantContext, sample: &str) -> Option<Genotype> {
    let genotype = record
        .genotypes
        .iter()
        .find(|genotype| genotype.sample_name == sample)?;
    let index_of = |allele: &htsjdk_vcf::allele::Allele| -> usize {
        record
            .alleles
            .iter()
            .position(|other| other == allele)
            .unwrap_or(0)
    };
    let alleles = if genotype.alleles.iter().any(|allele| allele.is_no_call()) {
        None
    } else {
        Some(genotype.alleles.iter().map(index_of).collect())
    };
    Some(Genotype {
        alleles,
        // `getGQ()` and `getDP()` answer -1 for a field that is not there, which is the value the
        // floors compare against.
        gq: genotype.gq.unwrap_or(-1),
        dp: genotype.dp.unwrap_or(-1),
        ad: genotype.ad.clone(),
        pl: genotype.pl.clone(),
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let trios = arg(&args, "TRIOS=").ok_or("TRIOS= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let number = |key: &str, default: i32| -> i32 {
        arg(&args, key)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let list = |key: &str, default: &[&str]| -> Vec<String> {
        let given = args_of(&args, key);
        if given.is_empty() {
            default.iter().map(|value| (*value).to_string()).collect()
        } else {
            // Picard's parser APPENDS to a list argument's default rather than replacing it.
            let mut all: Vec<String> = default.iter().map(|value| (*value).to_string()).collect();
            all.extend(given);
            all
        }
    };
    let options = Options {
        min_gq: number("MIN_GQ=", 30),
        min_dp: number("MIN_DP=", 0),
        min_het_fraction: arg(&args, "MIN_HET_FRACTION=")
            .and_then(|value| value.parse().ok())
            .unwrap_or(0.3),
        skip_chroms: list("SKIP_CHROMS=", &["MT", "chrM"]),
        male_chroms: list("MALE_CHROMS=", &["Y", "chrY"]),
        female_chroms: list("FEMALE_CHROMS=", &["X", "chrX"]),
        par_intervals: list(
            "PSEUDO_AUTOSOMAL_REGIONS=",
            &[
                "X:60001-2699520",
                "X:154931044-155260560",
                "chrX:10001-2781479",
                "chrX:155701383-156030895",
            ],
        )
        .iter()
        .filter_map(|text| parse_interval(text))
        .collect(),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let vcf = read_vcf(&std::fs::read_to_string(&input)?).map_err(|e| format!("{e:?}"))?;
    let pedigree = read_pedigree(&std::fs::read_to_string(&trios)?);
    let samples = &vcf.header.samples;

    let mut file = MetricsFile::new();
    file.add_header("FindMendelianViolations <command line>");
    file.add_header("Started on: <timestamp>");
    for trio in &pedigree {
        // "Removing trio due to the following missing samples in VCF": a pedigree row naming a
        // parent the file does not have is dropped, which is what the placeholder `0` rows are.
        if ![&trio.father, &trio.mother, &trio.offspring]
            .iter()
            .all(|name| samples.contains(name))
        {
            continue;
        }
        let sites: Vec<Site> = vcf
            .records
            .iter()
            .filter_map(|record| {
                Some(Site {
                    contig: record.contig.clone(),
                    position: record.start as i32,
                    filtered: record.filters.as_ref().is_some_and(|f| !f.is_empty()),
                    alleles: record
                        .alleles
                        .iter()
                        .map(|allele| allele.base_string())
                        .collect(),
                    father: genotype_of(record, &trio.father)?,
                    mother: genotype_of(record, &trio.mother)?,
                    child: genotype_of(record, &trio.offspring)?,
                })
            })
            .collect();
        let (metrics, _violations) = collect(&sites, trio, &options);
        file.add_metric(&Row(metrics));
    }
    std::fs::write(&output, file.write())?;
    Ok(())
}
