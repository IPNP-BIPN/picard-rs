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
use picard_analysis::args::arg;
use picard_analysis::find_mendelian_violations::{
    collect, Genotype, Interval, Metrics, Options, Sex, Site, Trio,
};

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
/// `PedFile.fromFile`: a line is split on tabs in `TAB_MODE` and on runs of whitespace otherwise,
/// any line without exactly six fields is skipped, and the trios are a `TreeMap` keyed by the
/// individual, so a later line replaces an earlier one and the order is the individuals' order.
fn read_pedigree(text: &str, tab_mode: bool) -> Vec<Trio> {
    let mut by_individual: std::collections::BTreeMap<String, Trio> =
        std::collections::BTreeMap::new();
    for line in text.lines() {
        // `String.split` keeps a leading empty field and drops the trailing ones.
        let mut fields: Vec<&str> = if tab_mode {
            line.split('\t').collect()
        } else {
            let mut fields = Vec::new();
            if line.starts_with(char::is_whitespace) {
                fields.push("");
            }
            fields.extend(line.split_whitespace());
            fields
        };
        while fields.len() > 1 && fields.last() == Some(&"") {
            fields.pop();
        }
        if fields.len() != 6 {
            continue;
        }
        by_individual.insert(
            fields[1].to_string(),
            Trio {
                family_id: fields[0].to_string(),
                father: fields[2].to_string(),
                mother: fields[3].to_string(),
                offspring: fields[1].to_string(),
                offspring_sex: match fields[4] {
                    "1" => Sex::Male,
                    "2" => Sex::Female,
                    _ => Sex::Unknown,
                },
            },
        );
    }
    by_individual.into_values().collect()
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

/// `FindMendelianViolations.writeAllViolations` for one family: the input's header lines and four
/// INFO lines over the trio's three samples, and each violating record cut to those samples with
/// its type and the input's AC, AF and AN beside it, in the dictionary's order with one record a
/// position (a `TreeSet` on the `VariantContextComparator`).
fn write_violations(
    dir: &str,
    vcf: &htsjdk_vcf::reader::VcfFile,
    trio: &picard_analysis::find_mendelian_violations::Trio,
    violations: &[(usize, &str)],
) -> Result<(), Box<dyn std::error::Error>> {
    use htsjdk_vcf::encoder::{MissingFields, VcfEncoder};
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
    let mut lines = vcf.header.lines.clone();
    lines.push(HeaderLine::info(
        "MV",
        Cardinality::Fixed(1),
        LineType::String,
        "Type of mendelian violation.",
    ));
    lines.push(HeaderLine::info(
        "AC_Orig",
        Cardinality::A,
        LineType::Integer,
        "Original AC",
    ));
    lines.push(HeaderLine::info(
        "AF_Orig",
        Cardinality::A,
        LineType::Float,
        "Original AF",
    ));
    lines.push(HeaderLine::info(
        "AN_Orig",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Original AN",
    ));
    let samples = vec![
        trio.mother.clone(),
        trio.father.clone(),
        trio.offspring.clone(),
    ];
    let header = VcfHeader {
        lines,
        samples: samples.clone(),
    };
    let contigs: Vec<String> = vcf
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Contig { fields, .. } => fields
                .iter()
                .find(|(k, _)| k == "ID")
                .map(|(_, v)| v.clone()),
            _ => None,
        })
        .collect();
    let order = |contig: &str| {
        contigs
            .iter()
            .position(|c| c == contig)
            .unwrap_or(usize::MAX)
    };
    let mut ordered: Vec<(usize, &str)> = violations.to_vec();
    ordered.sort_by_key(|(index, _)| {
        let record = &vcf.records[*index];
        (order(&record.contig), record.start)
    });
    ordered.dedup_by_key(|(index, _)| {
        let record = &vcf.records[*index];
        (record.contig.clone(), record.start)
    });
    let mut text = header.write();
    let encoder = VcfEncoder::new(&header).with_missing_fields(MissingFields::Allow);
    for (index, kind) in ordered {
        let mut record = vcf.records[index].clone();
        record.genotypes = record
            .genotypes
            .iter()
            .filter(|g| samples.contains(&g.sample_name))
            .cloned()
            .collect();
        let copy = |key: &str| {
            record
                .attributes
                .iter()
                .find(|(k, _)| k == key)
                .map(|(_, v)| v.clone())
        };
        let originals = [
            ("AC_Orig", copy("AC")),
            ("AF_Orig", copy("AF")),
            ("AN_Orig", copy("AN")),
        ];
        record.attributes.push((
            "MV".to_string(),
            htsjdk_vcf::variant::Value::Str(kind.to_string()),
        ));
        for (key, value) in originals {
            if let Some(value) = value {
                record.attributes.push((key.to_string(), value));
            }
        }
        text.push_str(&encoder.encode(&record).map_err(|e| format!("{e:?}"))?);
        text.push('\n');
    }
    let name: String = format!("{}.vcf", trio.family_id)
        .trim()
        .chars()
        .map(|c| {
            if c.is_whitespace() || "!\"#$%&'()*/:;<=>?@[]\\^`{|}~".contains(c) {
                '_'
            } else {
                c
            }
        })
        .collect();
    std::fs::write(std::path::Path::new(dir).join(name), text)?;
    Ok(())
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
    picard_analysis::argument_check::check("FindMendelianViolations");
    picard_analysis::output::install();
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
            // Picard's parser APPENDS to a list argument's default rather than replacing it, and
            // a `null` clears what the list holds so far.
            let mut all: Vec<String> = default.iter().map(|value| (*value).to_string()).collect();
            for value in given {
                if value == "null" {
                    all.clear();
                } else {
                    all.push(value);
                }
            }
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

    // `customCommandLineValidation`: a contig named as both sexes is refused before anything runs.
    {
        let mut both: Vec<&String> = options
            .female_chroms
            .iter()
            .filter(|name| options.male_chroms.contains(name))
            .collect();
        both.dedup();
        if !both.is_empty() {
            let names: Vec<&str> = both.iter().map(|name| name.as_str()).collect();
            eprintln!(
                "The following chromosomes were listed as both male and female sex chromosomes: [{}]",
                names.join(", ")
            );
            std::process::exit(1);
        }
    }
    let vcf_dir = arg(&args, "VCF_DIR=");
    // `IOUtil.assertDirectoryIsWritable(VCF_DIR)`, before anything is read.
    if let Some(dir) = &vcf_dir {
        if !std::path::Path::new(dir).is_dir() {
            eprintln!(
                "Exception in thread \"main\" htsjdk.samtools.SAMException: Directory does not exist: {}",
                picard_analysis::metrics_cli::absolute(dir)
            );
            std::process::exit(1);
        }
    }
    // `VariantProcessor` reads the input by interval, so it opens it requiring an index; a VCF
    // without one is refused before a record is read.
    if !std::path::Path::new(&format!("{input}.idx")).exists()
        && !std::path::Path::new(&format!("{input}.tbi")).exists()
    {
        eprintln!(
            "Exception in thread \"main\" htsjdk.tribble.TribbleException: An index is required, but none found with file ending .idx, for input source: {}",
            picard_analysis::fingerprint::file_uri(&input)
        );
        std::process::exit(1);
    }
    let vcf = picard_analysis::vcf_io::read_vcf_file(&input)?;
    let tab_mode =
        arg(&args, "TAB_MODE=").is_some_and(|value| picard_analysis::args::is_true(&value));
    let pedigree = read_pedigree(&std::fs::read_to_string(&trios)?, tab_mode);
    let samples = &vcf.header.samples;

    let mut by_family: std::collections::HashMap<String, Vec<(usize, &'static str)>> =
        std::collections::HashMap::new();
    let mut file = MetricsFile::new();
    for header in picard_analysis::command_line::default_headers("FindMendelianViolations") {
        file.add_header(&header);
    }
    for trio in &pedigree {
        // "Removing trio due to the following missing samples in VCF": a pedigree row naming a
        // parent the file does not have is dropped, which is what the placeholder `0` rows are.
        if ![&trio.father, &trio.mother, &trio.offspring]
            .iter()
            .all(|name| samples.contains(name))
        {
            continue;
        }
        let mut from_record: Vec<usize> = Vec::new();
        let sites: Vec<Site> = vcf
            .records
            .iter()
            .enumerate()
            .filter_map(|(index, record)| {
                let site = Some(Site {
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
                });
                from_record.push(index);
                site
            })
            .collect();
        let (metrics, violations) = collect(&sites, trio, &options);
        file.add_metric(&Row(metrics));
        by_family.entry(trio.family_id.clone()).or_default().extend(
            violations
                .iter()
                .map(|(site, violation)| (from_record[*site], violation.name())),
        );
    }
    std::fs::write(&output, file.write())?;
    // `writeAllViolations` walks the WHOLE pedigree, a `TreeMap` by individual, not the trios the
    // metrics kept: a founder's row, whose parents are both `0`, builds a header with a duplicate
    // sample and is refused there.
    if let Some(dir) = &vcf_dir {
        for trio in &pedigree {
            let names = [&trio.mother, &trio.father, &trio.offspring];
            if names[0] == names[1] || names[0] == names[2] || names[1] == names[2] {
                eprintln!(
                    "Exception in thread \"main\" htsjdk.tribble.TribbleException$InvalidHeader: Your input file has a malformed header: BUG: VCF header has duplicate sample names"
                );
                std::process::exit(1);
            }
            let records = by_family.get(&trio.family_id).cloned().unwrap_or_default();
            write_violations(dir, &vcf, trio, &records)?;
        }
    }
    Ok(())
}
