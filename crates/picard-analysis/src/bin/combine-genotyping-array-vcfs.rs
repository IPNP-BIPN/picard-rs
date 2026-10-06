//! `CombineGenotypingArrayVcfs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.CombineGenotypingArrayVcfs.doWork` and `merge` at tag 3.4.0. The checks
//! a step makes are `picard_analysis::combine_genotyping_array_vcfs`; this is the lockstep walk,
//! the merged header and the merged record.
//!
//! The headers are merged by `smartMergeHeaders`, the first line of a key winning, and every line
//! naming one sample's run is dropped: the list is 3.4.0's, which adds the control intensities,
//! the dates, the scanner, the sample alias and `fileDate` to the older one. A record's counts are
//! recomputed over the merged genotypes; its other attributes are the first file's.

use htsjdk_vcf::chromosome_counts::calculate_chromosome_counts;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use htsjdk_vcf::variant::{Value, VariantContext};
use picard_analysis::combine_genotyping_array_vcfs::{
    check_attributes, check_step, sample_list, Refusal, Site,
};
use picard_analysis::fingerprint::reference_path;
use picard_analysis::metrics_cli::{thrown, Args};
use picard_analysis::vcf_io::{header_dictionary, read_path, unroll_paths, write_output, Record};

/// `sampleSpecificHeaders` at 3.4.0.
const SAMPLE_SPECIFIC: &[&str] = &[
    "pipelineVersion",
    "analysisVersionNumber",
    "autocallDate",
    "autocallGender",
    "chipWellBarcode",
    "expectedGender",
    "fingerprintGender",
    "gtcCallRate",
    "imagingDate",
    "p95Green",
    "p95Red",
    "sampleAlias",
    "scannerName",
    "Biotin(Bgnd)",
    "Biotin(High)",
    "DNP(Bgnd)",
    "DNP(High)",
    "Extension(A)",
    "Extension(C)",
    "Extension(G)",
    "Extension(T)",
    "Hyb(High)",
    "Hyb(Low)",
    "Hyb(Medium)",
    "NP(A)",
    "NP(C)",
    "NP(G)",
    "NP(T)",
    "NSB(Bgnd)Blue",
    "NSB(Bgnd)Green",
    "NSB(Bgnd)Purple",
    "NSB(Bgnd)Red",
    "Restore",
    "String(MM)",
    "String(PM)",
    "TargetRemoval",
    "fileDate",
];

fn text(value: &Value) -> String {
    match value {
        Value::Str(s) => s.clone(),
        Value::Int(i) => i.to_string(),
        Value::Double(d) => d.to_string(),
        Value::Bool(b) => b.to_string(),
        Value::List(l) => l.iter().map(text).collect::<Vec<_>>().join(","),
        Value::Missing => ".".to_string(),
    }
}

fn site(vc: &VariantContext) -> Site {
    Site {
        contig: vc.contig.clone(),
        start: vc.start,
        id: vc.id.clone(),
        reference: vc.alleles[0].display_string(),
        alternates: vc.alleles[1..].iter().map(|a| a.display_string()).collect(),
        attributes: vc
            .attributes
            .iter()
            .map(|(k, v)| (k.clone(), text(v)))
            .collect(),
    }
}

fn refuse(r: &Refusal) -> String {
    format!("{}: {}", r.class(), r.message())
}

/// `merge`.
fn merge(step: &[&VariantContext]) -> Result<VariantContext, String> {
    let first = step[0];
    let first_site = site(first);
    let mut out = first.clone();
    out.genotypes = Vec::new();
    let mut filters: Vec<String> = Vec::new();
    let mut applied = false;
    let mut log10 = None;
    for vc in step {
        let other = site(vc);
        check_step(&first_site, &other).map_err(|r| refuse(&r))?;
        out.genotypes.extend(vc.genotypes.iter().cloned());
        if log10.is_none() && vc.has_log10_p_error() {
            log10 = Some(vc.log10_p_error);
        }
        if let Some(f) = &vc.filters {
            applied = true;
            for name in f {
                if !filters.contains(name) {
                    filters.push(name.clone());
                }
            }
        }
        check_attributes(&first_site, &other).map_err(|r| refuse(&r))?;
    }
    out.log10_p_error = log10.unwrap_or(first.log10_p_error);
    filters.sort();
    out.filters = if applied { Some(filters) } else { None };
    let counts = calculate_chromosome_counts(&out, false, &[]);
    for (key, value) in counts.attributes() {
        match out.attributes.iter_mut().find(|(k, _)| *k == key) {
            Some(slot) => slot.1 = value,
            None => out.attributes.push((key, value)),
        }
    }
    out.attributes.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(out)
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT")]);
    let inputs = unroll_paths(&args.collection("INPUT", &[])).unwrap_or_else(|e| thrown(&e));
    let output = args.required("OUTPUT");
    let create_index = args.bool("CREATE_INDEX", false);
    let files: Vec<_> = inputs
        .iter()
        .map(|p| read_path(p).unwrap_or_else(|e| thrown(&e)))
        .collect();
    let dictionary = files
        .first()
        .and_then(|f| header_dictionary(&f.file.header));
    let named: Vec<(String, Vec<String>)> = inputs
        .iter()
        .zip(&files)
        .map(|(p, f)| (reference_path(p), f.file.header.samples.clone()))
        .collect();
    let samples = sample_list(&named).unwrap_or_else(|r| thrown(&refuse(&r)));
    if create_index && dictionary.is_none() {
        thrown(&refuse(&Refusal::NoSequenceDictionary));
    }
    let sources: Vec<Source> = files
        .iter()
        .map(|f| Source {
            header: &f.file.header,
            version: Some("VCFv4.2"),
        })
        .collect();
    let (mut lines, _) = smart_merge_headers(&sources, false)
        .unwrap_or_else(|e| thrown(&format!("{}: {}", e.class(), e.message())));
    lines.retain(|l| match l {
        HeaderLine::Unstructured { key, .. } | HeaderLine::Structured { key, .. } => {
            !SAMPLE_SPECIFIC.contains(&key.as_str())
        }
        _ => true,
    });
    // `new VCFHeader(lines, samples)` puts the lines in a `TreeSet`, where two contig lines compare
    // by their index alone: a second file's contig at an index the first file's already holds is
    // equal to it, and is dropped.
    let mut indices: Vec<i32> = Vec::new();
    lines.retain(|l| match l {
        HeaderLine::Contig { index, .. } => {
            if indices.contains(index) {
                false
            } else {
                indices.push(*index);
                true
            }
        }
        _ => true,
    });
    let header = VcfHeader { lines, samples };

    let mut records: Vec<Record> = Vec::new();
    let mut failure: Option<String> = None;
    let mut step = 0;
    loop {
        let current: Vec<Option<&VariantContext>> = files
            .iter()
            .map(|f| f.records.get(step).map(|r| &r.variant))
            .collect();
        let closed = current.iter().filter(|c| c.is_none()).count();
        if closed > 0 {
            if closed != files.len() {
                failure = Some(refuse(&Refusal::VariantCount));
            }
            break;
        }
        let present: Vec<&VariantContext> = current.into_iter().flatten().collect();
        match merge(&present) {
            Ok(vc) => records.push(Record {
                variant: vc,
                lazy_genotypes: None,
            }),
            Err(e) => {
                failure = Some(e);
                break;
            }
        }
        step += 1;
    }
    let index = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    if let Err(e) = write_output(&output, &header, &records, index) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
    if let Some(e) = failure {
        thrown(&e);
    }
}
