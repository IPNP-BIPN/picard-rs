//! `CombineGenotypingArrayVcfs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.CombineGenotypingArrayVcfs.doWork` and `merge` at tag 3.4.0. The read and
//! write round trip is `picard_analysis::vcf_io`, the header merge htsjdk-rs's
//! `smartMergeHeaders`, and the chromosome counts its `calculateChromosomeCounts`.
//!
//! The merge is a lockstep, not a merge by position: the nth record of every input is merged with
//! the nth of every other, and the walk stops as soon as ANY input runs out; if not all of them
//! did, the variant counts disagreed. Within one step every input is checked against the first in
//! turn -- locus, ID, REF, ALT count, each ALT, then its attributes -- so the first input to
//! disagree decides the message.
//!
//! The attribute check runs over each input's OWN attributes, in `HashMap` order, against the
//! first's, and skips the seven that may differ. A key only the first input has is never looked
//! for. `DP` is summed, and a positive sum is written back into the map `getAttributes()`
//! returned, which is unmodifiable: the run ends in an `UnsupportedOperationException` with no
//! message, so no depth ever reaches an output.
//!
//! The samples are each file's names in SORTED order (`getSampleNamesInOrder`), the files in the
//! order given, and a name seen in an earlier file is refused. The writer's options always include
//! `INDEX_ON_THE_FLY`, so a run with no sequence dictionary is refused by the writer builder even
//! when `CREATE_INDEX` is false -- by `doWork` itself when it is true.
//!
//! The headers are merged from a `HashSet<VCFHeader>`, whose iteration order is the headers'
//! identity hashes; where two inputs disagree on a line the merge keeps whichever comes first. The
//! port merges them in input order, and the corpus only gives inputs whose kept lines agree.

use htsjdk_vcf::chromosome_counts::{
    calculate_chromosome_counts, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY, ALLELE_NUMBER_KEY,
};
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext, NO_LOG10_PERROR};
use picard_analysis::java_hash_map::JavaHashMap;
use picard_analysis::java_number::parse_int;
use picard_analysis::metrics_cli::Args;
use picard_analysis::vcf_io::{
    die, header_dictionary, read_path, sample_names_in_order, unroll_paths, write_output, LazyVcf,
    Record, Sequence,
};

/// The header lines the merged file drops: each one sample's own.
const SAMPLE_SPECIFIC_HEADERS: [&str; 37] = [
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

/// The attributes allowed to differ between inputs: the three recalculated after the merge, and
/// four with "minor allowable changes".
const EXEMPT_ATTRIBUTES: [&str; 7] = ["AC", "AF", "AN", "devX_AB", "devY_AB", "SOURCE", "refSNP"];

const DEPTH_KEY: &str = "DP";

fn picard(message: &str) -> String {
    format!("picard.PicardException: {message}")
}

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

/// `vc.getAttributes().entrySet()`, in the order the codec's `HashMap` hands it out.
fn attributes_in_hash_order(vc: &VariantContext) -> Vec<(String, Value)> {
    let mut map = JavaHashMap::new();
    for (key, value) in &vc.attributes {
        map.put(key, value.clone());
    }
    map.iter()
        .map(|(k, v)| (k.to_string(), v.clone()))
        .collect()
}

fn attribute<'a>(vc: &'a VariantContext, key: &str) -> Option<&'a Value> {
    vc.attributes.iter().find(|(k, _)| k == key).map(|(_, v)| v)
}

/// `getAttributeAsInt(key, 0)` on what the codec stored, which is a string.
fn attribute_as_int(value: &Value) -> Result<i32, String> {
    match value {
        Value::Int(i) => Ok(*i as i32),
        Value::Str(s) => {
            parse_int(s).map_err(|message| format!("java.lang.NumberFormatException: {message}"))
        }
        _ => Err("java.lang.ClassCastException".to_string()),
    }
}

/// `merge`: one step's records, already known to be one per input.
fn merge(contexts: &[&VariantContext]) -> Result<VariantContext, String> {
    let first = contexts[0];
    let mut filters: Vec<String> = Vec::new();
    let mut depth: i32 = 0;
    let mut log10_p_error = NO_LOG10_PERROR;
    let mut any_filters_applied = false;
    let mut genotypes: Vec<Genotype> = Vec::new();

    for vc in contexts {
        if vc.start != first.start || vc.contig != first.contig {
            return Err(picard("Mismatch in loci among input VCFs"));
        }
        if vc.id != first.id {
            return Err(picard("Mismatch in ID field among input VCFs"));
        }
        if vc.reference() != first.reference() {
            return Err(picard("Mismatch in REF allele among input VCFs"));
        }
        // `checkThatAllelesMatch(vc, first)`, which checks the REF a second time.
        if vc.alternate_alleles().len() != first.alternate_alleles().len() {
            return Err(picard("Mismatch in ALT allele count among input VCFs"));
        }
        if vc.alternate_alleles() != first.alternate_alleles() {
            return Err(picard(&format!(
                "Mismatch in ALT allele among input VCFs for {}.{}",
                vc.contig, vc.start
            )));
        }

        genotypes.extend(vc.genotypes.iter().cloned());

        if log10_p_error == NO_LOG10_PERROR {
            log10_p_error = vc.log10_p_error;
        }
        if let Some(own) = &vc.filters {
            for filter in own {
                if !filters.contains(filter) {
                    filters.push(filter.clone());
                }
            }
            any_filters_applied = true;
        }

        if let Some(value) = attribute(vc, DEPTH_KEY) {
            depth += attribute_as_int(value)?;
        }

        for (key, value) in attributes_in_hash_order(vc) {
            if EXEMPT_ATTRIBUTES.contains(&key.as_str()) {
                continue;
            }
            match attribute(first, &key) {
                None => return Err(picard(&format!("Attribute '{key}' not found in all VCFs"))),
                Some(extant) if *extant != value => {
                    return Err(picard(&format!(
                        "Values for attribute '{key}' disagrees among input VCFs"
                    )))
                }
                Some(_) => {}
            }
        }
    }

    // `firstAttributes.put(DP, ...)` on `Collections.unmodifiableMap`.
    if depth > 0 {
        return Err("java.lang.UnsupportedOperationException".to_string());
    }

    let mut merged = VariantContext::new(&first.contig, first.start, first.alleles.clone());
    merged.stop = first.stop;
    merged.id = first.id.clone();
    merged.genotypes = genotypes;
    merged.log10_p_error = log10_p_error;
    if any_filters_applied {
        filters.sort();
        merged.filters = Some(filters);
    }

    // `calculateChromosomeCounts(builder, false)`: AN always, AC and AF while there is an
    // alternate allele, and removed when there is none.
    let counts = calculate_chromosome_counts(&merged, false, &[]);
    let mut attributes: Vec<(String, Value)> = first.attributes.clone();
    if !merged.genotypes.is_empty() {
        attributes.retain(|(key, _)| {
            ![ALLELE_NUMBER_KEY, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY].contains(&key.as_str())
        });
        attributes.extend(counts.attributes());
    }
    merged.attributes = attributes;
    Ok(merged)
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT")]);
    let inputs = args.all("INPUT");
    let output = args.required("OUTPUT");
    let create_index = args.bool("CREATE_INDEX", true);

    let paths = unroll_paths(&inputs).unwrap_or_else(|exception| die(&exception));
    for path in &paths {
        if !std::path::Path::new(path).is_file() {
            die(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(path)
            ));
        }
    }

    let mut files: Vec<LazyVcf> = Vec::with_capacity(paths.len());
    for path in &paths {
        files.push(read_path(path).unwrap_or_else(|exception| die(&exception)));
    }
    // `VCFFileReader.getSequenceDictionary(UNROLLED_INPUT.get(0))`.
    let dictionary: Option<Vec<Sequence>> = files
        .first()
        .and_then(|first| header_dictionary(&first.file.header));

    let mut sample_list: Vec<String> = Vec::new();
    for (path, vcf) in paths.iter().zip(&files) {
        for sample in sample_names_in_order(&vcf.file.header) {
            if sample_list.contains(&sample) {
                die(&format!(
                    "java.lang.IllegalArgumentException: Input file {} contains a sample entry \
                     ({sample}) that appears in another input file.",
                    absolute(path)
                ));
            }
            sample_list.push(sample);
        }
    }

    if create_index && dictionary.is_none() {
        die(&picard(
            "A sequence dictionary must be available (either through the input file or by \
             setting it explicitly) when creating indexed output.",
        ));
    }
    // `VariantContextWriterBuilder.build`: INDEX_ON_THE_FLY is one of its default options.
    if dictionary.is_none() {
        die(
            "java.lang.IllegalArgumentException: A reference dictionary is required for creating \
             Tribble indices on the fly",
        );
    }

    let versions: Vec<Option<&str>> = files
        .iter()
        .map(|f| f.file.header_version.map(|v| v.version_string()))
        .collect();
    let sources: Vec<Source> = files
        .iter()
        .zip(&versions)
        .map(|(f, version)| Source {
            header: &f.file.header,
            version: *version,
        })
        .collect();
    let (mut lines, _) = smart_merge_headers(&sources, false)
        .unwrap_or_else(|e| die(&format!("{}: {}", e.class(), e.message())));
    lines.retain(|line: &HeaderLine| !SAMPLE_SPECIFIC_HEADERS.contains(&line.key()));
    let header = VcfHeader {
        lines,
        samples: sample_list,
    };

    // The lockstep. A refusal still leaves what was written before it, but the run has failed.
    let mut records: Vec<Record> = Vec::new();
    let mut failure: Option<String> = None;
    let mut step = 0usize;
    loop {
        let available: Vec<&VariantContext> = files
            .iter()
            .filter_map(|f| f.records.get(step).map(|r| &r.variant))
            .collect();
        if available.len() != files.len() {
            if !available.is_empty() {
                failure = Some(picard("Mismatch in number of variants among input VCFs"));
            }
            break;
        }
        match merge(&available) {
            Ok(merged) => records.push(Record {
                variant: merged,
                lazy_genotypes: None,
            }),
            Err(exception) => {
                failure = Some(exception);
                break;
            }
        }
        step += 1;
    }

    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    if let Err(e) = write_output(&output, &header, &records, index_dictionary) {
        die(&format!("htsjdk.tribble.TribbleException: {e}"));
    }
    if let Some(exception) = failure {
        die(&exception);
    }
}
