//! `GatherVcfs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.GatherVcfs.doWork` at tag 3.4.0 for a plain-text output, which is the
//! conventional gather: the block-copying one is taken only when every input AND the output are
//! block compressed. The read and write round trip is `picard_analysis::vcf_io`.
//!
//! Nothing is merged. The output header is the first file's own, `COMMENT` added through
//! `addMetaDataLine` (so only the first comment given is written), and every file's records follow
//! in file order -- each file's genotype blocks copied or re-encoded as its own sample order
//! decides. The files must carry the same samples in the same ORDER, unlike MergeVcfs, which
//! compares them sorted.
//!
//! # Most refusals are a log line and exit status 1, not an exception
//!
//! The checks and the gather run inside `try { ... } catch (RuntimeException e)`, which logs
//! `There was a problem with gathering the INPUT.` with the exception run on after it, deletes the
//! output, and returns 1. Two things escape it. The index check comes before it, and is an uncaught
//! `PicardException`. And a dictionary mismatch is an `AssertionError`, which is not a
//! `RuntimeException`: it logs the two files it compared and then escapes as the JVM's uncaught
//! exception. A first file without contig lines under `CREATE_INDEX=false` gets as far as
//! `new VariantContextComparator(null)`, whose `NullPointerException` is caught like the rest.

use htsjdk_vcf::comparator::VariantContextComparator;
use htsjdk_vcf::variant::VariantContext;
use picard_analysis::vcf_io::{
    add_other_meta_data_line, assert_same_dictionary, die, header_dictionary, java_compare,
    log_error, read_path, unroll_paths, write_output, LazyVcf, Record,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// A collection argument: every occurrence appends, and the value `null` empties it (Barclay).
fn collection(args: &[String], keys: &[&str]) -> Vec<String> {
    let mut values = Vec::new();
    for a in args {
        if let Some(value) = keys.iter().find_map(|key| a.strip_prefix(key)) {
            if value == "null" {
                values.clear();
            } else {
                values.push(value.to_string());
            }
        }
    }
    values
}

const NPE_UNKNOWN_CONTIG: &str = "java.lang.NullPointerException: Cannot invoke \
     \"java.lang.Integer.intValue()\" because the return value of \"java.util.Map.get(Object)\" \
     is null";

/// How the `try` block ended, when it did not end well.
enum Failure {
    /// A `RuntimeException`, as its `toString()`: caught, logged, exit status 1.
    Caught(String),
    /// An `AssertionError` from `assertSameDictionary`, after the two lines it logs.
    Escaped(String),
}

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

/// `TreeSet.toString()` of the names in `a` that are not in `b`.
fn unique_to(a: &[String], b: &[String]) -> String {
    let mut names: Vec<&String> = a.iter().filter(|name| !b.contains(name)).collect();
    names.sort_by(|x, y| java_compare(x, y));
    names.dedup();
    let names: Vec<&str> = names.iter().map(|name| name.as_str()).collect();
    format!("[{}]", names.join(", "))
}

fn compare(
    comparator: &VariantContextComparator,
    a: &VariantContext,
    b: &VariantContext,
) -> Result<i32, Failure> {
    comparator
        .compare(a, b)
        .map_err(|_| Failure::Caught(NPE_UNKNOWN_CONTIG.to_string()))
}

/// `assertSameSamplesAndValidOrdering`: the files, reordered when asked.
fn assert_same_samples_and_valid_ordering(
    mut files: Vec<(String, LazyVcf)>,
    reorder: bool,
) -> Result<Vec<(String, LazyVcf)>, Failure> {
    let header = files[0].1.file.header.clone();
    let Some(dict) = header_dictionary(&header) else {
        return Err(Failure::Caught(
            "java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \"dictionary\" is \
             null"
                .to_string(),
        ));
    };
    let names: Vec<String> = dict.iter().map(|s| s.name.clone()).collect();
    let comparator = VariantContextComparator::from_contigs(&names)
        .map_err(|e| Failure::Caught(format!("{}: {}", e.class(), e.message())))?;
    let samples = header.samples.clone();

    if reorder {
        // `List.sort` is a stable merge sort; a file with no records sorts after every other.
        let mut failure: Option<Failure> = None;
        files.sort_by(|a, b| {
            let (first_a, first_b) = (a.1.records.first(), b.1.records.first());
            match (first_a, first_b) {
                (None, None) => std::cmp::Ordering::Equal,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (Some(_), None) => std::cmp::Ordering::Less,
                (Some(x), Some(y)) => match compare(&comparator, &x.variant, &y.variant) {
                    Ok(c) => c.cmp(&0),
                    Err(e) => {
                        failure.get_or_insert(e);
                        std::cmp::Ordering::Equal
                    }
                },
            }
        });
        if let Some(failure) = failure {
            return Err(failure);
        }
    }

    let mut last: Option<(&String, &VariantContext)> = None;
    for (path, vcf) in &files {
        match header_dictionary(&vcf.file.header) {
            None => {
                return Err(Failure::Caught(
                    "java.lang.NullPointerException: Cannot read field \"mSequences\" because \
                     \"that\" is null"
                        .to_string(),
                ))
            }
            Some(that) => {
                if let Err(message) = assert_same_dictionary(&dict, &that) {
                    log_error("GatherVcfs", &format!("File #1: {}", files[0].0));
                    log_error("GatherVcfs", &format!("File #2: {path}"));
                    return Err(Failure::Escaped(format!(
                        "java.lang.AssertionError: {message}"
                    )));
                }
            }
        }
        let these = &vcf.file.header.samples;
        if &samples != these {
            return Err(Failure::Caught(format!(
                "java.lang.IllegalArgumentException: VCFs do not have identical sample lists. \
                 Samples unique to first file: {}. Samples unique to {}: {}.",
                unique_to(&samples, these),
                absolute(path),
                unique_to(these, &samples)
            )));
        }
        if let Some(current) = vcf.records.first() {
            if let Some((last_path, last_context)) = last {
                if compare(&comparator, last_context, &current.variant)? >= 0 {
                    return Err(Failure::Caught(format!(
                        "java.lang.IllegalArgumentException: First record in file {} is not after \
                         first record in previous file {}",
                        absolute(path),
                        absolute(last_path)
                    )));
                }
            }
            last = Some((path, &current.variant));
        }
    }
    Ok(files)
}

/// `gatherConventionally`, up to the writing: the header and the records.
fn gather_conventionally(
    files: Vec<(String, LazyVcf)>,
    comments: &[String],
) -> Result<(htsjdk_vcf::header::VcfHeader, Vec<Record>), Failure> {
    let mut header = files[0].1.file.header.clone();
    for comment in comments {
        add_other_meta_data_line(&mut header, "GatherVcfs.comment", comment);
    }
    let contigs: Vec<_> = header
        .lines
        .iter()
        .filter(|line| matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }))
        .cloned()
        .collect();
    let comparator = VariantContextComparator::from_header_lines(&contigs)
        .map_err(|e| Failure::Caught(format!("{}: {}", e.class(), e.message())))?;

    let mut records: Vec<Record> = Vec::new();
    let mut last_file: Option<String> = None;
    for (path, vcf) in files {
        if let (Some(last_context), Some(first)) = (records.last(), vcf.records.first()) {
            if compare(&comparator, &first.variant, &last_context.variant)? <= 0 {
                let (vc, lc) = (&first.variant, &last_context.variant);
                return Err(Failure::Caught(format!(
                    "java.lang.IllegalArgumentException: First variant in file {} is at {}:{} but \
                     last variant in earlier file {} is at {}:{}",
                    absolute(&path),
                    vc.contig,
                    vc.start,
                    absolute(last_file.as_deref().unwrap_or_default()),
                    lc.contig,
                    lc.start
                )));
            }
        }
        records.extend(vcf.records);
        // Set for every file, an empty one included.
        last_file = Some(path);
    }
    Ok((header, records))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let inputs = collection(&args, &["INPUT=", "I="]);
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let comments = collection(&args, &["COMMENT=", "CO="]);
    let reorder = arg(&args, "REORDER_INPUT_BY_FIRST_VARIANT=")
        .or_else(|| arg(&args, "RI="))
        .map(|value| value == "true")
        .unwrap_or(false);
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    if !output.ends_with(".vcf") {
        return Err("only a .vcf OUTPUT is ported".into());
    }

    let paths = unroll_paths(&inputs).unwrap_or_else(|exception| die(&exception));

    // `VCFFileReader.getSequenceDictionary(unrolledPaths.get(0))`, outside the `try`.
    let first = read_path(&paths[0]).unwrap_or_else(|exception| die(&exception));
    let sequence_dictionary = header_dictionary(&first.file.header);
    if create_index && sequence_dictionary.is_none() {
        die(
            "picard.PicardException: In order to index the resulting VCF input VCFs must contain \
             ##contig lines.",
        );
    }

    let gathered = (|| {
        let mut files = vec![(paths[0].clone(), first)];
        for path in &paths[1..] {
            let vcf = read_path(path).map_err(Failure::Caught)?;
            files.push((path.clone(), vcf));
        }
        let files = assert_same_samples_and_valid_ordering(files, reorder)?;
        gather_conventionally(files, &comments)
    })();

    match gathered {
        Ok((header, records)) => {
            let index_dictionary = if create_index {
                sequence_dictionary.as_deref()
            } else {
                None
            };
            write_output(&output, &header, &records, index_dictionary)?;
            Ok(())
        }
        Err(Failure::Caught(exception)) => {
            log_error(
                "GatherVcfs",
                &format!("There was a problem with gathering the INPUT.{exception}"),
            );
            let _ = std::fs::remove_file(&output);
            std::process::exit(1);
        }
        Err(Failure::Escaped(exception)) => die(&exception),
    }
}
