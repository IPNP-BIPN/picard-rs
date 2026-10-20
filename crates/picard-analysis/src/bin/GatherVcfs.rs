//! `GatherVcfs` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.GatherVcfs.doWork` at tag 3.4.0. The conventional gather reads every record
//! and writes them through `picard_analysis::vcf_io`, a `.vcf.gz` as BGZF with its `.tbi`. The
//! block-copying one is taken instead when every input AND the output have a block-compressed
//! name: no record is parsed, no `COMMENT` is added and no index is written, whatever
//! `CREATE_INDEX` says.
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
use htsjdk_vcf::encoder::MissingFields;
use htsjdk_vcf::reader::VcfFile;
use htsjdk_vcf::variant::VariantContext;
use picard_analysis::args::arg;
use picard_analysis::vcf_io::{
    add_other_meta_data_line, assert_same_dictionary, die, header_dictionary, is_compressed_vcf,
    java_compare, log_error, thrown, unroll_paths, Record, VcfReader, VcfWriter, WriteError,
};

/// One input as `assertSameSamplesAndValidOrdering` opens it: its header and its first record. The
/// records are read, and their failures met, as `gatherConventionally` writes them.
struct Scanned {
    file: VcfFile,
    first: Option<Record>,
}

/// `new VCFFileReader(path, false)` and its first record, or the exception that read gave.
fn scan(path: &str) -> Result<Scanned, String> {
    let mut reader = VcfReader::open(path)?;
    let first = reader.try_next().transpose()?;
    Ok(Scanned {
        file: reader.file,
        first,
    })
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
    mut files: Vec<(String, Scanned)>,
    reorder: bool,
) -> Result<Vec<(String, Scanned)>, Failure> {
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
            let (first_a, first_b) = (a.1.first.as_ref(), b.1.first.as_ref());
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
        if let Some(current) = &vcf.first {
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

/// `gatherConventionally`, up to the writing: the header, the comparator its order check uses, and
/// the files in the order their records are written.
fn gather_conventionally(
    files: Vec<(String, Scanned)>,
    comments: &[String],
) -> Result<
    (
        htsjdk_vcf::header::VcfHeader,
        VariantContextComparator,
        Vec<String>,
    ),
    Failure,
> {
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
    Ok((
        header,
        comparator,
        files.into_iter().map(|(path, _)| path).collect(),
    ))
}

/// `gatherConventionally`'s loop: each file's records written as they are read, the first checked
/// against the last one written, as Picard does on the stream. A failure is the caught one: the
/// caller logs it and deletes the output. A write error is the writer's, thrown as it was.
fn write_conventionally(
    writer: &mut VcfWriter,
    comparator: &VariantContextComparator,
    paths: &[String],
) -> Result<Result<(), WriteError>, Failure> {
    let mut last: Option<Record> = None;
    let mut last_file: Option<&String> = None;
    for path in paths {
        // Opened once already, so a failure now is the file changing underneath.
        let mut reader = VcfReader::open(path).unwrap_or_else(|exception| die(&exception));
        let mut next = reader.try_next().transpose().map_err(Failure::Caught)?;
        if let (Some(last_context), Some(first)) = (&last, &next) {
            if compare(comparator, &first.variant, &last_context.variant)? <= 0 {
                let (vc, lc) = (&first.variant, &last_context.variant);
                return Err(Failure::Caught(format!(
                    "java.lang.IllegalArgumentException: First variant in file {} is at {}:{} but \
                     last variant in earlier file {} is at {}:{}",
                    absolute(path),
                    vc.contig,
                    vc.start,
                    absolute(last_file.map(String::as_str).unwrap_or_default()),
                    lc.contig,
                    lc.start
                )));
            }
        }
        while let Some(record) = next {
            if let Err(error) = writer.add(&record).or_else(thrown) {
                return Ok(Err(error));
            }
            last = Some(record);
            next = reader.try_next().transpose().map_err(Failure::Caught)?;
        }
        // Set for every file, an empty one included.
        last_file = Some(path);
    }
    Ok(Ok(()))
}

/// `IOUtil.hasBlockCompressedExtension`, and not `VCFFileReader.isBCF`: what `areAllBlockCompressed`
/// asks of each name, without opening the file.
fn block_compressed(path: &str) -> bool {
    let lower = path.to_lowercase();
    !path.ends_with(".bcf")
        && [".gz", ".gzip", ".bgz", ".bgzf"]
            .iter()
            .any(|extension| lower.ends_with(extension))
}

/// `gatherWithBlockCopying`: the first file's blocks as they are, and of every later file the rest
/// of the block its header ends in, deflated again as one new block, then its remaining blocks as
/// they are. Each file's terminator is left out and one is written at the end.
fn gather_with_block_copying(paths: &[String]) -> Result<Vec<u8>, Failure> {
    use htsjdk_bgzf::{check_termination, BgzfReader, FileTermination, EMPTY_GZIP_BLOCK};
    use std::io::Write;

    let io = |e: std::io::Error| {
        Failure::Caught(format!("htsjdk.samtools.util.RuntimeIOException: {e}"))
    };
    let mut out = Vec::new();
    for (i, path) in paths.iter().enumerate() {
        let data = std::fs::read(path).map_err(io)?;
        let termination = check_termination(&data);
        if termination == FileTermination::Defective {
            return Err(Failure::Caught(format!(
                "picard.PicardException: {} does not have a valid GZIP block at the end of the \
                 file.",
                absolute(path)
            )));
        }
        let mut position = 0;
        if i > 0 {
            let mut reader = BgzfReader::new(&data[..]);
            let mut last_byte_newline = true;
            // `blockIn.available()` reads the next block, and is zero at an empty one, which ends
            // the scan with the stream positioned after it.
            while let Some(block) = reader
                .next_block()
                .map_err(|e| Failure::Caught(format!("htsjdk.samtools.SAMException: {e:?}")))?
            {
                position = block.block_address as usize + block.block_compressed_size;
                if block.data.is_empty() {
                    break;
                }
                let first = block.data.iter().position(|&b| {
                    let newline = b == b'\n' || b == b'\r';
                    let found = last_byte_newline && !newline && b != b'#';
                    last_byte_newline = newline;
                    found
                });
                if let Some(first) = first {
                    let mut block_out = picard_analysis::output::bgzf_writer(Vec::new());
                    block_out.write_all(&block.data[first..]).map_err(io)?;
                    out.extend(block_out.into_inner_without_terminator().map_err(io)?);
                    break;
                }
            }
        }
        let skip_last = if termination == FileTermination::HasTerminatorBlock {
            EMPTY_GZIP_BLOCK.len()
        } else {
            0
        };
        let end = data.len() - skip_last;
        if end > position {
            out.extend_from_slice(&data[position..end]);
        }
    }
    out.extend_from_slice(&EMPTY_GZIP_BLOCK);
    Ok(out)
}

/// What the `try` block produced: records to write, or the bytes the block copy assembled.
enum Gathered {
    Conventionally(
        htsjdk_vcf::header::VcfHeader,
        VariantContextComparator,
        Vec<String>,
    ),
    WithBlockCopying(Vec<u8>),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("GatherVcfs");
    picard_analysis::output::install();
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
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(false);
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| picard_analysis::args::is_true(&value))
        .unwrap_or(true);
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }
    let paths = unroll_paths(&inputs).unwrap_or_else(|exception| die(&exception));

    // `VCFFileReader.getSequenceDictionary(unrolledPaths.get(0))`, outside the `try`.
    let first = scan(&paths[0]).unwrap_or_else(|exception| die(&exception));
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
            let vcf = scan(path).map_err(Failure::Caught)?;
            files.push((path.clone(), vcf));
        }
        let files = assert_same_samples_and_valid_ordering(files, reorder)?;
        let paths: Vec<String> = files.iter().map(|(path, _)| path.clone()).collect();
        if paths.iter().all(|path| block_compressed(path)) && block_compressed(&output) {
            if create_index {
                eprintln!(
                    "WARN\tGatherVcfs\tIndex creation not currently supported when gathering \
                     block compressed VCFs."
                );
            }
            return gather_with_block_copying(&paths).map(Gathered::WithBlockCopying);
        }
        // `gatherConventionally` builds its writer before it reads a record, and the builder
        // refuses a name it cannot place inside the `try`.
        if output.ends_with(".bcf") {
            return Err(Failure::Escaped("a .bcf OUTPUT is not ported".to_string()));
        }
        if !output.ends_with(".vcf") && !is_compressed_vcf(&output) {
            return Err(Failure::Caught(
                "java.lang.IllegalArgumentException: Output format type is not set, or could not \
                 be inferred from the output path. If a path was used, does it have a valid VCF \
                 extension (.vcf, .vcf.gz, .vcf.bgz, .bcf)?"
                    .to_string(),
            ));
        }
        let (header, comparator, paths) = gather_conventionally(files, &comments)?;
        Ok(Gathered::Conventionally(header, comparator, paths))
    })();

    match gathered {
        Ok(Gathered::WithBlockCopying(bytes)) => {
            std::fs::write(&output, bytes)?;
            Ok(())
        }
        Ok(Gathered::Conventionally(header, comparator, paths)) => {
            let index_dictionary = if create_index {
                sequence_dictionary.as_deref()
            } else {
                None
            };
            let mut writer =
                VcfWriter::create(&output, &header, index_dictionary, MissingFields::default())
                    .or_else(thrown)?;
            match write_conventionally(&mut writer, &comparator, &paths) {
                Ok(Ok(())) => writer.finish().or_else(thrown)?,
                Ok(Err(error)) => return Err(error.into()),
                Err(failure) => {
                    drop(writer);
                    fail(failure, &output);
                }
            }
            Ok(())
        }
        Err(failure) => fail(failure, &output),
    }
}

/// The `catch`: a caught exception is logged and the output deleted, exit status 1; an escaped one
/// is thrown out of `doWork`.
fn fail(failure: Failure, output: &str) -> ! {
    match failure {
        Failure::Caught(exception) => {
            log_error(
                "GatherVcfs",
                &format!("There was a problem with gathering the INPUT.{exception}"),
            );
            let _ = std::fs::remove_file(output);
            std::process::exit(1);
        }
        Failure::Escaped(exception) => die(&exception),
    }
}
