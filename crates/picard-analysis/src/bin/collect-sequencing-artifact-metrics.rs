//! `CollectSequencingArtifactMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.artifacts.CollectSequencingArtifactMetrics` at tag 3.4.0, around the
//! counters in `picard_analysis::collect_sequencing_artifact_metrics`:
//!
//! * `customCommandLineValidation`, every message at once, the context lengths in the `HashSet`
//!   order `CONTEXTS_TO_PRINT` iterates in;
//! * `SinglePassSamProgram.makeItSo`: the sort check, the reference walker that refuses to
//!   rewind, the stop at the first unmapped record (`usesNoRefReads` is false) or at `STOP_AFTER`;
//! * `setup`'s record filters, interval mask and dbSNP mask, and one counter per library, keyed
//!   in a `HashMap` whose iteration order is the order the libraries are written in;
//! * `finish`: four files from the counters, then the error summary folded from the pre-adapter
//!   details that were printed.

use picard_analysis::collect_sequencing_artifact_metrics::{
    error_summaries, file_names, ArtifactCounter, CounterMetrics,
};
use picard_analysis::java_hash_map::JavaHashMap;
use picard_analysis::metrics_cli::{
    check_coordinate_sorted, original_qualities, read_group, read_input, read_interval_list,
    refuse_validation, thrown, Args, DbSnpMask, IntervalMask, ReferenceWalker,
};

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_metrics::file::MetricsFile;

const TOOL: &str = "CollectSequencingArtifactMetrics";
const UNKNOWN_LIBRARY: &str = "UnknownLibrary";
const UNKNOWN_SAMPLE: &str = "UnknownSample";

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const REVERSE: u16 = 0x10;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;

fn new_metrics_file() -> MetricsFile {
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    file
}

fn write(path: &str, file: &MetricsFile) {
    if let Err(e) = std::fs::write(path, file.write()) {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Could not write to file {path}: {e}"
        ));
    }
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("R", "REFERENCE_SEQUENCE"),
        ("Q", "MINIMUM_QUALITY_SCORE"),
        ("MQ", "MINIMUM_MAPPING_QUALITY"),
        ("MIN_INS", "MINIMUM_INSERT_SIZE"),
        ("MAX_INS", "MAXIMUM_INSERT_SIZE"),
        ("UNPAIRED", "INCLUDE_UNPAIRED"),
        ("DUPES", "INCLUDE_DUPLICATES"),
        ("NON_PF", "INCLUDE_NON_PF_READS"),
        ("TANDEM", "TANDEM_READS"),
        ("EXT", "FILE_EXTENSION"),
        ("AS", "ASSUME_SORTED"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let reference = args.required("REFERENCE_SEQUENCE");
    let intervals_path = args.get("INTERVALS").map(str::to_string);
    let db_snp_path = args.get("DB_SNP").map(str::to_string);
    let minimum_quality = args.int("MINIMUM_QUALITY_SCORE", 20);
    let minimum_mapping_quality = args.int("MINIMUM_MAPPING_QUALITY", 30);
    let minimum_insert = args.int("MINIMUM_INSERT_SIZE", 60);
    let maximum_insert = args.int("MAXIMUM_INSERT_SIZE", 600);
    let include_unpaired = args.bool("INCLUDE_UNPAIRED", false);
    let include_duplicates = args.bool("INCLUDE_DUPLICATES", false);
    let include_non_pf = args.bool("INCLUDE_NON_PF_READS", false);
    let tandem = args.bool("TANDEM_READS", false);
    let use_oq = args.bool("USE_OQ", true);
    let context_size = args.int("CONTEXT_SIZE", 1);
    let extension = args.get("FILE_EXTENSION").map(str::to_string);
    let assume_sorted = args.bool("ASSUME_SORTED", true);
    let stop_after = args.int("STOP_AFTER", 0);

    // `CONTEXTS_TO_PRINT` is a `HashSet<String>`: its validation messages and its membership
    // test both go through it, and only the first reaches an output in an order.
    let mut contexts_to_print: JavaHashMap<()> = JavaHashMap::new();
    for context in args.all("CONTEXTS_TO_PRINT") {
        contexts_to_print.put(&context, ());
    }

    // customCommandLineValidation.
    let mut messages = Vec::new();
    let full_length = 2 * context_size + 1;
    if context_size < 0 {
        messages.push("CONTEXT_SIZE cannot be negative".to_string());
    }
    for (context, _) in contexts_to_print.iter() {
        if context.chars().count() as i64 != full_length {
            messages.push(format!(
                "Context {context} is not the length implied by CONTEXT_SIZE: {full_length}"
            ));
        }
    }
    if minimum_insert < 0 {
        messages.push("MINIMUM_INSERT_SIZE cannot be negative".to_string());
    }
    if maximum_insert < 0 {
        messages.push("MAXIMUM_INSERT_SIZE cannot be negative".to_string());
    }
    if maximum_insert > 0 && maximum_insert < minimum_insert {
        messages.push(
            "MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE unless set to 0"
                .to_string(),
        );
    }
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }
    let context_size = context_size as usize;

    // makeItSo.
    let (header, records) = read_input(&input);
    let contigs = read_fasta_file(&reference)
        .unwrap_or_else(|e| picard_analysis::metrics_cli::fail(&format!("{e:?}")));
    check_coordinate_sorted(&input, &header, assume_sorted);

    // setup.
    let names = file_names(&output, extension.as_deref());
    let mut samples: JavaHashMap<()> = JavaHashMap::new();
    let mut libraries: JavaHashMap<()> = JavaHashMap::new();
    for group in &header.read_groups {
        samples.put(group.attributes.get("SM").unwrap_or(UNKNOWN_SAMPLE), ());
        libraries.put(group.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY), ());
    }
    let intervals = intervals_path.as_deref().map(read_interval_list);
    let interval_mask = intervals.as_ref().map(IntervalMask::new);
    let db_snp = db_snp_path
        .as_deref()
        .map(|path| DbSnpMask::load(path, intervals.as_ref()));
    let maximum_insert = if maximum_insert == 0 {
        i64::from(i32::MAX)
    } else {
        maximum_insert
    };
    let sample_alias = samples
        .iter()
        .map(|(s, _)| s.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let mut counters: JavaHashMap<ArtifactCounter> = JavaHashMap::new();
    for (library, _) in libraries.iter() {
        counters.put(
            library,
            ArtifactCounter::new(&sample_alias, library, context_size, tandem),
        );
    }
    let mut counters: Vec<(String, ArtifactCounter)> = {
        let order: Vec<String> = counters.iter().map(|(k, _)| k.to_string()).collect();
        let mut out = Vec::new();
        for key in order {
            let counter = counters.remove(&key).expect("counter");
            out.push((key, counter));
        }
        out
    };

    let filtered = |flags: u16, mapq: u8, insert: i32| -> bool {
        if !include_non_pf && flags & QC_FAIL != 0 {
            return true;
        }
        if flags & SECONDARY != 0 {
            return true;
        }
        if !include_duplicates && flags & DUPLICATE != 0 {
            return true;
        }
        if flags & UNMAPPED != 0 {
            return true;
        }
        if i64::from(mapq) < minimum_mapping_quality {
            return true;
        }
        if !include_unpaired {
            if flags & PAIRED == 0 {
                return true;
            }
            let ins = i64::from(insert.unsigned_abs());
            if ins < minimum_insert || ins > maximum_insert {
                return true;
            }
        }
        false
    };

    let mut walker = ReferenceWalker::default();
    let mut upper: Option<(i32, Vec<u8>)> = None;
    let mut count = 0i64;
    for record in &records {
        if record.reference_index != -1 {
            walker.get(record.reference_index);
        }
        if !filtered(
            record.flags,
            record.mapping_quality,
            record.inferred_insert_size,
        ) {
            let library = read_group(&header, record)
                .map(|g| g.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY))
                .unwrap_or(UNKNOWN_LIBRARY);
            if !libraries.contains_key(library) {
                thrown(&format!(
                    "picard.PicardException: Record contains library that is missing from header: {library}"
                ));
            }
            let index = record.reference_index;
            if upper.as_ref().is_none_or(|(i, _)| *i != index) {
                let bases = contigs[index as usize].bases.to_ascii_uppercase();
                upper = Some((index, bases));
            }
            let reference_bases = &upper.as_ref().expect("reference").1;
            let contig_name = &contigs[index as usize].name;
            let quals = if use_oq {
                original_qualities(record).unwrap_or_else(|| record.base_qualities.clone())
            } else {
                record.base_qualities.clone()
            };
            let counter = &mut counters
                .iter_mut()
                .find(|(l, _)| l == library)
                .expect("library")
                .1;
            let negative = record.flags & REVERSE != 0;
            let paired = record.flags & PAIRED != 0;
            let second = record.flags & SECOND != 0;
            let full_length = 2 * context_size + 1;
            for block in
                htsjdk_bam::alignment_block::alignment_blocks(&record.cigar, record.alignment_start)
            {
                for offset in 0..block.length {
                    let read_pos = (block.read_start + offset) as usize;
                    let ref_pos = block.reference_start + offset;
                    let qual = quals[read_pos - 1] as i8;
                    if i64::from(qual) < minimum_quality {
                        continue;
                    }
                    let read_base = record.read_bases[read_pos - 1].to_ascii_uppercase();
                    if read_base == b'N' {
                        continue;
                    }
                    if let Some(mask) = &interval_mask {
                        if !mask.get(index, ref_pos) {
                            continue;
                        }
                    }
                    if let Some(db) = &db_snp {
                        if db.is_db_snp_site(contig_name, ref_pos) {
                            continue;
                        }
                    }
                    let start = i64::from(ref_pos) - context_size as i64 - 1;
                    if start < 0 || start as usize + full_length > reference_bases.len() {
                        continue;
                    }
                    let start = start as usize;
                    let context = &reference_bases[start..start + full_length];
                    if context.contains(&b'N') {
                        continue;
                    }
                    if !matches!(read_base, b'A' | b'C' | b'G' | b'T' | b'N') {
                        continue;
                    }
                    let context = std::str::from_utf8(context).unwrap_or("");
                    counter.count(context, read_base, negative, paired, second);
                }
            }
        }
        count += 1;
        if stop_after > 0 && count >= stop_after {
            break;
        }
        if record.reference_index == -1 {
            break;
        }
    }

    // finish.
    let mut pre_summary = new_metrics_file();
    let mut pre_detail = new_metrics_file();
    let mut bait_summary = new_metrics_file();
    let mut bait_detail = new_metrics_file();
    let mut error_summary = new_metrics_file();
    let printed =
        |context: &str| contexts_to_print.is_empty() || contexts_to_print.contains_key(context);
    let mut printed_pre = Vec::new();
    for (_, counter) in counters {
        let CounterMetrics {
            pre_adapter_summary,
            pre_adapter_detail,
            bait_bias_summary,
            bait_bias_detail,
        } = counter.finish();
        for m in &pre_adapter_summary {
            pre_summary.add_metric(m);
        }
        for m in &bait_bias_summary {
            bait_summary.add_metric(m);
        }
        for m in pre_adapter_detail {
            if printed(&m.context) {
                pre_detail.add_metric(&m);
                printed_pre.push(m);
            }
        }
        for m in &bait_bias_detail {
            if printed(&m.context) {
                bait_detail.add_metric(m);
            }
        }
    }
    write(&names[1], &pre_detail);
    write(&names[0], &pre_summary);
    write(&names[3], &bait_detail);
    write(&names[2], &bait_summary);
    for m in error_summaries(&printed_pre) {
        error_summary.add_metric(&m);
    }
    write(&names[4], &error_summary);
}
