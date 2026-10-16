//! `CollectOxoGMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectOxoGMetrics.doWork` at tag 3.4.0 around the counters in
//! `picard_analysis::collect_oxo_g_metrics`:
//!
//! * `customCommandLineValidation`, every message at once, the contexts in `HashSet` order;
//! * the read-group requirement, then dbSNP, then the locus iterator, whose constructor refuses a
//!   header sorted any way but by coordinate;
//! * htsjdk 4.2.0's `SamLocusIterator` as this tool configures it: no quality cutoff, the
//!   mapping-quality and PF cutoffs, the tool's own record filters in place of the default ones
//!   (so supplementary records are counted), only covered loci, and with `INTERVALS` only the
//!   loci inside them. An emitted locus carries every aligned base of every passing record;
//! * the rows in the iteration order of the `ListMap` keyed by context, libraries in the order of
//!   the `HashSet` they were read into.

use std::collections::BTreeMap;

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::header::SamHeader;
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::collect_oxo_g_metrics::{finish, Counts, Metrics};
use picard_analysis::java_hash_map::JavaHashMap;
use picard_analysis::metrics_cli::{
    fail, original_qualities, read_group, read_input, read_interval_list, refuse_validation,
    sort_order, thrown, Args, DbSnpMask, IntervalMask,
};

const TOOL: &str = "CollectOxoGMetrics";
const UNKNOWN_LIBRARY: &str = "UnknownLibrary";
const UNKNOWN_SAMPLE: &str = "UnknownSample";

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const REVERSE: u16 = 0x10;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;

struct Row {
    sample_alias: String,
    library: String,
    context: String,
    sites: i64,
    metrics: Metrics,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        "picard.analysis.CollectOxoGMetrics$CpcgMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "SAMPLE_ALIAS",
            "LIBRARY",
            "CONTEXT",
            "TOTAL_SITES",
            "TOTAL_BASES",
            "REF_NONOXO_BASES",
            "REF_OXO_BASES",
            "REF_TOTAL_BASES",
            "ALT_NONOXO_BASES",
            "ALT_OXO_BASES",
            "OXIDATION_ERROR_RATE",
            "OXIDATION_Q",
            "C_REF_REF_BASES",
            "G_REF_REF_BASES",
            "C_REF_ALT_BASES",
            "G_REF_ALT_BASES",
            "C_REF_OXO_ERROR_RATE",
            "C_REF_OXO_Q",
            "G_REF_OXO_ERROR_RATE",
            "G_REF_OXO_Q",
        ]
    }
    fn values(&self) -> Vec<Value> {
        let m = &self.metrics;
        vec![
            Value::Str(self.sample_alias.clone()),
            Value::Str(self.library.clone()),
            Value::Str(self.context.clone()),
            Value::Long(self.sites),
            Value::Long(m.total_bases),
            Value::Long(m.ref_nonoxo_bases),
            Value::Long(m.ref_oxo_bases),
            Value::Long(m.ref_total_bases),
            Value::Long(m.alt_nonoxo_bases),
            Value::Long(m.alt_oxo_bases),
            Value::Double(m.oxidation_error_rate),
            Value::Double(m.oxidation_q),
            Value::Long(m.c_ref_ref_bases),
            Value::Long(m.g_ref_ref_bases),
            Value::Long(m.c_ref_alt_bases),
            Value::Long(m.g_ref_alt_bases),
            Value::Double(m.c_ref_oxo_error_rate),
            Value::Double(m.c_ref_oxo_q),
            Value::Double(m.g_ref_oxo_error_rate),
            Value::Double(m.g_ref_oxo_q),
        ]
    }
}

/// One `Calculator`: a library, a context, and what it has counted.
struct Calculator {
    library: String,
    context: String,
    sites: i64,
    counts: Counts,
}

/// `SequenceUtil.reverseComplement`, which upper-cases as it complements.
fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|b| htsjdk_bam::sequence::complement(b.to_ascii_uppercase()))
        .collect()
}

/// The library a record's read group names, or the NullPointerException a record without one
/// raises in `computeAlleleFraction`.
fn library_of<'a>(header: &'a SamHeader, record: &htsjdk_bam::record::BamRecord) -> &'a str {
    match read_group(header, record) {
        Some(group) => group.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY),
        None => thrown(
            "java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.samtools.SAMReadGroupRecord.getLibrary()\" because the return value of \
             \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null",
        ),
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
        ("NON_PF", "INCLUDE_NON_PF_READS"),
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
    let include_non_pf = args.bool("INCLUDE_NON_PF_READS", true);
    let use_oq = args.bool("USE_OQ", true);
    let context_size = args.int("CONTEXT_SIZE", 1);
    let stop_after = args.int("STOP_AFTER", i64::from(i32::MAX));

    let mut requested: JavaHashMap<()> = JavaHashMap::new();
    for context in args.all("CONTEXTS") {
        requested.put(&context, ());
    }

    // customCommandLineValidation.
    let size = 1 + 2 * context_size;
    let mut messages = Vec::new();
    for (context, _) in requested.iter() {
        if context.chars().count() as i64 != size {
            messages.push(format!(
                "Context {context} is not {size} long as implied by CONTEXT_SIZE={context_size}"
            ));
        } else if context.as_bytes()[context.len() / 2] != b'C' {
            messages.push(format!(
                "Middle base of context sequence {context} must be C"
            ));
        }
    }
    if minimum_insert < 0 {
        messages.push("MINIMUM_INSERT_SIZE cannot be negative".to_string());
    }
    if maximum_insert < 0 {
        messages.push("MAXIMUM_INSERT_SIZE cannot be negative".to_string());
    }
    if maximum_insert < minimum_insert {
        messages.push("MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE".to_string());
    }
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }
    let context_size = context_size.max(0) as usize;

    // doWork.
    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| fail(&format!("{e:?}")));
    let (header, records) = read_input(&input);
    if header.read_groups.is_empty() {
        thrown(
            "picard.PicardException: This analysis requires a read group entry in the alignment \
             file header",
        );
    }
    let mut samples: JavaHashMap<()> = JavaHashMap::new();
    let mut libraries: JavaHashMap<()> = JavaHashMap::new();
    for group in &header.read_groups {
        samples.put(group.attributes.get("SM").unwrap_or(UNKNOWN_SAMPLE), ());
        libraries.put(group.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY), ());
    }

    // The contexts: the requested set, or makeContextStrings' set of every kmer centred on C,
    // inserted in `generateAllKmers` order.
    let contexts: JavaHashMap<()> = if requested.is_empty() {
        let mut all: Vec<String> = vec![String::new()];
        for _ in 0..size {
            all = all
                .iter()
                .flat_map(|p| b"ACGT".iter().map(move |&b| format!("{p}{}", b as char)))
                .collect();
        }
        let mut set = JavaHashMap::new();
        for kmer in all {
            if kmer.as_bytes()[context_size] == b'C' {
                set.put(&kmer, ());
            }
        }
        set
    } else {
        requested
    };
    // `ListMap<String, Calculator>`: a HashMap keyed by context, each list in library order.
    let mut calculators: JavaHashMap<Vec<Calculator>> = JavaHashMap::new();
    for (context, _) in contexts.iter() {
        let list = libraries
            .iter()
            .map(|(library, _)| Calculator {
                library: library.to_string(),
                context: context.to_string(),
                sites: 0,
                counts: Counts::default(),
            })
            .collect();
        calculators.put(context, list);
    }
    let mut calculators: Vec<(String, Vec<Calculator>)> = {
        let order: Vec<String> = calculators.iter().map(|(k, _)| k.to_string()).collect();
        order
            .into_iter()
            .map(|k| {
                let v = calculators.remove(&k).expect("context");
                (k, v)
            })
            .collect()
    };

    let db_snp = db_snp_path.as_deref().map(|p| DbSnpMask::load(p, None));

    // The locus iterator's constructor.
    let order = sort_order(&header);
    if order != "unsorted" && order != "coordinate" {
        thrown(
            "htsjdk.samtools.SAMException: SamLocusIterator cannot operate on a SAM file that is \
             not coordinate sorted.",
        );
    }
    let intervals = intervals_path.as_deref().map(read_interval_list);
    let mask = intervals.as_ref().map(IntervalMask::new);

    // The pileups: every covered locus, with the record and read offset of each aligned base.
    let insert_filter = minimum_insert > 0 || maximum_insert > 0;
    let mut pileups: BTreeMap<(i32, i32), Vec<(usize, usize)>> = BTreeMap::new();
    for (index, record) in records.iter().enumerate() {
        let flags = record.flags;
        // The tool's filters, which replace the iterator's defaults.
        if flags & SECONDARY != 0 || flags & DUPLICATE != 0 {
            continue;
        }
        if insert_filter {
            if flags & PAIRED == 0 {
                continue;
            }
            let ins = i64::from(record.inferred_insert_size.unsigned_abs());
            if ins < minimum_insert || ins > maximum_insert {
                continue;
            }
        }
        if record.reference_index == -1 {
            break;
        }
        if flags & UNMAPPED != 0
            || i64::from(record.mapping_quality) < minimum_mapping_quality
            || (!include_non_pf && flags & QC_FAIL != 0)
        {
            continue;
        }
        for block in alignment_blocks(&record.cigar, record.alignment_start) {
            for i in 0..block.length {
                let offset = (block.read_start + i - 1) as usize;
                pileups
                    .entry((record.reference_index, block.reference_start + i))
                    .or_default()
                    .push((index, offset));
            }
        }
    }

    let mut sites = 0i64;
    for (&(sequence, position), pile) in &pileups {
        if let Some(m) = &mask {
            if !m.get(sequence, position) {
                continue;
            }
        }
        let contig = &contigs[sequence as usize];
        if let Some(db) = &db_snp {
            if db.is_db_snp_site(&contig.name, position) {
                continue;
            }
        }
        let bases = &contig.bases;
        let pos = position as usize;
        if pos <= context_size || pos > bases.len() - context_size {
            continue;
        }
        let index = pos - 1;
        let base = bases[index].to_ascii_uppercase();
        if base != b'C' && base != b'G' {
            continue;
        }
        let window: Vec<u8> = bases[index - context_size..=index + context_size]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let context = if base == b'C' {
            window
        } else {
            reverse_complement(&window)
        };
        let context = String::from_utf8_lossy(&context).to_string();
        let Some((_, list)) = calculators.iter_mut().find(|(c, _)| *c == context) else {
            continue;
        };
        for calculator in list.iter_mut() {
            // computeAlleleFraction.
            let mut counts = [0i64; 4]; // controlA, oxidatedA, controlC, oxidatedC
            let alternate = if base == b'C' { b'A' } else { b'T' };
            for &(r, offset) in pile {
                let record = &records[r];
                let quality = if use_oq {
                    original_qualities(record)
                        .map(|q| q[offset])
                        .unwrap_or(record.base_qualities[offset])
                } else {
                    record.base_qualities[offset]
                };
                if i64::from(quality as i8) < minimum_quality {
                    continue;
                }
                if calculator.library != library_of(&header, record) {
                    continue;
                }
                let read_base = record.read_bases[offset];
                let as_read = if record.flags & REVERSE != 0 {
                    htsjdk_bam::sequence::complement(read_base)
                } else {
                    read_base
                };
                let read = if record.flags & PAIRED != 0 && record.flags & SECOND != 0 {
                    2
                } else {
                    1
                };
                if read_base == base {
                    match (as_read, read) {
                        (b'G', 1) | (b'C', 2) => counts[3] += 1,
                        (b'G', 2) | (b'C', 1) => counts[2] += 1,
                        _ => {}
                    }
                } else if read_base == alternate {
                    match (as_read, read) {
                        (b'T', 1) | (b'A', 2) => counts[1] += 1,
                        (b'T', 2) | (b'A', 1) => counts[0] += 1,
                        _ => {}
                    }
                }
            }
            if counts.iter().sum::<i64>() > 0 {
                calculator.sites += 1;
                let c = &mut calculator.counts;
                if base == b'C' {
                    c.ref_c_control_a += counts[0];
                    c.ref_c_oxidated_a += counts[1];
                    c.ref_c_control_c += counts[2];
                    c.ref_c_oxidated_c += counts[3];
                } else {
                    c.ref_g_control_a += counts[0];
                    c.ref_g_oxidated_a += counts[1];
                    c.ref_g_control_c += counts[2];
                    c.ref_g_oxidated_c += counts[3];
                }
            }
        }
        sites += 1;
        if sites >= stop_after {
            break;
        }
    }

    let sample_alias = samples
        .iter()
        .map(|(s, _)| s.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for (_, list) in calculators.drain(..) {
        for calculator in list {
            file.add_metric(&Row {
                sample_alias: sample_alias.clone(),
                library: calculator.library,
                context: calculator.context,
                sites: calculator.sites,
                metrics: finish(&calculator.counts),
            });
        }
    }
    if let Err(e) = std::fs::write(&output, file.write()) {
        fail(&format!("{e}"));
    }
}
