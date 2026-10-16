//! `CollectSamErrorMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.SamErrorMetric.CollectSamErrorMetrics.doWork` at tag 3.4.0 around the counters
//! in `picard_analysis::collect_sam_error_metrics`:
//!
//! * `--ERROR_METRICS` is a collection with a default, so each value is APPENDED to the twenty-seven
//!   the tool carries and `null` empties the list first; the aggregators are built before any file
//!   is read, and two with the same suffix are refused;
//! * `customCommandLineValidation`, every message at once;
//! * the VCF must be indexed, in both of the tool's modes, and `INTERVAL_ITERATOR` without a VCF
//!   dereferences a reader that was never opened at the first locus it looks at;
//! * the sort-order refusal, and `SamLocusIterator` as the tool configures it: indels included,
//!   only covered loci, the two cutoffs, and with `INTERVALS` only the loci inside every list;
//! * the loop: a draw from `Random(42)` per locus against `PROBABILITY`, the known sites taken out,
//!   `MAX_LOCI` counted over what is left, then one file per aggregator, named
//!   `<OUTPUT>.<suffix><FILE_EXTENSION>`.
//!
//! The reference and the reads are placed on one concatenated coordinate (every contig after the
//! one before it, with one position between them), which is what lets the single-contig counters
//! of the library walk a file of several contigs without knowing.

use htsjdk_bam::cigar::Op;
use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::collect_sam_error_metrics::{
    aggregation_suffix, collect_joint, pileup, processed_loci, suffixes, Calculator, Options, Read,
    Stratifier, Table, DEFAULT_ERROR_METRICS,
};
use picard_analysis::metrics_cli::{
    fail, read_input, read_interval_list, refuse_validation, sort_order, thrown, Args, IntervalMask,
};

const TOOL: &str = "CollectSamErrorMetrics";

/// One row of any of the three tables.
struct Row {
    class: &'static str,
    columns: &'static [&'static str],
    values: Vec<Value>,
}

impl MetricBean for Row {
    fn class_name(&self) -> &str {
        self.class
    }
    fn columns(&self) -> &[&'static str] {
        self.columns
    }
    fn values(&self) -> Vec<Value> {
        self.values.clone()
    }
}

const BASE_COLUMNS: &[&str] = &["ERROR_BASES", "Q_SCORE", "COVARIATE", "TOTAL_BASES"];
const OVERLAPPING_COLUMNS: &[&str] = &[
    "NUM_BASES_WITH_OVERLAPPING_READS",
    "NUM_DISAGREES_WITH_REFERENCE_ONLY",
    "DISAGREES_WITH_REFERENCE_ONLY_Q",
    "NUM_DISAGREES_WITH_REF_AND_MATE",
    "DISAGREES_WITH_REF_AND_MATE_ONLY_Q",
    "NUM_THREE_WAYS_DISAGREEMENT",
    "THREE_WAYS_DISAGREEMENT_ONLY_Q",
    "COVARIATE",
    "TOTAL_BASES",
];
const INDEL_COLUMNS: &[&str] = &[
    "NUM_INSERTIONS",
    "NUM_INSERTED_BASES",
    "INSERTIONS_Q",
    "NUM_DELETIONS",
    "NUM_DELETED_BASES",
    "DELETIONS_Q",
    "ERROR_BASES",
    "Q_SCORE",
    "COVARIATE",
    "TOTAL_BASES",
];

fn long(value: u64) -> Value {
    Value::Long(value as i64)
}

fn int(value: i32) -> Value {
    Value::Long(i64::from(value))
}

/// The metrics file one table is written as.
fn file_of(table: &Table) -> MetricsFile {
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    match table {
        Table::Base(rows) => {
            for row in rows {
                file.add_metric(&Row {
                    class: "picard.sam.SamErrorMetric.BaseErrorMetric",
                    columns: BASE_COLUMNS,
                    values: vec![
                        long(row.error_bases),
                        int(row.q_score),
                        Value::Str(row.covariate.clone()),
                        long(row.total_bases),
                    ],
                });
            }
        }
        Table::Overlapping(rows) => {
            for row in rows {
                file.add_metric(&Row {
                    class: "picard.sam.SamErrorMetric.OverlappingErrorMetric",
                    columns: OVERLAPPING_COLUMNS,
                    values: vec![
                        long(row.bases_with_overlapping_reads),
                        long(row.disagrees_with_reference_only),
                        int(row.disagrees_with_reference_only_q),
                        long(row.disagrees_with_ref_and_mate),
                        int(row.disagrees_with_ref_and_mate_q),
                        long(row.three_ways_disagreement),
                        int(row.three_ways_disagreement_q),
                        Value::Str(row.covariate.clone()),
                        long(row.total_bases),
                    ],
                });
            }
        }
        Table::Indel(rows) => {
            for row in rows {
                file.add_metric(&Row {
                    class: "picard.sam.SamErrorMetric.IndelErrorMetric",
                    columns: INDEL_COLUMNS,
                    values: vec![
                        long(row.insertions),
                        long(row.inserted_bases),
                        int(row.insertions_q),
                        long(row.deletions),
                        long(row.deleted_bases),
                        int(row.deletions_q),
                        long(row.error_bases),
                        int(row.q_score),
                        Value::Str(row.covariate.clone()),
                        long(row.total_bases),
                    ],
                });
            }
        }
    }
    file
}

/// One directive, `ERROR(:STRATIFIER)*`, as the calculator and the stratifiers it names.
fn parse_directive(directive: &str) -> (Calculator, Vec<Stratifier>) {
    let mut terms = directive.split(':').map(str::trim);
    let calculator = terms.next().unwrap_or("");
    let stratifiers = terms
        .map(|name| {
            Stratifier::parse(name).unwrap_or_else(|| {
                fail(&format!(
                    "{TOOL}: the stratifier {name} is not ported, so a run that names it cannot \
                     be answered"
                ))
            })
        })
        .collect();
    let calculator = Calculator::parse(calculator).unwrap_or_else(|| {
        thrown(&format!(
            "java.lang.IllegalArgumentException: No enum constant \
             picard.sam.SamErrorMetric.ErrorType.{calculator}"
        ))
    });
    (calculator, stratifiers)
}

/// The index file `VCFFileReader(path, true)` demands next to a VCF.
fn has_index(vcf: &str) -> bool {
    let index = if vcf.ends_with(".gz") { ".tbi" } else { ".idx" };
    std::path::Path::new(&format!("{vcf}{index}")).exists()
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("R", "REFERENCE_SEQUENCE"),
        ("V", "VCF"),
        ("L", "INTERVALS"),
        ("MQ", "MIN_MAPPING_Q"),
        ("BQ", "MIN_BASE_Q"),
        ("PE", "PRIOR_Q"),
        ("MAX", "MAX_LOCI"),
        ("LH", "LONG_HOMOPOLYMER"),
        ("LBS", "LOCATION_BIN_SIZE"),
        ("P", "PROBABILITY"),
        ("EXT", "FILE_EXTENSION"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let reference = args.required("REFERENCE_SEQUENCE");
    let vcf = args.get("VCF").map(str::to_string);
    let interval_iterator = args.bool("INTERVAL_ITERATOR", false);
    let extension = args.get("FILE_EXTENSION").unwrap_or("").to_string();
    let min_mapping_q = args.int("MIN_MAPPING_Q", 20);
    let min_base_q = args.int("MIN_BASE_Q", 20);
    let prior_q = args.int("PRIOR_Q", 30);
    let max_loci = args.int("MAX_LOCI", 0);
    let long_homopolymer = args.int("LONG_HOMOPOLYMER", 6);
    let probability = args.double("PROBABILITY", 1.0);

    // customCommandLineValidation.
    let mut messages = Vec::new();
    if args.get("ERROR_VALUE").is_some() {
        messages.push(
            "ERROR_VALUE is a fake argument that is only there to show what are the different \
             Error aggregation options. Please use it within the ERROR_METRICS argument."
                .to_string(),
        );
    }
    if args.get("STRATIFIER_VALUE").is_some() {
        messages.push(
            "STRATIFIER_VALUE is a fake argument that is only there to show what are the \
             different Stratification options. Please use it within the STRATIFIER_VALUE \
             argument."
                .to_string(),
        );
    }
    if min_mapping_q < 0 {
        messages.push(format!(
            "MIN_MAPPING_Q must be non-negative. found value: {min_mapping_q}"
        ));
    }
    if min_base_q < 0 {
        messages.push(format!(
            "MIN_BASE_Q must be non-negative. found value: {min_base_q}"
        ));
    }
    if prior_q < 0 {
        messages.push(format!("PRIOR_Q must be 2 or more. found value: {prior_q}"));
    }
    if max_loci < 0 {
        messages.push(format!(
            "MAX_LOCI must be non-negative. found value: {max_loci}"
        ));
    }
    if long_homopolymer < 0 {
        messages.push(format!(
            "LONG_HOMOPOLYMER must be non-negative. found value: {long_homopolymer}"
        ));
    }
    if !(0.0..=1.0).contains(&probability) {
        messages.push(format!(
            "PROBABILITY must be between 0 and 1. found value: {probability:?}"
        ));
    }
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }

    // initializeAggregationState: the aggregators, before a file is read.
    let directives = args.collection("ERROR_METRICS", &DEFAULT_ERROR_METRICS);
    // Every directive is parsed and its suffix claimed in order, so a duplicate is found before
    // any stratifier this port does not bin is asked for; only a run that survives that and still
    // names one cannot be answered.
    for directive in &directives {
        if aggregation_suffix(directive).is_none() {
            parse_directive(directive);
        }
    }
    let written = match suffixes(&directives) {
        Ok(written) => written,
        Err(refusal) => thrown(&format!(
            "java.lang.IllegalArgumentException: {}",
            refusal.message()
        )),
    };
    let aggregators: Vec<(Calculator, Vec<Stratifier>)> =
        directives.iter().map(|d| parse_directive(d)).collect();

    // processData: the reader, the reference walker, then the variant source.
    let (header, records) = read_input(&input);
    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| fail(&format!("{e:?}")));
    if let Some(path) = &vcf {
        if !has_index(path) {
            if interval_iterator {
                // The reader that failed to open is closed in `finally`, and that NPE replaces
                // the Tribble exception it hides.
                thrown(
                    "java.lang.NullPointerException: Cannot invoke \
                     \"htsjdk.variant.vcf.VCFFileReader.close()\" because \
                     \"this.vcfFileReader\" is null",
                );
            }
            thrown(&format!(
                "htsjdk.tribble.TribbleException: An index is required, but none found with \
                 file ending .idx, for input source: file://{path}"
            ));
        }
    }

    // createSamLocusAndReferenceIterator.
    if sort_order(&header) != "coordinate" {
        thrown("picard.PicardException: Input BAM must be sorted by coordinate");
    }
    let mut masks: Vec<IntervalMask> = Vec::new();
    for path in args.all("INTERVALS") {
        masks.push(IntervalMask::new(&read_interval_list(&path)));
    }

    // One coordinate for every contig: `offsets[i] + position`, a position between contigs so a
    // position of zero never lands on the last base of the one before.
    let mut offsets: Vec<i32> = Vec::new();
    let mut joined: Vec<u8> = Vec::new();
    for contig in &contigs {
        offsets.push(joined.len() as i32);
        joined.extend_from_slice(&contig.bases);
        joined.push(b'N');
    }
    let contig_of = |position: i32| -> (i32, i32) {
        let index = offsets.partition_point(|&offset| offset < position) - 1;
        (index as i32, position - offsets[index])
    };

    let reads: Vec<Read> = records
        .iter()
        .filter(|record| record.reference_index >= 0)
        .map(|record| {
            let group = match record.tags.get(Tag::new(b"RG")) {
                Some(TagValue::Str(id)) => id.to_string(),
                _ => String::new(),
            };
            let mate_start = if record.mate_reference_index >= 0 {
                offsets[record.mate_reference_index as usize] + record.mate_alignment_start
            } else {
                0
            };
            Read {
                name: record.read_name.clone(),
                start: offsets[record.reference_index as usize] + record.alignment_start,
                bases: record.read_bases.clone(),
                qualities: record.base_qualities.clone(),
                flags: record.flags,
                mate_start,
                cigar: record
                    .cigar
                    .elements
                    .iter()
                    .map(|e| {
                        let operator = match e.op {
                            Op::Eq => '=',
                            op => char::from(op.to_char()),
                        };
                        (e.length as usize, operator)
                    })
                    .collect(),
                mapping_quality: record.mapping_quality,
                read_group: group,
                insert_size: record.inferred_insert_size,
            }
        })
        .collect();

    // The sites the VCF takes out: every unfiltered record's whole span.
    let mut known_sites: Vec<i32> = Vec::new();
    if let Some(path) = &vcf {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 7 || !(f[6] == "." || f[6] == "PASS") {
                continue;
            }
            let Some(index) = contigs.iter().position(|c| c.name == f[0]) else {
                continue;
            };
            let start: i32 = f[1].parse().unwrap_or(0);
            for position in start..start + f[3].len() as i32 {
                known_sites.push(offsets[index] + position);
            }
        }
    }

    let options = Options {
        min_mapping_q: min_mapping_q.clamp(0, 255) as u8,
        min_base_q: min_base_q.clamp(0, 255) as u8,
        prior_q: prior_q as i32,
        max_loci: max_loci as u64,
        known_sites,
        probability,
    };
    let mut loci = pileup(&reads, &options);
    if !masks.is_empty() {
        loci.retain(|locus| {
            let (index, position) = contig_of(locus.position);
            masks.iter().all(|mask| mask.get(index, position))
        });
    }
    if interval_iterator && vcf.is_none() && !loci.is_empty() {
        // `checkLocus` reads the reader the run never opened, at the first locus that survives
        // the draw; with `PROBABILITY` at one that is the first locus.
        let mut random = picard_analysis::theoretical_sensitivity::JavaRandom::new(42);
        if loci.iter().any(|_| random.next_double() <= probability) {
            thrown(
                "java.lang.NullPointerException: Cannot invoke \
                 \"htsjdk.variant.vcf.VCFFileReader.query(htsjdk.samtools.util.Locatable)\" \
                 because \"vcfFileReader\" is null",
            );
        }
    }
    let loci = processed_loci(loci, &options);

    for ((calculator, stratifiers), suffix) in aggregators.iter().zip(&written) {
        let table = collect_joint(
            &reads,
            &joined,
            1,
            &loci,
            *calculator,
            stratifiers,
            &options,
        );
        let path = format!("{output}.{suffix}{extension}");
        if let Err(e) = std::fs::write(&path, file_of(&table).write()) {
            fail(&format!("{e}"));
        }
    }
}
