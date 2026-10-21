//! `CollectMultipleMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectMultipleMetrics` at tag 3.4.0. The tool is a dispatcher: every
//! `PROGRAM` becomes a `SinglePassSamProgram` built by `makeInstance` with its own output names
//! under the `OUTPUT` prefix and a handful of fixed arguments, and one `makeItSo` pass runs them
//! all. Here the five default programs are that pass's collectors, the library's own, built as
//! `makeInstance` builds them and each fed the records on a thread of its own; any other program,
//! and one given an `EXTRA_ARGUMENT`, is the repository's port of that collector, run as its binary
//! with exactly the arguments `makeInstance` gives it, beside the pass:
//!
//! * `PROGRAM` is a set appended to its default five and emptied by `null`; an empty set is
//!   refused by `customCommandLineValidation`;
//! * in program order, a program that needs a reference or a refFlat without one is refused
//!   before anything is read;
//! * the pass is `makeItSo`'s, shared: the sort check is made once with this tool's
//!   `ASSUME_SORTED`, and the reference walker sees every record up to `STOP_AFTER` and, unless
//!   every program ignores them, past the first unmapped one. Both refusals are made here, and a
//!   refusal leaves no program's files, which is what the reference leaves when no program's
//!   `setup` fails first;
//! * the metrics files carry this tool's command line, not the program's: `setDefaultHeaders`.
//!
//! The charts are R's and are not drawn. `EXTRA_ARGUMENT` hands a program arguments of its own.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::sync_channel;
use std::sync::Arc;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::alignment_summary::{LevelCollectors, Options as AlignmentSummaryOptions};
use picard_analysis::cycle::{CollectBaseDistributionByCycle, MeanQualityByCycle};
use picard_analysis::insert_size::{Distributors, Options as InsertSizeOptions};
use picard_analysis::metrics_cli::{
    check_coordinate_sorted, fail, open_input, refuse_validation, thrown, Args, ReferenceWalker,
};
use picard_analysis::quality_score_distribution::{
    Options as QualityOptions, QualityScoreDistribution,
};

const TOOL: &str = "CollectMultipleMetrics";

/// Records handed to the collectors at a time.
const BATCH: usize = 1024;

/// `CollectMultipleMetrics.Program`, in declaration order.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Program {
    CollectAlignmentSummaryMetrics,
    CollectInsertSizeMetrics,
    QualityScoreDistribution,
    MeanQualityByCycle,
    CollectBaseDistributionByCycle,
    CollectGcBiasMetrics,
    RnaSeqMetrics,
    CollectSequencingArtifactMetrics,
    CollectQualityYieldMetrics,
}

impl Program {
    fn parse(name: &str) -> Option<Program> {
        Some(match name {
            "CollectAlignmentSummaryMetrics" => Program::CollectAlignmentSummaryMetrics,
            "CollectInsertSizeMetrics" => Program::CollectInsertSizeMetrics,
            "QualityScoreDistribution" => Program::QualityScoreDistribution,
            "MeanQualityByCycle" => Program::MeanQualityByCycle,
            "CollectBaseDistributionByCycle" => Program::CollectBaseDistributionByCycle,
            "CollectGcBiasMetrics" => Program::CollectGcBiasMetrics,
            "RnaSeqMetrics" => Program::RnaSeqMetrics,
            "CollectSequencingArtifactMetrics" => Program::CollectSequencingArtifactMetrics,
            "CollectQualityYieldMetrics" => Program::CollectQualityYieldMetrics,
            _ => return None,
        })
    }

    fn name(self) -> &'static str {
        match self {
            Program::CollectAlignmentSummaryMetrics => "CollectAlignmentSummaryMetrics",
            Program::CollectInsertSizeMetrics => "CollectInsertSizeMetrics",
            Program::QualityScoreDistribution => "QualityScoreDistribution",
            Program::MeanQualityByCycle => "MeanQualityByCycle",
            Program::CollectBaseDistributionByCycle => "CollectBaseDistributionByCycle",
            Program::CollectGcBiasMetrics => "CollectGcBiasMetrics",
            Program::RnaSeqMetrics => "RnaSeqMetrics",
            Program::CollectSequencingArtifactMetrics => "CollectSequencingArtifactMetrics",
            Program::CollectQualityYieldMetrics => "CollectQualityYieldMetrics",
        }
    }

    /// The port that runs it, whose name is also the tool name it writes into its metrics headers.
    fn port(self) -> &'static str {
        match self {
            Program::RnaSeqMetrics => "CollectRnaSeqMetrics",
            other => other.name(),
        }
    }

    fn needs_reference_sequence(self) -> bool {
        matches!(
            self,
            Program::CollectGcBiasMetrics | Program::CollectSequencingArtifactMetrics
        )
    }

    fn needs_refflat_file(self) -> bool {
        self == Program::RnaSeqMetrics
    }

    /// `usesNoRefReads`: only these two stop at the unmapped reads on their own.
    fn uses_no_ref_reads(self) -> bool {
        !matches!(
            self,
            Program::CollectInsertSizeMetrics | Program::CollectSequencingArtifactMetrics
        )
    }
}

/// Everything `makeInstance` is handed.
struct Shared<'a> {
    outbase: &'a str,
    outext: &'a str,
    input: &'a str,
    reference: Option<&'a str>,
    levels: &'a [String],
    db_snp: Option<&'a str>,
    intervals: Option<&'a str>,
    refflat: Option<&'a str>,
    ignore_sequence: &'a [String],
    include_unpaired: bool,
    stop_after: i64,
    skip_plotting: bool,
}

/// `makeInstance`: the program's own arguments, and the metrics files it will write.
fn instance(program: Program, s: &Shared) -> (Vec<String>, Vec<String>) {
    let (ob, ext) = (s.outbase, s.outext);
    let mut args = vec![format!("INPUT={}", s.input)];
    let mut outputs = Vec::new();
    // `setReferenceSequence`, which every program but the quality yield is given.
    if let Some(r) = s.reference {
        if program != Program::CollectQualityYieldMetrics {
            args.push(format!("REFERENCE_SEQUENCE={r}"));
        }
    }
    // The shared pass made the sort check already; the programs are told the input is sorted.
    let sorted = "ASSUME_SORTED=true".to_string();
    let stop = format!("STOP_AFTER={}", s.stop_after);
    let levels = |args: &mut Vec<String>| {
        args.push("METRIC_ACCUMULATION_LEVEL=null".to_string());
        for level in s.levels {
            args.push(format!("METRIC_ACCUMULATION_LEVEL={level}"));
        }
    };
    let mut output = |args: &mut Vec<String>, name: &str, extension: &str| {
        let path = format!("{ob}{extension}{ext}");
        args.push(format!("{name}={path}"));
        outputs.push(path);
    };
    match program {
        Program::CollectAlignmentSummaryMetrics => {
            output(&mut args, "OUTPUT", ".alignment_summary_metrics");
            if !s.skip_plotting {
                args.push(format!("HISTOGRAM_FILE={ob}.read_length_histogram.pdf"));
            }
            // This port appends to its ALL_READS default and has no `null`; the set it is given
            // always holds ALL_READS here, which is the same set.
            for level in s.levels {
                args.push(format!("METRIC_ACCUMULATION_LEVEL={level}"));
            }
            // The shared pass stops every program at once, this one included.
            args.push(stop);
        }
        Program::CollectInsertSizeMetrics => {
            output(&mut args, "OUTPUT", ".insert_size_metrics");
            if !s.skip_plotting {
                args.push(format!("Histogram_FILE={ob}.insert_size_histogram.pdf"));
            }
            args.extend([sorted, stop]);
        }
        Program::QualityScoreDistribution => {
            output(&mut args, "OUTPUT", ".quality_distribution_metrics");
            if !s.skip_plotting {
                args.push(format!("CHART_OUTPUT={ob}.quality_distribution.pdf"));
            }
            args.extend([sorted, stop]);
        }
        Program::MeanQualityByCycle => {
            output(&mut args, "OUTPUT", ".quality_by_cycle_metrics");
            if !s.skip_plotting {
                args.push(format!("CHART_OUTPUT={ob}.quality_by_cycle.pdf"));
            }
            args.extend([sorted, stop]);
        }
        Program::CollectBaseDistributionByCycle => {
            output(&mut args, "OUTPUT", ".base_distribution_by_cycle_metrics");
            if !s.skip_plotting {
                args.push(format!("CHART_OUTPUT={ob}.base_distribution_by_cycle.pdf"));
            }
            args.extend([sorted, stop]);
        }
        Program::CollectGcBiasMetrics => {
            output(&mut args, "OUTPUT", ".gc_bias.detail_metrics");
            output(&mut args, "SUMMARY_OUTPUT", ".gc_bias.summary_metrics");
            if !s.skip_plotting {
                args.push(format!("CHART_OUTPUT={ob}.gc_bias.pdf"));
            }
            levels(&mut args);
            args.extend([
                "SCAN_WINDOW_SIZE=100".to_string(),
                "MINIMUM_GENOME_FRACTION=1.0E-5".to_string(),
                "IS_BISULFITE_SEQUENCED=false".to_string(),
                "ALSO_IGNORE_DUPLICATES=false".to_string(),
                sorted,
                stop,
            ]);
        }
        Program::RnaSeqMetrics => {
            output(&mut args, "OUTPUT", ".rna_metrics");
            if !s.skip_plotting {
                args.push(format!("CHART_OUTPUT={ob}.rna_coverage.pdf"));
            }
            levels(&mut args);
            if let Some(i) = s.intervals {
                args.push(format!("RIBOSOMAL_INTERVALS={i}"));
            }
            for sequence in s.ignore_sequence {
                args.push(format!("IGNORE_SEQUENCE={sequence}"));
            }
            if let Some(r) = s.refflat {
                args.push(format!("REF_FLAT={r}"));
            }
            args.extend([
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND".to_string(),
                sorted,
                stop,
            ]);
        }
        Program::CollectSequencingArtifactMetrics => {
            args.push(format!("OUTPUT={ob}"));
            args.push(format!("FILE_EXTENSION={ext}"));
            for e in [
                ".pre_adapter_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".bait_bias_summary_metrics",
                ".bait_bias_detail_metrics",
                ".error_summary_metrics",
            ] {
                outputs.push(format!("{ob}{e}{ext}"));
            }
            if let Some(d) = s.db_snp {
                args.push(format!("DB_SNP={d}"));
            }
            if let Some(i) = s.intervals {
                args.push(format!("INTERVALS={i}"));
            }
            args.push(format!("INCLUDE_UNPAIRED={}", s.include_unpaired));
            args.extend([sorted, stop]);
        }
        Program::CollectQualityYieldMetrics => {
            output(&mut args, "OUTPUT", ".quality_yield_metrics");
            args.push(stop);
        }
    }
    (args, outputs)
}

fn main() {
    picard_analysis::argument_check::check("CollectMultipleMetrics");
    picard_analysis::output::install();
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("R", "REFERENCE_SEQUENCE"),
        ("AS", "ASSUME_SORTED"),
        ("LEVEL", "METRIC_ACCUMULATION_LEVEL"),
        ("EXT", "FILE_EXTENSION"),
        ("UNPAIRED", "INCLUDE_UNPAIRED"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let reference = args.get("REFERENCE_SEQUENCE").map(str::to_string);
    let assume_sorted = args.bool("ASSUME_SORTED", true);
    let stop_after = args.int("STOP_AFTER", 0);
    let levels = {
        let mut out: Vec<String> = Vec::new();
        for level in args.collection("METRIC_ACCUMULATION_LEVEL", &["ALL_READS"]) {
            if !matches!(
                level.as_str(),
                "ALL_READS" | "SAMPLE" | "LIBRARY" | "READ_GROUP"
            ) {
                fail(&format!(
                    "Argument 'METRIC_ACCUMULATION_LEVEL' cannot be set to '{level}'"
                ));
            }
            if !out.contains(&level) {
                out.push(level);
            }
        }
        out
    };
    let file_extension = args.get("FILE_EXTENSION").map(str::to_string);
    let mut programs: Vec<Program> = Vec::new();
    for name in args.collection(
        "PROGRAM",
        &[
            "CollectAlignmentSummaryMetrics",
            "CollectBaseDistributionByCycle",
            "CollectInsertSizeMetrics",
            "MeanQualityByCycle",
            "QualityScoreDistribution",
        ],
    ) {
        let program = Program::parse(&name)
            .unwrap_or_else(|| fail(&format!("Argument 'PROGRAM' cannot be set to '{name}'")));
        if !programs.contains(&program) {
            programs.push(program);
        }
    }
    let intervals = args.get("INTERVALS").map(str::to_string);
    let db_snp = args.get("DB_SNP").map(str::to_string);
    let refflat = args.get("REF_FLAT").map(str::to_string);
    let ignore_sequence = args.all("IGNORE_SEQUENCE");
    let include_unpaired = args.bool("INCLUDE_UNPAIRED", false);
    // `processAdditionalArguments`: `<PROGRAM>::<ARGUMENT_AND_VALUE>`, the value optionally after
    // a space, grouped by program in the order given.
    let pattern = regex::Regex::new(r"^(?s)(.*)::(.+?)( +(.+))?$").expect("a valid pattern");
    let mut extra: Vec<(Program, Vec<String>)> = Vec::new();
    for value in args.all("EXTRA_ARGUMENT") {
        let Some(groups) = pattern.captures(&value) else {
            thrown(&format!(
                "org.broadinstitute.barclay.argparser.CommandLineException: couldn't understand EXTRA_ARGUMENT {value} it doesn't conform to the form '<PROGRAM>::<ARGUMENT_AND_VALUE>'."
            ));
        };
        let name = &groups[1];
        let Some(program) = Program::parse(name) else {
            thrown(&format!(
                "org.broadinstitute.barclay.argparser.CommandLineException: Couldn't find program with value {name}"
            ));
        };
        let mut words = vec![groups[2].to_string()];
        if let Some(optional) = groups.get(4).filter(|m| !m.as_str().is_empty()) {
            words.push(optional.as_str().to_string());
        }
        match extra.iter_mut().find(|(p, _)| *p == program) {
            Some((_, list)) => list.extend(words),
            None => extra.push((program, words)),
        }
    }

    // DO_NOT_CREATE_PLOTS, new in Picard 3.5.0: the programs are not handed a chart file.
    let skip_plotting = args.bool("DO_NOT_CREATE_PLOTS", false);

    // customCommandLineValidation.
    let mut messages = Vec::new();
    if programs.is_empty() {
        messages.push("No programs specified with PROGRAM".to_string());
    }
    if !skip_plotting && picard_analysis::metrics_cli::running_in_gatk_lite_docker() {
        messages.push("Histogram files cannot be written because that requires R, which is not available in the GATK Lite Docker image.".to_string());
    }
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }

    // doWork, before the pass.
    for program in &programs {
        if program.needs_reference_sequence() && reference.is_none() {
            thrown(&format!(
                "picard.PicardException: The {} program needs a REF Sequence, please set \
                 REFERENCE_SEQUENCE in the command line",
                program.name()
            ));
        }
        if program.needs_refflat_file() && refflat.is_none() {
            thrown(&format!(
                "picard.PicardException: The {} program needs a gene annotations file, please \
                 set REF_FLAT in the command line",
                program.name()
            ));
        }
    }

    // Every program's extra arguments are taken by its instance; any left over named a program
    // that was not asked for.
    let leftover: Vec<&(Program, Vec<String>)> = extra
        .iter()
        .filter(|(p, _)| !programs.contains(p))
        .collect();
    if !leftover.is_empty() {
        thrown(&format!(
            "org.broadinstitute.barclay.argparser.CommandLineException: EXTRA_ARGUMENT values were provided, but corresponding PROGRAM wasn't requested:{}",
            leftover
                .iter()
                .map(|(p, list)| format!("{}::[{}]", p.name(), list.join(", ")))
                .collect::<String>()
        ));
    }

    let here: PathBuf = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(PathBuf::from))
        .unwrap_or_default();
    let shared = Shared {
        skip_plotting,
        outbase: &output,
        outext: file_extension.as_deref().unwrap_or(""),
        input: &input,
        reference: reference.as_deref(),
        levels: &levels,
        db_snp: db_snp.as_deref(),
        intervals: intervals.as_deref(),
        refflat: refflat.as_deref(),
        ignore_sequence: &ignore_sequence,
        include_unpaired,
        stop_after,
    };
    let mut runs = Vec::new();
    for program in &programs {
        let binary = program.port();
        let (mut program_args, outputs) = instance(*program, &shared);
        // The program's own parser takes its extra arguments over what this tool set, except the
        // four `makeItSo` reads from this tool rather than from the program.
        if let Some((_, words)) = extra.iter().find(|(p, _)| p == program) {
            let mut given: Vec<String> = Vec::new();
            let mut words = words.iter().peekable();
            while let Some(word) = words.next() {
                let pair = match word.strip_prefix("--") {
                    Some(name) if !name.contains('=') => {
                        format!("{name}={}", words.next().cloned().unwrap_or_default())
                    }
                    Some(name) => name.to_string(),
                    None => word.clone(),
                };
                let name = pair.split('=').next().unwrap_or_default();
                if !matches!(
                    name,
                    "INPUT"
                        | "I"
                        | "REFERENCE_SEQUENCE"
                        | "R"
                        | "ASSUME_SORTED"
                        | "AS"
                        | "STOP_AFTER"
                ) {
                    given.push(pair);
                }
            }
            given.append(&mut program_args);
            program_args = given;
        }
        runs.push((*program, binary, program_args, outputs));
    }
    let existed: Vec<Vec<bool>> = runs
        .iter()
        .map(|(_, _, _, outputs)| outputs.iter().map(|p| Path::new(p).exists()).collect())
        .collect();

    // makeItSo's shared refusal: the sort check.
    let (header, records) = open_input(&input).into_parts();
    check_coordinate_sorted(&input, &header, assume_sorted);

    // The programs with a collector here take the records of this one pass; the others, and any
    // given arguments of its own, are their port's process over the input, run side by side.
    let mut slots: Vec<Slot> = Vec::new();
    for (program, binary, program_args, _) in &runs {
        let in_process = !extra.iter().any(|(p, _)| p == program);
        slots.push(
            match (in_process, Collector::new(*program, &levels, &header)) {
                (true, Some(collector)) => Slot::Here(Ok(collector)),
                _ => {
                    let mut child = Command::new(here.join(binary))
                        .args(program_args)
                        // The reference the shared pass would have read the input with, whichever
                        // program this is.
                        .env(
                            picard_analysis::cram_input::READER_REFERENCE_ENV,
                            reference.as_deref().unwrap_or(""),
                        )
                        .stdin(Stdio::null())
                        .stdout(Stdio::null())
                        .stderr(Stdio::piped())
                        .spawn()
                        .unwrap_or_else(|e| fail(&format!("cannot run {binary}: {e}")));
                    let mut stderr = child.stderr.take().expect("a piped stderr");
                    let drain = std::thread::spawn(move || {
                        let mut text = Vec::new();
                        let _ = stderr.read_to_end(&mut text);
                        text
                    });
                    Slot::Process(child, drain)
                }
            },
        );
    }
    // The pass: the reference walker sees every record up to `STOP_AFTER` and, unless every
    // program ignores them, past the first unmapped one.
    let contigs = match (&reference, slots.iter().any(Slot::wants_bases)) {
        (Some(path), true) => match read_fasta_file(path) {
            Ok(contigs) => contigs,
            Err(e) => {
                let message = format!("Error: {:?}\n", format!("{e:?}"));
                for slot in slots.iter_mut() {
                    if let Slot::Here(state @ Ok(Collector::AlignmentSummary(_))) = slot {
                        *state = Err(message.clone());
                    }
                }
                Vec::new()
            }
        },
        _ => Vec::new(),
    };
    let mut children: Vec<&mut Child> = Vec::new();
    let mut states: Vec<&mut Result<Collector, String>> = Vec::new();
    for slot in slots.iter_mut() {
        match slot {
            Slot::Here(state) => states.push(state),
            Slot::Process(child, _) => children.push(child),
        }
    }
    let any_use_no_ref_reads = programs.iter().any(|p| p.uses_no_ref_reads());
    // Each collector takes the records on a thread of its own, in batches, so that the pass lasts
    // as long as its slowest collector rather than as long as all of them.
    std::thread::scope(|scope| {
        let (header, contigs) = (&header, &contigs);
        let mut senders = Vec::new();
        for state in states {
            let (sender, batches) = sync_channel::<Arc<Vec<BamRecord>>>(4);
            senders.push(sender);
            scope.spawn(move || {
                for batch in batches {
                    for record in batch.iter() {
                        let bases = usize::try_from(record.reference_index)
                            .ok()
                            .and_then(|i| contigs.get(i))
                            .map(|c| c.bases.as_slice());
                        if let Some(Err(failure)) =
                            state.as_mut().ok().map(|c| c.accept(header, record, bases))
                        {
                            *state = Err(failure);
                        }
                    }
                }
            });
        }
        let send = |batch: Vec<BamRecord>| {
            let batch = Arc::new(batch);
            for sender in &senders {
                let _ = sender.send(Arc::clone(&batch));
            }
        };
        let mut batch = Vec::with_capacity(BATCH);
        let mut walker = ReferenceWalker::default();
        let mut count = 0i64;
        for record in records {
            if reference.is_some() && record.reference_index != -1 {
                // A refusal of the pass leaves no program's files, as if none had run.
                if walker.rewinds(record.reference_index) {
                    for child in children.iter_mut() {
                        let _ = child.kill();
                        let _ = child.wait();
                    }
                    for ((_, _, _, outputs), existed) in runs.iter().zip(&existed) {
                        for (path, existed) in outputs.iter().zip(existed) {
                            if !existed {
                                let _ = std::fs::remove_file(path);
                            }
                        }
                    }
                }
                walker.get(record.reference_index);
            } else if reference.is_none() && senders.is_empty() {
                break;
            }
            count += 1;
            let last = (stop_after > 0 && count >= stop_after)
                || (!any_use_no_ref_reads && record.reference_index == -1);
            batch.push(record);
            if batch.len() == BATCH {
                send(std::mem::replace(&mut batch, Vec::with_capacity(BATCH)));
            }
            if last {
                break;
            }
        }
        send(batch);
        drop(senders);
    });

    // The results are taken in program order. A file a later program created is removed when an
    // earlier one fails, as if it had never run.
    // `setDefaultHeaders(getDefaultHeaders())`: this tool's command line and start head every
    // file, in place of the ones the program would have written for itself.
    let headers = picard_analysis::command_line::default_headers(TOOL);
    for (k, ((_, binary, _, outputs), slot)) in runs.iter().zip(slots).enumerate() {
        let failure = match slot {
            Slot::Here(state) => state
                .and_then(|collector| collector.write(&outputs[0], &headers))
                .err()
                .map(|text| (text, 1)),
            Slot::Process(child, drain) => {
                let status = child
                    .wait_with_output()
                    .unwrap_or_else(|e| fail(&format!("cannot run {binary}: {e}")))
                    .status;
                let stderr = drain.join().expect("stderr thread");
                if status.success() {
                    rewrite_headers(binary, outputs, &headers);
                    None
                } else {
                    Some((
                        String::from_utf8_lossy(&stderr).into_owned(),
                        status.code().unwrap_or(1),
                    ))
                }
            }
        };
        if let Some((text, code)) = failure {
            for ((_, _, _, later), existed) in runs.iter().zip(&existed).skip(k + 1) {
                for (path, existed) in later.iter().zip(existed) {
                    if !existed {
                        let _ = std::fs::remove_file(path);
                    }
                }
            }
            eprint!("{text}");
            std::process::exit(code);
        }
    }
}

/// One program of the pass: a collector here, which holds the failure it ends in once it has one,
/// or its port's process.
// One per program: the size is never paid more than nine times.
#[allow(clippy::large_enum_variant)]
enum Slot {
    Here(Result<Collector, String>),
    Process(Child, std::thread::JoinHandle<Vec<u8>>),
}

impl Slot {
    fn wants_bases(&self) -> bool {
        matches!(self, Slot::Here(Ok(Collector::AlignmentSummary(_))))
    }
}

/// The collectors of the programs this tool runs itself, built as their port builds them from
/// the arguments `makeInstance` hands it.
enum Collector {
    AlignmentSummary(LevelCollectors),
    /// `true` once it has taken an unmapped record: `usesNoRefReads()` is false.
    InsertSize(Distributors, bool),
    QualityScoreDistribution(Box<QualityScoreDistribution>),
    MeanQualityByCycle(MeanQualityByCycle),
    BaseDistributionByCycle(CollectBaseDistributionByCycle),
}

impl Collector {
    fn new(program: Program, levels: &[String], header: &SamHeader) -> Option<Collector> {
        Some(match program {
            Program::CollectAlignmentSummaryMetrics => Collector::AlignmentSummary(
                LevelCollectors::new(&AlignmentSummaryOptions::default(), levels, header),
            ),
            // `makeInstance` does not hand it the levels: it runs at its default, all reads.
            Program::CollectInsertSizeMetrics => Collector::InsertSize(
                Distributors::new(
                    &["ALL_READS".to_string()],
                    header,
                    InsertSizeOptions::default(),
                ),
                false,
            ),
            Program::QualityScoreDistribution => Collector::QualityScoreDistribution(Box::new(
                QualityScoreDistribution::new(QualityOptions::default()),
            )),
            Program::MeanQualityByCycle => Collector::MeanQualityByCycle(Default::default()),
            Program::CollectBaseDistributionByCycle => {
                Collector::BaseDistributionByCycle(Default::default())
            }
            _ => return None,
        })
    }

    /// Takes one record, or says what the program fails with.
    fn accept(
        &mut self,
        header: &SamHeader,
        record: &BamRecord,
        bases: Option<&[u8]>,
    ) -> Result<(), String> {
        match self {
            Collector::AlignmentSummary(c) => c.accept(record, bases),
            Collector::InsertSize(_, true) => {}
            Collector::InsertSize(c, done) => {
                c.accept(header, record)
                    .map_err(|thrown| format!("Exception in thread \"main\" {thrown}\n"))?;
                *done = record.reference_index == -1;
            }
            Collector::QualityScoreDistribution(c) => c.accept(record),
            Collector::MeanQualityByCycle(c) => c.accept(record),
            Collector::BaseDistributionByCycle(c) => c.accept(record),
        }
        Ok(())
    }

    /// Writes the program's metrics file under `headers`, or says what the program fails with.
    fn write(self, output: &str, headers: &[String]) -> Result<(), String> {
        let io = |e: std::io::Error| format!("Error: {e:?}\n");
        let file = match self {
            Collector::AlignmentSummary(c) => c
                .metrics_file(headers)
                .map_err(|e| format!("Error: {e:?}\n"))?,
            Collector::InsertSize(c, _) => {
                return match c.write(output, headers, None, None).map_err(io)? {
                    Some(failure) => Err(failure),
                    None => Ok(()),
                }
            }
            Collector::QualityScoreDistribution(c) => c.metrics_file(headers),
            Collector::MeanQualityByCycle(c) => {
                let mut file = MetricsFile::new();
                for header in headers {
                    file.add_header(header);
                }
                file.histograms = c.finish();
                file
            }
            Collector::BaseDistributionByCycle(c) => {
                let mut file = MetricsFile::new();
                for header in headers {
                    file.add_header(header);
                }
                for metric in c.finish() {
                    file.add_metric(&metric);
                }
                file
            }
        };
        std::fs::write(output, file.write()).map_err(io)
    }
}

/// `setDefaultHeaders` on a file a program's process wrote: its own command line and start give
/// way to this tool's.
fn rewrite_headers(binary: &str, outputs: &[String], headers: &[String; 2]) {
    let [command_line, started_on] = headers;
    // `port` is the tool name the program writes into its own headers.
    let own = format!("# {binary} ");
    for path in outputs {
        if let Ok(text) = std::fs::read_to_string(path) {
            let mut replaced_command = false;
            let mut replaced_start = false;
            let lines: Vec<String> = text
                .split('\n')
                .map(|line| {
                    if !replaced_command && line.starts_with(&own) {
                        replaced_command = true;
                        format!("# {command_line}")
                    } else if !replaced_start && line.starts_with("# Started on: ") {
                        replaced_start = true;
                        format!("# {started_on}")
                    } else {
                        line.to_string()
                    }
                })
                .collect();
            let _ = std::fs::write(path, lines.join("\n"));
        }
    }
}
