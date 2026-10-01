//! `CollectMultipleMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.CollectMultipleMetrics` at tag 3.4.0. The tool is a dispatcher: every
//! `PROGRAM` becomes a `SinglePassSamProgram` built by `makeInstance` with its own output names
//! under the `OUTPUT` prefix and a handful of fixed arguments, and one `makeItSo` pass runs them
//! all. Each program here is the repository's own port of that collector, run as its binary with
//! exactly the arguments `makeInstance` gives it:
//!
//! * `PROGRAM` is a set appended to its default five and emptied by `null`; an empty set is
//!   refused by `customCommandLineValidation`;
//! * in program order, a program that needs a reference or a refFlat without one is refused
//!   before anything is read;
//! * the pass is `makeItSo`'s, shared: the sort check is made once with this tool's
//!   `ASSUME_SORTED`, and the reference walker sees every record up to `STOP_AFTER` and, unless
//!   every program ignores them, past the first unmapped one. Both refusals are made here, before
//!   any program runs, which is where the reference reports them when no program's `setup`
//!   fails first;
//! * the metrics files carry this tool's command line, not the program's: `setDefaultHeaders`.
//!
//! The charts are R's and are not drawn. `EXTRA_ARGUMENT` is not ported.

use std::path::PathBuf;
use std::process::Command;

use picard_analysis::metrics_cli::{
    check_coordinate_sorted, fail, read_input, refuse_validation, thrown, Args, ReferenceWalker,
};

const TOOL: &str = "CollectMultipleMetrics";

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

    /// The port that runs it, and the tool name that port writes into its metrics headers.
    fn port(self) -> (&'static str, &'static str) {
        match self {
            Program::CollectAlignmentSummaryMetrics => (
                "collect-alignment-summary-metrics",
                "CollectAlignmentSummaryMetrics",
            ),
            Program::CollectInsertSizeMetrics => {
                ("collect-insert-size-metrics", "CollectInsertSizeMetrics")
            }
            Program::QualityScoreDistribution => {
                ("quality-score-distribution", "QualityScoreDistribution")
            }
            Program::MeanQualityByCycle => ("mean-quality-by-cycle", "MeanQualityByCycle"),
            Program::CollectBaseDistributionByCycle => (
                "collect-base-distribution-by-cycle",
                "CollectBaseDistributionByCycle",
            ),
            Program::CollectGcBiasMetrics => ("collect-gc-bias-metrics", "CollectGcBiasMetrics"),
            Program::RnaSeqMetrics => ("collect-rna-seq-metrics", "CollectRnaSeqMetrics"),
            Program::CollectSequencingArtifactMetrics => (
                "collect-sequencing-artifact-metrics",
                "CollectSequencingArtifactMetrics",
            ),
            Program::CollectQualityYieldMetrics => (
                "collect-quality-yield-metrics",
                "CollectQualityYieldMetrics",
            ),
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
            args.push(format!("HISTOGRAM_FILE={ob}.read_length_histogram.pdf"));
            // This port appends to its ALL_READS default and has no `null`; the set it is given
            // always holds ALL_READS here, which is the same set.
            for level in s.levels {
                args.push(format!("METRIC_ACCUMULATION_LEVEL={level}"));
            }
        }
        Program::CollectInsertSizeMetrics => {
            output(&mut args, "OUTPUT", ".insert_size_metrics");
            args.push(format!("Histogram_FILE={ob}.insert_size_histogram.pdf"));
            args.extend([sorted, stop]);
        }
        Program::QualityScoreDistribution => {
            output(&mut args, "OUTPUT", ".quality_distribution_metrics");
            args.push(format!("CHART_OUTPUT={ob}.quality_distribution.pdf"));
            args.extend([sorted, stop]);
        }
        Program::MeanQualityByCycle => {
            output(&mut args, "OUTPUT", ".quality_by_cycle_metrics");
            args.push(format!("CHART_OUTPUT={ob}.quality_by_cycle.pdf"));
            args.extend([sorted, stop]);
        }
        Program::CollectBaseDistributionByCycle => {
            output(&mut args, "OUTPUT", ".base_distribution_by_cycle_metrics");
            args.push(format!("CHART_OUTPUT={ob}.base_distribution_by_cycle.pdf"));
            args.extend([sorted, stop]);
        }
        Program::CollectGcBiasMetrics => {
            output(&mut args, "OUTPUT", ".gc_bias.detail_metrics");
            output(&mut args, "SUMMARY_OUTPUT", ".gc_bias.summary_metrics");
            args.push(format!("CHART_OUTPUT={ob}.gc_bias.pdf"));
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
            args.push(format!("CHART_OUTPUT={ob}.rna_coverage.pdf"));
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
        }
    }
    (args, outputs)
}

fn main() {
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
    if !args.all("EXTRA_ARGUMENT").is_empty() {
        fail("EXTRA_ARGUMENT is not ported");
    }

    // customCommandLineValidation.
    if programs.is_empty() {
        refuse_validation(TOOL, &["No programs specified with PROGRAM".to_string()]);
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

    // makeItSo's shared refusals: the sort check, then the walker over the records the pass
    // would reach.
    let (header, records) = read_input(&input);
    check_coordinate_sorted(&input, &header, assume_sorted);
    if reference.is_some() {
        let any_use_no_ref_reads = programs.iter().any(|p| p.uses_no_ref_reads());
        let mut walker = ReferenceWalker::default();
        let mut count = 0i64;
        for record in &records {
            if record.reference_index != -1 {
                walker.get(record.reference_index);
            }
            count += 1;
            if stop_after > 0 && count >= stop_after {
                break;
            }
            if !any_use_no_ref_reads && record.reference_index == -1 {
                break;
            }
        }
    }

    let here: PathBuf = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(PathBuf::from))
        .unwrap_or_default();
    let shared = Shared {
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
    for program in &programs {
        let (binary, tool) = program.port();
        let (program_args, outputs) = instance(*program, &shared);
        let result = Command::new(here.join(binary))
            .args(&program_args)
            .output()
            .unwrap_or_else(|e| fail(&format!("cannot run {binary}: {e}")));
        if !result.status.success() {
            eprint!("{}", String::from_utf8_lossy(&result.stderr));
            std::process::exit(result.status.code().unwrap_or(1));
        }
        // `setDefaultHeaders(getDefaultHeaders())`: this tool's command line heads every file.
        let from = format!("\n# {tool} <command line>\n");
        let to = format!("\n# {TOOL} <command line>\n");
        for path in outputs {
            if let Ok(text) = std::fs::read_to_string(&path) {
                let _ = std::fs::write(&path, text.replacen(&from, &to, 1));
            }
        }
    }
}
