//! `CollectMultipleMetrics`: which programs one pass runs, and which files they land on.
//!
//! The tool is a dispatcher rather than an arithmetic. It builds one instance per `PROGRAM`, hands
//! every one of them the same records, and lets each write its own files, which is why the numbers
//! it produces are the standalone tools' and are measured with them rather than here.
//!
//! What is here is the dispatch: the eight programs the enum declares, the five the default set
//! holds, the files each one lands on, what each needs before a record is read, and how one
//! `EXTRA_ARGUMENT` reaches one program and no other.
//!
//! Ported from `picard.analysis.CollectMultipleMetrics`.

/// One of the programs the enum declares.
///
/// The order is the enum's own, which is not the default set's: `PROGRAM` is a `LinkedHashSet`
/// built from five of these in a different order again.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
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

/// Every program, in the enum's own declaration order.
pub const PROGRAMS: [Program; 9] = [
    Program::CollectAlignmentSummaryMetrics,
    Program::CollectInsertSizeMetrics,
    Program::QualityScoreDistribution,
    Program::MeanQualityByCycle,
    Program::CollectBaseDistributionByCycle,
    Program::CollectGcBiasMetrics,
    Program::RnaSeqMetrics,
    Program::CollectSequencingArtifactMetrics,
    Program::CollectQualityYieldMetrics,
];

/// The default `PROGRAM` set, in the order the `LinkedHashSet` is built.
///
/// Five of the nine, and the three that need more than the reads are not among them: a plain run
/// needs no reference at all.
pub const DEFAULT_PROGRAMS: [Program; 5] = [
    Program::CollectAlignmentSummaryMetrics,
    Program::CollectBaseDistributionByCycle,
    Program::CollectInsertSizeMetrics,
    Program::MeanQualityByCycle,
    Program::QualityScoreDistribution,
];

impl Program {
    /// The name a command line writes, which is the enum constant's own.
    pub fn name(self) -> &'static str {
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

    /// The extensions `getExtensions()` reports, metrics first and chart second.
    ///
    /// They are not all pairs: two of the programs write one file and one writes three, and the
    /// list is what the tool's own help prints under `PROGRAM`.
    pub fn extensions(self) -> &'static [&'static str] {
        match self {
            Program::CollectAlignmentSummaryMetrics => {
                &[".alignment_summary_metrics", ".read_length_histogram.pdf"]
            }
            Program::CollectInsertSizeMetrics => {
                &[".insert_size_metrics", ".insert_size_histogram.pdf"]
            }
            Program::QualityScoreDistribution => {
                &[".quality_distribution_metrics", ".quality_distribution.pdf"]
            }
            Program::MeanQualityByCycle => &[".quality_by_cycle_metrics", ".quality_by_cycle.pdf"],
            Program::CollectBaseDistributionByCycle => &[
                ".base_distribution_by_cycle_metrics",
                ".base_distribution_by_cycle.pdf",
            ],
            Program::CollectGcBiasMetrics => &[
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            ],
            Program::RnaSeqMetrics => &[".rna_metrics", ".rna_coverage.pdf"],
            Program::CollectSequencingArtifactMetrics => &[
                ".bait_bias_detail_metrics",
                ".bait_bias_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".pre_adapter_summary_metrics",
                ".error_summary_metrics",
            ],
            Program::CollectQualityYieldMetrics => &[".quality_yield_metrics"],
        }
    }

    /// Whether the run is refused without a `REFERENCE_SEQUENCE`, before any record is read.
    pub fn needs_reference_sequence(self) -> bool {
        matches!(
            self,
            Program::CollectGcBiasMetrics | Program::CollectSequencingArtifactMetrics
        )
    }

    /// Whether the run is refused without a `REF_FLAT`, which is a second and separate check.
    pub fn needs_refflat_file(self) -> bool {
        matches!(self, Program::RnaSeqMetrics)
    }

    /// Whether a `METRIC_ACCUMULATION_LEVEL` other than the default reaches this program.
    ///
    /// The three that do not are still RUN with a level that was overridden; the level is ignored
    /// and a warning is logged, which is why this is not a refusal.
    pub fn supports_metric_accumulation_level(self) -> bool {
        matches!(
            self,
            Program::CollectAlignmentSummaryMetrics
                | Program::CollectInsertSizeMetrics
                | Program::CollectGcBiasMetrics
                | Program::RnaSeqMetrics
        )
    }

    /// The name a `PROGRAM` value resolves to, or nothing.
    pub fn parse(name: &str) -> Option<Program> {
        PROGRAMS.into_iter().find(|program| program.name() == name)
    }
}

/// `FILE_EXTENSION` lands on the metrics files and NOT on the charts.
///
/// `makeInstance` appends `outext` to the metrics file's name and hands the chart the base name
/// and its own extension, so `EXT=.txt` renames one half of a program's output.
pub fn file_names(basename: &str, program: Program, file_extension: Option<&str>) -> Vec<String> {
    program
        .extensions()
        .iter()
        .map(|extension| {
            if extension.ends_with(".pdf") {
                format!("{basename}{extension}")
            } else {
                format!("{basename}{extension}{}", file_extension.unwrap_or(""))
            }
        })
        .collect()
}

/// What `customCommandLineValidation` and the loop in `doWork` refuse, in that order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// `PROGRAM` was emptied and nothing put back.
    NoPrograms,
    /// A program that needs the reference was named without one.
    NeedsReferenceSequence(Program),
    /// A program that needs the annotations was named without them.
    NeedsRefflatFile(Program),
    /// An `EXTRA_ARGUMENT` that does not name a program and an argument.
    ExtraArgumentMalformed(String),
    /// An `EXTRA_ARGUMENT` whose program name resolves to nothing.
    ExtraArgumentUnknownProgram(String),
    /// An `EXTRA_ARGUMENT` for a program that is not in `PROGRAM`.
    ExtraArgumentNotRequested(Program),
}

/// `No programs specified with PROGRAM`, which is a validation message and not an exception.
pub const NO_PROGRAMS_MESSAGE: &str = "No programs specified with PROGRAM";

impl Refusal {
    /// The message the reference writes, which for the first is printed under the usage and for
    /// the rest is an exception's own.
    pub fn message(&self) -> String {
        match self {
            Refusal::NoPrograms => NO_PROGRAMS_MESSAGE.to_string(),
            Refusal::NeedsReferenceSequence(program) => format!(
                "The {} program needs a REF Sequence, please set REFERENCE_SEQUENCE in the \
                 command line",
                program.name()
            ),
            Refusal::NeedsRefflatFile(program) => format!(
                "The {} program needs a gene annotations file, please set REF_FLAT in the \
                 command line",
                program.name()
            ),
            Refusal::ExtraArgumentMalformed(value) => format!(
                "couldn't understand EXTRA_ARGUMENT {value} it doesn't conform to the form \
                 '<PROGRAM>::<ARGUMENT_AND_VALUE>'."
            ),
            Refusal::ExtraArgumentUnknownProgram(name) => {
                format!("Couldn't find program with value {name}")
            }
            Refusal::ExtraArgumentNotRequested(program) => format!(
                "EXTRA_ARGUMENT values were provided, but corresponding PROGRAM wasn't requested:{}",
                program.name()
            ),
        }
    }
}

/// One `EXTRA_ARGUMENT`, taken apart the way the pattern does.
///
/// The pattern is `(?<program>.*)::(?<argumentAndValue>.+?)( +(?<optionalValue>.+))?`, which is
/// reluctant in the middle: `A::--B C` is the argument `--B` and the value `C`, two entries, while
/// `A::B=C` is one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtraArgument {
    pub program: Program,
    pub values: Vec<String>,
}

/// `processAdditionalArguments`, one value at a time.
pub fn extra_argument(value: &str) -> Result<ExtraArgument, Refusal> {
    let Some((name, rest)) = value.split_once("::") else {
        return Err(Refusal::ExtraArgumentMalformed(value.to_string()));
    };
    if rest.is_empty() {
        return Err(Refusal::ExtraArgumentMalformed(value.to_string()));
    }
    let Some(program) = Program::parse(name) else {
        return Err(Refusal::ExtraArgumentUnknownProgram(name.to_string()));
    };
    // The reluctant group takes as little as it can, so the split is at the FIRST run of spaces
    // and the rest of the line is one value however many spaces are in it.
    let mut values = Vec::new();
    match rest.find(' ') {
        None => values.push(rest.to_string()),
        Some(position) => {
            values.push(rest[..position].to_string());
            let tail = rest[position..].trim_start_matches(' ');
            if !tail.is_empty() {
                values.push(tail.to_string());
            }
        }
    }
    Ok(ExtraArgument { program, values })
}

/// `doWork`'s loop, as far as it decides anything: which programs run, and what refuses first.
///
/// The checks are in the reference's own order, which is why a program that needs a reference is
/// refused before an `EXTRA_ARGUMENT` that names a program nobody asked for: the loop runs first,
/// and the leftover check runs after it.
pub fn plan(
    programs: &[Program],
    reference_sequence: bool,
    refflat: bool,
    extra_arguments: &[&str],
) -> Result<Vec<(Program, Vec<String>)>, Refusal> {
    if programs.is_empty() {
        return Err(Refusal::NoPrograms);
    }
    let mut extras: Vec<ExtraArgument> = Vec::new();
    for value in extra_arguments {
        extras.push(extra_argument(value)?);
    }
    let mut planned = Vec::new();
    for program in programs {
        if program.needs_reference_sequence() && !reference_sequence {
            return Err(Refusal::NeedsReferenceSequence(*program));
        }
        if program.needs_refflat_file() && !refflat {
            return Err(Refusal::NeedsRefflatFile(*program));
        }
        let mut values = Vec::new();
        for extra in &extras {
            if extra.program == *program {
                values.extend(extra.values.iter().cloned());
            }
        }
        planned.push((*program, values));
    }
    for extra in &extras {
        if !programs.contains(&extra.program) {
            return Err(Refusal::ExtraArgumentNotRequested(extra.program));
        }
    }
    Ok(planned)
}

/// Whether the bytes of the charts are a claim this port makes, which they are not.
///
/// Every chart is drawn by R, and the same fixture drawn twice gives files that differ byte for
/// byte. The golden records their names and that they are not empty, and says nothing else about
/// them; a port that produced identical bytes would be reproducing a coincidence.
pub const CHART_BYTES_ARE_REPRODUCIBLE: bool = false;
