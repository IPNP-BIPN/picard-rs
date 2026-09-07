//! `GenotypeConcordance`: which pair of genotype states a site is filed under, and what that pair
//! contributes to the contingency table.
//!
//! Reading the VCFs is not ported. What is ported is the state a genotype resolves to, the scheme
//! that maps a pair of states to a contingency, and the three files the basename stands for.
//!
//! Ported from `picard.vcf.GenotypeConcordance`, `picard.vcf.GenotypeConcordanceStates`,
//! `picard.vcf.GA4GHScheme` and `picard.vcf.GA4GHSchemeWithMissingAsHomRef` in Picard 3.4.0.

/// The three files `--OUTPUT` stands for.
pub const SUMMARY_METRICS_FILE_EXTENSION: &str = ".genotype_concordance_summary_metrics";
pub const DETAILED_METRICS_FILE_EXTENSION: &str = ".genotype_concordance_detail_metrics";
pub const CONTINGENCY_METRICS_FILE_EXTENSION: &str = ".genotype_concordance_contingency_metrics";

/// The names the three files take from one basename, in the order the tool assigns them.
pub fn file_names(basename: &str) -> [String; 3] {
    [
        format!("{basename}{SUMMARY_METRICS_FILE_EXTENSION}"),
        format!("{basename}{DETAILED_METRICS_FILE_EXTENSION}"),
        format!("{basename}{CONTINGENCY_METRICS_FILE_EXTENSION}"),
    ]
}

/// What one truth genotype resolves to.
///
/// A missing site, a no-call, a filter, a low quality and a low depth are STATES here and not
/// exclusions: the site is still counted, under another name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TruthState {
    Missing,
    HomRef,
    HetRefVar1,
    HetVar1Var2,
    HomVar1,
    NoCall,
    LowGq,
    LowDp,
    VcFiltered,
    GtFiltered,
    IsMixed,
}

/// What one call genotype resolves to.
///
/// The call side has six states the truth side has no name for, because a call may carry an allele
/// the truth never mentioned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallState {
    Missing,
    HomRef,
    HetRefVar1,
    HetRefVar2,
    HetRefVar3,
    HetVar1Var2,
    HetVar1Var3,
    HetVar3Var4,
    HomVar1,
    HomVar2,
    HomVar3,
    NoCall,
    VcFiltered,
    GtFiltered,
    LowGq,
    LowDp,
    IsMixed,
}

impl TruthState {
    /// The name the metrics file writes, which is the enum's own and not Rust's spelling.
    pub fn name(self) -> &'static str {
        TRUTH_NAMES[TRUTH_ORDER
            .iter()
            .position(|state| *state == self)
            .expect("a state")]
    }
}

impl CallState {
    /// The name the metrics file writes.
    pub fn name(self) -> &'static str {
        CALL_NAMES[CALL_ORDER
            .iter()
            .position(|state| *state == self)
            .expect("a state")]
    }
}

/// The four counters the contingency table holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContingencyState {
    Tp,
    Fp,
    Tn,
    Fn,
}

impl ContingencyState {
    /// The name the detail file writes, which is what `CONTINGENCY_VALUES` holds.
    pub fn name(self) -> &'static str {
        match self {
            ContingencyState::Tp => "TP",
            ContingencyState::Fp => "FP",
            ContingencyState::Tn => "TN",
            ContingencyState::Fn => "FN",
        }
    }
}

/// One cell of a scheme: what the pair contributes, or that the pair cannot happen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cell {
    /// The contingency values, which may be none at all.
    Values(&'static [ContingencyState]),
    /// `NA`: a pair the reference says its own code cannot reach.
    Unreachable,
}

/// The GA4GH scheme: what each (call, truth) pair contributes.
pub const GA4GH: [(CallState, [Cell; 11]); 17] = [
    (
        CallState::Missing,
        [
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HomRef,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HetRefVar1,
        [
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HetRefVar2,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[
                ContingencyState::Fp,
                ContingencyState::Tn,
                ContingencyState::Fn,
            ]),
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HetRefVar3,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HetVar1Var2,
        [
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp]),
            Cell::Values(&[
                ContingencyState::Tp,
                ContingencyState::Fp,
                ContingencyState::Fn,
            ]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HetVar1Var3,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[
                ContingencyState::Tp,
                ContingencyState::Fp,
                ContingencyState::Fn,
            ]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HetVar3Var4,
        [
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HomVar1,
        [
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Tp]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HomVar2,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HomVar3,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::NoCall,
        [
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::VcFiltered,
        [
            Cell::Values(&[]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::GtFiltered,
        [
            Cell::Values(&[]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::LowGq,
        [
            Cell::Values(&[]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::LowDp,
        [
            Cell::Values(&[]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::IsMixed,
        [
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
];

/// The same scheme with a missing site read as homozygous reference.
pub const GA4GH_MISSING_AS_HOM_REF: [(CallState, [Cell; 11]); 17] = [
    (
        CallState::Missing,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HomRef,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HetRefVar1,
        [
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HetRefVar2,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[
                ContingencyState::Fp,
                ContingencyState::Tn,
                ContingencyState::Fn,
            ]),
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HetRefVar3,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HetVar1Var2,
        [
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp]),
            Cell::Values(&[
                ContingencyState::Tp,
                ContingencyState::Fp,
                ContingencyState::Fn,
            ]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HetVar1Var3,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[
                ContingencyState::Tp,
                ContingencyState::Fp,
                ContingencyState::Fn,
            ]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HetVar3Var4,
        [
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HomVar1,
        [
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fp]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Tp]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::HomVar2,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Tp, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::HomVar3,
        [
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Values(&[ContingencyState::Fp, ContingencyState::Fn]),
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
            Cell::Unreachable,
        ],
    ),
    (
        CallState::NoCall,
        [
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::VcFiltered,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::GtFiltered,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::LowGq,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::LowDp,
        [
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn]),
            Cell::Values(&[ContingencyState::Tn, ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[ContingencyState::Fn]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
    (
        CallState::IsMixed,
        [
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
            Cell::Values(&[]),
        ],
    ),
];

/// The truth states, in the order the scheme's columns are written.
pub const TRUTH_ORDER: [TruthState; 11] = [
    TruthState::Missing,
    TruthState::HomRef,
    TruthState::HetRefVar1,
    TruthState::HetVar1Var2,
    TruthState::HomVar1,
    TruthState::NoCall,
    TruthState::LowGq,
    TruthState::LowDp,
    TruthState::VcFiltered,
    TruthState::GtFiltered,
    TruthState::IsMixed,
];

/// The names the truth column is written under, in the order of [`TRUTH_ORDER`].
pub const TRUTH_NAMES: [&str; 11] = [
    "MISSING",
    "HOM_REF",
    "HET_REF_VAR1",
    "HET_VAR1_VAR2",
    "HOM_VAR1",
    "NO_CALL",
    "LOW_GQ",
    "LOW_DP",
    "VC_FILTERED",
    "GT_FILTERED",
    "IS_MIXED",
];

/// The call states, in the order the scheme's rows are written.
pub const CALL_ORDER: [CallState; 17] = [
    CallState::Missing,
    CallState::HomRef,
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
    CallState::HomVar1,
    CallState::HomVar2,
    CallState::HomVar3,
    CallState::NoCall,
    CallState::VcFiltered,
    CallState::GtFiltered,
    CallState::LowGq,
    CallState::LowDp,
    CallState::IsMixed,
];

/// The names the call column is written under, in the order of [`CALL_ORDER`].
pub const CALL_NAMES: [&str; 17] = [
    "MISSING",
    "HOM_REF",
    "HET_REF_VAR1",
    "HET_REF_VAR2",
    "HET_REF_VAR3",
    "HET_VAR1_VAR2",
    "HET_VAR1_VAR3",
    "HET_VAR3_VAR4",
    "HOM_VAR1",
    "HOM_VAR2",
    "HOM_VAR3",
    "NO_CALL",
    "VC_FILTERED",
    "GT_FILTERED",
    "LOW_GQ",
    "LOW_DP",
    "IS_MIXED",
];

/// The state a name in the metrics file stands for.
pub fn truth_state(name: &str) -> Option<TruthState> {
    TRUTH_NAMES
        .iter()
        .position(|written| *written == name)
        .map(|index| TRUTH_ORDER[index])
}

/// The call state a name in the metrics file stands for.
pub fn call_state(name: &str) -> Option<CallState> {
    CALL_NAMES
        .iter()
        .position(|written| *written == name)
        .map(|index| CALL_ORDER[index])
}

/// What one pair contributes under a scheme.
pub fn contingency(
    scheme: &[(CallState, [Cell; 11])],
    call: CallState,
    truth: TruthState,
) -> Option<Cell> {
    let column = TRUTH_ORDER.iter().position(|state| *state == truth)?;
    scheme
        .iter()
        .find(|(state, _)| *state == call)
        .map(|(_, row)| row[column])
}

/// The `CONTINGENCY_VALUES` column, as the detail file writes it.
///
/// An empty cell writes `EMPTY`, not an empty column. The reference reaches that string twice
/// over: `getContingencyStateString` answers `EMPTY` for an array of length zero, and the scheme's
/// `EMPTY` is itself a one-element array holding the `EMPTY` state, whose join is the same word.
pub fn contingency_values(cell: Cell) -> Option<String> {
    match cell {
        Cell::Unreachable => None,
        Cell::Values([]) => Some("EMPTY".to_string()),
        Cell::Values(values) => Some(
            values
                .iter()
                .map(|value| value.name())
                .collect::<Vec<_>>()
                .join(","),
        ),
    }
}

/// The same column for a cell the scheme says is unreachable, which the detail file writes as
/// `NA` when `OUTPUT_ALL_ROWS` asks for every row rather than the ones that happened.
pub fn contingency_string(cell: Cell) -> String {
    contingency_values(cell).unwrap_or_else(|| "NA".to_string())
}

/// The truth states in the order the ENUM declares them, which is the order the detail file's rows
/// take. It is the same order the scheme's columns are written in.
pub const TRUTH_DECLARATION_ORDER: [TruthState; 11] = TRUTH_ORDER;

/// The call states in the order the ENUM declares them, which is NOT the order the scheme's rows
/// are written in: the scheme groups `HET_REF_VAR2` and `HET_REF_VAR3` next to `HET_REF_VAR1`,
/// while the enum puts the five comparable states first and the six incomparable ones after them.
/// The detail file's rows follow the ENUM, so the two orders have to be kept apart.
pub const CALL_DECLARATION_ORDER: [CallState; 17] = [
    CallState::Missing,
    CallState::HomRef,
    CallState::HetRefVar1,
    CallState::HetVar1Var2,
    CallState::HomVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
    CallState::HomVar2,
    CallState::HomVar3,
    CallState::NoCall,
    CallState::LowGq,
    CallState::LowDp,
    CallState::VcFiltered,
    CallState::GtFiltered,
    CallState::IsMixed,
];

/// `GenotypeConcordanceStateCodes.ordinal()`: the code two states share when they say the same
/// thing about the site.
///
/// The diagonal of the concordance calculation is defined on these codes and not on the states, so
/// a truth `HOM_VAR1` and a call `HOM_VAR1` agree while a call `HOM_VAR2` -- which carries the
/// INCOMPARABLE code, as all six call-only states do -- agrees with nothing, itself included.
pub fn truth_code(state: TruthState) -> u8 {
    match state {
        TruthState::Missing => 0,
        TruthState::HomRef => 1,
        TruthState::HetRefVar1 => 2,
        TruthState::HetVar1Var2 => 3,
        TruthState::HomVar1 => 4,
        TruthState::NoCall => 5,
        TruthState::LowGq => 6,
        TruthState::LowDp => 7,
        TruthState::VcFiltered => 8,
        TruthState::GtFiltered => 9,
        TruthState::IsMixed => 10,
    }
}

/// The same codes on the call side. `INCOMPARABLE_CODE` is eleven, and six states carry it.
pub fn call_code(state: CallState) -> u8 {
    match state {
        CallState::Missing => 0,
        CallState::HomRef => 1,
        CallState::HetRefVar1 => 2,
        CallState::HetVar1Var2 => 3,
        CallState::HomVar1 => 4,
        CallState::NoCall => 5,
        CallState::LowGq => 6,
        CallState::LowDp => 7,
        CallState::VcFiltered => 8,
        CallState::GtFiltered => 9,
        CallState::IsMixed => 10,
        CallState::HetRefVar2
        | CallState::HetRefVar3
        | CallState::HetVar1Var3
        | CallState::HetVar3Var4
        | CallState::HomVar2
        | CallState::HomVar3 => 11,
    }
}

/// `GenotypeConcordanceCounts.HET_TRUTH_STATES` and its four companions, which are what the
/// summary metrics' six ratios are taken over.
pub const HET_TRUTH_STATES: [TruthState; 2] = [TruthState::HetRefVar1, TruthState::HetVar1Var2];
pub const HOM_VAR_TRUTH_STATES: [TruthState; 1] = [TruthState::HomVar1];
/// The VAR set carries `HOM_REF` and `MISSING` as well, because specificity needs the sites where
/// there was nothing to find.
pub const VAR_TRUTH_STATES: [TruthState; 5] = [
    TruthState::HetRefVar1,
    TruthState::HetVar1Var2,
    TruthState::HomVar1,
    TruthState::HomRef,
    TruthState::Missing,
];
pub const HET_CALL_STATES: [CallState; 6] = [
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
];
pub const HOM_VAR_CALL_STATES: [CallState; 3] =
    [CallState::HomVar1, CallState::HomVar2, CallState::HomVar3];
pub const VAR_CALL_STATES: [CallState; 9] = [
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
    CallState::HomVar1,
    CallState::HomVar2,
    CallState::HomVar3,
];

/// The contingency table, as counts of (truth, call) pairs.
///
/// The counts are `double` in the reference and not integers, because `MISSING_SITES_HOM_REF` adds
/// the whole uncovered genome to one cell as a computed difference. `count` casts to a long the
/// way `getCount` does; nothing else rounds.
#[derive(Debug, Clone, Default)]
pub struct Counts {
    cells: Vec<((TruthState, CallState), f64)>,
}

impl Counts {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn increment(&mut self, truth: TruthState, call: CallState) {
        self.increment_by(truth, call, 1.0);
    }

    pub fn increment_by(&mut self, truth: TruthState, call: CallState, by: f64) {
        match self.cells.iter_mut().find(|(key, _)| *key == (truth, call)) {
            Some((_, value)) => *value += by,
            None => self.cells.push(((truth, call), by)),
        }
    }

    /// `getCount`, which is the bin's value CAST to a long: a cell holding 2.5 counts as two.
    pub fn count(&self, truth: TruthState, call: CallState) -> i64 {
        self.raw(truth, call) as i64
    }

    pub fn raw(&self, truth: TruthState, call: CallState) -> f64 {
        self.cells
            .iter()
            .find(|(key, _)| *key == (truth, call))
            .map_or(0.0, |(_, value)| *value)
    }

    /// `getCounterSize`, the sum of every bin, which is what the missing-missing cell is computed
    /// against.
    pub fn size(&self) -> f64 {
        self.cells.iter().map(|(_, value)| *value).sum()
    }

    /// `validateCountsAgainstScheme`: a pair the scheme calls unreachable must not have happened.
    ///
    /// The refusal names the two states, and it is a `PicardException` rather than a wrong number:
    /// a count in an NA cell means the state machine reached a pair its own scheme says it cannot.
    pub fn validate_against(
        &self,
        scheme: &[(CallState, [Cell; 11])],
    ) -> Result<(), (TruthState, CallState)> {
        for truth in TRUTH_DECLARATION_ORDER {
            for call in CALL_DECLARATION_ORDER {
                if self.count(truth, call) > 0
                    && contingency(scheme, call, truth) == Some(Cell::Unreachable)
                {
                    return Err((truth, call));
                }
            }
        }
        Ok(())
    }

    fn cell_values(
        scheme: &[(CallState, [Cell; 11])],
        truth: TruthState,
        call: CallState,
    ) -> &'static [ContingencyState] {
        match contingency(scheme, call, truth) {
            Some(Cell::Values(values)) => values,
            // An unreachable cell contributes nothing to a ratio: its array holds `NA`, which is
            // neither a TP nor an FP nor an FN nor a TN.
            _ => &[],
        }
    }

    /// `getSensitivity`: TP / (TP + FN) over a subset of truth states.
    ///
    /// A cell that contributes several contingency values is added ONCE PER VALUE, so a pair that
    /// is both a TP and an FN puts its count in the numerator once and in the denominator twice.
    /// That is the reference's arithmetic and not a rounding of it.
    pub fn sensitivity(&self, scheme: &[(CallState, [Cell; 11])], truths: &[TruthState]) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for truth in truths {
            for call in CALL_DECLARATION_ORDER {
                let count = self.count(*truth, call) as f64;
                for value in Self::cell_values(scheme, *truth, call) {
                    match value {
                        ContingencyState::Tp => {
                            numerator += count;
                            denominator += count;
                        }
                        ContingencyState::Fn => denominator += count,
                        _ => {}
                    }
                }
            }
        }
        numerator / denominator
    }

    /// `Ppv`: TP / (TP + FP) over a subset of CALL states.
    pub fn ppv(&self, scheme: &[(CallState, [Cell; 11])], calls: &[CallState]) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for call in calls {
            for truth in TRUTH_DECLARATION_ORDER {
                let count = self.count(truth, *call) as f64;
                for value in Self::cell_values(scheme, truth, *call) {
                    match value {
                        ContingencyState::Tp => {
                            numerator += count;
                            denominator += count;
                        }
                        ContingencyState::Fp => denominator += count,
                        _ => {}
                    }
                }
            }
        }
        numerator / denominator
    }

    /// `getSpecificity`: TN / (FP + TN) over a subset of truth states.
    pub fn specificity(&self, scheme: &[(CallState, [Cell; 11])], truths: &[TruthState]) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for truth in truths {
            for call in CALL_DECLARATION_ORDER {
                let count = self.count(*truth, call) as f64;
                for value in Self::cell_values(scheme, *truth, call) {
                    match value {
                        ContingencyState::Tn => {
                            numerator += count;
                            denominator += count;
                        }
                        ContingencyState::Fp => denominator += count,
                        _ => {}
                    }
                }
            }
        }
        numerator / denominator
    }

    /// `calculateGenotypeConcordanceUtil`: the diagonal over everything counted.
    ///
    /// Three things decide what is counted. `missing_sites` false drops every pair where either
    /// side is MISSING; `include_hom_ref` false drops every pair where neither side is a variant,
    /// which is what makes the non-reference concordance a different number; and the diagonal is
    /// the pair whose two CODES are equal, so the six call-only states are never on it.
    ///
    /// An empty denominator is `NaN`, which the metrics file writes as `?`.
    pub fn genotype_concordance(&self, missing_sites: bool, include_hom_ref: bool) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for truth in TRUTH_DECLARATION_ORDER {
            for call in CALL_DECLARATION_ORDER {
                if !missing_sites && (truth == TruthState::Missing || call == CallState::Missing) {
                    continue;
                }
                if include_hom_ref || is_var(truth, call) {
                    let count = self.count(truth, call) as f64;
                    if truth_code(truth) == call_code(call) {
                        numerator += count;
                    }
                    denominator += count;
                }
            }
        }
        if denominator > 0.0 {
            numerator / denominator
        } else {
            f64::NAN
        }
    }

    /// `getContingencyStateCounts`, in the order the reference's enum declares: TP, FP, TN, FN,
    /// and then EMPTY. `NA` is counted too and never written.
    pub fn contingency_counts(&self, scheme: &[(CallState, [Cell; 11])]) -> ContingencyCounts {
        let mut counts = ContingencyCounts::default();
        for truth in TRUTH_DECLARATION_ORDER {
            for call in CALL_DECLARATION_ORDER {
                let count = self.count(truth, call);
                match contingency(scheme, call, truth) {
                    // `EMPTY` is a one-element array holding the EMPTY state, so an empty cell
                    // adds its count to the EMPTY counter exactly once.
                    Some(Cell::Values([])) => counts.empty += count,
                    Some(Cell::Values(values)) => {
                        for value in values {
                            match value {
                                ContingencyState::Tp => counts.tp += count,
                                ContingencyState::Fp => counts.fp += count,
                                ContingencyState::Tn => counts.tn += count,
                                ContingencyState::Fn => counts.fn_ += count,
                            }
                        }
                    }
                    Some(Cell::Unreachable) | None => {}
                }
            }
        }
        counts
    }
}

/// The five counters the contingency metrics file writes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ContingencyCounts {
    pub tp: i64,
    pub tn: i64,
    pub fp: i64,
    pub fn_: i64,
    pub empty: i64,
}

/// `isVar`: whether EITHER side of the pair says there is a variant.
pub fn is_var(truth: TruthState, call: CallState) -> bool {
    matches!(
        truth,
        TruthState::HomVar1 | TruthState::HetRefVar1 | TruthState::HetVar1Var2
    ) || matches!(
        call,
        CallState::HetRefVar1
            | CallState::HetRefVar2
            | CallState::HetRefVar3
            | CallState::HetVar1Var2
            | CallState::HetVar1Var3
            | CallState::HetVar3Var4
            | CallState::HomVar1
            | CallState::HomVar2
            | CallState::HomVar3
    )
}

/// One sample's genotype at a site, reduced to what the state machine reads.
///
/// A no-call allele is the string `.`, which is what `Allele.getBaseString` answers for one, so
/// the splice below can compare against it the way the reference does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypeView {
    pub alleles: Vec<String>,
    pub is_filtered: bool,
    /// `-1` when the field is absent, which is the value the thresholds skip on.
    pub gq: i32,
    pub dp: i32,
}

impl GenotypeView {
    /// `Genotype.isNoCall()`: every allele is a no-call.
    pub fn is_no_call(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(|allele| allele == NO_CALL_STRING)
    }

    /// `Genotype.isMixed()`: called on one chromosome and not on the other.
    pub fn is_mixed(&self) -> bool {
        !self.is_no_call() && self.alleles.iter().any(|allele| allele == NO_CALL_STRING)
    }
}

/// A site as one side of the comparison sees it, already subset to its own sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteView {
    /// `VariantContext.isMixed()`, which is a property of the ALLELES: a site carrying both a SNP
    /// and an indel alternate. It is not the genotype's mixedness.
    pub is_mixed: bool,
    pub is_filtered: bool,
    pub reference: String,
    pub genotype: GenotypeView,
}

/// `Allele.NO_CALL_STRING`.
pub const NO_CALL_STRING: &str = ".";
/// `Allele.SPAN_DEL_STRING`, which a splice leaves alone.
pub const SPAN_DEL_STRING: &str = "*";

/// `GenotypeConcordanceStateCodes`, the answer of the checks that run before any allele is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StateCode {
    Missing,
    HomRef,
    HetRefVar1,
    HetVar1Var2,
    HomVar1,
    NoCall,
    LowGq,
    LowDp,
    VcFiltered,
    GtFiltered,
    IsMixed,
}

/// `getStateCode`: the state a site takes before its alleles are looked at, or none.
///
/// The order of the checks is the answer when several hold at once, and the reference says so:
/// "if a variant context has BOTH GQ and DP less than the specified threshold, then it will be of
/// Truth/Call State LOW_GQ". A genotype called on one chromosome only comes back as `NO_CALL`
/// rather than as a state of its own.
pub fn state_code(site: Option<&SiteView>, min_gq: i32, min_dp: i32) -> Option<StateCode> {
    let site = site?;
    if site.is_mixed {
        return Some(StateCode::IsMixed);
    }
    if site.is_filtered {
        return Some(StateCode::VcFiltered);
    }
    let genotype = &site.genotype;
    if genotype.is_no_call() {
        return Some(StateCode::NoCall);
    }
    if genotype.is_filtered {
        return Some(StateCode::GtFiltered);
    }
    if genotype.gq != -1 && genotype.gq < min_gq {
        return Some(StateCode::LowGq);
    }
    if genotype.dp != -1 && genotype.dp < min_dp {
        return Some(StateCode::LowDp);
    }
    if genotype.is_mixed() {
        return Some(StateCode::NoCall);
    }
    None
}

fn truth_of_code(code: StateCode) -> TruthState {
    match code {
        StateCode::Missing => TruthState::Missing,
        StateCode::HomRef => TruthState::HomRef,
        StateCode::HetRefVar1 => TruthState::HetRefVar1,
        StateCode::HetVar1Var2 => TruthState::HetVar1Var2,
        StateCode::HomVar1 => TruthState::HomVar1,
        StateCode::NoCall => TruthState::NoCall,
        StateCode::LowGq => TruthState::LowGq,
        StateCode::LowDp => TruthState::LowDp,
        StateCode::VcFiltered => TruthState::VcFiltered,
        StateCode::GtFiltered => TruthState::GtFiltered,
        StateCode::IsMixed => TruthState::IsMixed,
    }
}

fn call_of_code(code: StateCode) -> CallState {
    match code {
        StateCode::Missing => CallState::Missing,
        StateCode::HomRef => CallState::HomRef,
        StateCode::HetRefVar1 => CallState::HetRefVar1,
        StateCode::HetVar1Var2 => CallState::HetVar1Var2,
        StateCode::HomVar1 => CallState::HomVar1,
        StateCode::NoCall => CallState::NoCall,
        StateCode::LowGq => CallState::LowGq,
        StateCode::LowDp => CallState::LowDp,
        StateCode::VcFiltered => CallState::VcFiltered,
        StateCode::GtFiltered => CallState::GtFiltered,
        StateCode::IsMixed => CallState::IsMixed,
    }
}

/// `spliceOrAppendString`: put the reference's extra bases where they belong in an allele.
fn splice_or_append(destination: &str, to_insert: &str, insert_at: usize) -> String {
    if destination == SPAN_DEL_STRING {
        return destination.to_string();
    }
    if insert_at <= destination.len() {
        format!(
            "{}{}{}",
            &destination[..insert_at],
            to_insert,
            &destination[insert_at..]
        )
    } else {
        format!("{destination}{to_insert}")
    }
}

/// What `normalizeAlleles` produced: the ordered set the states are indices into, and the four
/// alleles as they read after the reference lengths were reconciled.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Alleles {
    pub all: Vec<String>,
    pub truth1: Option<String>,
    pub truth2: Option<String>,
    pub call1: Option<String>,
    pub call2: Option<String>,
}

/// `normalizeAlleles`: one list of alleles both sides index into.
///
/// Two things happen here and neither is obvious.
///
/// The references are reconciled first. Two files can spell the same indel against different
/// reference alleles (`TCAA/T` against `TCAACAA/TCAA`), so the shorter reference must be a prefix
/// of the longer one and its extra bases are spliced into the shorter side's alleles. A pair whose
/// references are the same length and different is an error rather than a mismatch.
///
/// Then the ORDER is decided by the call and not by the truth. The truth's two alleles go in
/// first, but if either call allele landed past index one, the two variant slots are emptied and
/// refilled with the truth's alleles REVERSED, so that the allele the call shares becomes var1.
/// That is why `A -> C/G` against `C/A` is HET_VAR1_VAR2 against HET_REF_VAR1 whichever way the
/// truth wrote its two alleles.
pub fn normalize_alleles(
    truth: Option<&SiteView>,
    call: Option<&SiteView>,
    ignore_filtered: bool,
) -> Result<Alleles, String> {
    let truth_genotype = truth.filter(|site| !site.is_mixed && !site.is_filtered);
    let call_genotype =
        call.filter(|site| !site.is_mixed && (ignore_filtered || !site.is_filtered));

    let mut truth_ref = truth_genotype.map(|site| site.reference.clone());
    let mut call_ref = call_genotype.map(|site| site.reference.clone());

    let alleles_of = |site: Option<&SiteView>| -> Result<Option<(String, String)>, String> {
        match site {
            None => Ok(None),
            Some(site) => {
                if site.genotype.alleles.len() != 2 {
                    return Err("does not have exactly 2 alleles".to_string());
                }
                Ok(Some((
                    site.genotype.alleles[0].clone(),
                    site.genotype.alleles[1].clone(),
                )))
            }
        }
    };
    let (mut truth1, mut truth2) = match alleles_of(truth_genotype)? {
        Some((a, b)) => (Some(a), Some(b)),
        None => (None, None),
    };
    let (mut call1, mut call2) = match alleles_of(call_genotype)? {
        Some((a, b)) => (Some(a), Some(b)),
        None => (None, None),
    };

    if let (Some(t), Some(c)) = (truth_ref.clone(), call_ref.clone()) {
        if t != c {
            let splice = |allele: &mut Option<String>, suffix: &str, at: usize| {
                if let Some(value) = allele {
                    if value != NO_CALL_STRING {
                        *value = splice_or_append(value, suffix, at);
                    }
                }
            };
            if t.len() < c.len() {
                let suffix = c
                    .strip_prefix(t.as_str())
                    .ok_or_else(|| "Ref alleles mismatch".to_string())?
                    .to_string();
                let at = t.len();
                splice(&mut truth1, &suffix, at);
                splice(&mut truth2, &suffix, at);
                truth_ref = Some(format!("{t}{suffix}"));
            } else if t.len() > c.len() {
                let suffix = t
                    .strip_prefix(c.as_str())
                    .ok_or_else(|| "Ref alleles mismatch".to_string())?
                    .to_string();
                let at = c.len();
                splice(&mut call1, &suffix, at);
                splice(&mut call2, &suffix, at);
                call_ref = Some(format!("{c}{suffix}"));
            } else {
                return Err("Ref alleles mismatch".to_string());
            }
        }
    }

    let mut all: Vec<String> = Vec::new();
    let smart_add = |all: &mut Vec<String>, allele: &Option<String>| {
        if let Some(value) = allele {
            if !all.contains(value) {
                all.push(value.clone());
            }
        }
    };
    if truth_genotype.is_some() || call_genotype.is_some() {
        let zeroth = if truth_genotype.is_none() {
            call_ref.clone()
        } else {
            truth_ref.clone()
        };
        smart_add(&mut all, &zeroth);
    }
    if truth_genotype.is_some() {
        smart_add(&mut all, &truth1);
        smart_add(&mut all, &truth2);
    }
    if call_genotype.is_some() {
        let index_of = |all: &[String], allele: &Option<String>| -> i32 {
            allele
                .as_ref()
                .and_then(|value| all.iter().position(|other| other == value))
                .map_or(-1, |position| position as i32)
        };
        if index_of(&all, &call1) > 1 || index_of(&all, &call2) > 1 {
            all.remove(2);
            all.remove(1);
            smart_add(&mut all, &truth2);
            smart_add(&mut all, &truth1);
        }
        smart_add(&mut all, &call1);
        smart_add(&mut all, &call2);
    }

    Ok(Alleles {
        all,
        truth1,
        truth2,
        call1,
        call2,
    })
}

/// `TruthState.getHom` / `getVar`, which refuse an index they have no name for.
fn truth_hom(index: i32) -> Result<TruthState, String> {
    match index {
        0 => Ok(TruthState::HomRef),
        1 => Ok(TruthState::HomVar1),
        _ => Err("Shouldn't be here.".to_string()),
    }
}

fn truth_var(first: i32, second: i32) -> Result<TruthState, String> {
    match (first, second) {
        (0, 1) | (1, 0) => Ok(TruthState::HetRefVar1),
        (1, 2) | (2, 1) => Ok(TruthState::HetVar1Var2),
        _ => Err("Shouldn't be here.".to_string()),
    }
}

fn call_hom(index: i32) -> Result<CallState, String> {
    match index {
        0 => Ok(CallState::HomRef),
        1 => Ok(CallState::HomVar1),
        2 => Ok(CallState::HomVar2),
        3 => Ok(CallState::HomVar3),
        _ => Err("Shouldn't be here.".to_string()),
    }
}

/// `CallState.getHet`, which sorts the two indices first, so `2/1` and `1/2` are one state.
fn call_het(first: i32, second: i32) -> Result<CallState, String> {
    let (low, high) = if first > second {
        (second, first)
    } else {
        (first, second)
    };
    match (low, high) {
        (0, 1) => Ok(CallState::HetRefVar1),
        (0, 2) => Ok(CallState::HetRefVar2),
        (0, 3) => Ok(CallState::HetRefVar3),
        (1, 2) => Ok(CallState::HetVar1Var2),
        (1, 3) => Ok(CallState::HetVar1Var3),
        // Not a mistake: VAR2/VAR3 is symbolic, so it is folded into VAR3/VAR4.
        (2, 3) | (3, 4) => Ok(CallState::HetVar3Var4),
        _ => Err("Shouldn't be here.".to_string()),
    }
}

/// `determineState`: the pair of states one site contributes.
///
/// `ignore_filtered` does not skip a filtered call: it CLEARS the state the filter would have
/// given, so the call is read from its alleles like any other. A filtered truth site has no such
/// escape.
pub fn determine_state(
    truth: Option<&SiteView>,
    call: Option<&SiteView>,
    min_gq: i32,
    min_dp: i32,
    ignore_filtered: bool,
) -> Result<(TruthState, CallState), String> {
    let truth_code_value = state_code(truth, min_gq, min_dp);
    let mut call_code_value = state_code(call, min_gq, min_dp);
    if ignore_filtered && call_code_value == Some(StateCode::VcFiltered) {
        call_code_value = None;
    }
    let mut truth_state = truth_code_value.map(truth_of_code);
    let mut call_state = call_code_value.map(call_of_code);

    let alleles = normalize_alleles(
        if truth_state.is_none() { truth } else { None },
        if call_state.is_none() { call } else { None },
        ignore_filtered,
    )?;
    let index_of = |allele: &Option<String>| -> i32 {
        allele
            .as_ref()
            .and_then(|value| alleles.all.iter().position(|other| other == value))
            .map_or(-1, |position| position as i32)
    };

    if truth_state.is_none() {
        let first = index_of(&alleles.truth1);
        let second = index_of(&alleles.truth2);
        truth_state = Some(if first == second {
            truth_hom(first)?
        } else {
            truth_var(first, second)?
        });
    }
    if call_state.is_none() {
        let first = index_of(&alleles.call1);
        let second = index_of(&alleles.call2);
        call_state = Some(if first == second {
            call_hom(first)?
        } else {
            call_het(first, second)?
        });
    }
    Ok((
        truth_state.expect("a truth state"),
        call_state.expect("a call state"),
    ))
}
