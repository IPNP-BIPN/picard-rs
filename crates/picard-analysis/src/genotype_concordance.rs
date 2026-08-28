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
pub fn contingency_values(cell: Cell) -> Option<String> {
    match cell {
        Cell::Unreachable => None,
        Cell::Values(values) => Some(
            values
                .iter()
                .map(|value| value.name())
                .collect::<Vec<_>>()
                .join(","),
        ),
    }
}
