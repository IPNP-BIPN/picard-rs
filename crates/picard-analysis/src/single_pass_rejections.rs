//! The refusals `SinglePassSamProgram` reaches before a single record is counted, ported from
//! Picard 3.4.0.
//!
//! `CollectQualityYieldMetrics` and `CollectAlignmentSummaryMetrics` share a driver, and the
//! driver decides three things before either of them sees a read: whether the input's sort order
//! is acceptable, whether an obsolete argument was passed, and, once the sort check is bypassed,
//! what happens when the reference walker is asked to go backwards.
//!
//! The covering arrays are generated over combinations the tool ACCEPTS (decision 0009), so a row
//! spent being rejected is a row not spent on the tool. That leaves these paths untested by
//! construction, and they are behaviour: a port that happily processes a queryname-sorted BAM
//! where Picard refuses is not a byte-identical port, it is a different tool with the same name.
//!
//! # The sort check names the file and offers the escape in the same breath
//!
//! ```java
//! throw new PicardException("File " + input.getAbsolutePath() + " should be coordinate sorted but "
//!         + "the header says the sort order is " + sortOrder.name() + ". If you believe the file "
//!         + "to be coordinate sorted you may pass ASSUME_SORTED=true");
//! ```
//!
//! The message carries the ABSOLUTE path, the sort order that was found, and the argument that
//! would bypass the check. All three are in the golden.
//!
//! # And taking the escape moves the failure rather than removing it
//!
//! With `ASSUME_SORTED=true` the driver proceeds, and the reference walker is then asked for a
//! contig it has already passed. That refusal comes from htsjdk rather than from Picard, and it
//! names the two indices rather than the contigs:
//!
//! ```text
//! htsjdk.samtools.SAMException: Requesting earlier reference sequence: 0 < 1
//! ```
//!
//! So the two arguments do not choose between success and failure. They choose which of two
//! exceptions a queryname-sorted input produces, from two different libraries.

/// `htsjdk.samtools.SAMFileHeader.SortOrder`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
    Unknown,
}

impl SortOrder {
    /// `SortOrder.name()`, which is the enum constant as written.
    pub fn name(self) -> &'static str {
        match self {
            SortOrder::Unsorted => "unsorted",
            SortOrder::Queryname => "queryname",
            SortOrder::Coordinate => "coordinate",
            SortOrder::Duplicate => "duplicate",
            SortOrder::Unknown => "unknown",
        }
    }
}

/// What the driver refuses, and which library refuses it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Rejection {
    /// `FLOW_MODE` on `CollectQualityYieldMetrics`, which moved to another tool.
    FlowModeObsolete,
    /// A sort order the program cannot use, without `ASSUME_SORTED`.
    NotCoordinateSorted {
        /// The input's absolute path, as `getAbsolutePath` gives it.
        path: String,
        found: SortOrder,
    },
    /// The reference walker asked to go backwards, which is what `ASSUME_SORTED=true` leads to on
    /// a queryname-sorted input.
    EarlierReferenceSequence { requested: i32, current: i32 },
}

impl Rejection {
    pub fn java_class(&self) -> &'static str {
        match self {
            Rejection::FlowModeObsolete | Rejection::NotCoordinateSorted { .. } => {
                "picard.PicardException"
            }
            Rejection::EarlierReferenceSequence { .. } => "htsjdk.samtools.SAMException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            Rejection::FlowModeObsolete => "FLOW_MODE is obsolete. Flow support now provided by \
                 CollectQualityYieldMetricsFlow"
                .to_string(),
            Rejection::NotCoordinateSorted { path, found } => format!(
                "File {path} should be coordinate sorted but the header says the sort order is \
                 {}. If you believe the file to be coordinate sorted you may pass \
                 ASSUME_SORTED=true",
                found.name()
            ),
            Rejection::EarlierReferenceSequence { requested, current } => {
                format!("Requesting earlier reference sequence: {requested} < {current}")
            }
        }
    }

    /// `<class>: <message>`, which is how the harness records a throwable.
    pub fn thrown(&self) -> String {
        format!("{}: {}", self.java_class(), self.message())
    }
}

/// `CollectQualityYieldMetrics.customCommandLineValidation`, in effect: the obsolete argument is
/// refused before anything else is looked at.
pub fn check_flow_mode(flow_mode: bool) -> Result<(), Rejection> {
    if flow_mode {
        return Err(Rejection::FlowModeObsolete);
    }
    Ok(())
}

/// `SinglePassSamProgram.makeItSo`'s sort-order assertion.
///
/// `assume_sorted` does not make the input sorted; it makes the driver stop asking.
pub fn check_sort_order(
    path: &str,
    found: SortOrder,
    assume_sorted: bool,
) -> Result<(), Rejection> {
    if assume_sorted || found == SortOrder::Coordinate {
        return Ok(());
    }
    Err(Rejection::NotCoordinateSorted {
        path: path.to_string(),
        found,
    })
}

/// `ReferenceSequenceFileWalker.get`, which refuses to rewind.
///
/// The walker keeps the index it last served and compares the request against it, so a
/// queryname-sorted input reaches this as soon as two consecutive records are on contigs in
/// descending order.
pub fn walk_reference(current: Option<i32>, requested: i32) -> Result<i32, Rejection> {
    if let Some(current) = current {
        if requested < current {
            return Err(Rejection::EarlierReferenceSequence { requested, current });
        }
    }
    Ok(requested)
}
