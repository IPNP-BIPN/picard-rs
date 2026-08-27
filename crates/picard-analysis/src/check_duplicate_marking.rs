//! `CheckDuplicateMarking`: whether every record of a query name agrees about its duplicate flag.
//!
//! Reading the file and sorting it by query name are not ported: the tool sorts the input itself
//! when it is not already query-name sorted, and then walks it once. What is ported is that walk,
//! which is the whole verdict.
//!
//! Ported from `picard.sam.markduplicates.CheckDuplicateMarking` in Picard 3.4.0.

/// `CheckDuplicateMarking.Mode`, which decides what is skipped before anything is compared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    All,
    PrimaryOnly,
    PrimaryMappedOnly,
    PrimaryProperPairOnly,
}

/// One record, reduced to what the walk reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub name: String,
    pub duplicate: bool,
    pub secondary_or_supplementary: bool,
    pub unmapped: bool,
    pub proper_pair: bool,
}

/// The three `continue`s at the top of `checkDuplicateMarkingsInIterator`, in their own order.
///
/// Every mode but `ALL` skips the secondary and supplementary records, `PRIMARY_MAPPED_ONLY`
/// skips the unmapped ones as well, and `PRIMARY_PROPER_PAIR_ONLY` skips everything that is not a
/// proper pair, which includes every unpaired read. `PRIMARY_PROPER_PAIR_ONLY` does NOT test the
/// unmapped flag: an unmapped record that is somehow flagged a proper pair is kept.
pub fn is_skipped(record: &Record, mode: Mode) -> bool {
    if mode != Mode::All && record.secondary_or_supplementary {
        return true;
    }
    if mode == Mode::PrimaryMappedOnly && record.unmapped {
        return true;
    }
    if mode == Mode::PrimaryProperPairOnly && !record.proper_pair {
        return true;
    }
    false
}

/// What one walk found: the bad records' names in the order they were met, one line each.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Verdict {
    pub bad_names: Vec<String>,
}

impl Verdict {
    /// `numBadRecords > 0 ? 1 : 0`, which is the count's SIGN and never the count.
    pub fn exit_code(&self) -> i32 {
        i32::from(!self.bad_names.is_empty())
    }
}

/// `checkDuplicateMarkingsInIterator` over records already in query-name order.
///
/// The comparison is against the FIRST record of each name and not against the record before it,
/// because a disagreeing record returns early WITHOUT updating the remembered flag. A name whose
/// flags go true, false, false therefore reports two bad records and writes its name twice.
///
/// The mode's skip happens before the tally, so a skipped record does not become the name's first
/// one either: what is remembered is the first record the mode KEPT.
pub fn check(records: &[Record], mode: Mode) -> Verdict {
    let mut verdict = Verdict::default();
    let mut current_name: Option<String> = None;
    let mut current_duplicate = false;
    for record in records {
        if is_skipped(record, mode) {
            continue;
        }
        if current_name.as_deref() != Some(record.name.as_str()) {
            current_name = Some(record.name.clone());
            current_duplicate = record.duplicate;
        } else if record.duplicate != current_duplicate {
            verdict.bad_names.push(record.name.clone());
        }
    }
    verdict
}
