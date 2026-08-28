//! `CollectDuplicateMetrics`: the marks a duplicate-marked file already carries, tallied.
//!
//! The tool never marks anything itself. Reading the file is not ported; the tally is, along with
//! the derived fields and the halving that happens once the walk is over.
//!
//! Ported from `picard.sam.markduplicates.CollectDuplicateMetrics`,
//! `picard.sam.DuplicationMetrics` and
//! `picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram` in Picard 3.4.0.

use std::collections::BTreeMap;

use crate::jumping_library::estimate_library_size;

/// `LibraryIdGenerator.getLibraryName`, where the read group names no library.
pub const UNKNOWN_LIBRARY: &str = "Unknown Library";

/// One record, reduced to what the two tallies read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub library: String,
    pub duplicate: bool,
    pub secondary_or_supplementary: bool,
    pub unmapped: bool,
    pub paired: bool,
    pub mate_unmapped: bool,
}

/// `DuplicationMetrics`, as the walk fills it in.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DuplicationMetrics {
    pub library: String,
    pub unpaired_reads_examined: i64,
    /// Counts READS while the walk runs and pairs once it is over.
    pub read_pairs_examined: i64,
    pub secondary_or_supplementary_reads: i64,
    pub unmapped_reads: i64,
    pub unpaired_read_duplicates: i64,
    /// Likewise.
    pub read_pair_duplicates: i64,
    /// Always zero here: this tool never looks for an optical duplicate.
    pub read_pair_optical_duplicates: i64,
}

impl DuplicationMetrics {
    /// `addReadToLibraryMetrics`: a CHAIN of four, so a read reaches exactly one counter.
    ///
    /// An unmapped read is counted as unmapped and nothing else, and a secondary or supplementary
    /// one as secondary and nothing else. A paired read whose mate is unmapped falls to the
    /// unpaired counter, which is why a half-mapped pair contributes one of each.
    pub fn add_read(&mut self, record: &Record) {
        if record.unmapped {
            self.unmapped_reads += 1;
        } else if record.secondary_or_supplementary {
            self.secondary_or_supplementary_reads += 1;
        } else if !record.paired || record.mate_unmapped {
            self.unpaired_reads_examined += 1;
        } else {
            self.read_pairs_examined += 1;
        }
    }

    /// `addDuplicateReadToMetrics`, whose guard is its OWN and not the one above: an unmapped or
    /// secondary duplicate is counted nowhere.
    pub fn add_duplicate(&mut self, record: &Record) {
        if record.secondary_or_supplementary || record.unmapped {
            return;
        }
        if !record.paired || record.mate_unmapped {
            self.unpaired_read_duplicates += 1;
        } else {
            self.read_pair_duplicates += 1;
        }
    }

    /// `finalizeAndWriteMetrics`' halving, which is an INTEGER division: a lone paired read
    /// reports no pairs examined at all.
    pub fn halve_the_pairs(&mut self) {
        self.read_pairs_examined /= 2;
        self.read_pair_duplicates /= 2;
    }

    /// `calculateDerivedFields`' library size, computed from the optical count this tool leaves
    /// at zero, which is what makes it slightly wrong by the tool's own admission.
    pub fn estimated_library_size(&self) -> Option<i64> {
        estimate_library_size(
            self.read_pairs_examined - self.read_pair_optical_duplicates,
            self.read_pairs_examined - self.read_pair_duplicates,
        )
    }

    /// `calculateDerivedFields`' percentage, which weights the PAIRS by two on both sides.
    pub fn percent_duplication(&self) -> f64 {
        let examined = self.unpaired_reads_examined + self.read_pairs_examined * 2;
        if examined == 0 {
            return 0.0;
        }
        (self.unpaired_read_duplicates + self.read_pair_duplicates * 2) as f64 / examined as f64
    }
}

/// The whole run: one row per library, in the order the libraries are named.
///
/// The rows come from `LibraryIdGenerator`, which is built from the HEADER, so every library a
/// read group names gets a row whether a read ever used it or not. A file with no reads at all
/// still writes a row of zeros.
pub fn collect(header_libraries: &[String], records: &[Record]) -> Vec<DuplicationMetrics> {
    let mut by_library: BTreeMap<&str, DuplicationMetrics> = BTreeMap::new();
    let mut order: Vec<&str> = Vec::new();
    for library in header_libraries {
        if by_library
            .insert(
                library.as_str(),
                DuplicationMetrics {
                    library: library.clone(),
                    ..Default::default()
                },
            )
            .is_none()
        {
            order.push(library.as_str());
        }
    }
    for record in records {
        let metrics = by_library
            .entry(record.library.as_str())
            .or_insert_with(|| {
                order.push(record.library.as_str());
                DuplicationMetrics {
                    library: record.library.clone(),
                    ..Default::default()
                }
            });
        metrics.add_read(record);
        if record.duplicate {
            metrics.add_duplicate(record);
        }
    }
    let mut rows: Vec<DuplicationMetrics> = order
        .into_iter()
        .map(|library| by_library[library].clone())
        .collect();
    for row in &mut rows {
        row.halve_the_pairs();
    }
    rows
}

/// `finalizeAndWriteMetrics`' last step: the coverage histogram is written only for a file of ONE
/// library, and only when that library has an estimated size.
pub fn writes_a_histogram(rows: &[DuplicationMetrics]) -> bool {
    rows.len() == 1 && rows[0].estimated_library_size().is_some()
}
