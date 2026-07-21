//! The gene annotation model.
//!
//! Ported from `picard.annotation.Gene`, its inner `Transcript` and `Exon`, and
//! `picard.annotation.LocusFunction`, tag 3.4.0.
//!
//! Everything here is **1-based inclusive**. The refFlat file that feeds it is 0-based half-open,
//! and the conversion happens in the reader (`refflat.rs`), not here, so a coordinate that
//! reaches this module has already been rebased. Stated because the two conventions differ by
//! one at the start and agree at the end (a half-open end equals an inclusive end), and mixing
//! them is the classic way to be off by one at exactly one boundary.
//!
//! The single load-bearing subtlety is `LocusFunction`'s **ordinal order**, which htsjdk's own
//! comment calls out: a base that is coding in one transcript and UTR in another is counted as
//! coding, by taking the higher ordinal. htsjdk-rs decision 0020 rests on this being a max
//! reduction, which is why the order of the enum is not cosmetic and is reproduced exactly.

/// `LocusFunction`, in its declared order, which **is** its meaning: a higher ordinal wins when
/// two transcripts disagree about a base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LocusFunction {
    Intergenic = 0,
    Intronic = 1,
    Utr = 2,
    Coding = 3,
    Ribosomal = 4,
}

impl LocusFunction {
    pub fn name(self) -> &'static str {
        match self {
            LocusFunction::Intergenic => "INTERGENIC",
            LocusFunction::Intronic => "INTRONIC",
            LocusFunction::Utr => "UTR",
            LocusFunction::Coding => "CODING",
            LocusFunction::Ribosomal => "RIBOSOMAL",
        }
    }
}

/// `Gene.Transcript.Exon`: a 1-based inclusive interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Exon {
    pub start: i32,
    pub end: i32,
}

/// `Gene.Transcript`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub name: String,
    pub transcription_start: i32,
    pub transcription_end: i32,
    pub coding_start: i32,
    pub coding_end: i32,
    pub exons: Vec<Exon>,
    /// The number of bases across all exons, accumulated as exons are added.
    length: i32,
}

impl Transcript {
    pub fn new(
        name: &str,
        transcription_start: i32,
        transcription_end: i32,
        coding_start: i32,
        coding_end: i32,
    ) -> Self {
        Transcript {
            name: name.to_string(),
            transcription_start,
            transcription_end,
            coding_start,
            coding_end,
            exons: Vec::new(),
            length: 0,
        }
    }

    /// `addExon(start, end)`: appends and grows `length` by the exon's base count.
    ///
    /// `CoordMath.getLength(start, end)` is `end - start + 1`, the inclusive count.
    pub fn add_exon(&mut self, start: i32, end: i32) {
        self.exons.push(Exon { start, end });
        self.length += end - start + 1;
    }

    pub fn length(&self) -> i32 {
        self.length
    }

    /// `utr(locus)`: outside the coding range on either side.
    fn is_utr(&self, locus: i32) -> bool {
        locus < self.coding_start || locus > self.coding_end
    }

    /// `inExon(locus)`.
    ///
    /// The early return on `exon.start > locus` assumes the exons are sorted ascending, which
    /// the reader guarantees (and rejects a refFlat where they are not). It is not a mere
    /// optimisation: it is the reason a locus before the first exon is intronic rather than
    /// scanned against every exon.
    fn in_exon(&self, locus: i32) -> bool {
        for exon in &self.exons {
            if exon.start > locus {
                return false;
            }
            if locus >= exon.start && locus <= exon.end {
                return true;
            }
        }
        false
    }

    /// `assignLocusFunctionForRange(start, locusFunctions)`.
    ///
    /// Writes the per-base `LocusFunction` for the window beginning at genome position `start`,
    /// taking the **max** by ordinal so a stronger classification from another transcript is not
    /// overwritten. The `> CODING` guard skips positions already at the ceiling (only RIBOSOMAL
    /// is higher, and that comes from the ribosomal interval list, not from a transcript).
    pub fn assign_locus_function_for_range(
        &self,
        start: i32,
        locus_functions: &mut [LocusFunction],
    ) {
        let window_end = start + locus_functions.len() as i32 - 1;
        let from = start.max(self.transcription_start);
        let to = self.transcription_end.min(window_end);
        let mut i = from;
        while i <= to {
            let idx = (i - start) as usize;
            if locus_functions[idx] > LocusFunction::Coding {
                i += 1;
                continue;
            }
            let lf = if self.in_exon(i) {
                if self.is_utr(i) {
                    LocusFunction::Utr
                } else {
                    LocusFunction::Coding
                }
            } else {
                LocusFunction::Intronic
            };
            if lf > locus_functions[idx] {
                locus_functions[idx] = lf;
            }
            i += 1;
        }
    }

    /// `getTranscriptCoordinate(genomeCoordinate)`: 1-based position within the spliced
    /// transcript, or `-1` when the genome coordinate falls in an intron.
    pub fn transcript_coordinate(&self, genome_coordinate: i32) -> i32 {
        let mut exon_offset = 0;
        for exon in &self.exons {
            if genome_coordinate >= exon.start && genome_coordinate <= exon.end {
                return (genome_coordinate - exon.start + 1) + exon_offset;
            }
            exon_offset += exon.end - exon.start + 1;
        }
        -1
    }

    /// `addCoverageCounts(genomeStart, genomeEnd, coverage)`.
    ///
    /// Note the loop is `for i in genomeStart..genomeEnd` — **half-open**, so the base at
    /// `genomeEnd` is not counted. htsjdk passes `CoordMath.getEnd(start, length)` as the end,
    /// which is inclusive, so the effect is that the last base of the block is dropped from
    /// coverage. Reproduced exactly; it is the kind of boundary a reimplementation silently
    /// "fixes".
    pub fn add_coverage_counts(&self, genome_start: i32, genome_end: i32, coverage: &mut [i32]) {
        let mut i = genome_start;
        while i < genome_end {
            let tx_base = self.transcript_coordinate(i);
            if tx_base > 0 {
                coverage[(tx_base - 1) as usize] += 1;
            }
            i += 1;
        }
    }
}

/// `Gene`: a named, stranded interval carrying its transcripts.
///
/// `Gene extends Interval`, and the extent (`start`, `end`) is the min transcription start and
/// max transcription end across its transcripts, computed by the reader.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    pub name: String,
    pub transcripts: Vec<Transcript>,
}

impl Gene {
    pub fn new(contig: &str, start: i32, end: i32, negative_strand: bool, name: &str) -> Self {
        Gene {
            contig: contig.to_string(),
            start,
            end,
            negative_strand,
            name: name.to_string(),
            transcripts: Vec::new(),
        }
    }

    /// The gene's extent, `end - start + 1`.
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn transcript() -> Transcript {
        // A two-exon transcript, coding in the middle, so all four functions are reachable.
        let mut t = Transcript::new("tx", 100, 300, 150, 250);
        t.add_exon(100, 180);
        t.add_exon(220, 300);
        t
    }

    #[test]
    fn the_ordinal_order_is_the_priority_order() {
        assert!(LocusFunction::Coding > LocusFunction::Utr);
        assert!(LocusFunction::Utr > LocusFunction::Intronic);
        assert!(LocusFunction::Ribosomal > LocusFunction::Coding);
    }

    #[test]
    fn exon_length_is_inclusive() {
        let t = transcript();
        // 100..180 is 81 bases, 220..300 is 81 bases.
        assert_eq!(t.length(), 162);
    }

    #[test]
    fn a_base_is_coding_utr_or_intronic_by_position() {
        let t = transcript();
        let mut lf = vec![LocusFunction::Intergenic; 201];
        t.assign_locus_function_for_range(100, &mut lf);
        // 160 is in exon 1 and inside the coding range: CODING.
        assert_eq!(lf[160 - 100], LocusFunction::Coding);
        // 120 is in exon 1 but before codingStart 150: UTR.
        assert_eq!(lf[120 - 100], LocusFunction::Utr);
        // 200 is between the exons: INTRONIC.
        assert_eq!(lf[200 - 100], LocusFunction::Intronic);
    }

    /// A stronger classification already present is not overwritten by a weaker one.
    #[test]
    fn assignment_takes_the_max_not_the_last() {
        let t = transcript();
        let mut lf = vec![LocusFunction::Intergenic; 201];
        lf[160 - 100] = LocusFunction::Ribosomal; // higher than anything the transcript assigns
        t.assign_locus_function_for_range(100, &mut lf);
        assert_eq!(
            lf[160 - 100],
            LocusFunction::Ribosomal,
            "ribosomal survives"
        );
    }

    #[test]
    fn transcript_coordinate_splices_out_the_intron() {
        let t = transcript();
        assert_eq!(t.transcript_coordinate(100), 1, "first base of exon 1");
        assert_eq!(t.transcript_coordinate(180), 81, "last base of exon 1");
        assert_eq!(
            t.transcript_coordinate(220),
            82,
            "first base of exon 2 follows directly"
        );
        assert_eq!(t.transcript_coordinate(200), -1, "in the intron");
    }

    /// The coverage loop is half-open, so the last base of the range is not counted.
    #[test]
    fn coverage_is_half_open_at_the_end() {
        let t = transcript();
        let mut coverage = vec![0; t.length() as usize];
        // Cover genome 100..=102 by passing end = 103 (exclusive), as htsjdk does with the
        // inclusive getEnd. Bases 100, 101, 102 count; 103 would be exon position 4 but is not
        // reached.
        t.add_coverage_counts(100, 103, &mut coverage);
        assert_eq!(&coverage[..4], &[1, 1, 1, 0]);
    }

    #[test]
    fn a_gene_extent_is_inclusive() {
        let g = Gene::new("chr1", 100, 300, false, "GENEA");
        assert_eq!(g.length(), 201);
    }
}
