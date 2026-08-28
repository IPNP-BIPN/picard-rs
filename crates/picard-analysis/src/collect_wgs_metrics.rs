//! `CollectWgsMetrics`: how deeply each base of the reference is covered, and why the rest of the
//! reads did not count.
//!
//! Walking the file and building the pileup are not ported. What is ported is the accounting: the
//! seven exclusion counters, which are a partition of the bases that did not reach the histogram,
//! and the metrics derived from the depths that did.
//!
//! Ported from `picard.analysis.CollectWgsMetrics` and `picard.analysis.WgsMetrics` in
//! Picard 3.4.0.

/// `CollectWgsMetrics.MINIMUM_MAPPING_QUALITY`.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 20;
/// `CollectWgsMetrics.MINIMUM_BASE_QUALITY`.
pub const DEFAULT_MINIMUM_BASE_QUALITY: i32 = 20;
/// `CollectWgsMetrics.COVERAGE_CAP`.
pub const DEFAULT_COVERAGE_CAP: i32 = 250;
/// `CollectWgsMetrics.LOCUS_ACCUMULATION_CAP`.
pub const DEFAULT_LOCUS_ACCUMULATION_CAP: i32 = 100_000;

/// Why one base of one read did not reach the depth histogram, or that it did.
///
/// The order is the order the tool tests in, and it is a CHAIN: a base that fails an earlier test
/// never reaches a later one, which is what makes the seven counters a partition rather than
/// seven independent tallies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// It counted, and the base's depth went up by one.
    Counted,
    /// The read was flagged an adapter read.
    Adapter,
    /// The read's mapping quality was under the floor. A WHOLE read goes here.
    MappingQuality,
    /// The read was a duplicate. A whole read again.
    Duplicate,
    /// The read was unpaired, or paired with its mate unmapped, and `COUNT_UNPAIRED` was not set.
    Unpaired,
    /// The base's quality was under the floor. An `N` base is here whatever its quality says.
    BaseQuality,
    /// The base was already counted from this pair's other end.
    Overlap,
    /// The base's depth had already reached the cap.
    Capped,
}

/// One read, reduced to what the chain reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Read {
    pub adapter: bool,
    pub mapping_quality: i32,
    pub duplicate: bool,
    pub paired: bool,
    pub mate_unmapped: bool,
}

/// The arguments the chain consults.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arguments {
    pub minimum_mapping_quality: i32,
    pub minimum_base_quality: i32,
    pub coverage_cap: i32,
    pub count_unpaired: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            minimum_mapping_quality: DEFAULT_MINIMUM_MAPPING_QUALITY,
            minimum_base_quality: DEFAULT_MINIMUM_BASE_QUALITY,
            coverage_cap: DEFAULT_COVERAGE_CAP,
            count_unpaired: false,
        }
    }
}

/// `CollectRawWgsMetrics`' four defaults, which are the whole of that tool.
///
/// It is `CollectWgsMetrics` with the mapping-quality floor dropped to nothing, the base-quality
/// floor dropped to three, and the two caps raised past anything a real file reaches. Every other
/// rule is the same one, which is why one port serves both.
pub const RAW_MINIMUM_MAPPING_QUALITY: i32 = 0;
pub const RAW_MINIMUM_BASE_QUALITY: i32 = 3;
pub const RAW_COVERAGE_CAP: i32 = 100_000;
pub const RAW_LOCUS_ACCUMULATION_CAP: i32 = 200_000;

/// `CollectRawWgsMetrics`' arguments, which differ from the default in four values.
pub fn raw_arguments() -> Arguments {
    Arguments {
        minimum_mapping_quality: RAW_MINIMUM_MAPPING_QUALITY,
        minimum_base_quality: RAW_MINIMUM_BASE_QUALITY,
        coverage_cap: RAW_COVERAGE_CAP,
        count_unpaired: false,
    }
}

/// The tests that take a WHOLE read, in the order the tool applies them.
pub fn read_fate(read: &Read, arguments: &Arguments) -> Option<Fate> {
    if read.adapter {
        return Some(Fate::Adapter);
    }
    if read.mapping_quality < arguments.minimum_mapping_quality {
        return Some(Fate::MappingQuality);
    }
    if read.duplicate {
        return Some(Fate::Duplicate);
    }
    // A pair with one end unmapped is unpaired for this purpose, which is what makes an
    // unmapped mate reach the same counter as a read that was never paired.
    if !arguments.count_unpaired && (!read.paired || read.mate_unmapped) {
        return Some(Fate::Unpaired);
    }
    None
}

/// The tests that take a single BASE, once its read has survived.
///
/// An `N` is excluded by QUALITY and not by a rule of its own: the tool reads it as quality zero,
/// so it lands under any floor above nought.
pub fn base_fate(
    base: u8,
    quality: i32,
    already_counted_here: bool,
    depth_here: i32,
    arguments: &Arguments,
) -> Fate {
    let quality = if base.eq_ignore_ascii_case(&b'N') {
        0
    } else {
        quality
    };
    if quality < arguments.minimum_base_quality {
        return Fate::BaseQuality;
    }
    if already_counted_here {
        return Fate::Overlap;
    }
    if depth_here >= arguments.coverage_cap {
        return Fate::Capped;
    }
    Fate::Counted
}

/// What one run counted, in bases.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Counts {
    pub counted: i64,
    pub adapter: i64,
    pub mapping_quality: i64,
    pub duplicate: i64,
    pub unpaired: i64,
    pub base_quality: i64,
    pub overlap: i64,
    pub capped: i64,
}

impl Counts {
    pub fn add(&mut self, fate: Fate) {
        let counter = match fate {
            Fate::Counted => &mut self.counted,
            Fate::Adapter => &mut self.adapter,
            Fate::MappingQuality => &mut self.mapping_quality,
            Fate::Duplicate => &mut self.duplicate,
            Fate::Unpaired => &mut self.unpaired,
            Fate::BaseQuality => &mut self.base_quality,
            Fate::Overlap => &mut self.overlap,
            Fate::Capped => &mut self.capped,
        };
        *counter += 1;
    }

    /// Every base the walk looked at, counted or not, which is what the percentages divide by.
    pub fn total(&self) -> i64 {
        self.counted
            + self.adapter
            + self.mapping_quality
            + self.duplicate
            + self.unpaired
            + self.base_quality
            + self.overlap
            + self.capped
    }

    /// The excluded bases, which is the total less the counted ones.
    pub fn excluded(&self) -> i64 {
        self.total() - self.counted
    }

    /// One counter as a fraction of the total, which is what a `PCT_EXC_` column holds. It is a
    /// fraction and not a percentage, whatever the name says.
    pub fn fraction(&self, count: i64) -> f64 {
        if self.total() == 0 {
            0.0
        } else {
            count as f64 / self.total() as f64
        }
    }

    /// `PCT_EXC_TOTAL`, which is the seven summed and therefore never over one.
    pub fn excluded_fraction(&self) -> f64 {
        self.fraction(self.excluded())
    }
}

/// `MEAN_COVERAGE`: the counted bases over the TERRITORY, not over the covered bases, so an
/// uncovered half halves it.
pub fn mean_coverage(counted: i64, genome_territory: i64) -> f64 {
    if genome_territory == 0 {
        0.0
    } else {
        counted as f64 / genome_territory as f64
    }
}

/// `GENOME_TERRITORY`: the reference's non-`N` bases, which is smaller than its length.
pub fn genome_territory(reference: &[u8]) -> i64 {
    reference
        .iter()
        .filter(|base| !base.eq_ignore_ascii_case(&b'N'))
        .count() as i64
}
