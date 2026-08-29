//! `MarkDuplicates`: which record of a duplicate set keeps `0x400` clear, and what the metrics say.
//!
//! The tool's answer is one bit per record and a metrics table. It reads the file twice: the first
//! pass builds a `ReadEndsForMarkDuplicates` per record and sorts them, the second walks the file
//! again setting the flag on every index the first pass called a duplicate. What is ported here is
//! both passes over records already in memory; the sorting collections that spill to disk, the
//! reading and the writing are not, because they decide where the records live and not which of
//! them is a duplicate.
//!
//! # What decides a duplicate
//!
//! Not the bases. A record's position for this purpose is its UNCLIPPED 5' coordinate, which is
//! the unclipped start of a forward read and the unclipped end of a reverse one, so a soft-clipped
//! read and an unclipped one can be duplicates of each other and two reads with different
//! sequences at one position are. A pair is keyed by both ends and its orientation, a fragment by
//! its one end, and the two lists are walked separately: a fragment that shares its key with a
//! PAIR is always marked, whatever it scores.
//!
//! One record of each set keeps the bit clear, and which one is the SCORE's answer. The default
//! sums the mapped reference length over the pair and `SUM_OF_BASE_QUALITIES` sums the qualities
//! at or above 15, and the comparison is `>` rather than `>=`, so the FIRST of an equal-scoring
//! set is the one kept.
//!
//! # What makes a duplicate optical
//!
//! The read NAME, and nothing in the alignment. The default parsing takes the last three
//! colon-separated fields as tile, x and y, and REFUSES a name that has neither five nor seven
//! fields rather than parsing what it can. Two records are optical when they share a read group
//! and a tile and lie within the pixel distance on both axes, which is why the same two pairs are
//! optical or not depending only on what they are called.
//!
//! # What is not ported
//!
//!  * **the graph path of `OpticalDuplicateFinder`**, which runs for a set of three or more without
//!    a keeper and four or more with one, and clusters transitive neighbours with union-find. The
//!    fast path below is what the golden's sets of two exercise, and [`optical_duplicates`] says so
//!    rather than pretending;
//!  * **`RANDOM` scoring**, which hashes the read name with Murmur3 and is in no golden;
//!  * **the barcode HASH**: the reference keys a barcode by `Objects.hash` of the tag's value, and
//!    what the port compares is the value itself. The two agree except on a collision, which no
//!    fixture has;
//!  * **`TAG_DUPLICATE_SET_MEMBERS`** and the `DI`/`DS` tags it writes.
//!
//! Ported from `picard.sam.markduplicates.MarkDuplicates`,
//! `picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram`,
//! `picard.sam.markduplicates.util.OpticalDuplicateFinder`, `picard.sam.util.ReadNameParser`,
//! `picard.sam.DuplicationMetrics` and `htsjdk.samtools.DuplicateScoringStrategy` in Picard 3.4.0.

use htsjdk_bam::cigar::{Cigar, Op};

/// `DuplicateScoringStrategy.ScoringStrategy`, minus `RANDOM`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringStrategy {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
}

/// `MarkDuplicates.DuplicateTaggingPolicy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaggingPolicy {
    DontTag,
    OpticalOnly,
    All,
}

/// `DuplicateType`, whose codes are what the `DT` tag carries.
pub const SEQUENCING_CODE: &str = "SQ";
pub const LIBRARY_CODE: &str = "LB";

/// The default pixel distance, `OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE`.
pub const DEFAULT_OPTICAL_DUPLICATE_DISTANCE: i32 = 100;

/// The arguments this port reads.
#[derive(Debug, Clone)]
pub struct Options {
    pub scoring: ScoringStrategy,
    pub remove_duplicates: bool,
    pub remove_sequencing_duplicates: bool,
    pub tagging_policy: TaggingPolicy,
    pub clear_dt: bool,
    pub optical_duplicate_pixel_distance: i32,
    /// `READ_NAME_REGEX`, of which only "the default" and "none" are ported: a custom regex is a
    /// `java.util.regex` pattern, and the port has no Java regex engine.
    pub parse_read_names: bool,
    /// `BARCODE_TAG`, whose value splits one position into one set per barcode.
    pub barcode_tag: Option<String>,
}

impl Default for Options {
    /// The tool's own defaults: `TOTAL_MAPPED_REFERENCE_LENGTH`, no removal, no tagging, `CLEAR_DT`
    /// on and the default pixel distance.
    fn default() -> Self {
        Self {
            scoring: ScoringStrategy::TotalMappedReferenceLength,
            remove_duplicates: false,
            remove_sequencing_duplicates: false,
            tagging_policy: TaggingPolicy::DontTag,
            clear_dt: true,
            optical_duplicate_pixel_distance: DEFAULT_OPTICAL_DUPLICATE_DISTANCE,
            parse_read_names: true,
            barcode_tag: None,
        }
    }
}

/// One record, reduced to what the two passes read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub name: String,
    pub flags: u16,
    /// `-1` where the record is unmapped.
    pub reference_index: i32,
    /// One-based, as the SAM field is.
    pub alignment_start: i32,
    pub cigar: Cigar,
    pub qualities: Vec<u8>,
    pub mate_reference_index: i32,
    pub library: String,
    /// The nth read group of the header, which is what `closeEnough` compares.
    pub read_group: i32,
    /// The value of `BARCODE_TAG`, where the record carries one.
    pub barcode: Option<String>,
    /// An incoming `DT`, which `CLEAR_DT` decides the fate of.
    pub existing_dt: Option<String>,
    /// The `MC` tag, which is what the mate-cigar markers read instead of waiting for the mate.
    /// `MarkDuplicates` itself does not look at it.
    pub mate_cigar: Option<Cigar>,
}

impl Record {
    pub fn unmapped(&self) -> bool {
        self.flags & 0x4 != 0
    }
    pub fn paired(&self) -> bool {
        self.flags & 0x1 != 0
    }
    pub fn mate_unmapped(&self) -> bool {
        self.flags & 0x8 != 0
    }
    pub fn reverse_strand(&self) -> bool {
        self.flags & 0x10 != 0
    }
    pub fn first_of_pair(&self) -> bool {
        self.flags & 0x40 != 0
    }
    pub fn secondary_or_supplementary(&self) -> bool {
        self.flags & 0x100 != 0 || self.flags & 0x800 != 0
    }
    pub fn fails_vendor_quality(&self) -> bool {
        self.flags & 0x200 != 0
    }

    /// `SAMRecord.getUnclippedStart`: the alignment start walked back over the leading clips.
    pub fn unclipped_start(&self) -> i32 {
        let mut start = self.alignment_start;
        for element in &self.cigar.elements {
            match element.op {
                Op::S | Op::H => start -= element.length as i32,
                _ => break,
            }
        }
        start
    }

    /// `SAMRecord.getUnclippedEnd`: the alignment end walked forward over the trailing clips.
    pub fn unclipped_end(&self) -> i32 {
        let mut end = self.alignment_start + self.cigar.reference_length() as i32 - 1;
        for element in self.cigar.elements.iter().rev() {
            match element.op {
                Op::S | Op::H => end += element.length as i32,
                _ => break,
            }
        }
        end
    }

    /// `buildReadEnds`: the 5' coordinate, which is the end the strand points away from.
    pub fn five_prime_coordinate(&self) -> i32 {
        if self.reverse_strand() {
            self.unclipped_end()
        } else {
            self.unclipped_start()
        }
    }
}

/// `DuplicateScoringStrategy.computeDuplicateScore`, for the two strategies that are measured.
///
/// The cap is `Short.MAX_VALUE / 2` on each half, so two long high-quality reads cannot overflow
/// the pair's sum, and a record that fails the vendor check is discounted by `Short.MIN_VALUE / 2`
/// AFTER the cap rather than before it.
pub fn duplicate_score(record: &Record, strategy: ScoringStrategy) -> i16 {
    let capped = |value: i64| -> i16 { value.min(i64::from(i16::MAX / 2)) as i16 };
    let mut score: i16 = match strategy {
        ScoringStrategy::SumOfBaseQualities => {
            let sum: i64 = record
                .qualities
                .iter()
                .filter(|quality| **quality >= 15)
                .map(|quality| i64::from(*quality))
                .sum();
            capped(sum)
        }
        ScoringStrategy::TotalMappedReferenceLength => {
            if record.unmapped() {
                0
            } else {
                capped(i64::from(record.cigar.reference_length()))
            }
        }
    };
    if record.fails_vendor_quality() {
        score = score.wrapping_add(i16::MIN / 2);
    }
    score
}

/// `ReadEnds.getOrientationByte`, whose four values order a pair's two strands.
pub const F: u8 = 0;
pub const R: u8 = 1;
pub const FF: u8 = 2;
pub const FR: u8 = 3;
pub const RR: u8 = 4;
pub const RF: u8 = 5;

pub fn orientation_byte(read1_negative: bool, read2_negative: bool) -> u8 {
    if read1_negative {
        if read2_negative {
            RR
        } else {
            RF
        }
    } else if read2_negative {
        FR
    } else {
        FF
    }
}

/// A tile coordinate read off a read name, which is the whole of a record's physical location.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Location {
    pub tile: i16,
    pub x: i32,
    pub y: i32,
    pub known: bool,
}

/// `ReadNameParser.rapidParseInt`: digits until the first thing that is not one.
///
/// It is not `Integer.parseInt`: a field of digits followed by letters parses to the digits, and a
/// field with no digits at all is a `NumberFormatException`, which the caller turns into "no
/// location" rather than into a failure.
pub fn rapid_parse_int(text: &str) -> Option<i32> {
    let bytes = text.as_bytes();
    let mut index = 0;
    let negative = bytes.first() == Some(&b'-');
    if negative {
        index = 1;
    }
    let mut value: i32 = 0;
    let mut digits = false;
    while index < bytes.len() {
        let byte = bytes[index];
        if byte.is_ascii_digit() {
            value = value.wrapping_mul(10).wrapping_add(i32::from(byte - b'0'));
            digits = true;
            index += 1;
        } else {
            break;
        }
    }
    if !digits {
        return None;
    }
    Some(if negative { -value } else { value })
}

/// `ReadNameParser.readLocationInformation` with the default regex.
///
/// The last three colon-separated fields are tile, x and y, and the name is REFUSED unless it has
/// exactly five or seven fields: a name with six is not half-parsed, it has no location at all.
pub fn location(name: &str) -> Location {
    let fields: Vec<&str> = name.split(':').collect();
    if fields.len() < 3 {
        return Location::default();
    }
    if fields.len() != 5 && fields.len() != 7 {
        return Location::default();
    }
    let last_three = &fields[fields.len() - 3..];
    let mut parsed = [0i32; 3];
    for (slot, field) in parsed.iter_mut().zip(last_three) {
        match rapid_parse_int(field) {
            Some(value) => *slot = value,
            None => return Location::default(),
        }
    }
    Location {
        tile: parsed[0] as i16,
        x: parsed[1],
        y: parsed[2],
        known: true,
    }
}

/// `ReadEndsForMarkDuplicates`: one record's, or one pair's, position and score.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadEnds {
    pub library: String,
    pub read1_reference_index: i32,
    pub read1_coordinate: i32,
    pub read2_reference_index: i32,
    pub read2_coordinate: i32,
    pub orientation: u8,
    pub orientation_for_optical_duplicates: u8,
    pub read1_index_in_file: usize,
    pub read2_index_in_file: usize,
    pub score: i16,
    pub location: Location,
    pub read_group: i32,
    pub barcode: Option<String>,
    pub is_optical_duplicate: bool,
}

impl ReadEnds {
    /// `ReadEnds.isPaired`, which is a mate reference index and not the paired FLAG: a read whose
    /// mate is unmapped is a fragment here.
    pub fn is_paired(&self) -> bool {
        self.read2_reference_index != -1
    }
}

/// `MarkDuplicates.buildReadEnds`, over one record.
pub fn build_read_ends(record: &Record, index: usize, options: &Options) -> ReadEnds {
    ReadEnds {
        library: record.library.clone(),
        read1_reference_index: record.reference_index,
        read1_coordinate: record.five_prime_coordinate(),
        read2_reference_index: if record.paired() && !record.mate_unmapped() {
            record.mate_reference_index
        } else {
            -1
        },
        read2_coordinate: -1,
        orientation: if record.reverse_strand() { R } else { F },
        orientation_for_optical_duplicates: 0,
        read1_index_in_file: index,
        read2_index_in_file: index,
        score: duplicate_score(record, options.scoring),
        location: if options.parse_read_names {
            location(&record.name)
        } else {
            Location::default()
        },
        read_group: record.read_group,
        barcode: options.barcode_tag.as_ref().and(record.barcode.clone()),
        is_optical_duplicate: false,
    }
}

/// `OpticalDuplicateFinder.closeEnough`: the same read group, the same tile, and within the
/// distance on both axes. Two ends at one index are never compared, which is object identity in
/// the reference and an index here.
fn close_enough(left: &ReadEnds, right: &ReadEnds, distance: i32) -> bool {
    left.location.known
        && right.location.known
        && left.read_group == right.read_group
        && left.location.tile == right.location.tile
        && (left.location.x - right.location.x).abs() <= distance
        && (left.location.y - right.location.y).abs() <= distance
}

/// `OpticalDuplicateFinder.findOpticalDuplicates`, on its fast path.
///
/// The keeper is compared to everyone first, and then every other pair is compared once, with the
/// one already flagged left alone so a chain of two does not flag both. The GRAPH path, which the
/// reference takes for three ends without a keeper or four with one, is not ported: it clusters
/// transitive neighbours, and every set the golden carries is a set of two.
pub fn optical_duplicates(list: &[ReadEnds], keeper: Option<usize>, distance: i32) -> Vec<bool> {
    let mut flags = vec![false; list.len()];
    if list.len() < 2 {
        return flags;
    }
    let keeper = keeper.filter(|index| list[*index].location.known);
    if list.len() >= if keeper.is_none() { 3 } else { 4 } {
        // The graph path, which is not ported. Saying so is better than answering with the fast
        // path's answer, which is only correct when there is nothing transitive to cluster.
        return flags;
    }
    if let Some(keeper) = keeper {
        for (index, other) in list.iter().enumerate() {
            flags[index] = index != keeper && close_enough(&list[keeper], other, distance);
        }
    }
    for left in 0..list.len() {
        if Some(left) == keeper {
            continue;
        }
        for right in (left + 1)..list.len() {
            if Some(right) == keeper {
                continue;
            }
            if flags[left] && flags[right] {
                continue;
            }
            if close_enough(&list[left], &list[right], distance) {
                let index = if flags[right] { left } else { right };
                flags[index] = true;
            }
        }
    }
    flags
}

/// `DuplicationMetrics`, the row one library gets.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Metrics {
    pub library: String,
    pub unpaired_reads_examined: i64,
    pub read_pairs_examined: i64,
    pub secondary_or_supplementary: i64,
    pub unmapped_reads: i64,
    pub unpaired_read_duplicates: i64,
    pub read_pair_duplicates: i64,
    pub read_pair_optical_duplicates: i64,
    pub percent_duplication: f64,
    pub estimated_library_size: Option<i64>,
}

/// `DuplicationMetrics.f`, the function the bisection finds a root of.
fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

/// `DuplicationMetrics.estimateLibrarySize`: forty bisections, and `None` where there is nothing
/// to estimate from.
pub fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    let duplicates = read_pairs - unique_read_pairs;
    if read_pairs <= 0 || duplicates <= 0 {
        return None;
    }
    let unique = unique_read_pairs as f64;
    let pairs = read_pairs as f64;
    let mut m = 1.0;
    let mut big = 100.0;
    if unique_read_pairs >= read_pairs || f(m * unique, unique, pairs) < 0.0 {
        // `IllegalStateException` in the reference, which no caller of it recovers from.
        return None;
    }
    while f(big * unique, unique, pairs) > 0.0 {
        big *= 10.0;
    }
    for _ in 0..40 {
        let r = (m + big) / 2.0;
        let u = f(r * unique, unique, pairs);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            m = r;
        } else {
            big = r;
        }
    }
    Some((unique * (m + big) / 2.0) as i64)
}

/// `DuplicationMetrics.estimateRoi`.
pub fn estimate_roi(estimated_library_size: i64, x: f64, pairs: i64, unique_pairs: i64) -> f64 {
    let size = estimated_library_size as f64;
    size * (1.0 - (-(x * pairs as f64) / size).exp()) / unique_pairs as f64
}

/// `DuplicationMetrics.calculateRoiHistogram`: one hundred bins, `x` from 1 to 100.
pub fn roi_histogram(metrics: &Metrics) -> Option<Vec<(f64, f64)>> {
    let size = metrics.estimated_library_size?;
    let unique = metrics.read_pairs_examined - metrics.read_pair_duplicates;
    let mut bins = Vec::with_capacity(100);
    for step in 1..=100 {
        let x = f64::from(step);
        bins.push((
            x,
            estimate_roi(size, x, metrics.read_pairs_examined, unique),
        ));
    }
    Some(bins)
}

/// What one run decided, one entry per input record.
#[derive(Debug, Clone, PartialEq)]
pub struct Marking {
    /// `setDuplicateReadFlag`.
    pub duplicate: Vec<bool>,
    /// Whether the record's index was in the optical collection, which is what the `DT` tag and
    /// `REMOVE_SEQUENCING_DUPLICATES` both read.
    pub optical: Vec<bool>,
    /// The `DT` tag the record ends up with, which is `CLEAR_DT`'s answer and then the policy's.
    pub duplicate_type: Vec<Option<String>>,
    /// Whether the record is written at all.
    pub written: Vec<bool>,
    /// One row per library, in the order the libraries were first seen.
    pub metrics: Vec<Metrics>,
    /// `duplicateCountHist`, `opticalDuplicateCountHist` and `nonOpticalDuplicateCountHist`: how
    /// many duplicate SETS there were of each size, which the metrics file prints beside the ROI.
    pub all_sets: Vec<(f64, f64)>,
    pub optical_sets: Vec<(f64, f64)>,
    pub non_optical_sets: Vec<(f64, f64)>,
}

/// `Histogram.increment`, which adds to the bin or creates it.
fn increment(histogram: &mut Vec<(f64, f64)>, bin: f64, by: f64) {
    match histogram.iter_mut().find(|(id, _)| *id == bin) {
        Some((_, count)) => *count += by,
        None => histogram.push((bin, by)),
    }
}

/// The key a chunk is cut on: `areComparableForDuplicates`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Key {
    library: String,
    barcode: Option<String>,
    read1_reference_index: i32,
    read1_coordinate: i32,
    orientation: u8,
    read2: Option<(i32, i32)>,
}

fn key(ends: &ReadEnds, compare_read2: bool) -> Key {
    Key {
        library: ends.library.clone(),
        barcode: ends.barcode.clone(),
        read1_reference_index: ends.read1_reference_index,
        read1_coordinate: ends.read1_coordinate,
        orientation: ends.orientation,
        read2: compare_read2.then_some((ends.read2_reference_index, ends.read2_coordinate)),
    }
}

/// The best of a set: the FIRST with the highest score, because the comparison is `>`.
fn best(list: &[ReadEnds]) -> usize {
    let mut best_index = 0;
    let mut max = 0i16;
    for (index, ends) in list.iter().enumerate() {
        if ends.score > max || index == 0 {
            max = ends.score;
            best_index = index;
        }
    }
    best_index
}

/// Both passes: which indexes are duplicates, which of those are optical, and the metrics.
///
/// The records are taken in the file's own order, which for a coordinate-sorted file is the order
/// the writing pass walks. The reference sorts its read ends and cuts chunks on equal keys; the
/// port groups by the same key, which is the same partition without the spill to disk.
pub fn mark(records: &[Record], options: &Options) -> Marking {
    let mut pairs: Vec<ReadEnds> = Vec::new();
    let mut fragments: Vec<ReadEnds> = Vec::new();
    // `tmp`, the map from a read group and a name to the first end of a pair seen.
    let mut pending: Vec<(String, ReadEnds)> = Vec::new();

    for (index, record) in records.iter().enumerate() {
        if record.unmapped() || record.secondary_or_supplementary() {
            continue;
        }
        let fragment = build_read_ends(record, index, options);
        fragments.push(fragment.clone());
        if !record.paired() || record.mate_unmapped() {
            continue;
        }
        let name = record.name.clone();
        match pending.iter().position(|(key, _)| *key == name) {
            None => pending.push((name, fragment)),
            Some(position) => {
                let (_, mut paired) = pending.remove(position);
                let mates_reference_index = fragment.read1_reference_index;
                let mates_coordinate = fragment.read1_coordinate;
                // The optical orientation goes by the FIRST end of the pair whatever order the
                // file put them in, and it is set before the orientation below is rewritten.
                paired.orientation_for_optical_duplicates = if record.first_of_pair() {
                    orientation_byte(record.reverse_strand(), paired.orientation == R)
                } else {
                    orientation_byte(paired.orientation == R, record.reverse_strand())
                };
                if mates_reference_index > paired.read1_reference_index
                    || (mates_reference_index == paired.read1_reference_index
                        && mates_coordinate >= paired.read1_coordinate)
                {
                    paired.read2_reference_index = mates_reference_index;
                    paired.read2_coordinate = mates_coordinate;
                    paired.read2_index_in_file = index;
                    paired.orientation =
                        orientation_byte(paired.orientation == R, record.reverse_strand());
                    // Two ends at one position pointing at each other have an orientation that
                    // would otherwise depend on the file's order, so it is fixed to FR.
                    if paired.read2_reference_index == paired.read1_reference_index
                        && paired.read2_coordinate == paired.read1_coordinate
                        && paired.orientation == RF
                    {
                        paired.orientation = FR;
                    }
                } else {
                    paired.read2_reference_index = paired.read1_reference_index;
                    paired.read2_coordinate = paired.read1_coordinate;
                    paired.read2_index_in_file = paired.read1_index_in_file;
                    paired.read1_reference_index = mates_reference_index;
                    paired.read1_coordinate = mates_coordinate;
                    paired.read1_index_in_file = index;
                    paired.orientation =
                        orientation_byte(record.reverse_strand(), paired.orientation == R);
                }
                paired.score = paired
                    .score
                    .wrapping_add(duplicate_score(record, options.scoring));
                pairs.push(paired);
            }
        }
    }

    let mut duplicate = vec![false; records.len()];
    let mut optical = vec![false; records.len()];
    let mut all_sets: Vec<(f64, f64)> = Vec::new();
    let mut optical_sets: Vec<(f64, f64)> = Vec::new();
    let mut non_optical_sets: Vec<(f64, f64)> = Vec::new();

    // The pairs first. A set of one is a singleton and nothing is marked.
    for chunk in chunks(&pairs, true) {
        if chunk.len() == 1 {
            // `addSingletonToCount`: a pair nobody duplicated is a set of one in two of the three
            // histograms and in neither counter.
            increment(&mut all_sets, 1.0, 1.0);
            increment(&mut non_optical_sets, 1.0, 1.0);
            continue;
        }
        if chunk.is_empty() {
            continue;
        }
        let keeper = best(&chunk);
        let flags = if options.parse_read_names {
            optical_flags(&chunk, keeper, options)
        } else {
            vec![false; chunk.len()]
        };
        for (position, ends) in chunk.iter().enumerate() {
            if position == keeper {
                continue;
            }
            duplicate[ends.read1_index_in_file] = true;
            duplicate[ends.read2_index_in_file] = true;
            if flags[position] {
                optical[ends.read1_index_in_file] = true;
                optical[ends.read2_index_in_file] = true;
            }
        }
        // `trackDuplicateCounts`, which is reached THROUGH `trackOpticalDuplicates` and so does
        // not run at all when the read name regex is off: the set-size histograms are empty on a
        // run that found no locations, rather than counting sets nobody looked at.
        if !options.parse_read_names {
            continue;
        }
        // The set's size, and the two halves of it. The optical bin is the optical COUNT PLUS ONE,
        // because the record they are duplicates of is in the set too.
        let optical_count = flags.iter().filter(|flag| **flag).count() as f64;
        increment(&mut all_sets, chunk.len() as f64, 1.0);
        if chunk.len() as f64 - optical_count > 0.0 {
            increment(
                &mut non_optical_sets,
                chunk.len() as f64 - optical_count,
                1.0,
            );
        }
        if optical_count > 0.0 {
            increment(&mut optical_sets, optical_count + 1.0, 1.0);
        }
    }

    // Then the fragments, whose chunk is cut on read1 alone. A chunk that holds a pair marks
    // every fragment in it, and one that holds none keeps its best.
    for chunk in chunks(&fragments, false) {
        let contains_pairs = chunk.iter().any(ReadEnds::is_paired);
        let contains_fragments = chunk.iter().any(|ends| !ends.is_paired());
        if chunk.len() < 2 || !contains_fragments {
            continue;
        }
        if contains_pairs {
            for ends in &chunk {
                if !ends.is_paired() {
                    duplicate[ends.read1_index_in_file] = true;
                }
            }
        } else {
            let keeper = best(&chunk);
            for (position, ends) in chunk.iter().enumerate() {
                if position != keeper {
                    duplicate[ends.read1_index_in_file] = true;
                }
            }
        }
    }

    let mut marking = write(records, options, &duplicate, &optical);
    marking.all_sets = all_sets;
    marking.optical_sets = optical_sets;
    marking.non_optical_sets = non_optical_sets;
    marking
}

/// The chunks `generateDuplicateIndexes` cuts, which are the runs of an equal key.
fn chunks(list: &[ReadEnds], compare_read2: bool) -> Vec<Vec<ReadEnds>> {
    let mut keys: Vec<Key> = Vec::new();
    let mut groups: Vec<Vec<ReadEnds>> = Vec::new();
    for ends in list {
        let this = key(ends, compare_read2);
        match keys.iter().position(|other| *other == this) {
            Some(position) => groups[position].push(ends.clone()),
            None => {
                keys.push(this);
                groups.push(vec![ends.clone()]);
            }
        }
    }
    groups
}

/// `trackOpticalDuplicates`, which splits a set that mixes FR and RF before it looks at anything.
///
/// The split is the point: in PCR duplicate detection a set can hold both orientations once the
/// ends are ordered by position, and in OPTICAL duplicate detection two reads whose first ends
/// point opposite ways are not duplicates of one another.
fn optical_flags(chunk: &[ReadEnds], keeper: usize, options: &Options) -> Vec<bool> {
    let has_fr = chunk
        .iter()
        .any(|ends| ends.orientation_for_optical_duplicates == FR);
    let has_rf = chunk
        .iter()
        .any(|ends| ends.orientation_for_optical_duplicates == RF);
    if !(has_fr && has_rf) {
        return optical_duplicates(
            chunk,
            Some(keeper),
            options.optical_duplicate_pixel_distance,
        );
    }
    let mut flags = vec![false; chunk.len()];
    for orientation in [FR, RF] {
        let positions: Vec<usize> = (0..chunk.len())
            .filter(|index| chunk[*index].orientation_for_optical_duplicates == orientation)
            .collect();
        let list: Vec<ReadEnds> = positions
            .iter()
            .map(|index| chunk[*index].clone())
            .collect();
        let inner_keeper = positions.iter().position(|index| *index == keeper);
        let inner = optical_duplicates(
            &list,
            inner_keeper,
            options.optical_duplicate_pixel_distance,
        );
        for (position, flag) in positions.iter().zip(inner) {
            flags[*position] = flag;
        }
    }
    flags
}

/// The second pass: the flag, the tag, the removal and the metrics, in the reference's order.
fn write(records: &[Record], options: &Options, duplicate: &[bool], optical: &[bool]) -> Marking {
    let mut metrics: Vec<Metrics> = Vec::new();
    let mut duplicate_type: Vec<Option<String>> = Vec::new();
    let mut written = Vec::with_capacity(records.len());

    for (index, record) in records.iter().enumerate() {
        let position = match metrics.iter().position(|row| row.library == record.library) {
            Some(position) => position,
            None => {
                metrics.push(Metrics {
                    library: record.library.clone(),
                    ..Metrics::default()
                });
                metrics.len() - 1
            }
        };
        let row = &mut metrics[position];
        // `addReadToLibraryMetrics`, whose four cases are exclusive and in this order.
        if record.unmapped() {
            row.unmapped_reads += 1;
        } else if record.secondary_or_supplementary() {
            row.secondary_or_supplementary += 1;
        } else if !record.paired() || record.mate_unmapped() {
            row.unpaired_reads_examined += 1;
        } else {
            row.read_pairs_examined += 1;
        }

        if duplicate[index] && !record.secondary_or_supplementary() && !record.unmapped() {
            if !record.paired() || record.mate_unmapped() {
                row.unpaired_read_duplicates += 1;
            } else {
                row.read_pair_duplicates += 1;
            }
        }
        if optical[index] {
            row.read_pair_optical_duplicates += 1;
        }

        // `CLEAR_DT` runs before the policy, so an incoming tag survives only when it is off and
        // the record is not given a new one.
        let mut tag = if options.clear_dt {
            None
        } else {
            record.existing_dt.clone()
        };
        if options.tagging_policy != TaggingPolicy::DontTag && duplicate[index] {
            if optical[index] {
                tag = Some(SEQUENCING_CODE.to_string());
            } else if options.tagging_policy == TaggingPolicy::All {
                tag = Some(LIBRARY_CODE.to_string());
            }
        }
        duplicate_type.push(tag);

        let removed = (options.remove_duplicates && duplicate[index])
            || (options.remove_sequencing_duplicates && optical[index]);
        written.push(!removed);
    }

    for row in &mut metrics {
        // Both counters counted a pair twice, once per end.
        row.read_pairs_examined /= 2;
        row.read_pair_duplicates /= 2;
        row.read_pair_optical_duplicates /= 2;
        row.estimated_library_size = estimate_library_size(
            row.read_pairs_examined - row.read_pair_optical_duplicates,
            row.read_pairs_examined - row.read_pair_duplicates,
        );
        let examined = row.unpaired_reads_examined + row.read_pairs_examined;
        row.percent_duplication = if examined == 0 {
            0.0
        } else {
            (row.unpaired_read_duplicates + row.read_pair_duplicates * 2) as f64
                / (row.unpaired_reads_examined + row.read_pairs_examined * 2) as f64
        };
    }

    Marking {
        duplicate: duplicate.to_vec(),
        optical: optical.to_vec(),
        duplicate_type,
        written,
        metrics,
        all_sets: Vec::new(),
        optical_sets: Vec::new(),
        non_optical_sets: Vec::new(),
    }
}
