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
//!
//! Ported from `picard.sam.markduplicates.MarkDuplicates`,
//! `picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram`,
//! `picard.sam.markduplicates.util.OpticalDuplicateFinder`, `picard.sam.util.ReadNameParser`,
//! `picard.sam.DuplicationMetrics` and `htsjdk.samtools.DuplicateScoringStrategy` in Picard 3.4.0.

use foldhash::{HashMap, HashMapExt};
use std::hash::BuildHasher;
use std::sync::Arc;

use htsjdk_bam::cigar::{Cigar, Op};

/// `DuplicateScoringStrategy.ScoringStrategy`, minus `RANDOM`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringStrategy {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
    /// `RANDOM`, which is not random: the read NAME is hashed with Murmur3 seeded at 1 and the low
    /// fourteen bits become the score, so both ends of a pair score alike and two runs over one
    /// file keep the same reads.
    Random,
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
    /// `READ_NAME_REGEX`, which `parse_read_names` is false for when it is `null`.
    pub read_name_regex: ReadNameRegex,
    /// `MAX_OPTICAL_DUPLICATE_SET_SIZE`: a list handed to the optical finder that is longer than
    /// this gets no optical duplicates at all.
    pub max_optical_duplicate_set_size: i64,
    /// `BARCODE_TAG`, whose value splits one position into one set per barcode.
    pub barcode_tag: Option<String>,
    /// Not an argument of any tool: the two mate-cigar markers score a pair from ONE end, adding
    /// the mate's reference length off the `MC` tag (`DuplicateScoringStrategy.compare(..., true)`
    /// in `MarkQueue`), and `MarkDuplicates` scores each end alone. It decides which read of a set
    /// is kept, so the same file can come out marked differently under the two tools.
    pub assume_mate_cigar: bool,
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
            read_name_regex: ReadNameRegex::Default,
            max_optical_duplicate_set_size: 300_000,
            barcode_tag: None,
            assume_mate_cigar: false,
        }
    }
}

/// What the marking reads of a cigar: the reference length and the clips at either end, which
/// are all the unclipped coordinates and the reference-length score need. Twelve bytes in place of
/// the elements, because one is held per record of the file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Span {
    pub reference_length: u32,
    /// The soft and hard clips before the first element that is neither.
    pub leading_clip: i32,
    /// The soft and hard clips after the last element that is neither.
    pub trailing_clip: i32,
}

impl Span {
    pub fn of(cigar: &Cigar) -> Self {
        let clip = |elements: &mut dyn Iterator<Item = &htsjdk_bam::cigar::CigarElement>| -> i32 {
            let mut clipped = 0i32;
            for element in elements {
                match element.op {
                    Op::S | Op::H => clipped += element.length as i32,
                    _ => break,
                }
            }
            clipped
        };
        Self {
            reference_length: cigar.reference_length(),
            leading_clip: clip(&mut cigar.elements.iter()),
            trailing_clip: clip(&mut cigar.elements.iter().rev()),
        }
    }

    /// `getUnclippedStart`, for a record aligned at `alignment_start` (one-based).
    pub fn unclipped_start(&self, alignment_start: i32) -> i32 {
        alignment_start - self.leading_clip
    }

    /// `getUnclippedEnd`, for a record aligned at `alignment_start` (one-based).
    pub fn unclipped_end(&self, alignment_start: i32) -> i32 {
        alignment_start + self.reference_length as i32 - 1 + self.trailing_clip
    }
}

/// `SUM_OF_BASE_QUALITIES` before the vendor discount: the qualities at or above 15, capped at
/// `Short.MAX_VALUE / 2`. It is the only thing the score reads of the qualities, so a record keeps
/// this and not them.
pub fn base_quality_score(qualities: &[u8]) -> i16 {
    let sum: i64 = qualities
        .iter()
        .filter(|quality| **quality >= 15)
        .map(|quality| i64::from(*quality))
        .sum();
    sum.min(i64::from(i16::MAX / 2)) as i16
}

/// One record, reduced to what the two passes read. The cigar, the qualities and the mate cigar
/// are kept as what is read of them ([`Span`], [`base_quality_score`] and the read length), so a
/// file's worth of these holds no per-record allocation but the name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub name: String,
    pub flags: u16,
    /// `-1` where the record is unmapped.
    pub reference_index: i32,
    /// One-based, as the SAM field is.
    pub alignment_start: i32,
    /// The record's own cigar, as [`Span::of`] reads it.
    pub span: Span,
    /// The number of BASES, which is `SAMRecord.getReadLength` and what `toString` prints. Not the
    /// number of qualities: a read whose `QUAL` is `*` has none and still has its length.
    pub read_length: u32,
    /// [`base_quality_score`] of the record's qualities.
    pub base_quality_score: i16,
    pub mate_reference_index: i32,
    /// Shared between the records of one library, which is most of them.
    pub library: Arc<str>,
    /// The nth read group of the header, which is what `closeEnough` compares.
    pub read_group: i32,
    /// The value of `BARCODE_TAG`, where the record carries one.
    pub barcode: Option<String>,
    /// An incoming `DT`, which `CLEAR_DT` decides the fate of.
    pub existing_dt: Option<String>,
    /// The `MC` tag as [`Span::of`] reads it, which is what the mate-cigar markers read instead of
    /// waiting for the mate. `MarkDuplicates` itself does not look at it.
    pub mate_span: Option<Span>,
    /// The mate's alignment start, one-based, which the molecular identifier of a forward read is
    /// built from. `MarkDuplicates` reads the mate's position off the mate's own record instead.
    pub mate_alignment_start: i32,
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
    /// `SAMRecord.getMateNegativeStrandFlag`, which the mate-cigar path needs to orient a pair
    /// from one end.
    pub fn mate_reverse_strand(&self) -> bool {
        self.flags & 0x20 != 0
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
        self.span.unclipped_start(self.alignment_start)
    }

    /// `SAMRecord.getUnclippedEnd`: the alignment end walked forward over the trailing clips.
    pub fn unclipped_end(&self) -> i32 {
        self.span.unclipped_end(self.alignment_start)
    }

    /// `SAMUtils.getMateUnclippedStart`, off the `MC` tag, where the record carries one.
    pub fn mate_unclipped_start(&self) -> Option<i32> {
        Some(self.mate_span?.unclipped_start(self.mate_alignment_start))
    }

    /// `SAMUtils.getMateUnclippedEnd`, off the `MC` tag, where the record carries one.
    pub fn mate_unclipped_end(&self) -> Option<i32> {
        Some(self.mate_span?.unclipped_end(self.mate_alignment_start))
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
    duplicate_score_with(record, strategy, false)
}

/// `DuplicateScoringStrategy.computeDuplicateScore(record, strategy, assumeMateCigar)`.
///
/// `assume_mate_cigar` is what the two mate-cigar markers pass: under
/// `TOTAL_MAPPED_REFERENCE_LENGTH` a paired read with a mapped mate scores its OWN reference
/// length plus its MATE's, read off the `MC` tag, so a pair is scored as a pair from one end.
/// `MarkDuplicates` passes false and scores each end alone, which is why the same file can pick a
/// different read to keep under the two tools.
pub fn duplicate_score_with(
    record: &Record,
    strategy: ScoringStrategy,
    assume_mate_cigar: bool,
) -> i16 {
    let capped = |value: i64| -> i16 { value.min(i64::from(i16::MAX / 2)) as i16 };
    let mut score: i16 = match strategy {
        ScoringStrategy::SumOfBaseQualities => record.base_quality_score,
        ScoringStrategy::TotalMappedReferenceLength => {
            let own = if record.unmapped() {
                0
            } else {
                capped(i64::from(record.span.reference_length))
            };
            match (
                assume_mate_cigar,
                record.paired() && !record.mate_unmapped(),
            ) {
                (true, true) => own.wrapping_add(match &record.mate_span {
                    Some(span) => capped(i64::from(span.reference_length)),
                    // `SAMUtils.getMateCigar` returns null and the addition is of nothing; the
                    // caller has already refused a pair with no `MC` or dropped it.
                    None => 0,
                }),
                _ => own,
            }
        }
        // `score += hash & 0b11_1111_1111_1111; score -= Short.MIN_VALUE / 4;` -- the shift moves
        // the value off zero so that a vendor-failed read can be discounted twice, once per mate,
        // without underflowing.
        ScoringStrategy::Random => {
            let hashed = (crate::murmur3::hash_unencoded_chars(&record.name, 1)
                & 0b11_1111_1111_1111) as i16;
            hashed.wrapping_sub(i16::MIN / 4)
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

/// `READ_NAME_REGEX`: the optimised default, `null` (no locations at all), or a pattern.
#[derive(Debug, Clone, Default)]
pub enum ReadNameRegex {
    #[default]
    Default,
    Off,
    /// The pattern, anchored at both ends because `Matcher.matches` is.
    Pattern(regex::Regex),
}

impl ReadNameRegex {
    /// The argument as given: absent or the default's own text is the default, `null` is off.
    pub fn from_arg(value: Option<&str>) -> Result<Self, String> {
        match value {
            None
            | Some("<optimized capture of last three ':' separated fields as numeric values>") => {
                Ok(ReadNameRegex::Default)
            }
            Some("null") => Ok(ReadNameRegex::Off),
            Some(pattern) => regex::Regex::new(&format!("^(?:{pattern})$"))
                .map(ReadNameRegex::Pattern)
                .map_err(|e| format!("java.util.regex.PatternSyntaxException: {e}")),
        }
    }

    pub fn is_off(&self) -> bool {
        matches!(self, ReadNameRegex::Off)
    }
}

/// `Integer.parseInt`: an optional sign and decimal digits, nothing else, within `int`.
fn java_parse_int(text: &str) -> Option<i32> {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    text.parse::<i32>().ok()
}

/// `ReadNameParser.readLocationInformation` under any `READ_NAME_REGEX`: groups one to three of
/// a whole-name match are the tile (cast to `short`), x and y.
pub fn location_with(name: &str, regex: &ReadNameRegex) -> Location {
    match regex {
        ReadNameRegex::Default => location(name),
        ReadNameRegex::Off => Location::default(),
        ReadNameRegex::Pattern(pattern) => {
            let Some(groups) = pattern.captures(name) else {
                return Location::default();
            };
            let field = |i: usize| groups.get(i).and_then(|m| java_parse_int(m.as_str()));
            match (field(1), field(2), field(3)) {
                (Some(tile), Some(x), Some(y)) => Location {
                    tile: tile as i16,
                    x,
                    y,
                    known: true,
                },
                _ => Location::default(),
            }
        }
    }
}

/// `ReadEndsForMarkDuplicates`: one record's, or one pair's, position and score.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadEnds {
    pub library: Arc<str>,
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
        score: duplicate_score_with(record, options.scoring, options.assume_mate_cigar),
        location: if options.parse_read_names {
            location_with(&record.name, &options.read_name_regex)
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

/// `OpticalDuplicateFinder.findOpticalDuplicates`: the fast path for a small set, the graph for a
/// larger one.
///
/// The threshold is read off the keeper the CALLER passed, not off the one that survives the
/// location check: a keeper without a location still puts a set of three on the fast path.
///
/// On the fast path the keeper is compared to everyone first, and then every other pair is
/// compared once, with the one already flagged left alone so a chain of two does not flag both.
pub fn optical_duplicates(list: &[ReadEnds], keeper: Option<usize>, distance: i32) -> Vec<bool> {
    let mut flags = vec![false; list.len()];
    if list.len() < 2 {
        return flags;
    }
    let threshold = if keeper.is_none() { 3 } else { 4 };
    // `keeperOrNull`: a keeper with no location is no keeper.
    let keeper = keeper.filter(|index| list[*index].location.known);
    if list.len() >= threshold {
        return optical_duplicates_with_graph(list, keeper, distance);
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

/// `getOpticalDuplicatesFlagWithGraph`: cluster the ends that lie within the distance of each
/// other, then keep one read per cluster and flag the rest.
///
/// The clustering is what the fast path cannot do. Three ends in a row, each within the distance
/// of the next but the outer two further apart than that, are ONE optical cluster here and two
/// unrelated pairs there, so the count differs by one on exactly the shape a pairwise test gets
/// wrong.
///
/// Three details decide the answer.
///
/// The edges are only ever drawn inside a bucket of `(read group, tile)`, and the bucket key is
/// `(read group << 16) + tile` -- an integer, so two different pairs sharing that value share a
/// bucket. Inside one, the comparison drops the read-group and tile tests it no longer needs and
/// compares the two coordinates alone.
///
/// An end with no location joins no bucket and is a cluster of its own, which is why it is never
/// flagged.
///
/// The read kept in a cluster is the one with the SMALLEST x, ties broken by the smallest y --
/// except in the keeper's own cluster, where the keeper is kept whatever its coordinates are. The
/// walk over the ends is in index order, so a smaller read met later displaces the one already
/// held and flags it.
fn optical_duplicates_with_graph(
    list: &[ReadEnds],
    keeper: Option<usize>,
    distance: i32,
) -> Vec<bool> {
    let mut flags = vec![false; list.len()];

    // The graph's clusters, as a union-find over the ends' indices.
    let mut parent: Vec<usize> = (0..list.len()).collect();
    fn find(parent: &mut [usize], index: usize) -> usize {
        let mut root = index;
        while parent[root] != root {
            root = parent[root];
        }
        let mut walk = index;
        while parent[walk] != root {
            let next = parent[walk];
            parent[walk] = root;
            walk = next;
        }
        root
    }

    let mut buckets: Vec<(i32, Vec<usize>)> = Vec::new();
    for (index, end) in list.iter().enumerate() {
        if !end.location.known {
            continue;
        }
        let key = (end.read_group << 16) + i32::from(end.location.tile);
        match buckets.iter_mut().find(|(existing, _)| *existing == key) {
            Some((_, members)) => members.push(index),
            None => buckets.push((key, vec![index])),
        }
    }
    for (_, members) in &buckets {
        for (position, left) in members.iter().enumerate() {
            for right in &members[position + 1..] {
                // `closeEnoughShort`: the bucket has already settled the read group and the tile.
                if (list[*left].location.x - list[*right].location.x).abs() <= distance
                    && (list[*left].location.y - list[*right].location.y).abs() <= distance
                {
                    let a = find(&mut parent, *left);
                    let b = find(&mut parent, *right);
                    if a != b {
                        parent[a] = b;
                    }
                }
            }
        }
    }
    let cluster: Vec<usize> = (0..list.len())
        .map(|index| find(&mut parent, index))
        .collect();

    // `clusterToRepresentativeRead`, walked in index order. The reference's map is keyed on the
    // end's index, and every key is smaller than the table it lives in, so its iteration is that
    // order too.
    let mut representative: Vec<Option<usize>> = vec![None; list.len()];
    if let Some(keeper) = keeper {
        representative[cluster[keeper]] = Some(keeper);
    }
    let keeper_cluster = keeper.map(|index| cluster[index]);
    for index in 0..list.len() {
        let group = cluster[index];
        match representative[group] {
            Some(held) if Some(index) != keeper => {
                let smaller = list[index].location.x < list[held].location.x
                    || (list[index].location.x == list[held].location.x
                        && list[index].location.y < list[held].location.y);
                if keeper_cluster != Some(group) && smaller {
                    flags[held] = true;
                    representative[group] = Some(index);
                } else {
                    flags[index] = true;
                }
            }
            Some(_) => {}
            None => representative[group] = Some(index),
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
    /// `TAG_DUPLICATE_SET_MEMBERS`: for each record in a duplicate set of pairs, the index in the
    /// file of the set's representative read (its best end's first read) and the set's size,
    /// which the `DI` and `DS` tags carry.
    pub duplicate_set: Vec<Option<(usize, usize)>>,
}

/// `Histogram.increment`, which adds to the bin or creates it.
pub(crate) fn increment(histogram: &mut Vec<(f64, f64)>, bin: f64, by: f64) {
    match histogram.iter_mut().find(|(id, _)| *id == bin) {
        Some((_, count)) => *count += by,
        None => histogram.push((bin, by)),
    }
}

/// The key a chunk is cut on: `areComparableForDuplicates`. It borrows its strings, because
/// there is one per read end.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Key<'a> {
    library: &'a str,
    barcode: Option<&'a str>,
    read1_reference_index: i32,
    read1_coordinate: i32,
    orientation: u8,
    read2: Option<(i32, i32)>,
}

fn key(ends: &ReadEnds, compare_read2: bool) -> Key<'_> {
    key_with(&ends.library, ends.barcode.as_deref(), ends, compare_read2)
}

/// [`key`], with the strings borrowed from somewhere that outlives `ends`.
fn key_with<'a>(
    library: &'a str,
    barcode: Option<&'a str>,
    ends: &ReadEnds,
    compare_read2: bool,
) -> Key<'a> {
    Key {
        library,
        barcode,
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
    mark_sorted(records, options, false)
}

/// [`mark`], with the header's sort order. Under `queryname` every read end takes the index of the
/// FIRST record of its run of one read name rather than its own, and the second pass flags the
/// whole run (see [`queryname_pass`]).
pub fn mark_sorted(records: &[Record], options: &Options, queryname: bool) -> Marking {
    // Sized once: at most one fragment per record and one pair per two.
    let mut pairs: Vec<ReadEnds> = Vec::with_capacity(records.len() / 2);
    let mut fragments: Vec<ReadEnds> = Vec::with_capacity(records.len());
    // `tmp`, the map from a read group and a name to the first end of a pair seen. A key is in it
    // at most once: its second record takes it out. The end itself is the fragment already pushed,
    // so the map holds its place in `fragments`: the table is sized for half the file and touched
    // all over, and a whole end per slot was most of its memory.
    let mut pending: HashMap<(i32, &str), u32> = HashMap::with_capacity(records.len() / 2);
    // `duplicateIndex`, the index of the first record of the current run of one name.
    let mut run_start = 0;

    for (record_index, record) in records.iter().enumerate() {
        if queryname && (record_index == 0 || record.name != records[record_index - 1].name) {
            run_start = record_index;
        }
        // `indexForRead`.
        let index = if queryname { run_start } else { record_index };
        if record.unmapped() || record.secondary_or_supplementary() {
            continue;
        }
        let fragment = build_read_ends(record, index, options);
        if !record.paired() || record.mate_unmapped() {
            fragments.push(fragment);
            continue;
        }
        // The read group's id and the name: one name in two groups is two templates.
        let name = (record.read_group, record.name.as_str());
        match pending.remove(&name) {
            None => {
                pending.insert(name, fragments.len() as u32);
                fragments.push(fragment);
            }
            Some(first) => {
                let mut paired = fragments[first as usize].clone();
                let mates_reference_index = fragment.read1_reference_index;
                let mates_coordinate = fragment.read1_coordinate;
                fragments.push(fragment);
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

    // `duplicateIndexes`, `opticalDuplicateIndexes` and `representativeReadIndicesForDuplicates`:
    // under queryname both ends of a pair share an index, which then goes in once, and two runs
    // can still put one index in twice.
    let mut duplicate_indexes: Vec<usize> = Vec::new();
    let mut optical_indexes: Vec<usize> = Vec::new();
    let mut representatives: Vec<(usize, usize, usize)> = Vec::new();
    // `opticalDuplicatesByLibraryId`, counted per end in a set, not per record written.
    let mut optical_pairs: Vec<(String, i64)> = Vec::new();
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
        // `addRepresentativeReadIndex`, whose best is chosen by the same rule as the keeper.
        let representative = chunk[keeper].read1_index_in_file;
        for ends in &chunk {
            representatives.push((ends.read1_index_in_file, representative, chunk.len()));
            if ends.read2_index_in_file != ends.read1_index_in_file {
                representatives.push((ends.read2_index_in_file, representative, chunk.len()));
            }
        }
        let flags = if options.parse_read_names {
            optical_flags(&chunk, keeper, options)
        } else {
            vec![false; chunk.len()]
        };
        for (position, ends) in chunk.iter().enumerate() {
            if position == keeper {
                continue;
            }
            duplicate_indexes.push(ends.read1_index_in_file);
            if ends.read2_index_in_file != ends.read1_index_in_file {
                duplicate_indexes.push(ends.read2_index_in_file);
            }
            if flags[position] {
                optical_indexes.push(ends.read1_index_in_file);
                if ends.read2_index_in_file != ends.read1_index_in_file {
                    optical_indexes.push(ends.read2_index_in_file);
                }
            }
        }
        let optical_ends = flags.iter().filter(|flag| **flag).count() as i64;
        if optical_ends > 0 {
            match optical_pairs
                .iter_mut()
                .find(|(library, _)| library.as_str() == &*chunk[0].library)
            {
                Some((_, count)) => *count += optical_ends,
                None => optical_pairs.push((chunk[0].library.to_string(), optical_ends)),
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
                    duplicate_indexes.push(ends.read1_index_in_file);
                }
            }
        } else {
            let keeper = best(&chunk);
            for (position, ends) in chunk.iter().enumerate() {
                if position != keeper {
                    duplicate_indexes.push(ends.read1_index_in_file);
                }
            }
        }
    }

    let (duplicate, optical, duplicate_set) = if queryname {
        queryname_pass(records, duplicate_indexes, optical_indexes, representatives)
    } else {
        let mut duplicate = vec![false; records.len()];
        let mut optical = vec![false; records.len()];
        let mut duplicate_set = vec![None; records.len()];
        duplicate_indexes.iter().for_each(|&i| duplicate[i] = true);
        optical_indexes.iter().for_each(|&i| optical[i] = true);
        for (index, representative, size) in representatives {
            duplicate_set[index] = Some((representative, size));
        }
        (duplicate, optical, duplicate_set)
    };
    let mut marking = write(records, options, &duplicate, &optical, Some(&optical_pairs));
    marking.all_sets = all_sets;
    marking.optical_sets = optical_sets;
    marking.non_optical_sets = non_optical_sets;
    marking.duplicate_set = duplicate_set;
    marking
}

/// Per record: duplicate, optical, and the `DI`/`DS` pair.
type Flags = (Vec<bool>, Vec<bool>, Vec<Option<(usize, usize)>>);

/// The second pass under queryname, which walks the sorted index collections the way the
/// reference does: the next index is drawn only once the record is past it AND its name differs
/// from the last one marked, and a record past the index with that same name is marked too. So
/// the secondary, supplementary and unmapped records of a duplicate's run are flagged, and a run
/// of one name that a merge put beside another run of it is flagged with it.
fn queryname_pass(
    records: &[Record],
    mut duplicate_indexes: Vec<usize>,
    mut optical_indexes: Vec<usize>,
    mut representatives: Vec<(usize, usize, usize)>,
) -> Flags {
    // `NO_SUCH_INDEX`.
    const NONE: i64 = i64::MAX;
    duplicate_indexes.sort_unstable();
    optical_indexes.sort_unstable();
    // A stable sort, as the in-memory `SortingCollection` is, on the read's index alone.
    representatives.sort_by_key(|&(index, _, _)| index);
    let mut duplicates = duplicate_indexes.into_iter().map(|i| i as i64);
    let mut opticals = optical_indexes.into_iter().map(|i| i as i64);
    let mut sets = representatives.into_iter();
    let mut next_duplicate = duplicates.next().unwrap_or(NONE);
    let mut next_optical = opticals.next().unwrap_or(NONE);
    // An empty collection leaves the set index at -1, which no record equals.
    let (mut next_in_set, mut representative, mut size) =
        sets.next().map_or((-1, 0, 0), |(i, r, s)| (i as i64, r, s));
    let mut duplicate_name: Option<&str> = None;
    let mut representative_name: Option<&str> = None;

    let mut duplicate = vec![false; records.len()];
    let mut optical = vec![false; records.len()];
    let mut duplicate_set = vec![None; records.len()];
    for (index, record) in records.iter().enumerate() {
        let i = index as i64;
        let name = Some(record.name.as_str());
        // `nextIndexIfNeeded`.
        if i > next_duplicate && name != duplicate_name {
            next_duplicate = duplicates.next().unwrap_or(NONE);
        }
        duplicate[index] = i == next_duplicate || (i > next_duplicate && name == duplicate_name);
        if i > next_optical && name != duplicate_name {
            next_optical = opticals.next().unwrap_or(NONE);
        }
        optical[index] = i == next_optical || (i > next_optical && name == duplicate_name);
        if i > next_in_set && name != representative_name {
            if let Some((next, r, s)) = sets.next() {
                (next_in_set, representative, size) = (next as i64, r, s);
            }
        }
        let in_set = i == next_in_set || (i > next_in_set && name == representative_name);
        if in_set && !record.secondary_or_supplementary() && !record.unmapped() {
            duplicate_set[index] = Some((representative, size));
            representative_name = name;
        }
        // Moved only after both flags, as the reference does.
        if duplicate[index] {
            duplicate_name = name;
        }
    }
    (duplicate, optical, duplicate_set)
}

/// The records at one position, grouped the way a set is cut: by library, barcode, reference,
/// 5' coordinate and orientation, over the FRAGMENT list, which is one entry per record.
///
/// `MarkDuplicates` itself never needs this list; `UmiAwareMarkDuplicatesWithMateCigar` does,
/// because it splits each set by UMI before anything is marked.
pub fn sets_by_position(records: &[Record], options: &Options) -> Vec<Vec<usize>> {
    let mut keys: HashMap<Key, usize> = HashMap::with_capacity(records.len());
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (index, record) in records.iter().enumerate() {
        if record.unmapped() || record.secondary_or_supplementary() {
            continue;
        }
        let ends = build_read_ends(record, index, options);
        let next = groups.len();
        // The strings `build_read_ends` copied, borrowed from the record instead.
        let barcode = options.barcode_tag.as_ref().and(record.barcode.as_deref());
        let position = *keys
            .entry(key_with(&record.library, barcode, &ends, false))
            .or_insert(next);
        if position == next {
            groups.push(Vec::new());
        }
        groups[position].push(index);
    }
    groups
}

/// The chunks `generateDuplicateIndexes` cuts, which are the runs of an equal key: in the order
/// their keys first appear, each in the list's order.
///
/// The grouping is a group number per end, found in a table that holds each group's first end
/// rather than its key, and a counting sort of the indexes by it, so only the chunk being looked
/// at is ever copied out.
fn chunks(list: &[ReadEnds], compare_read2: bool) -> impl Iterator<Item = Vec<ReadEnds>> + '_ {
    let state = foldhash::fast::RandomState::default();
    let key_of = |index: u32| key(&list[index as usize], compare_read2);
    let mut group_of: Vec<u32> = Vec::with_capacity(list.len());
    let mut sizes: Vec<u32> = Vec::new();
    {
        // (the group's first end, the group)
        let mut table: hashbrown::HashTable<(u32, u32)> =
            hashbrown::HashTable::with_capacity(list.len());
        for (index, ends) in list.iter().enumerate() {
            let key = key(ends, compare_read2);
            let next = sizes.len() as u32;
            let group = match table.entry(
                state.hash_one(&key),
                |&(first, _)| key_of(first) == key,
                |&(first, _)| state.hash_one(key_of(first)),
            ) {
                hashbrown::hash_table::Entry::Occupied(entry) => entry.get().1,
                hashbrown::hash_table::Entry::Vacant(entry) => {
                    entry.insert((index as u32, next));
                    sizes.push(0);
                    next
                }
            };
            sizes[group as usize] += 1;
            group_of.push(group);
        }
    }
    // `bounds[g]..bounds[g + 1]` is group g's run in `order`.
    let mut bounds: Vec<u32> = Vec::with_capacity(sizes.len() + 1);
    let mut next = 0;
    bounds.push(0);
    for size in &sizes {
        next += size;
        bounds.push(next);
    }
    let mut fill = bounds[..sizes.len()].to_vec();
    drop(sizes);
    let mut order = vec![0u32; list.len()];
    for (index, group) in group_of.into_iter().enumerate() {
        order[fill[group as usize] as usize] = index as u32;
        fill[group as usize] += 1;
    }
    (0..bounds.len() - 1).map(move |group| {
        let mut chunk: Vec<ReadEnds> = order[bounds[group] as usize..bounds[group + 1] as usize]
            .iter()
            .map(|&index| list[index as usize].clone())
            .collect();
        // `ReadEndsMDComparator` past the key: the mate's position, then the physical location
        // (-1 where none was parsed), then the indexes in the file. Within one file the location
        // rarely ties apart from the index; two merged files with different name styles do.
        chunk.sort_by_key(|ends| {
            let l = ends.location;
            let location = if l.known {
                (i32::from(l.tile), l.x, l.y)
            } else {
                (-1, -1, -1)
            };
            (
                ends.read2_reference_index,
                ends.read2_coordinate,
                location,
                ends.read1_index_in_file,
                ends.read2_index_in_file,
            )
        });
        chunk
    })
}

/// `trackOpticalDuplicates`, which splits a set that mixes FR and RF before it looks at anything.
///
/// The split is the point: in PCR duplicate detection a set can hold both orientations once the
/// ends are ordered by position, and in OPTICAL duplicate detection two reads whose first ends
/// point opposite ways are not duplicates of one another.
pub(crate) fn optical_flags(chunk: &[ReadEnds], keeper: usize, options: &Options) -> Vec<bool> {
    let has_fr = chunk
        .iter()
        .any(|ends| ends.orientation_for_optical_duplicates == FR);
    let has_rf = chunk
        .iter()
        .any(|ends| ends.orientation_for_optical_duplicates == RF);
    let limited = |list: &[ReadEnds], keeper: Option<usize>| -> Vec<bool> {
        if list.len() as i64 > options.max_optical_duplicate_set_size {
            vec![false; list.len()]
        } else {
            optical_duplicates(list, keeper, options.optical_duplicate_pixel_distance)
        }
    };
    if !(has_fr && has_rf) {
        return limited(chunk, Some(keeper));
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
        let inner = limited(&list, inner_keeper);
        for (position, flag) in positions.iter().zip(inner) {
            flags[*position] = flag;
        }
    }
    flags
}

/// The second pass: the flag, the tag, the removal and the metrics, in the reference's order.
/// The writing pass over decisions someone else made: the flags, the `DT` tags, what is written
/// at all, and the metrics. `MarkDuplicates` reaches it through [`mark`]; the mate-cigar iterator
/// reaches it with its own decisions.
pub fn marking_from(
    records: &[Record],
    options: &Options,
    duplicate: &[bool],
    optical: &[bool],
) -> Marking {
    write(records, options, duplicate, optical, None)
}

/// `optical_pairs` is `MarkDuplicates`' own optical count per library, in pairs; without it the
/// count is the records flagged optical, halved.
fn write(
    records: &[Record],
    options: &Options,
    duplicate: &[bool],
    optical: &[bool],
    optical_pairs: Option<&[(String, i64)]>,
) -> Marking {
    let mut metrics: Vec<Metrics> = Vec::new();
    let mut duplicate_type: Vec<Option<String>> = Vec::new();
    let mut written = Vec::with_capacity(records.len());

    for (index, record) in records.iter().enumerate() {
        let position = match metrics
            .iter()
            .position(|row| *row.library == *record.library)
        {
            Some(position) => position,
            None => {
                metrics.push(Metrics {
                    library: record.library.to_string(),
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
        row.read_pair_optical_duplicates = match optical_pairs {
            Some(pairs) => pairs
                .iter()
                .find(|(library, _)| *library == row.library)
                .map_or(0, |(_, count)| *count),
            None => row.read_pair_optical_duplicates / 2,
        };
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
        duplicate_set: vec![None; records.len()],
    }
}

/// One `DuplicationMetrics` row, in its fields' declaration order.
struct DuplicationRow<'a>(&'a Metrics);

impl htsjdk_metrics::file::MetricBean for DuplicationRow<'_> {
    fn class_name(&self) -> &str {
        "picard.sam.DuplicationMetrics"
    }
    fn columns(&self) -> &[&'static str] {
        &[
            "LIBRARY",
            "UNPAIRED_READS_EXAMINED",
            "READ_PAIRS_EXAMINED",
            "SECONDARY_OR_SUPPLEMENTARY_RDS",
            "UNMAPPED_READS",
            "UNPAIRED_READ_DUPLICATES",
            "READ_PAIR_DUPLICATES",
            "READ_PAIR_OPTICAL_DUPLICATES",
            "PERCENT_DUPLICATION",
            "ESTIMATED_LIBRARY_SIZE",
        ]
    }
    fn values(&self) -> Vec<htsjdk_metrics::file::Value> {
        use htsjdk_metrics::file::Value;
        let m = self.0;
        vec![
            Value::Str(m.library.clone()),
            Value::Long(m.unpaired_reads_examined),
            Value::Long(m.read_pairs_examined),
            Value::Long(m.secondary_or_supplementary),
            Value::Long(m.unmapped_reads),
            Value::Long(m.unpaired_read_duplicates),
            Value::Long(m.read_pair_duplicates),
            Value::Long(m.read_pair_optical_duplicates),
            Value::Double(m.percent_duplication),
            m.estimated_library_size.map_or(Value::Null, Value::Long),
        ]
    }
}

/// `finalizeAndWriteMetrics`: a row per library, the return-on-investment histogram when there is
/// one library, and the three set-size histograms beside it (an empty one is not printed).
pub fn duplication_metrics_text(
    tool: &str,
    metrics: &[Metrics],
    header_libraries: &[String],
    all_sets: &[(f64, f64)],
    optical_sets: &[(f64, f64)],
    non_optical_sets: &[(f64, f64)],
) -> String {
    use htsjdk_metrics::file::{Histogram, MetricsFile};
    // `Histogram<Double>` prints its keys with `Double.toString`.
    let key = |k: f64| {
        if k.fract() == 0.0 && k.abs() < 1e7 {
            format!("{k:.1}")
        } else {
            format!("{k}")
        }
    };
    let histogram = |value_label: &str, bin_label: &str, bins: &[(f64, f64)]| Histogram {
        bin_label: bin_label.to_string(),
        value_label: value_label.to_string(),
        key_class: "java.lang.Double".to_string(),
        bins: bins.iter().map(|&(k, v)| (key(k), v)).collect(),
    };
    let mut file = MetricsFile::new();
    for header in crate::command_line::default_headers(tool) {
        file.add_header(&header);
    }
    // `new LibraryIdGenerator(header)` makes a row for every library the header's read groups
    // name, so a library with no record left is a row of zeroes.
    let mut all: Vec<Metrics> = metrics.to_vec();
    for library in header_libraries {
        if !all.iter().any(|m| &m.library == library) {
            all.push(Metrics {
                library: library.clone(),
                ..Metrics::default()
            });
        }
    }
    let metrics = &all[..];
    // `LibraryIdGenerator` keeps its metrics in a `TreeMap`: the rows come out by library name.
    let mut sorted_metrics: Vec<&Metrics> = metrics.iter().collect();
    sorted_metrics.sort_by(|a, b| a.library.cmp(&b.library));
    for m in sorted_metrics {
        file.add_metric(&DuplicationRow(m));
    }
    if metrics.len() == 1 {
        if let Some(roi) = roi_histogram(&metrics[0]) {
            file.histograms.push(histogram("CoverageMult", "BIN", &roi));
        }
    }
    let sorted = |bins: &[(f64, f64)]| {
        let mut bins = bins.to_vec();
        bins.sort_by(|a, b| a.0.total_cmp(&b.0));
        bins
    };
    let (all, optical, non_optical) = (
        sorted(all_sets),
        sorted(optical_sets),
        sorted(non_optical_sets),
    );
    file.histograms
        .push(histogram("all_sets", "set_size", &all));
    file.histograms
        .push(histogram("optical_sets", "set_size", &optical));
    file.histograms
        .push(histogram("non_optical_sets", "set_size", &non_optical));
    file.write()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn named(name: &str, flags: u16) -> Record {
        Record {
            name: name.to_string(),
            flags,
            reference_index: 0,
            alignment_start: 1,
            span: Span::of(&Default::default()),
            read_length: 0,
            base_quality_score: 0,
            mate_reference_index: -1,
            library: Arc::from(""),
            read_group: 0,
            barcode: None,
            existing_dt: None,
            mate_span: None,
            mate_alignment_start: 0,
        }
    }

    #[test]
    fn queryname_pass_flags_the_run_and_lags_on_a_repeated_index() {
        // A secondary of a duplicate's run is flagged and a later run of the same name is not,
        // but the set tag, whose collection has run dry, still matches that later run by name.
        let records = [
            named("a", 0),
            named("a", 0x100),
            named("b", 0),
            named("a", 0),
        ];
        let (duplicate, _, set) = queryname_pass(&records, vec![0], vec![], vec![(0, 0, 2)]);
        assert_eq!(duplicate, [true, true, false, false]);
        assert_eq!(set, [Some((0, 2)), None, None, Some((0, 2))]);

        // One index in twice holds the walk back a record, so the reference misses index 1.
        let records = [named("a", 0), named("b", 0), named("c", 0)];
        let (duplicate, optical, _) = queryname_pass(&records, vec![1, 0, 0], vec![0], vec![]);
        assert_eq!(duplicate, [true, false, false]);
        assert_eq!(optical, [true, false, false]);
    }
}
