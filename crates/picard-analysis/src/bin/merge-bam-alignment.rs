//! `MergeBamAlignment` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.MergeBamAlignment`, `SamAlignmentMerger`, `AbstractAlignmentMerger`,
//! `MultiHitAlignedReadIterator`, `HitsForInsert` and the four primary-alignment strategies at tag
//! 3.4.0, with the htsjdk 4.2.0 pieces they call: `SamPairUtil`, `CigarUtil`, `SAMUtils`,
//! `SAMRecord.reverseComplement`, `SAMSequenceDictionary.mergeDictionaries`,
//! `OverclippedReadFilter` and `SAMRecordCoordinateComparator`.
//!
//! The merge walks the unmapped BAM, which is in queryname order, beside the aligned reads grouped
//! by name. A read with no alignment goes through as it was (unless `ALIGNED_READS_ONLY`); a read
//! with alignments gets one output record per hit, cloned when there is more than one hit or a
//! supplementary one. Each record keeps the unmapped BAM's bases, qualities and tags and takes the
//! aligner's position, cigar, mapping quality, strand and non-reserved tags.
//!
//! Three things happen where a reader might not look for them:
//!
//! * **The aligned iterator runs one group ahead,** and a group with more than one hit is resolved
//!   by the strategy when it is READ, not when it is merged. The strategies draw from a
//!   `java.util.Random(1)` that lives as long as the run, so the draws happen in read order.
//! * **An `IllegalStateException` restarts the merge.** `SamAlignmentMerger.mergeAlignment` catches
//!   one -- the aligned reads out of queryname order, or the merge finding them behind the unmapped
//!   reads -- sorts the aligned reads by queryname and runs the whole merge again. The strategy and
//!   its generator are not rebuilt, so the second run continues the first run's draws.
//! * **Only coordinate output gets NM, MD and UQ.** They are recomputed against the reference as
//!   the sorted records are written; queryname and unsorted output is written as the records are
//!   merged, with whatever tags the aligner gave them.

use std::cmp::Ordering;
use std::collections::HashMap;

use htsjdk_bam::cigar::{clip_end_of_read, Cigar, CigarElement, Op};
use htsjdk_bam::header::{ProgramRecord, SamHeader, SequenceRecord};
use htsjdk_bam::pair::set_mate_info;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{reverse_complement, reverse_qualities};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::java_hash_map::string_hash_code;
use picard_analysis::metrics_cli::{absolute, read_input, refuse_validation, thrown};
use picard_analysis::set_nm_md_and_uq_tags::{fix_record, Options as NmOptions};

const TOOL: &str = "MergeBamAlignment";

const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const NEGATIVE: u16 = 0x10;
const MATE_NEGATIVE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

const SEED: i64 = 0x5DEECE66D;

/// `java.util.Random`.
struct JavaRandom {
    seed: i64,
}

impl JavaRandom {
    fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ SEED) & ((1 << 48) - 1),
        }
    }

    fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(SEED).wrapping_add(0xB)) & ((1 << 48) - 1);
        (self.seed >> (48 - bits)) as i32
    }

    fn next_int(&mut self, bound: i32) -> i32 {
        let mut r = self.next(31);
        let m = bound - 1;
        if bound & m == 0 {
            return ((i64::from(bound) * i64::from(r)) >> 31) as i32;
        }
        let mut u = r;
        loop {
            r = u % bound;
            if u.wrapping_sub(r).wrapping_add(m) >= 0 {
                return r;
            }
            u = self.next(31);
        }
    }
}

/// An exception the merge can throw, by what catches it.
enum Failure {
    /// `IllegalStateException`, which `SamAlignmentMerger.mergeAlignment` catches once.
    IllegalState(String),
    /// Anything else, which ends the run: the full `class: message`.
    Other(String),
}

fn picard(message: String) -> Failure {
    Failure::Other(format!("picard.PicardException: {message}"))
}

// ---------------------------------------------------------------------------------------------
// Record helpers: the `SAMRecord` getters and setters the merge uses.

fn flag(rec: &BamRecord, bit: u16) -> bool {
    rec.flags & bit != 0
}

fn set_flag(rec: &mut BamRecord, bit: u16, on: bool) {
    if on {
        rec.flags |= bit;
    } else {
        rec.flags &= !bit;
    }
}

fn is_unmapped(rec: &BamRecord) -> bool {
    flag(rec, UNMAPPED)
}

fn is_secondary_or_supplementary(rec: &BamRecord) -> bool {
    flag(rec, SECONDARY) || flag(rec, SUPPLEMENTARY)
}

fn int_tag(rec: &BamRecord, name: &[u8; 2]) -> Option<i64> {
    match rec.tags.get(Tag::new(name)) {
        Some(TagValue::Int(v)) => Some(*v),
        _ => None,
    }
}

fn set_tag(rec: &mut BamRecord, name: &[u8; 2], value: Option<TagValue>) {
    match value {
        Some(value) => rec.tags.insert(Tag::new(name), value),
        None => rec.tags.remove(Tag::new(name)),
    }
}

fn tag_name(tag: &Tag) -> String {
    tag.to_string()
}

/// `getAlignmentStart()` through the end, `getAlignmentEnd()`.
fn alignment_end(rec: &BamRecord) -> i32 {
    rec.alignment_end()
}

fn unclipped_start(rec: &BamRecord) -> i32 {
    let mut start = rec.alignment_start;
    for element in &rec.cigar.elements {
        match element.op {
            Op::S | Op::H => start -= element.length as i32,
            _ => break,
        }
    }
    start
}

fn unclipped_end(rec: &BamRecord) -> i32 {
    let mut end = alignment_end(rec);
    for element in rec.cigar.elements.iter().rev() {
        match element.op {
            Op::S | Op::H => end += element.length as i32,
            _ => break,
        }
    }
    end
}

/// `SAMRecord.toString()`.
fn describe(rec: &BamRecord, header: &SamHeader) -> String {
    let mut out = rec.read_name.clone();
    if flag(rec, PAIRED) {
        out.push_str(if flag(rec, FIRST_OF_PAIR) {
            " 1/2"
        } else {
            " 2/2"
        });
    }
    out.push_str(&format!(" {}b", rec.read_bases.len()));
    if is_unmapped(rec) {
        out.push_str(" unmapped read.");
    } else {
        let contig = header
            .sequences
            .get(rec.reference_index.max(0) as usize)
            .map_or("*", |s| s.name.as_str());
        out.push_str(&format!(
            " aligned to {contig}:{}-{}.",
            rec.alignment_start,
            alignment_end(rec)
        ));
    }
    out
}

/// One `AlignmentBlock`: 1-based read start, 1-based reference start, length.
fn alignment_blocks(cigar: &Cigar, alignment_start: i32) -> Vec<(i32, i32, i32)> {
    let mut blocks = Vec::new();
    let mut read_base = 1;
    let mut ref_base = alignment_start;
    for element in &cigar.elements {
        let length = element.length as i32;
        match element.op {
            Op::H | Op::P => {}
            Op::S | Op::I => read_base += length,
            Op::N | Op::D => ref_base += length,
            Op::M | Op::Eq | Op::X => {
                blocks.push((read_base, ref_base, length));
                read_base += length;
                ref_base += length;
            }
        }
    }
    blocks
}

fn cigar_maps_no_bases_to_ref(cigar: &Cigar) -> bool {
    !cigar
        .elements
        .iter()
        .any(|e| e.op.consumes_read_bases() && e.op.consumes_reference_bases())
}

/// `SAMUtils.makeReadUnmapped`, whose reverse-complement is in place and uses the DEFAULT tag sets.
fn make_read_unmapped(rec: &mut BamRecord) {
    if flag(rec, NEGATIVE) {
        reverse_complement_record(
            rec,
            &["E2".to_string(), "SQ".to_string()],
            &["OQ".to_string(), "U2".to_string()],
        )
        .ok();
        set_flag(rec, NEGATIVE, false);
    }
    set_flag(rec, DUPLICATE, false);
    rec.reference_index = -1;
    rec.alignment_start = 0;
    rec.cigar = Cigar::new(Vec::new());
    rec.mapping_quality = 0;
    rec.inferred_insert_size = 0;
    set_flag(rec, SECONDARY, false);
    set_flag(rec, SUPPLEMENTARY, false);
    set_flag(rec, PROPER_PAIR, false);
    set_flag(rec, UNMAPPED, true);
}

/// `SAMRecord.reverseComplement(tagsToRevcomp, tagsToReverse, inplace)`.
fn reverse_complement_record(
    rec: &mut BamRecord,
    to_reverse_complement: &[String],
    to_reverse: &[String],
) -> Result<(), Failure> {
    reverse_complement(&mut rec.read_bases);
    reverse_qualities(&mut rec.base_qualities);
    for name in to_reverse_complement {
        let tag = Tag::new(&two(name));
        let value = match rec.tags.get(tag) {
            None => continue,
            Some(TagValue::Str(s)) => {
                let mut bytes = s.clone().into_bytes();
                reverse_complement(&mut bytes);
                TagValue::Str(String::from_utf8_lossy(&bytes).into_owned())
            }
            Some(TagValue::ByteArray { values, unsigned }) => {
                let mut bytes: Vec<u8> = values.iter().map(|b| *b as u8).collect();
                reverse_complement(&mut bytes);
                TagValue::ByteArray {
                    values: bytes.into_iter().map(|b| b as i8).collect(),
                    unsigned: *unsigned,
                }
            }
            Some(other) => {
                return Err(Failure::Other(format!(
                    "java.lang.UnsupportedOperationException: Don't know how to reverse \
                     complement: {other:?}"
                )))
            }
        };
        rec.tags.insert(tag, value);
    }
    for name in to_reverse {
        let tag = Tag::new(&two(name));
        let value = match rec.tags.get(tag) {
            None => continue,
            Some(TagValue::Str(s)) => TagValue::Str(s.chars().rev().collect()),
            Some(TagValue::ByteArray { values, unsigned }) => TagValue::ByteArray {
                values: values.iter().rev().copied().collect(),
                unsigned: *unsigned,
            },
            Some(TagValue::ShortArray { values, unsigned }) => TagValue::ShortArray {
                values: values.iter().rev().copied().collect(),
                unsigned: *unsigned,
            },
            Some(TagValue::IntArray { values, unsigned }) => TagValue::IntArray {
                values: values.iter().rev().copied().collect(),
                unsigned: *unsigned,
            },
            Some(TagValue::FloatArray(values)) => {
                TagValue::FloatArray(values.iter().rev().copied().collect())
            }
            Some(other) => {
                return Err(Failure::Other(format!(
                    "java.lang.UnsupportedOperationException: Don't know how to reverse: \
                     {other:?}"
                )))
            }
        };
        rec.tags.insert(tag, value);
    }
    Ok(())
}

/// A two-letter tag name as bytes; anything else is padded, which no valid tag needs.
fn two(name: &str) -> [u8; 2] {
    let bytes = name.as_bytes();
    [
        bytes.first().copied().unwrap_or(b' '),
        bytes.get(1).copied().unwrap_or(b' '),
    ]
}

// ---------------------------------------------------------------------------------------------
// CigarUtil.

/// `Cigar.isValid`: whether it has any validation error.
fn cigar_is_valid(cigar: &Cigar) -> bool {
    let elements = &cigar.elements;
    if elements.is_empty() {
        return true;
    }
    let real = |op: Op| matches!(op, Op::M | Op::I | Op::D | Op::N | Op::Eq | Op::X);
    let indel = |op: Op| matches!(op, Op::I | Op::D);
    let mut seen_real = false;
    let mut valid = true;
    for (i, element) in elements.iter().enumerate() {
        if element.length == 0 {
            valid = false;
        }
        let op = element.op;
        let last = elements.len() - 1;
        match op {
            Op::H => {
                if i != 0 && i != last {
                    valid = false;
                }
            }
            Op::S => {
                if i == 0 || i == last {
                } else if i == 1 {
                    if elements.len() == 3 && elements[2].op == Op::H {
                    } else if elements[0].op != Op::H {
                        valid = false;
                    }
                } else if i == last - 1 {
                    if elements[last].op != Op::H {
                        valid = false;
                    }
                } else {
                    valid = false;
                }
            }
            Op::P => {
                if i == 0 {
                } else if i == last || !real(elements[i - 1].op) || !real(elements[i + 1].op) {
                    valid = false;
                }
            }
            _ => {
                seen_real = true;
                if indel(op) {
                    for next in &elements[i + 1..] {
                        if (real(next.op) && !indel(next.op)) || next.op == Op::P {
                            break;
                        }
                        if indel(next.op) && next.op == op {
                            valid = false;
                        }
                    }
                }
            }
        }
    }
    valid && seen_real
}

/// `CigarUtil.isValidCigar`.
fn is_valid_cigar(rec: &BamRecord, cigar: &Cigar) -> bool {
    if cigar.elements.is_empty() {
        return false;
    }
    if !cigar_is_valid(cigar) {
        return false;
    }
    rec.read_bases.len() == cigar.read_length() as usize
}

/// `CigarUtil.clip3PrimeEndOfRead(rec, clipFrom, operator)`.
fn clip_3prime_end_of_read(rec: &mut BamRecord, clip_from: i32, hard: bool) -> Result<(), Failure> {
    let cigar = rec.cigar.clone();
    let negative = flag(rec, NEGATIVE);
    if !is_valid_cigar(rec, &cigar) {
        return Ok(());
    }
    let original_read_length = rec.read_bases.len() as i32;
    let original_reference_length = cigar.reference_length() as i32;
    let mut old = cigar.elements.clone();
    if negative {
        old.reverse();
    }
    let operator = if hard { Op::H } else { Op::S };
    let mut new_elements = clip_end_of_read(clip_from, &old, operator);
    if negative {
        new_elements.reverse();
    }
    let new_cigar = Cigar::new(new_elements);
    if negative {
        let size_change = original_reference_length - new_cigar.reference_length() as i32;
        if size_change > 0 {
            rec.alignment_start += size_change;
        } else if size_change < 0 {
            return Err(Failure::Other(format!(
                "htsjdk.samtools.SAMException: The clipped length {} is longer than the old \
                 unclipped length {original_reference_length}",
                new_cigar.reference_length()
            )));
        }
    }
    rec.cigar = new_cigar.clone();
    if hard {
        let keep = |v: &Vec<u8>| -> Vec<u8> {
            if negative {
                v[(original_read_length - clip_from + 1) as usize..].to_vec()
            } else {
                v[..(clip_from - 1) as usize].to_vec()
            }
        };
        rec.read_bases = keep(&rec.read_bases);
        rec.base_qualities = keep(&rec.base_qualities);
    }
    let has_mapped_bases = new_cigar
        .elements
        .iter()
        .any(|e| e.op.consumes_reference_bases() && e.op.consumes_read_bases());
    if new_cigar.reference_length() as i32 != original_reference_length {
        for name in [b"NM", b"MD", b"UQ"] {
            rec.tags.remove(Tag::new(name));
        }
    }
    if !has_mapped_bases {
        set_flag(rec, UNMAPPED, true);
        rec.cigar = Cigar::new(Vec::new());
        rec.reference_index = -1;
        rec.alignment_start = 0;
        rec.mapping_quality = 0;
        rec.inferred_insert_size = 0;
    } else if !is_valid_cigar(rec, &new_cigar) {
        return Err(Failure::IllegalState(format!(
            "Invalid new Cigar: {} ({}) for {}",
            new_cigar.to_text(),
            cigar_list(&old),
            rec.read_name
        )));
    }
    Ok(())
}

/// `List<CigarElement>.toString()`.
fn cigar_list(elements: &[CigarElement]) -> String {
    let parts: Vec<String> = elements
        .iter()
        .map(|e| format!("{}{}", e.length, e.op.to_char() as char))
        .collect();
    format!("[{}]", parts.join(", "))
}

/// `CigarUtil.addSoftClippedBasesToEndsOfCigar`.
fn add_soft_clipped_bases_to_ends_of_cigar(
    cigar: &Cigar,
    negative: bool,
    three_prime: i32,
    five_prime: i32,
) -> Cigar {
    let mut elements = cigar.elements.clone();
    if negative {
        elements.reverse();
    }
    if three_prime > 0 {
        let mut bases = three_prime;
        if elements.last().is_some_and(|e| e.op == Op::S) {
            bases += elements.pop().map_or(0, |e| e.length as i32);
        }
        elements.push(CigarElement {
            length: bases as u32,
            op: Op::S,
        });
    }
    if five_prime > 0 {
        let mut bases = five_prime;
        if elements
            .first()
            .is_some_and(|e| matches!(e.op, Op::S | Op::H))
        {
            bases += elements.remove(0).length as i32;
        }
        elements.insert(
            0,
            CigarElement {
                length: bases as u32,
                op: Op::S,
            },
        );
    }
    if negative {
        elements.reverse();
    }
    Cigar::new(elements)
}

// ---------------------------------------------------------------------------------------------
// SamPairUtil.

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Orientation {
    Fr,
    Rf,
    Tandem,
}

fn pair_orientation(rec: &BamRecord) -> Orientation {
    let reverse = flag(rec, NEGATIVE);
    if reverse == flag(rec, MATE_NEGATIVE) {
        return Orientation::Tandem;
    }
    let positive_five_prime = if reverse {
        i64::from(rec.mate_alignment_start)
    } else {
        i64::from(rec.alignment_start)
    };
    let negative_five_prime = if reverse {
        i64::from(alignment_end(rec))
    } else {
        i64::from(rec.alignment_start) + i64::from(rec.inferred_insert_size)
    };
    if positive_five_prime < negative_five_prime {
        Orientation::Fr
    } else {
        Orientation::Rf
    }
}

fn is_proper_pair(first: &BamRecord, second: &BamRecord, expected: &[Orientation]) -> bool {
    if is_unmapped(first) || is_unmapped(second) {
        return false;
    }
    if first.reference_index < 0 || first.reference_index != second.reference_index {
        return false;
    }
    expected.contains(&pair_orientation(first))
}

fn set_proper_pair_flags(rec1: &mut BamRecord, rec2: &mut BamRecord, expected: &[Orientation]) {
    let proper = !is_unmapped(rec1) && !is_unmapped(rec2) && is_proper_pair(rec1, rec2, expected);
    set_flag(rec1, PROPER_PAIR, proper);
    set_flag(rec2, PROPER_PAIR, proper);
}

fn set_mate_information_on_supplemental(
    supplemental: &mut BamRecord,
    mate: &BamRecord,
    set_mate_cigar: bool,
) {
    supplemental.mate_reference_index = mate.reference_index;
    supplemental.mate_alignment_start = mate.alignment_start;
    set_flag(supplemental, MATE_NEGATIVE, flag(mate, NEGATIVE));
    set_flag(supplemental, MATE_UNMAPPED, is_unmapped(mate));
    supplemental.inferred_insert_size = -mate.inferred_insert_size;
    if set_mate_cigar && !is_unmapped(mate) {
        set_tag(
            supplemental,
            b"MC",
            Some(TagValue::Str(mate.cigar.to_text())),
        );
    } else {
        set_tag(supplemental, b"MC", None);
    }
    set_tag(
        supplemental,
        b"MQ",
        Some(TagValue::Int(i64::from(mate.mapping_quality))),
    );
}

// ---------------------------------------------------------------------------------------------
// The aligned side: HitsForInsert, the strategies, MultiHitAlignedReadIterator.

#[derive(Clone, Copy, PartialEq, Eq)]
enum Strategy {
    BestMapq,
    EarliestFragment,
    BestEndMapq,
    MostDistant,
}

#[derive(Default, Clone)]
struct Hits {
    first: Vec<Option<BamRecord>>,
    second: Vec<Option<BamRecord>>,
    supplemental_first: Vec<BamRecord>,
    supplemental_second: Vec<BamRecord>,
}

impl Hits {
    fn num_hits(&self) -> usize {
        self.first.len().max(self.second.len())
    }

    fn first_of_pair(&self, i: usize) -> Option<&BamRecord> {
        self.first.get(i).and_then(|r| r.as_ref())
    }

    fn second_of_pair(&self, i: usize) -> Option<&BamRecord> {
        self.second.get(i).and_then(|r| r.as_ref())
    }

    fn representative(&self) -> &BamRecord {
        self.first
            .iter()
            .flatten()
            .next()
            .or_else(|| self.second.iter().flatten().next())
            .expect("a group with hits has a read")
    }

    fn read_name(&self) -> &str {
        &self.representative().read_name
    }

    fn has_supplemental_hits(&self) -> bool {
        !(self.supplemental_first.is_empty() && self.supplemental_second.is_empty())
    }

    fn set_primary_alignment(&mut self, primary: usize) {
        for i in 0..self.num_hits() {
            let not_primary = i != primary;
            if let Some(Some(r)) = self.first.get_mut(i) {
                set_flag(r, SECONDARY, not_primary);
            }
            if let Some(Some(r)) = self.second.get_mut(i) {
                set_flag(r, SECONDARY, not_primary);
            }
        }
    }

    fn index_of_earliest_primary(&self) -> Option<usize> {
        (0..self.num_hits()).find(|&i| {
            self.first_of_pair(i)
                .is_some_and(|r| !is_secondary_or_supplementary(r))
                || self
                    .second_of_pair(i)
                    .is_some_and(|r| !is_secondary_or_supplementary(r))
        })
    }

    /// `coordinateByHitIndex`.
    fn coordinate_by_hit_index(&mut self) {
        let by_hi = |a: &Option<BamRecord>, b: &Option<BamRecord>| -> Ordering {
            let hi = |r: &Option<BamRecord>| r.as_ref().and_then(|r| int_tag(r, b"HI"));
            match (hi(a), hi(b)) {
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Greater,
                (Some(_), None) => Ordering::Less,
                (Some(x), Some(y)) => x.cmp(&y),
            }
        };
        self.first.sort_by(by_hi);
        self.second.sort_by(by_hi);
        let mut i = 0;
        while i < self.first.len().min(self.second.len()) {
            let left = self.first[i].as_ref().and_then(|r| int_tag(r, b"HI"));
            let right = self.second[i].as_ref().and_then(|r| int_tag(r, b"HI"));
            match (left, right) {
                (Some(l), Some(r)) => {
                    if l < r {
                        self.second.insert(i, None);
                    } else if r < l {
                        self.first.insert(i, None);
                    }
                }
                (None, Some(_)) => self.first.insert(i, None),
                (Some(_), None) => {}
                (None, None) => self.second.insert(i, None),
            }
            i += 1;
        }
        for i in 0..self.num_hits() {
            let both = self.first_of_pair(i).is_some() && self.second_of_pair(i).is_some();
            if both {
                for list in [&mut self.first, &mut self.second] {
                    if let Some(Some(r)) = list.get_mut(i) {
                        set_tag(r, b"HI", Some(TagValue::Int(i as i64)));
                    }
                }
            } else if let Some(Some(r)) = self.first.get_mut(i) {
                set_tag(r, b"HI", None);
            } else if let Some(Some(r)) = self.second.get_mut(i) {
                set_tag(r, b"HI", None);
            }
        }
    }

    fn tally_primary(list: &[Option<BamRecord>]) -> u8 {
        let mut seen = false;
        for r in list.iter().flatten() {
            if !is_secondary_or_supplementary(r) {
                if seen {
                    return 2;
                }
                seen = true;
            }
        }
        u8::from(seen)
    }
}

fn combine_mapqs(m1: i32, m2: i32) -> i32 {
    let m1 = if m1 == 255 { 1 } else { m1 * 100 };
    let m2 = if m2 == 255 { 1 } else { m2 * 100 };
    m1 + m2
}

fn compare_mapqs(mapq1: i32, mapq2: i32) -> i32 {
    if mapq1 == mapq2 {
        0
    } else if mapq1 == 0 {
        -1
    } else if mapq2 == 0 {
        1
    } else if mapq1 == 255 {
        -1
    } else if mapq2 == 255 {
        1
    } else {
        mapq1 - mapq2
    }
}

/// The strategy, and the generator it draws ties from.
struct Picker {
    strategy: Strategy,
    random: JavaRandom,
}

impl Picker {
    fn pick(&mut self, hits: &mut Hits, header: &SamHeader) -> Result<(), Failure> {
        match self.strategy {
            Strategy::BestMapq => self.best_mapq(hits),
            Strategy::EarliestFragment => self.earliest_fragment(hits, header),
            Strategy::BestEndMapq => self.best_end_mapq(hits),
            Strategy::MostDistant => self.most_distant(hits),
        }
        .map(|_| ())
    }

    fn best_mapq(&mut self, hits: &mut Hits) -> Result<(), Failure> {
        hits.coordinate_by_hit_index();
        let first = Hits::tally_primary(&hits.first);
        let second = Hits::tally_primary(&hits.second);
        if (first == 0 && second == 0) || first == 2 || second == 2 {
            let mut indices = Vec::new();
            let mut best = -1;
            for i in 0..hits.num_hits() {
                let a = hits
                    .first_of_pair(i)
                    .map_or(0, |r| i32::from(r.mapping_quality));
                let b = hits
                    .second_of_pair(i)
                    .map_or(0, |r| i32::from(r.mapping_quality));
                let this = combine_mapqs(a, b);
                if this > best {
                    best = this;
                    indices.clear();
                }
                if this == best {
                    indices.push(i);
                }
            }
            let primary = if indices.len() == 1 {
                indices[0]
            } else {
                indices[self.random.next_int(indices.len() as i32) as usize]
            };
            hits.set_primary_alignment(primary);
        }
        Ok(())
    }

    fn earliest_fragment(&mut self, hits: &mut Hits, header: &SamHeader) -> Result<(), Failure> {
        let mut earliest = Vec::new();
        let mut earliest_base = i32::MAX;
        let mut best_mapq = -1;
        for i in 0..hits.num_hits() {
            let rec = hits.first[i].as_ref().expect("a fragment hit");
            if flag(rec, PAIRED) {
                return Err(Failure::Other(format!(
                    "java.lang.UnsupportedOperationException: getFragment called for paired \
                     read: {}",
                    describe(rec, header)
                )));
            }
            if is_unmapped(rec) {
                continue;
            }
            let blocks = alignment_blocks(&rec.cigar, rec.alignment_start);
            let first_base = if flag(rec, NEGATIVE) {
                let last = blocks[blocks.len() - 1];
                rec.read_bases.len() as i32 - (last.0 + last.2 - 1) + 1
            } else {
                blocks[0].0
            };
            let mapq = i32::from(rec.mapping_quality);
            if first_base < earliest_base || (first_base == earliest_base && mapq > best_mapq) {
                earliest.clear();
                earliest.push(i);
                earliest_base = first_base;
                best_mapq = mapq;
            } else if first_base == earliest_base && mapq == best_mapq {
                earliest.push(i);
            }
        }
        let primary = if earliest.len() == 1 {
            earliest[0]
        } else {
            earliest[self.random.next_int(earliest.len() as i32) as usize]
        };
        hits.set_primary_alignment(primary);
        Ok(())
    }

    fn best_end_mapq(&mut self, hits: &mut Hits) -> Result<(), Failure> {
        let order = |a: &Option<BamRecord>, b: &Option<BamRecord>| -> Ordering {
            let (a, b) = (a.as_ref().expect("hit"), b.as_ref().expect("hit"));
            match (is_unmapped(a), is_unmapped(b)) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                _ => (-compare_mapqs(i32::from(a.mapping_quality), i32::from(b.mapping_quality)))
                    .cmp(&0),
            }
        };
        hits.first.sort_by(order);
        hits.second.sort_by(order);
        for list in [&mut hits.first, &mut hits.second] {
            if list.is_empty() {
                continue;
            }
            let best = list[0].as_ref().map_or(0, |r| r.mapping_quality);
            let mut i = 1;
            while i < list.len() && list[i].as_ref().map_or(0, |r| r.mapping_quality) == best {
                i += 1;
            }
            let chosen = self.random.next_int(i as i32) as usize;
            if chosen != 0 {
                list.swap(0, chosen);
            }
        }
        hits.set_primary_alignment(0);
        if !flag(hits.representative(), PAIRED) {
            return Ok(());
        }
        if hits.first.len() <= 1 || hits.second.len() <= 1 {
            return Ok(());
        }
        for _ in 0..hits.first.len() - 1 {
            hits.second.insert(1, None);
        }
        Ok(())
    }

    fn most_distant(&mut self, hits: &mut Hits) -> Result<(), Failure> {
        // Records are identified by their index in their list.
        let mut first_best = (-1, Vec::<usize>::new());
        let mut second_best = (-1, Vec::<usize>::new());
        let consider = |best: &mut (i32, Vec<usize>), index: usize, mapq: i32| {
            if best.0 == -1 {
                best.0 = mapq;
                best.1.push(index);
            } else {
                let cmp = compare_mapqs(best.0, mapq);
                if cmp < 0 {
                    best.0 = mapq;
                    best.1.clear();
                    best.1.push(index);
                } else if cmp == 0 {
                    best.1.push(index);
                }
            }
        };
        let mut by_sequence: HashMap<i32, Vec<usize>> = HashMap::new();
        for (i, rec) in hits.first.iter().enumerate() {
            let rec = rec.as_ref().expect("hit");
            if is_unmapped(rec) {
                return Err(Failure::IllegalState(String::new()));
            }
            consider(&mut first_best, i, i32::from(rec.mapping_quality));
            by_sequence.entry(rec.reference_index).or_default().push(i);
        }
        let mut pair_best: (i32, i32, Vec<(usize, usize)>) = (-1, -1, Vec::new());
        for (j, second) in hits.second.iter().enumerate() {
            let second = second.as_ref().expect("hit");
            if is_unmapped(second) {
                return Err(Failure::IllegalState(String::new()));
            }
            consider(&mut second_best, j, i32::from(second.mapping_quality));
            if let Some(firsts) = by_sequence.get(&second.reference_index) {
                for &i in firsts {
                    let first = hits.first[i].as_ref().expect("hit");
                    let mapq = combine_mapqs(
                        i32::from(first.mapping_quality),
                        i32::from(second.mapping_quality),
                    );
                    let distance = first.alignment_end().max(alignment_end(second))
                        - first.alignment_start.min(second.alignment_start)
                        + 1;
                    if distance > pair_best.0 || (distance == pair_best.0 && mapq > pair_best.1) {
                        pair_best = (distance, mapq, vec![(i, j)]);
                    } else if distance == pair_best.0 && mapq == pair_best.1 {
                        pair_best.2.push((i, j));
                    }
                }
            }
        }
        let (best_first, best_second) = if pair_best.0 != -1 {
            let pick = self.random.next_int(pair_best.2.len() as i32) as usize;
            let (i, j) = pair_best.2[pick];
            (Some(i), Some(j))
        } else {
            let a = if first_best.0 != -1 {
                Some(first_best.1[self.random.next_int(first_best.1.len() as i32) as usize])
            } else {
                None
            };
            let b = if second_best.0 != -1 {
                Some(second_best.1[self.random.next_int(second_best.1.len() as i32) as usize])
            } else {
                None
            };
            (a, b)
        };
        if let Some(i) = best_first {
            let rec = hits.first.remove(i);
            hits.first.insert(0, rec);
        }
        if let Some(j) = best_second {
            let rec = hits.second.remove(j);
            hits.second.insert(0, rec);
        }
        hits.set_primary_alignment(0);
        if hits.first.len() <= 1 || hits.second.len() <= 1 {
            return Ok(());
        }
        for _ in 0..hits.first.len() - 1 {
            hits.second.insert(1, None);
        }
        Ok(())
    }
}

/// `MultiHitAlignedReadIterator`'s `replaceHardWithSoftClips`.
fn replace_hard_with_soft_clips(rec: &mut BamRecord) {
    if is_unmapped(rec) || rec.cigar.elements.is_empty() {
        return;
    }
    let elements = rec.cigar.elements.clone();
    let first = elements[0];
    let last = if elements.len() == 1 {
        None
    } else {
        Some(elements[elements.len() - 1])
    };
    let start_hard = if first.op == Op::H { first.length } else { 0 } as usize;
    let end_hard = match last {
        Some(e) if e.op == Op::H => e.length as usize,
        _ => 0,
    };
    if start_hard + end_hard > 0 {
        let length = rec.read_bases.len() + start_hard + end_hard;
        let mut bases = vec![b'N'; length];
        bases[start_hard..start_hard + rec.read_bases.len()].copy_from_slice(&rec.read_bases);
        let mut quals = vec![2u8; length];
        quals[start_hard..start_hard + rec.base_qualities.len()]
            .copy_from_slice(&rec.base_qualities);
        let mut new_elements = elements.clone();
        if start_hard > 0 {
            new_elements[0] = CigarElement {
                length: first.length,
                op: Op::S,
            };
        }
        if end_hard > 0 {
            let n = new_elements.len();
            new_elements[n - 1] = CigarElement {
                length: last.map_or(0, |e| e.length),
                op: Op::S,
            };
        }
        rec.read_bases = bases;
        rec.base_qualities = quals;
        rec.cigar = Cigar::new(new_elements);
    }
}

/// `MultiHitAlignedReadIterator`, over the aligned reads that survive both filters.
struct AlignedIterator<'a> {
    records: Vec<BamRecord>,
    at: usize,
    the_next: Option<Hits>,
    picker: &'a mut Picker,
}

impl<'a> AlignedIterator<'a> {
    fn new(
        records: Vec<BamRecord>,
        picker: &'a mut Picker,
        header: &SamHeader,
    ) -> Result<Self, Failure> {
        let mut iterator = AlignedIterator {
            records,
            at: 0,
            the_next: None,
            picker,
        };
        iterator.advance(header)?;
        Ok(iterator)
    }

    fn has_next(&self) -> bool {
        self.the_next.is_some()
    }

    fn next(&mut self, header: &SamHeader) -> Result<Option<Hits>, Failure> {
        let current = self.the_next.take();
        if current.is_some() {
            self.advance(header)?;
        }
        Ok(current)
    }

    fn advance(&mut self, header: &SamHeader) -> Result<(), Failure> {
        while self.at < self.records.len() {
            let hits = self.next_maybe_empty(header)?;
            if hits.num_hits() > 0 {
                self.the_next = Some(hits);
                return Ok(());
            }
        }
        self.the_next = None;
        Ok(())
    }

    fn next_maybe_empty(&mut self, header: &SamHeader) -> Result<Hits, Failure> {
        let name = self.records[self.at].read_name.clone();
        let mut hits = Hits::default();
        let mut is_paired: Option<bool> = None;
        loop {
            let mut rec = self.records[self.at].clone();
            self.at += 1;
            replace_hard_with_soft_clips(&mut rec);
            if let Some(peek) = self.records.get(self.at) {
                if rec.read_name.as_str() > peek.read_name.as_str() {
                    return Err(Failure::IllegalState(
                        "Underlying iterator is not queryname sorted".to_string(),
                    ));
                }
            }
            let paired = flag(&rec, PAIRED);
            match is_paired {
                None => is_paired = Some(paired),
                Some(p) if p != paired => {
                    return Err(picard(format!(
                        "Got a mix of paired and unpaired alignments for read {name}"
                    )))
                }
                _ => {}
            }
            if !paired || flag(&rec, FIRST_OF_PAIR) {
                if flag(&rec, SUPPLEMENTARY) {
                    hits.supplemental_first.push(rec);
                } else {
                    hits.first.push(Some(rec));
                }
            } else if flag(&rec, SECOND_OF_PAIR) {
                if flag(&rec, SUPPLEMENTARY) {
                    hits.supplemental_second.push(rec);
                } else {
                    hits.second.push(Some(rec));
                }
            } else {
                return Err(picard(format!(
                    "Read is marked as pair but neither first or second: {name}"
                )));
            }
            match self.records.get(self.at) {
                Some(peek) if peek.read_name == name => {}
                _ => break,
            }
        }
        if hits.num_hits() <= 1 {
            for list in [&mut hits.first, &mut hits.second] {
                if let Some(Some(r)) = list.get_mut(0) {
                    set_tag(r, b"HI", None);
                    set_flag(r, SECONDARY, false);
                }
            }
        } else {
            self.picker.pick(&mut hits, header)?;
        }
        Ok(hits)
    }
}

// ---------------------------------------------------------------------------------------------
// The merge.

#[derive(Clone, Copy, PartialEq, Eq)]
enum UnmappingStrategy {
    CopyToTag,
    DoNotChange,
    DoNotChangeInvalid,
    MoveToTag,
}

impl UnmappingStrategy {
    fn reset_mapping_information(self) -> bool {
        self == UnmappingStrategy::MoveToTag
    }
    fn populate_oa_tag(self) -> bool {
        matches!(
            self,
            UnmappingStrategy::CopyToTag | UnmappingStrategy::MoveToTag
        )
    }
    fn keep_valid(self) -> bool {
        self != UnmappingStrategy::DoNotChangeInvalid
    }
}

struct Settings {
    clip_adapters: bool,
    bisulfite: bool,
    aligned_reads_only: bool,
    max_gaps: i32,
    retain: Vec<String>,
    remove: Vec<String>,
    reverse: Vec<String>,
    reverse_complement: Vec<String>,
    read1_trim: i32,
    read2_trim: i32,
    expected: Vec<Orientation>,
    sort_order: String,
    add_mate_cigar: bool,
    unmap_contaminants: bool,
    min_unclipped_bases: i32,
    unmapping: UnmappingStrategy,
    clip_overlapping: bool,
    hard_clip_overlapping: bool,
    keep_aligner_proper_pair: bool,
    include_secondary: bool,
    add_pg_tag: bool,
}

/// `isReservedTag`.
fn is_reserved_tag(tag: &str) -> bool {
    let first = tag.chars().next().unwrap_or(' ');
    first.is_lowercase() || matches!(first, 'X' | 'Y' | 'Z')
}

/// `OverclippedReadFilter.filterOut` with `filterSingleEndClips = false`.
fn overclipped(rec: &BamRecord, threshold: i32) -> bool {
    let mut aligned = 0;
    let mut soft_clip_blocks = 0;
    let mut last: Option<Op> = None;
    for element in &rec.cigar.elements {
        if element.op == Op::S {
            if last != Some(Op::S) {
                soft_clip_blocks += 1;
            }
        } else if element.op.consumes_read_bases() {
            aligned += element.length as i32;
        }
        last = Some(element.op);
    }
    aligned < threshold && soft_clip_blocks >= 2
}

struct Merger<'s> {
    settings: &'s Settings,
    header: SamHeader,
    program_id: Option<String>,
}

impl Merger<'_> {
    /// `setValuesFromAlignment`.
    fn set_values_from_alignment(
        &self,
        rec: &mut BamRecord,
        alignment: &BamRecord,
        aligned_header: &SamHeader,
    ) -> Result<(), Failure> {
        if !is_unmapped(rec) {
            return Err(picard(
                "UNMAPPED_BAM contains mapped reads.  If you would like to use this file as the \
                 UNMAPPED_BAM, first revert it using RevertSam."
                    .to_string(),
            ));
        }
        for (tag, value) in alignment.tags.iter() {
            let name = tag_name(tag);
            if (!is_reserved_tag(&name) || self.settings.retain.contains(&name))
                && !self.settings.remove.contains(&name)
            {
                rec.tags.insert(*tag, value.clone());
            }
        }
        set_flag(rec, UNMAPPED, is_unmapped(alignment));
        rec.reference_index = if alignment.reference_index < 0 {
            -1
        } else {
            let name = &aligned_header.sequences[alignment.reference_index as usize].name;
            self.header
                .sequences
                .iter()
                .position(|s| &s.name == name)
                .map_or(-1, |i| i as i32)
        };
        rec.alignment_start = alignment.alignment_start;
        set_flag(rec, NEGATIVE, flag(alignment, NEGATIVE));
        set_flag(rec, SECONDARY, flag(alignment, SECONDARY));
        set_flag(rec, SUPPLEMENTARY, flag(alignment, SUPPLEMENTARY));
        if !is_unmapped(alignment) {
            rec.cigar = alignment.cigar.clone();
            rec.mapping_quality = alignment.mapping_quality;
        }
        if flag(rec, PAIRED) {
            set_flag(rec, PROPER_PAIR, flag(alignment, PROPER_PAIR));
        }
        if flag(rec, NEGATIVE) {
            reverse_complement_record(
                rec,
                &self.settings.reverse_complement,
                &self.settings.reverse,
            )?;
        }
        Ok(())
    }

    /// `createNewCigarIfMapsOffEndOfReference`.
    fn clipped_off_end(
        &self,
        reference_index: i32,
        alignment_end: i32,
        read_length: i32,
        cigar: &Cigar,
    ) -> Option<Cigar> {
        let sequence = self.header.sequences.get(reference_index.max(0) as usize)?;
        let overhang = alignment_end - sequence.length;
        if overhang <= 0 {
            return None;
        }
        let mut clip_from = read_length - overhang + 1;
        if let Some(last) = cigar.elements.last() {
            if last.op == Op::S {
                clip_from -= last.length as i32;
            }
        }
        Some(Cigar::new(clip_end_of_read(
            clip_from,
            &cigar.elements,
            Op::S,
        )))
    }

    /// `updateCigarForTrimmedOrClippedBases`.
    fn update_cigar_for_trimmed_or_clipped_bases(
        &self,
        rec: &mut BamRecord,
        alignment: &BamRecord,
    ) -> Result<(), Failure> {
        let alignment_read_length = alignment.read_bases.len() as i32;
        let original_read_length = rec.read_bases.len() as i32;
        let trimmed = if !flag(rec, PAIRED) || flag(rec, FIRST_OF_PAIR) {
            self.settings.read1_trim
        } else {
            self.settings.read2_trim
        };
        let not_written = original_read_length - (alignment_read_length + trimmed);
        // `createNewCigarsIfMapsOffEndOfReference`.
        if !is_unmapped(rec) {
            if let Some(cigar) = self.clipped_off_end(
                rec.reference_index,
                alignment_end(rec),
                rec.read_bases.len() as i32,
                &rec.cigar,
            ) {
                rec.cigar = cigar;
            }
        }
        if flag(rec, PAIRED) && !flag(rec, MATE_UNMAPPED) {
            if let Some(TagValue::Str(mate_cigar)) = rec.tags.get(Tag::new(b"MC")).cloned() {
                if let Ok(parsed) = htsjdk_bam::text_parse::parse_cigar(&mate_cigar) {
                    let mate_end = rec.mate_alignment_start + parsed.reference_length() as i32 - 1;
                    if let Some(clipped) = self.clipped_off_end(
                        rec.mate_reference_index,
                        mate_end,
                        parsed.read_length() as i32,
                        &parsed,
                    ) {
                        set_tag(rec, b"MC", Some(TagValue::Str(clipped.to_text())));
                    }
                }
            }
        }
        rec.cigar = add_soft_clipped_bases_to_ends_of_cigar(
            &rec.cigar,
            flag(rec, NEGATIVE),
            not_written,
            trimmed,
        );
        if self.settings.clip_adapters {
            if let Some(xt) = int_tag(rec, b"XT") {
                clip_3prime_end_of_read(rec, xt as i32, false)?;
            }
        }
        Ok(())
    }

    /// `transferAlignmentInfoToFragment`.
    fn transfer_to_fragment(
        &self,
        unaligned: &mut BamRecord,
        aligned: &BamRecord,
        aligned_header: &SamHeader,
        contaminant: bool,
    ) -> Result<(), Failure> {
        self.set_values_from_alignment(unaligned, aligned, aligned_header)?;
        self.update_cigar_for_trimmed_or_clipped_bases(unaligned, aligned)?;
        let beyond_end = || {
            let length = aligned_header
                .sequences
                .get(aligned.reference_index.max(0) as usize)
                .map_or(0, |s| s.length);
            length < aligned.alignment_start
        };
        if cigar_maps_no_bases_to_ref(&unaligned.cigar) || beyond_end() {
            make_read_unmapped(unaligned);
        } else if contaminant {
            let strategy = self.settings.unmapping;
            if strategy.populate_oa_tag() {
                let contig = aligned_header
                    .sequences
                    .get(aligned.reference_index.max(0) as usize)
                    .map_or("*".to_string(), |s| s.name.clone());
                let nm = int_tag(aligned, b"NM").map_or(String::new(), |v| v.to_string());
                let oa = format!(
                    "{contig},{},{},{},{nm};",
                    aligned.alignment_start,
                    if aligned.cigar.elements.is_empty() {
                        "*".to_string()
                    } else {
                        aligned.cigar.to_text()
                    },
                    aligned.mapping_quality
                );
                set_tag(unaligned, b"OA", Some(TagValue::Str(oa)));
            }
            if strategy.reset_mapping_information() {
                unaligned.reference_index = -1;
                unaligned.alignment_start = 0;
                set_tag(unaligned, b"NM", None);
            }
            set_flag(unaligned, UNMAPPED, true);
            if strategy.keep_valid() {
                unaligned.mapping_quality = 0;
                unaligned.cigar = Cigar::new(Vec::new());
            }
            let comment = match unaligned.tags.get(Tag::new(b"CO")) {
                Some(TagValue::Str(s)) => format!("{s} | "),
                _ => String::new(),
            };
            set_tag(
                unaligned,
                b"CO",
                Some(TagValue::Str(format!(
                    "{comment}Cross-species contamination"
                ))),
            );
        }
        Ok(())
    }

    /// `transferAlignmentInfoToPairedRead`.
    fn transfer_to_pair(
        &self,
        first: &mut BamRecord,
        second: &mut BamRecord,
        first_aligned: Option<&BamRecord>,
        second_aligned: Option<&BamRecord>,
        aligned_header: &SamHeader,
        contaminant: bool,
    ) -> Result<(), Failure> {
        if let Some(a) = first_aligned {
            self.transfer_to_fragment(first, a, aligned_header, contaminant)?;
        }
        if let Some(a) = second_aligned {
            self.transfer_to_fragment(second, a, aligned_header, contaminant)?;
        }
        if self.settings.clip_overlapping {
            clip_for_overlapping_reads(first, second, self.settings.hard_clip_overlapping)?;
        }
        set_mate_info(second, first, self.settings.add_mate_cigar);
        if !self.settings.keep_aligner_proper_pair {
            set_proper_pair_flags(second, first, &self.settings.expected);
        }
        Ok(())
    }

    fn maybe_set_pg_tag(&self, rec: &mut BamRecord) {
        if let Some(id) = &self.program_id {
            if self.settings.add_pg_tag {
                set_tag(rec, b"PG", Some(TagValue::Str(id.clone())));
            }
        }
    }

    fn is_contaminant(&self, hits: &Hits) -> Result<bool, Failure> {
        if hits.num_hits() == 0 {
            return Ok(false);
        }
        let Some(primary) = hits.index_of_earliest_primary() else {
            return Err(Failure::IllegalState(
                "No primary alignment was found, despite having nonzero hits.".to_string(),
            ));
        };
        let threshold = self.settings.min_unclipped_bases;
        Ok(
            match (hits.first_of_pair(primary), hits.second_of_pair(primary)) {
                (Some(a), Some(b)) => overclipped(a, threshold) || overclipped(b, threshold),
                (Some(a), None) => overclipped(a, threshold),
                (None, Some(b)) => overclipped(b, threshold),
                (None, None) => {
                    return Err(Failure::IllegalState(
                        "Neither read1 or read2 exist for chosen primary alignment".to_string(),
                    ))
                }
            },
        )
    }
}

/// `clipForOverlappingReads`.
fn clip_for_overlapping_reads(
    read1: &mut BamRecord,
    read2: &mut BamRecord,
    hard: bool,
) -> Result<(), Failure> {
    let overlaps = read1.reference_index == read2.reference_index
        && read1.alignment_start <= alignment_end(read2)
        && read2.alignment_start <= alignment_end(read1);
    if !is_unmapped(read1)
        && !is_unmapped(read2)
        && flag(read1, NEGATIVE) != flag(read2, NEGATIVE)
        && overlaps
    {
        let (pos, neg) = if flag(read1, NEGATIVE) {
            (read2, read1)
        } else {
            (read1, read2)
        };
        clip_3prime_ends_to_5prime_ends(pos, neg, false, false)?;
        if hard {
            clip_3prime_ends_to_5prime_ends(pos, neg, true, true)?;
        }
    }
    Ok(())
}

fn clip_3prime_ends_to_5prime_ends(
    pos: &mut BamRecord,
    neg: &mut BamRecord,
    hard: bool,
    use_unclipped_ends: bool,
) -> Result<(), Failure> {
    let neg_end = if use_unclipped_ends {
        unclipped_end(neg)
    } else {
        alignment_end(neg)
    };
    let pos_start = if use_unclipped_ends {
        unclipped_start(pos)
    } else {
        pos.alignment_start
    };
    let pos_3prime = read_position_ignoring_soft_clips(pos, neg_end);
    if pos_3prime > 0 && pos_3prime < pos.read_bases.len() as i32 {
        clip_3prime_end_of_read_maybe_hard(pos, pos_3prime + 1, hard)?;
    }
    let neg_5prime = read_position_ignoring_soft_clips(neg, pos_start - 1);
    let neg_first = if neg_5prime > 0 {
        (neg.read_bases.len() as i32 + 1) - neg_5prime
    } else {
        0
    };
    if neg_first > 0 {
        clip_3prime_end_of_read_maybe_hard(neg, neg_first, hard)?;
    }
    Ok(())
}

/// `getReadPositionAtReferencePositionIgnoreSoftClips`.
fn read_position_ignoring_soft_clips(rec: &BamRecord, position: i32) -> i32 {
    let mut shift = 0;
    let mut found_non_clip = false;
    let mut elements = Vec::new();
    for element in &rec.cigar.elements {
        if element.op == Op::S {
            elements.push(CigarElement {
                length: element.length,
                op: Op::M,
            });
            if !found_non_clip {
                shift += element.length as i32;
            }
        } else {
            if !matches!(element.op, Op::S | Op::H) {
                found_non_clip = true;
            }
            elements.push(*element);
        }
    }
    read_position_at_reference_position(
        &Cigar::new(elements),
        rec.alignment_start,
        position + shift,
    )
}

/// `SAMRecord.getReadPositionAtReferencePosition(rec, pos, true)`.
fn read_position_at_reference_position(cigar: &Cigar, start: i32, position: i32) -> i32 {
    if position <= 0 {
        return 0;
    }
    let mut last_offset = 0;
    for (read_start, ref_start, length) in alignment_blocks(cigar, start) {
        // `CoordMath.getEnd(refStart, length) >= pos`.
        if ref_start + length > position {
            if position < ref_start {
                return last_offset;
            }
            return position - ref_start + read_start;
        }
        last_offset = read_start + length - 1;
    }
    0
}

/// `AbstractAlignmentMerger.clip3PrimeEndOfRead`.
fn clip_3prime_end_of_read_maybe_hard(
    rec: &mut BamRecord,
    clip_from: i32,
    hard: bool,
) -> Result<(), Failure> {
    if hard {
        if rec.tags.get(Tag::new(b"XB")).is_some() || rec.tags.get(Tag::new(b"XQ")).is_some() {
            return Err(picard(format!(
                "Record {} already contains tags for restoring hard-clipped bases.  This \
                 operation will permanently erase information if it proceeds.",
                rec.read_name
            )));
        }
        let length = rec.read_bases.len() as i32;
        let (from, to) = if flag(rec, NEGATIVE) {
            (0, length - clip_from + 1)
        } else {
            (clip_from - 1, length)
        };
        let mut bases = rec.read_bases[from as usize..to as usize].to_vec();
        let mut quals: Vec<u8> = rec.base_qualities[from as usize..to as usize]
            .iter()
            .map(|q| q + 33)
            .collect();
        if flag(rec, NEGATIVE) {
            reverse_complement(&mut bases);
            quals.reverse();
        }
        set_tag(
            rec,
            b"XB",
            Some(TagValue::Str(String::from_utf8_lossy(&bases).into_owned())),
        );
        set_tag(
            rec,
            b"XQ",
            Some(TagValue::Str(String::from_utf8_lossy(&quals).into_owned())),
        );
    }
    clip_3prime_end_of_read(rec, clip_from, hard)
}

/// Where a merge pass ended up.
struct PassOutput {
    records: Vec<BamRecord>,
}

/// `mergeAlignment`, once. `force_sort` is the second attempt's queryname sort of the aligned
/// reads.
#[allow(clippy::too_many_arguments)]
fn merge_pass(
    merger: &mut Merger,
    picker: &mut Picker,
    unmapped_header: &SamHeader,
    unmapped: &[BamRecord],
    aligned_header: &SamHeader,
    aligned: &[BamRecord],
    reference_dictionary: &[SequenceRecord],
    force_sort: bool,
) -> Result<PassOutput, Failure> {
    let settings = merger.settings;
    merger.header.read_groups = unmapped_header.read_groups.clone();

    // `getQuerynameSortedAlignedRecords`, then the gap filter.
    let mut stream: Vec<BamRecord> = aligned.to_vec();
    if force_sort {
        stream.sort_by(htsjdk_bam::query_name::compare);
    }
    let stream: Vec<BamRecord> = stream
        .into_iter()
        .filter(|rec| {
            if settings.max_gaps == -1 {
                return true;
            }
            let gaps = rec
                .cigar
                .elements
                .iter()
                .filter(|e| matches!(e.op, Op::I | Op::D))
                .count() as i32;
            gaps <= settings.max_gaps
        })
        .filter(|rec| !is_unmapped(rec) && !cigar_maps_no_bases_to_ref(&rec.cigar))
        .collect();

    let mut iterator = AlignedIterator::new(stream, picker, aligned_header)?;

    // `getDictionaryForMergedBam`.
    merger.header.sequences = merge_dictionaries(&aligned_header.sequences, reference_dictionary)?;

    let mut next_aligned = iterator.next(aligned_header)?;

    if let Some(id) = &merger.program_id {
        if unmapped_header.programs.iter().any(|pg| &pg.id == id) {
            return Err(picard(
                "Program Record ID already in use in unmapped BAM file.".to_string(),
            ));
        }
    }

    let mut sink: Vec<BamRecord> = Vec::new();
    let add_if_not_filtered = |sink: &mut Vec<BamRecord>, rec: BamRecord| {
        if settings.include_secondary || !flag(&rec, SECONDARY) {
            sink.push(rec);
        }
    };

    let mut index = 0;
    while index < unmapped.len() {
        let mut rec = unmapped[index].clone();
        index += 1;
        merger.maybe_set_pg_tag(&mut rec);
        let mut second_of_pair: Option<BamRecord> = None;
        if flag(&rec, PAIRED) {
            let Some(next) = unmapped.get(index) else {
                return Err(Failure::Other(
                    "java.util.NoSuchElementException".to_string(),
                ));
            };
            let mut second = next.clone();
            index += 1;
            merger.maybe_set_pg_tag(&mut second);
            if rec.read_name != second.read_name {
                return Err(picard(format!(
                    "Second read from pair not found in unmapped bam: {}, {}",
                    rec.read_name, second.read_name
                )));
            }
            if !flag(&rec, FIRST_OF_PAIR) {
                return Err(picard(format!(
                    "First record in unmapped bam is not first of pair: {}",
                    rec.read_name
                )));
            }
            if !flag(&second, PAIRED) {
                return Err(picard(format!(
                    "Second record in unmapped bam is not marked as paired: {}",
                    second.read_name
                )));
            }
            if !flag(&second, SECOND_OF_PAIR) {
                return Err(picard(format!(
                    "Second record in unmapped bam is not second of pair: {}",
                    second.read_name
                )));
            }
            second_of_pair = Some(second);
        }

        let matched = next_aligned
            .as_ref()
            .is_some_and(|hits| hits.read_name() == rec.read_name);
        if matched {
            let hits = next_aligned.take().expect("matched");
            let clone = hits.num_hits() > 1 || hits.has_supplemental_hits();
            let contaminant = settings.unmap_contaminants && merger.is_contaminant(&hits)?;
            if flag(&rec, PAIRED) {
                let second_original = second_of_pair.clone().expect("paired");
                let mut r1_primary: Option<BamRecord> = None;
                let mut r2_primary: Option<BamRecord> = None;
                let mut reused_first = Some(rec.clone());
                let mut reused_second = Some(second_original.clone());
                for i in 0..hits.num_hits() {
                    let first_aligned = hits.first_of_pair(i);
                    let second_aligned = hits.second_of_pair(i);
                    let is_primary = first_aligned
                        .is_some_and(|r| !is_secondary_or_supplementary(r))
                        || second_aligned.is_some_and(|r| !is_secondary_or_supplementary(r));
                    let (mut first_out, mut second_out) = if clone {
                        (rec.clone(), second_original.clone())
                    } else {
                        (
                            reused_first.take().unwrap_or_else(|| rec.clone()),
                            reused_second
                                .take()
                                .unwrap_or_else(|| second_original.clone()),
                        )
                    };
                    merger.transfer_to_pair(
                        &mut first_out,
                        &mut second_out,
                        first_aligned,
                        second_aligned,
                        aligned_header,
                        contaminant,
                    )?;
                    if is_primary {
                        r1_primary = Some(first_out.clone());
                        r2_primary = Some(second_out.clone());
                    }
                    if !is_unmapped(&first_out) || is_primary {
                        add_if_not_filtered(&mut sink, first_out);
                    }
                    if !is_unmapped(&second_out) || is_primary {
                        add_if_not_filtered(&mut sink, second_out);
                    }
                }
                for is_read1 in [true, false] {
                    let supplementals = if is_read1 {
                        &hits.supplemental_first
                    } else {
                        &hits.supplemental_second
                    };
                    let source = if is_read1 { &rec } else { &second_original };
                    let mate_primary = if is_read1 { &r2_primary } else { &r1_primary };
                    for supplemental in supplementals {
                        let mut out = source.clone();
                        merger.transfer_to_fragment(
                            &mut out,
                            supplemental,
                            aligned_header,
                            contaminant,
                        )?;
                        if let Some(mate) = mate_primary {
                            set_mate_information_on_supplemental(
                                &mut out,
                                mate,
                                settings.add_mate_cigar,
                            );
                        }
                        if !is_unmapped(&out) {
                            add_if_not_filtered(&mut sink, out);
                        }
                    }
                }
            } else {
                for i in 0..hits.num_hits() {
                    let fragment = hits.first[i].as_ref().expect("a fragment hit");
                    if flag(fragment, PAIRED) {
                        return Err(Failure::Other(format!(
                            "java.lang.UnsupportedOperationException: getFragment called for \
                             paired read: {}",
                            describe(fragment, aligned_header)
                        )));
                    }
                    let mut out = rec.clone();
                    let is_primary = !is_secondary_or_supplementary(fragment);
                    merger.transfer_to_fragment(&mut out, fragment, aligned_header, contaminant)?;
                    if !is_unmapped(&out) || is_primary {
                        add_if_not_filtered(&mut sink, out);
                    }
                }
                for supplemental in &hits.supplemental_first {
                    let mut out = rec.clone();
                    merger.transfer_to_fragment(
                        &mut out,
                        supplemental,
                        aligned_header,
                        contaminant,
                    )?;
                    if !is_unmapped(&out) {
                        add_if_not_filtered(&mut sink, out);
                    }
                }
            }
            next_aligned = iterator.next(aligned_header)?;
        } else {
            if let Some(hits) = &next_aligned {
                if rec.read_name.as_str() > hits.read_name() {
                    return Err(Failure::IllegalState(format!(
                        "Aligned record iterator ({}) is behind the unmapped reads ({})",
                        hits.read_name(),
                        rec.read_name
                    )));
                }
            }
            if !settings.aligned_reads_only {
                sink.push(rec);
                if let Some(second) = second_of_pair {
                    sink.push(second);
                }
            }
        }
    }
    if iterator.has_next() {
        let name = iterator
            .next(aligned_header)?
            .map(|hits| hits.read_name().to_string())
            .unwrap_or_default();
        return Err(Failure::IllegalState(format!(
            "Reads remaining on alignment iterator: {name}!"
        )));
    }
    Ok(PassOutput { records: sink })
}

/// `SAMSequenceDictionary.mergeDictionaries(aligned, reference, [M5, LN])`.
fn merge_dictionaries(
    aligned: &[SequenceRecord],
    reference: &[SequenceRecord],
) -> Result<Vec<SequenceRecord>, Failure> {
    let names =
        |d: &[SequenceRecord]| -> Vec<String> { d.iter().map(|s| s.name.clone()).collect() };
    if names(aligned) != names(reference) {
        return Err(Failure::Other(format!(
            "java.lang.IllegalArgumentException: Do not use this function to merge dictionaries \
             with different sequences in them. Sequences must be in the same order as well. \
             Found [{}] and [{}].",
            names(aligned).join(", "),
            names(reference).join(", ")
        )));
    }
    let mut merged = Vec::new();
    for (s1, s2) in aligned.iter().zip(reference) {
        let mut record = SequenceRecord::new(&s1.name, 0);
        // `allTags` is a HashSet: its iteration order is the order the merged record gets them.
        let mut tags: Vec<String> = Vec::new();
        for (key, _) in s1.attributes.iter().chain(s2.attributes.iter()) {
            if !tags.iter().any(|t| t == key) {
                tags.push(key.to_string());
            }
        }
        tags.sort_by_key(|t| {
            let h = string_hash_code(t) as u32;
            ((h ^ (h >> 16)) & 15, 0)
        });
        for tag in tags {
            let v1 = s1.attributes.get(&tag);
            let v2 = s2.attributes.get(&tag);
            if let (Some(a), Some(b)) = (v1, v2) {
                if a != b && (tag == "M5" || tag == "LN") {
                    return Err(Failure::Other(format!(
                        "java.lang.IllegalArgumentException: Cannot merge dictionaries. Found \
                         sequence entry for which tags differ: {} and tag {tag} has the two \
                         values: {a} and {b}.",
                        s1.name
                    )));
                }
            }
            record
                .attributes
                .set(&tag, v1.or(v2).expect("one of them has it"));
        }
        if s1.length != 0 && s2.length != 0 && s1.length != s2.length {
            return Err(Failure::Other(format!(
                "java.lang.IllegalArgumentException: Cannot merge the two dictionaries. Found \
                 sequence entry for which lengths differ: {} has lengths {} and {}",
                s1.name, s1.length, s2.length
            )));
        }
        record.length = if s1.length == 0 { s2.length } else { s1.length };
        merged.push(record);
    }
    Ok(merged)
}

// ---------------------------------------------------------------------------------------------
// The command line.

struct CommandLine {
    pairs: Vec<(String, String)>,
}

impl CommandLine {
    fn from_env() -> Self {
        const ALIASES: [(&str, &str); 17] = [
            ("UNMAPPED", "UNMAPPED_BAM"),
            ("ALIGNED", "ALIGNED_BAM"),
            ("R1_ALIGNED", "READ1_ALIGNED_BAM"),
            ("R2_ALIGNED", "READ2_ALIGNED_BAM"),
            ("O", "OUTPUT"),
            ("PG", "PROGRAM_RECORD_ID"),
            ("PG_VERSION", "PROGRAM_GROUP_VERSION"),
            ("PG_COMMAND", "PROGRAM_GROUP_COMMAND_LINE"),
            ("PG_NAME", "PROGRAM_GROUP_NAME"),
            ("PE", "PAIRED_RUN"),
            ("JUMP", "JUMP_SIZE"),
            ("MAX_GAPS", "MAX_INSERTIONS_OR_DELETIONS"),
            ("RV", "ATTRIBUTES_TO_REVERSE"),
            ("RC", "ATTRIBUTES_TO_REVERSE_COMPLEMENT"),
            ("ORIENTATIONS", "EXPECTED_ORIENTATIONS"),
            ("SO", "SORT_ORDER"),
            ("UNMAP_CONTAM", "UNMAP_CONTAMINANT_READS"),
        ];
        let mut pairs = Vec::new();
        for raw in std::env::args().skip(1) {
            let raw = raw.trim_start_matches('-');
            if let Some((name, value)) = raw.split_once('=') {
                let name = match name {
                    "R1_TRIM" => "READ1_TRIM",
                    "R2_TRIM" => "READ2_TRIM",
                    "MC" => "ADD_MATE_CIGAR",
                    "R" => "REFERENCE_SEQUENCE",
                    other => ALIASES
                        .iter()
                        .find(|(short, _)| *short == other)
                        .map_or(other, |(_, long)| long),
                };
                pairs.push((name.to_string(), value.to_string()));
            }
        }
        CommandLine { pairs }
    }

    fn has_been_set(&self, name: &str) -> bool {
        self.pairs.iter().any(|(n, _)| n == name)
    }

    fn get(&self, name: &str) -> Option<&str> {
        self.pairs
            .iter()
            .rev()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_str())
            .filter(|v| *v != "null")
    }

    /// A collection: each value appends to the default, `null` empties it.
    fn collection(&self, name: &str, default: &[&str]) -> Vec<String> {
        let mut out: Vec<String> = default.iter().map(|s| s.to_string()).collect();
        for (n, v) in &self.pairs {
            if n == name {
                if v == "null" {
                    out.clear();
                } else {
                    out.push(v.clone());
                }
            }
        }
        out
    }

    fn boolean(&self, name: &str, default: bool) -> bool {
        self.get(name)
            .map_or(default, |v| v.eq_ignore_ascii_case("true"))
    }

    fn int(&self, name: &str, default: i32) -> i32 {
        match self.get(name) {
            None => default,
            Some(v) => v.parse().unwrap_or_else(|_| {
                refuse_validation(
                    TOOL,
                    &[format!(
                        "Argument '{name}' cannot be set to '{v}': it is not a number"
                    )],
                )
            }),
        }
    }
}

fn main() {
    let args = CommandLine::from_env();

    // Barclay's mutex checks, in field order.
    for (name, mutex) in [
        (
            "ALIGNED_BAM",
            &["READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM"][..],
        ),
        ("READ1_ALIGNED_BAM", &["ALIGNED_BAM"][..]),
        ("READ2_ALIGNED_BAM", &["ALIGNED_BAM"][..]),
        ("JUMP_SIZE", &["EXPECTED_ORIENTATIONS"][..]),
        ("EXPECTED_ORIENTATIONS", &["JUMP_SIZE"][..]),
    ] {
        let provided: Vec<&str> = mutex
            .iter()
            .copied()
            .filter(|m| args.has_been_set(m))
            .collect();
        if args.has_been_set(name) && !provided.is_empty() {
            refuse_validation(
                TOOL,
                &[format!(
                    "Argument '{name}' cannot be used in conjunction with argument(s) {}",
                    provided.join(" ")
                )],
            );
        }
    }
    for name in ["UNMAPPED_BAM", "OUTPUT"] {
        if !args.has_been_set(name) {
            refuse_validation(
                TOOL,
                &[format!(
                    "Argument {name} was missing: Argument '{name}' is required"
                )],
            );
        }
    }

    let program_record_id = args.get("PROGRAM_RECORD_ID").map(str::to_string);
    let program_version = args.get("PROGRAM_GROUP_VERSION").map(str::to_string);
    let program_command_line = args.get("PROGRAM_GROUP_COMMAND_LINE").map(str::to_string);
    let program_name = args.get("PROGRAM_GROUP_NAME").map(str::to_string);
    let aligned_bams = args.collection("ALIGNED_BAM", &[]);
    let read1 = args.collection("READ1_ALIGNED_BAM", &[]);
    let read2 = args.collection("READ2_ALIGNED_BAM", &[]);

    // `customCommandLineValidation`.
    let any_pg =
        program_record_id.is_some() || program_version.is_some() || program_command_line.is_some();
    let all_pg =
        program_record_id.is_some() && program_version.is_some() && program_command_line.is_some();
    if any_pg && !all_pg {
        refuse_validation(
            TOOL,
            &["PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, and PROGRAM_GROUP_COMMAND_LINE must all \
               be supplied or none should be included."
                .to_string()],
        );
    }
    if read1.is_empty() != read2.is_empty() {
        refuse_validation(
            TOOL,
            &[
                "READ1_ALIGNED_BAM and READ2_ALIGNED_BAM must both be supplied or neither should \
               be included.  For single-end read use ALIGNED_BAM."
                    .to_string(),
            ],
        );
    }
    if aligned_bams.is_empty() && !(!read1.is_empty() && !read2.is_empty()) {
        refuse_validation(
            TOOL,
            &[
                "Either ALIGNED_BAM or the combination of READ1_ALIGNED_BAM and READ2_ALIGNED_BAM \
               must be supplied."
                    .to_string(),
            ],
        );
    }
    if aligned_bams.is_empty() {
        thrown(
            "java.lang.UnsupportedOperationException: READ1_ALIGNED_BAM and READ2_ALIGNED_BAM \
             are not supported by this port",
        );
    }

    let unmapped_path = args.get("UNMAPPED_BAM").unwrap_or_default().to_string();
    let output = args.get("OUTPUT").unwrap_or_default().to_string();
    let reference =
        args.get("REFERENCE_SEQUENCE")
            .unwrap_or_else(|| {
                refuse_validation(
                TOOL,
                &["Argument REFERENCE_SEQUENCE was missing: Argument 'REFERENCE_SEQUENCE' is \
                   required"
                    .to_string()],
            )
            })
            .to_string();

    let expected: Vec<Orientation> = if args.get("JUMP_SIZE").is_some() {
        vec![Orientation::Rf]
    } else {
        let given = args.collection("EXPECTED_ORIENTATIONS", &[]);
        if given.is_empty() {
            vec![Orientation::Fr]
        } else {
            given
                .iter()
                .map(|o| match o.as_str() {
                    "FR" => Orientation::Fr,
                    "RF" => Orientation::Rf,
                    "TANDEM" => Orientation::Tandem,
                    other => refuse_validation(
                        TOOL,
                        &[format!(
                            "Argument 'EXPECTED_ORIENTATIONS' cannot be set to '{other}'"
                        )],
                    ),
                })
                .collect()
        }
    };
    let strategy = match args.get("PRIMARY_ALIGNMENT_STRATEGY").unwrap_or("BestMapq") {
        "BestMapq" => Strategy::BestMapq,
        "EarliestFragment" => Strategy::EarliestFragment,
        "BestEndMapq" => Strategy::BestEndMapq,
        "MostDistant" => Strategy::MostDistant,
        other => refuse_validation(
            TOOL,
            &[format!(
                "Argument 'PRIMARY_ALIGNMENT_STRATEGY' cannot be set to '{other}'"
            )],
        ),
    };
    let unmapping = match args
        .get("UNMAPPED_READ_STRATEGY")
        .unwrap_or("DO_NOT_CHANGE")
    {
        "COPY_TO_TAG" => UnmappingStrategy::CopyToTag,
        "DO_NOT_CHANGE" => UnmappingStrategy::DoNotChange,
        "DO_NOT_CHANGE_INVALID" => UnmappingStrategy::DoNotChangeInvalid,
        "MOVE_TO_TAG" => UnmappingStrategy::MoveToTag,
        other => refuse_validation(
            TOOL,
            &[format!(
                "Argument 'UNMAPPED_READ_STRATEGY' cannot be set to '{other}'"
            )],
        ),
    };
    let settings = Settings {
        clip_adapters: args.boolean("CLIP_ADAPTERS", true),
        bisulfite: args.boolean("IS_BISULFITE_SEQUENCE", false),
        aligned_reads_only: args.boolean("ALIGNED_READS_ONLY", false),
        max_gaps: args.int("MAX_INSERTIONS_OR_DELETIONS", 1),
        retain: Vec::new(),
        remove: Vec::new(),
        reverse: dedup_sorted(args.collection("ATTRIBUTES_TO_REVERSE", &["OQ", "U2"])),
        reverse_complement: dedup_sorted(
            args.collection("ATTRIBUTES_TO_REVERSE_COMPLEMENT", &["E2", "SQ"]),
        ),
        read1_trim: args.int("READ1_TRIM", 0),
        read2_trim: args.int("READ2_TRIM", 0),
        expected,
        sort_order: args.get("SORT_ORDER").unwrap_or("coordinate").to_string(),
        add_mate_cigar: args.boolean("ADD_MATE_CIGAR", true),
        unmap_contaminants: args.boolean("UNMAP_CONTAMINANT_READS", false),
        min_unclipped_bases: args.int("MIN_UNCLIPPED_BASES", 32),
        unmapping,
        clip_overlapping: args.boolean("CLIP_OVERLAPPING_READS", true),
        hard_clip_overlapping: args.boolean("HARD_CLIP_OVERLAPPING_READS", false),
        keep_aligner_proper_pair: args.boolean("ALIGNER_PROPER_PAIR_FLAGS", false),
        include_secondary: args.boolean("INCLUDE_SECONDARY_ALIGNMENTS", true),
        add_pg_tag: args.boolean("ADD_PG_TAG_TO_READS", true),
    };
    // The constructor: `remove` overrides `retain`.
    let retain_given = args.collection("ATTRIBUTES_TO_RETAIN", &[]);
    let remove = args.collection("ATTRIBUTES_TO_REMOVE", &[]);
    let retain: Vec<String> = retain_given
        .into_iter()
        .filter(|a| !remove.contains(a))
        .collect();
    let settings = Settings {
        retain,
        remove,
        ..settings
    };

    // The constructors' checks.
    if settings.min_unclipped_bases < 0 {
        thrown("htsjdk.samtools.SAMException: unclippedBasesThreshold must be non-negative");
    }
    for path in std::iter::once(&unmapped_path)
        .chain(aligned_bams.iter())
        .chain(std::iter::once(&reference))
    {
        if !std::path::Path::new(path).is_file() {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(path)
            ));
        }
    }

    let (unmapped_header, unmapped) = read_input(&unmapped_path);
    if aligned_bams.len() != 1 {
        thrown(
            "java.lang.UnsupportedOperationException: more than one ALIGNED_BAM is not \
             supported by this port",
        );
    }
    let (aligned_header, aligned) = read_input(&aligned_bams[0]);
    let dictionary_text =
        std::fs::read_to_string(std::path::Path::new(&reference).with_extension("dict"))
            .unwrap_or_else(|_| {
                thrown(&format!(
                    "picard.PicardException: No sequence dictionary found for {}.  Use Picard's \
             CreateSequenceDictionary to create a sequence dictionary.",
                    absolute(&reference)
                ))
            });
    let reference_dictionary = htsjdk_bam::sam_file::read_sam(&dictionary_text)
        .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMFormatException: {e:?}")))
        .0
        .sequences;
    let contigs = htsjdk_bam::fasta::read_fasta_file(&reference)
        .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMException: {e:?}")));

    // The program record: the command line's, or the aligned file's when it has exactly one.
    let mut header = SamHeader::new();
    header.set_sort_order("coordinate");
    let mut program_id = None;
    if let Some(id) = &program_record_id {
        let mut pg = ProgramRecord::new(id);
        if let Some(v) = &program_version {
            pg.attributes.set("VN", v);
        }
        if let Some(v) = &program_command_line {
            pg.attributes.set("CL", v);
        }
        if let Some(v) = &program_name {
            pg.attributes.set("PN", v);
        }
        header.programs.push(pg);
        program_id = Some(id.clone());
    } else if aligned_header.programs.len() == 1 {
        header.programs.push(aligned_header.programs[0].clone());
        program_id = Some(aligned_header.programs[0].id.clone());
    }

    let mut picker = Picker {
        strategy,
        random: JavaRandom::new(1),
    };
    let mut merger = Merger {
        settings: &settings,
        header,
        program_id,
    };
    let mut result = merge_pass(
        &mut merger,
        &mut picker,
        &unmapped_header,
        &unmapped,
        &aligned_header,
        &aligned,
        &reference_dictionary,
        false,
    );
    if let Err(Failure::IllegalState(_)) = result {
        result = merge_pass(
            &mut merger,
            &mut picker,
            &unmapped_header,
            &unmapped,
            &aligned_header,
            &aligned,
            &reference_dictionary,
            true,
        );
    }
    let mut records = match result {
        Ok(pass) => pass.records,
        Err(Failure::IllegalState(message)) => {
            thrown(&format!("java.lang.IllegalStateException: {message}"))
        }
        Err(Failure::Other(message)) => thrown(&message),
    };

    let mut header = merger.header.clone();
    if settings.sort_order == "coordinate" {
        records.sort_by(htsjdk_bam::coordinate::compare);
        let bases: HashMap<&str, &[u8]> = contigs
            .iter()
            .map(|c| (c.name.as_str(), c.bases.as_slice()))
            .collect();
        for rec in &mut records {
            if is_unmapped(rec) {
                continue;
            }
            let name = header.sequences[rec.reference_index as usize].name.clone();
            let reference_bases = bases.get(name.as_str()).copied().unwrap_or(&[]);
            fix_record(
                rec,
                reference_bases,
                NmOptions {
                    is_bisulfite_sequence: settings.bisulfite,
                    set_only_uq: false,
                },
            );
        }
    }
    header.set_sort_order(&settings.sort_order);

    let mut writer = BamWriter::new(Vec::new(), &header)
        .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}")));
    for rec in &records {
        writer
            .write(rec)
            .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMException: {e:?}")));
    }
    let bytes = writer
        .finish()
        .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}")));
    if let Err(e) = std::fs::write(&output, bytes) {
        thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}"));
    }
}

/// A `TreeSet<String>`'s contents: sorted, without duplicates.
fn dedup_sorted(mut values: Vec<String>) -> Vec<String> {
    values.sort();
    values.dedup();
    values
}
