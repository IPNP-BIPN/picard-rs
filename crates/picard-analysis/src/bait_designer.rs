//! `BaitDesigner`: turning a list of targets into a list of probes to order.
//!
//! Each target is tiled with fixed-length baits, and what decides the answer is the strategy, the
//! bait's size and offset, and what is done with a target too small to tile: the default lays at
//! least two baits over it whatever its length, which means baits that hang off both ends.
//!
//! The design is written as interval lists, a FASTA, and the pool file that is the actual order.
//! The pool is not the design once: it is the design repeated until the plate is full, every
//! second copy reverse complemented, with the bait numbering restarting at one each time.
//!
//! Ported from `picard.util.BaitDesigner` in Picard 3.4.0.

/// One bait: an interval with the bases under it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bait {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    pub name: String,
    pub bases: Vec<u8>,
}

impl Bait {
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// A target to be covered.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Target {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    pub name: String,
}

impl Target {
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// The three ways of laying baits along a target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    CenteredConstrained,
    FixedOffset,
    Simple,
}

/// What a run was asked for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    pub strategy: Strategy,
    pub bait_size: i32,
    pub bait_offset: i32,
    pub minimum_baits_per_target: i32,
    pub padding: i32,
    pub merge_nearby_targets: bool,
    pub design_on_target_strand: bool,
    pub left_primer: String,
    pub right_primer: String,
    pub pool_size: usize,
    pub fill_pools: bool,
    pub repeat_tolerance: i32,
    pub design_name: String,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            strategy: Strategy::FixedOffset,
            bait_size: 120,
            bait_offset: 80,
            minimum_baits_per_target: 2,
            padding: 0,
            merge_nearby_targets: true,
            design_on_target_strand: false,
            left_primer: "ATCGCACCAGCGTGT".to_string(),
            right_primer: "CACTGCGGCTCCTCA".to_string(),
            pool_size: 55_000,
            fill_pools: true,
            repeat_tolerance: 50,
            design_name: "design".to_string(),
        }
    }
}

/// A bait's name: the target's, and the index padded to the width of the count.
pub fn make_bait_name(target: &str, index: i32, total: i32) -> String {
    let total = total.to_string();
    let mut bait = index.to_string();
    while bait.len() < total.len() {
        bait.insert(0, '0');
    }
    format!("{target}_bait#{bait}")
}

/// How many baits a stretch is expected to take.
///
/// The expression reads like a ceiling and is not one: the rounding is applied to the length
/// before the division, and the division's own result is truncated, so this is the floor of
/// `(length - bait size) / offset` plus one.
pub fn estimate_baits(start: i32, end: i32, options: &Options) -> i32 {
    let length = end - start + 1;
    let tiled = ((length - options.bait_size) as f64 / f64::from(options.bait_offset)) as i32 + 1;
    options.minimum_baits_per_target.max(tiled)
}

/// The targets a design is actually laid over: padded, then merged where two of them would take
/// no more baits together than apart.
pub fn prepare_targets(
    targets: &[Target],
    reference_length: i32,
    options: &Options,
) -> Vec<Target> {
    let padded: Vec<Target> = targets
        .iter()
        .map(|target| Target {
            start: (target.start - options.padding).max(1),
            end: (target.end + options.padding).min(reference_length),
            ..target.clone()
        })
        .collect();
    if !options.merge_nearby_targets {
        return padded;
    }
    let mut merged: Vec<Target> = Vec::new();
    let mut iterator = padded.into_iter();
    let Some(mut previous) = iterator.next() else {
        return merged;
    };
    for next in iterator {
        let apart = estimate_baits(previous.start, previous.end, options)
            + estimate_baits(next.start, next.end, options);
        let together = estimate_baits(previous.start, next.end, options);
        if previous.contig == next.contig && apart >= together {
            previous = Target {
                end: previous.end.max(next.end),
                ..previous
            };
        } else {
            merged.push(previous);
            previous = next;
        }
    }
    merged.push(previous);
    merged
}

/// The bases under a bait, complemented where the design follows the target's strand.
fn bases_for(
    reference: &[u8],
    start: i32,
    end: i32,
    negative_strand: bool,
    options: &Options,
) -> Vec<u8> {
    let bases = reference[(start - 1) as usize..end as usize].to_vec();
    if options.design_on_target_strand && negative_strand {
        return reverse_complement(&bases);
    }
    bases
}

fn bait(
    target: &Target,
    start: i32,
    end: i32,
    index: i32,
    total: i32,
    reference: &[u8],
    options: &Options,
) -> Bait {
    Bait {
        contig: target.contig.clone(),
        start,
        end,
        negative_strand: target.negative_strand,
        name: make_bait_name(&target.name, index, total),
        bases: bases_for(reference, start, end, target.negative_strand, options),
    }
}

/// Baits kept inside the target where that is possible: a target no longer than a bait gets one
/// bait centred on it rather than a tiling that runs past both ends.
pub fn design_centered_constrained(
    target: &Target,
    reference: &[u8],
    options: &Options,
) -> Vec<Bait> {
    let bait_size = options.bait_size;
    if target.length() <= bait_size {
        let midpoint = target.start + target.length() / 2;
        let start = midpoint - bait_size / 2;
        return vec![bait(
            target,
            start,
            start + bait_size - 1,
            1,
            1,
            reference,
            options,
        )];
    }
    let count = 1
        + ((f64::from(target.length() - bait_size) / f64::from(options.bait_offset)).ceil() as i32);
    let first_start = target.start;
    let last_start = target.end - bait_size + 1;
    let shift = f64::from(last_start - first_start) / f64::from(count - 1);
    let mut baits = Vec::new();
    let mut index = 1;
    let mut start = first_start;
    while start <= last_start {
        baits.push(bait(
            target,
            start,
            start + bait_size - 1,
            index,
            count,
            reference,
            options,
        ));
        // Recomputed from the shift each time rather than accumulated, so the rounding does not
        // compound.
        start = first_start + (shift * f64::from(index)).round() as i32;
        index += 1;
    }
    baits
}

/// Baits at a fixed offset, allowed to hang off the ends.
///
/// A target shorter than the minimum number of baits would tile is widened, evenly and with the
/// odd base going to the right, before anything is laid over it.
pub fn design_fixed_offset(target: &Target, reference: &[u8], options: &Options) -> Vec<Bait> {
    let bait_size = options.bait_size;
    let bait_offset = options.bait_offset;
    let minimum = bait_size + bait_offset * (options.minimum_baits_per_target - 1);
    let reference_length = reference.len() as i32;

    let widened = if target.length() < minimum {
        let addon = minimum - target.length();
        let left = addon / 2;
        let right = addon - left;
        Target {
            start: (target.start - left).max(1),
            end: (target.end + right).min(reference_length),
            ..target.clone()
        }
    } else {
        target.clone()
    };

    let count =
        1 + ((f64::from(widened.length() - bait_size) / f64::from(bait_offset)).ceil() as i32);
    let baited_bases = bait_size + bait_offset * (count - 1);
    let first_start = (widened.start - (baited_bases - widened.length()) / 2).max(1);

    let mut baits = Vec::new();
    for index in 1..=count {
        let start = first_start + bait_offset * (index - 1);
        let end = start + bait_size - 1;
        if end > reference_length {
            break;
        }
        baits.push(bait(&widened, start, end, index, count, reference, options));
    }
    baits
}

/// Baits from the target's start until one would begin past its end.
pub fn design_simple(target: &Target, reference: &[u8], options: &Options) -> Vec<Bait> {
    let bait_size = options.bait_size;
    let reference_length = reference.len() as i32;
    let last_possible_start = target.end.min(reference_length - bait_size);
    let count = 1
        + ((f64::from(last_possible_start - target.start) / f64::from(options.bait_offset)).floor()
            as i32);
    let mut baits = Vec::new();
    let mut index = 0;
    let mut start = target.start;
    while start < last_possible_start {
        index += 1;
        baits.push(bait(
            target,
            start,
            start + bait_size - 1,
            index,
            count,
            reference,
            options,
        ));
        start += options.bait_offset;
    }
    baits
}

/// Lay one target, by whichever strategy was asked for.
pub fn design_target(target: &Target, reference: &[u8], options: &Options) -> Vec<Bait> {
    match options.strategy {
        Strategy::CenteredConstrained => design_centered_constrained(target, reference, options),
        Strategy::FixedOffset => design_fixed_offset(target, reference, options),
        Strategy::Simple => design_simple(target, reference, options),
    }
}

/// How many bases of a bait are neither a called base nor upper case.
pub fn masked_base_count(bases: &[u8]) -> i32 {
    bases
        .iter()
        .filter(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
        .count() as i32
}

/// Every bait of a whole design, the too-repetitive ones discarded.
pub fn design(targets: &[Target], reference: &[u8], options: &Options) -> Vec<Bait> {
    let prepared = prepare_targets(targets, reference.len() as i32, options);
    let mut baits = Vec::new();
    for target in &prepared {
        for bait in design_target(target, reference, options) {
            if masked_base_count(&bait.bases) <= options.repeat_tolerance {
                baits.push(bait);
            }
        }
    }
    baits
}

/// The reverse complement of a sequence.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
}

/// A bait's sequence as it is ordered: the primers on each end, and the whole thing complemented
/// where the copy calls for it.
pub fn bait_sequence(bait: &Bait, options: &Options, reverse: bool) -> String {
    let mut sequence = options.left_primer.clone().into_bytes();
    sequence.extend_from_slice(&bait.bases);
    sequence.extend_from_slice(options.right_primer.as_bytes());
    if reverse {
        sequence = reverse_complement(&sequence);
    }
    String::from_utf8(sequence)
        .expect("bases are ASCII")
        .to_uppercase()
}

/// How many times the design is repeated to fill a plate.
pub fn copies(baits: usize, options: &Options) -> usize {
    if options.fill_pools && baits < options.pool_size {
        return options.pool_size / baits;
    }
    1
}

/// The pool file: one row per bait per copy, named by the design and a number that restarts at one
/// with every copy, and every second copy reverse complemented.
pub fn pool_rows(baits: &[Bait], options: &Options) -> Vec<(String, String)> {
    let prefix = format!(
        "{}_",
        &options.design_name[..options.design_name.len().min(8)]
    );
    let mut rows = Vec::new();
    for copy in 0..copies(baits.len(), options) {
        let reverse = copy % 2 == 1;
        for (index, bait) in baits.iter().enumerate() {
            rows.push((
                format!("{prefix}{:06}", index + 1),
                bait_sequence(bait, options, reverse),
            ));
        }
    }
    rows
}

/// The pool file as it is written: the name, a tab, the sequence.
pub fn render_pool(rows: &[(String, String)]) -> String {
    rows.iter()
        .map(|(name, sequence)| format!("{name}\t{sequence}"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// One line of the bait interval list.
pub fn render_bait(bait: &Bait) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}",
        bait.contig,
        bait.start,
        bait.end,
        if bait.negative_strand { "-" } else { "+" },
        bait.name
    )
}

/// The files a design writes, by name, sorted the way a directory listing sorts them.
pub fn written_files(options: &Options) -> Vec<String> {
    let name = &options.design_name;
    let mut files = vec![
        format!("{name}.baits.interval_list"),
        format!("{name}.design.fasta"),
        format!("{name}.design_parameters.txt"),
        format!("{name}.targets.interval_list"),
    ];
    if options.pool_size > 0 {
        files.push(format!("{name}.pool0.design.fasta"));
        files.push(format!("{name}.pool0.design.txt"));
    }
    files.sort();
    files
}
