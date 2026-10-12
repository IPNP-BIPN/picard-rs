//! `DownsampleSam`, the default `ConstantMemory` strategy.
//!
//! Ports `picard.sam.DownsampleSam` driving `htsjdk.samtools.ConstantMemoryDownsamplingIterator` at
//! Picard 3.4.0 / htsjdk 4.2.0, for the default strategy: keep each read whose read-name
//! [`Murmur3`](htsjdk_bam::murmur3::Murmur3) hash falls at or below a probability-derived threshold,
//! and drop the rest, in input order. Because the hash is over the **read name**, both mates of a
//! template share a decision, so a template is kept or dropped as a unit.
//!
//! The keep predicate is a pure function of the read name and the seed, so the filter runs on all
//! cores; rayon's ordered `collect` keeps the survivors in input order, so the parallel output is the
//! same bytes as a serial pass (decision 0006). Unlike the metrics collectors' serial float folds and
//! unlike the thin `AddOATag` string stamp, the per-read work here is a real hash computation, so this
//! is where a multicore transform actually earns its keep.
//!
//! `doWork` adds a `@PG` provenance record whose `CL:` is the command line (temp paths, every option);
//! that record is **canonicalized away** in comparison, exactly as the metrics tools' command-line
//! header is, so this port does not reproduce it and the claim is over the surviving records and the
//! rest of the header. `ConstantMemory` is the default `STRATEGY`; [`keep_flags`] answers all three,
//! `HighAccuracy` and `Chained` included.

use htsjdk_bam::murmur3::Murmur3;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};
use rayon::prelude::*;

/// `RANDOM_SEED`'s default.
pub const DEFAULT_SEED: i32 = 1;

/// `ConstantMemoryDownsamplingIterator`'s threshold: a read is kept when its hash is `<=` this.
///
/// `maxHashValue = Integer.MIN_VALUE + (int) Math.round(range * proportion)`, where
/// `range = (long) Integer.MAX_VALUE - (long) Integer.MIN_VALUE`, in Java's 32-bit **wrapping**
/// arithmetic. The cast to `int` truncates the long to 32 bits.
///
/// `Math.round` is [`jmath::math::round`] and not `floor(x + 0.5)`, which is what this line used to
/// say: the two disagree on the double immediately below a half, where the addition rounds up to
/// exactly 1.0 and the arithmetic version rounds twice. Reaching that input needs a `PROBABILITY`
/// of about 1.16e-10, which the existing corpus does not carry, so this is a latent difference
/// rather than a measured one; it is fixed because it is the same function either way.
fn max_hash_value(probability: f64) -> i32 {
    let range = (i32::MAX as i64) - (i32::MIN as i64); // 4_294_967_295
    let rounded = jmath::math::round((range as f64) * probability);
    i32::MIN.wrapping_add(rounded as i32) // (int) cast truncates, the + wraps
}

/// `DownsampleSam` with the `ConstantMemory` strategy, for SAM input and output.
///
/// `probability` is `PROBABILITY` (0..=1); `seed` is `RANDOM_SEED` ([`DEFAULT_SEED`]).
pub fn downsample_sam(input_sam: &str, probability: f64, seed: i32) -> Result<String, ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;
    let hasher = Murmur3::new(seed);
    let max_hash = max_hash_value(probability);

    // Keep a read iff its name hashes to at or below the threshold (the iterator discards those
    // strictly above it). Independent per record, so the filter is parallel and order-preserving.
    let kept: Vec<BamRecord> = records
        .into_par_iter()
        .filter(|rec| hasher.hash_unencoded_chars(&rec.read_name) <= max_hash)
        .collect();

    Ok(write_sam(&header, &kept).expect("records that parsed re-encode as SAM text"))
}

/// `DownsamplingIteratorFactory.Strategy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    ConstantMemory,
    HighAccuracy,
    Chained,
}

/// `Math.min(a, b)` for doubles, which returns NaN when either side is NaN (Rust's `f64::min`
/// returns the other side).
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// `Math.max(a, b)` for doubles, NaN-propagating as [`java_min`] is.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `String.hashCode`, over UTF-16 code units.
fn java_string_hash(s: &str) -> i32 {
    s.encode_utf16()
        .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(c as i32))
}

/// The iteration order of a `new HashSet<String>()` the names were added to in this order:
/// buckets by `(h ^ h >>> 16) & (n - 1)` at the table size the inserts grew it to, insertion order
/// within a bucket (a resize splits a bucket without reordering it). A bucket that reaches the
/// treeification threshold is reordered (or the table resized early, below 64), which this does
/// not model and refuses rather than answer wrongly.
fn hash_set_order(names: &[String]) -> Result<Vec<String>, String> {
    let mut capacity = 16usize;
    let mut buckets: Vec<Vec<(i32, usize)>> = vec![Vec::new(); capacity];
    let spread = |h: i32| {
        let h = h as u32;
        (h ^ (h >> 16)) as usize
    };
    for (index, name) in names.iter().enumerate() {
        let hash = java_string_hash(name);
        let bucket = &mut buckets[spread(hash) & (capacity - 1)];
        bucket.push((hash, index));
        if bucket.len() > 8 {
            return Err("a HashSet bucket reached the treeification threshold".to_string());
        }
        if index + 1 > capacity * 3 / 4 {
            capacity *= 2;
            let mut grown: Vec<Vec<(i32, usize)>> = vec![Vec::new(); capacity];
            for entry in buckets.into_iter().flatten() {
                grown[spread(entry.0) & (capacity - 1)].push(entry);
            }
            buckets = grown;
        }
    }
    Ok(buckets
        .into_iter()
        .flatten()
        .map(|(_, index)| names[index].clone())
        .collect())
}

/// `ConstantMemoryDownsamplingIterator` as `ChainedDownsamplingIterator` drives it: a peeking
/// iterator over the records whose name hashes at or below the threshold, which skips the rest as
/// soon as it is built and after every `next`, counting both into statistics `reset` clears.
struct ConstantMemory {
    keep: Vec<bool>,
    position: usize,
    seen: i64,
    accepted: i64,
}

impl ConstantMemory {
    fn new(keep: Vec<bool>) -> Self {
        let mut it = ConstantMemory {
            keep,
            position: 0,
            seen: 0,
            accepted: 0,
        };
        it.advance();
        it
    }

    fn advance(&mut self) {
        while self.position < self.keep.len() && !self.keep[self.position] {
            self.position += 1;
            self.seen += 1;
        }
    }

    fn has_next(&self) -> bool {
        self.position < self.keep.len()
    }

    fn next(&mut self) -> usize {
        let index = self.position;
        self.position += 1;
        self.seen += 1;
        self.accepted += 1;
        self.advance();
        index
    }
}

/// Which records `DownsamplingIteratorFactory.make(strategy, proportion, accuracy, seed)` returns,
/// by index into `names` (the records' read names, in input order).
///
/// `HighAccuracyDownsamplingIterator` reads `ceil(1 / accuracy)` new templates at a time, keeps
/// `ceil` or `floor` of `proportion` of them depending on whether it is behind or ahead so far, and
/// picks the discards by shuffling the chunk's `HashSet` of names with `java.util.Random(seed)`;
/// a name decided in one chunk keeps its decision in every later one. `Chained` runs the same over
/// a `ConstantMemory` pass at a proportion widened by a 99.9% interval, and asks of each chunk only
/// what is left to reach the target given what that pass let through.
pub fn keep_flags(
    names: &[String],
    strategy: Strategy,
    proportion: f64,
    accuracy: f64,
    seed: i32,
) -> Result<Vec<bool>, String> {
    let hashed = |p: f64| {
        let hasher = Murmur3::new(seed);
        let max_hash = max_hash_value(p);
        names
            .iter()
            .map(|n| hasher.hash_unencoded_chars(n) <= max_hash)
            .collect::<Vec<bool>>()
    };
    if strategy == Strategy::ConstantMemory {
        return Ok(hashed(proportion));
    }
    if accuracy >= 1.0 || accuracy <= 1.0 / i32::MAX as f64 {
        return Err(
            "java.lang.IllegalArgumentException: Illegal value. Must be 1/MAX_INT < accuracy < 1"
                .to_string(),
        );
    }
    let mut underlying = match strategy {
        Strategy::Chained => {
            let widened = java_min(1.0, proportion + 3.3 * (proportion / 50000.0).sqrt());
            ConstantMemory::new(hashed(widened))
        }
        _ => ConstantMemory::new(vec![true; names.len()]),
    };
    let chained = strategy == Strategy::Chained;
    let templates_to_read = (1.0 / accuracy).ceil() as i32 as usize;
    let mut random = crate::theoretical_sensitivity::JavaRandom::new(seed as i64);
    let mut decisions: std::collections::HashMap<&str, bool> = std::collections::HashMap::new();
    let (mut total_templates, mut kept_templates) = (0i64, 0i64);
    let mut keep = vec![false; names.len()];
    loop {
        if chained {
            underlying.seen = 0;
            underlying.accepted = 0;
        }
        let mut records = Vec::new();
        let mut chunk: Vec<String> = Vec::new();
        let mut in_chunk = std::collections::HashSet::new();
        while underlying.has_next() && chunk.len() < templates_to_read {
            let index = underlying.next();
            records.push(index);
            let name = names[index].as_str();
            if decisions.contains_key(name) {
                continue;
            }
            if in_chunk.insert(name) {
                chunk.push(name.to_string());
            }
        }
        if records.is_empty() {
            break;
        }
        let read = chunk.len() as i64;
        let p = if chained {
            let prior = underlying.accepted as f64 / underlying.seen as f64;
            java_max(0.0, java_min(1.0, proportion / prior))
        } else {
            proportion
        };
        let raw = read as f64 * p;
        let to_keep = if (kept_templates as f64 / total_templates as f64) < p {
            raw.ceil() as i32 as i64
        } else {
            raw.floor() as i32 as i64
        };
        let mut shuffled = hash_set_order(&chunk)?;
        for i in (1..shuffled.len()).rev() {
            let j = random.next_int(i as i32 + 1) as usize;
            shuffled.swap(i, j);
        }
        let discarded: std::collections::HashSet<&str> = shuffled
            .iter()
            .take((read - to_keep).max(0) as usize)
            .map(String::as_str)
            .collect();
        total_templates += read;
        kept_templates += read - discarded.len() as i64;
        for index in records {
            let name = names[index].as_str();
            let decision = *decisions
                .entry(name)
                .or_insert_with(|| in_chunk.contains(name) && !discarded.contains(name));
            keep[index] = decision;
        }
    }
    Ok(keep)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input() -> String {
        let mut s = String::from(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:rg1\tSM:s\n",
        );
        for (i, pos) in (0..12).map(|i| (i, 100 + i * 100)) {
            s.push_str(&format!(
                "read{i}\t0\tchr1\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            ));
        }
        s
    }

    fn kept_names(sam: &str) -> Vec<String> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    }

    #[test]
    fn a_probability_of_one_keeps_every_read() {
        let out = downsample_sam(&input(), 1.0, DEFAULT_SEED).unwrap();
        assert_eq!(kept_names(&out).len(), 12);
    }

    #[test]
    fn a_probability_of_zero_keeps_nothing() {
        let out = downsample_sam(&input(), 0.0, DEFAULT_SEED).unwrap();
        assert!(kept_names(&out).is_empty());
    }

    #[test]
    fn half_keeps_the_reads_whose_name_hashes_at_or_below_zero() {
        // seed=1, P=0.5 -> maxHashValue wraps to 0. The set is fixed by htsjdk's Murmur3 and matches
        // what DownsampleSam emits (verified against the tool).
        let out = downsample_sam(&input(), 0.5, DEFAULT_SEED).unwrap();
        assert_eq!(
            kept_names(&out),
            ["read1", "read3", "read6", "read7", "read8", "read9"]
        );
    }

    #[test]
    fn the_threshold_arithmetic_matches_htsjdk() {
        assert_eq!(max_hash_value(0.5), 0);
        assert_eq!(max_hash_value(1.0), i32::MAX);
        assert_eq!(max_hash_value(0.0), i32::MIN);
    }

    #[test]
    fn the_surviving_reads_pass_through_verbatim_in_input_order() {
        let out = downsample_sam(&input(), 0.5, DEFAULT_SEED).unwrap();
        // read1 is unchanged from the input line.
        assert!(out.contains("read1\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1"));
    }
}
