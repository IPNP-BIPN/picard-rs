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
//! rest of the header. `ConstantMemory` is the default `STRATEGY`; `HighAccuracy` and `Chained`
//! (which buffer and use a second pass) are separate surfaces.

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
/// arithmetic. `Math.round` is `floor(x + 0.5)`, and the cast to `int` truncates the long to 32 bits.
fn max_hash_value(probability: f64) -> i32 {
    let range = (i32::MAX as i64) - (i32::MIN as i64); // 4_294_967_295
    let rounded = ((range as f64) * probability + 0.5).floor() as i64; // Math.round(range * proportion)
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
