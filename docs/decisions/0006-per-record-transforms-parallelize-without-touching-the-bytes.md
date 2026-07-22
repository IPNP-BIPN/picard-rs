# 0006. The per-record transforms parallelize without touching the bytes

**Status:** accepted
**Date:** 2026-07-22

## Why this exists

Single-threaded execution is the reference implementation's main weakness: Picard and GATK
process a SAM/BAM largely on one core. A byte-identical Rust port that also uses the other fifteen
cores is strictly better at the same output, and that multicore edge is a stated goal of the
program, not a nice-to-have. This records where the port takes it and, more importantly, why doing
so does not weaken the byte-identity claim.

## The rule

A transform may run on multiple cores **only where the parallel output is provably the same bytes
as the serial output.** Parallelism is never allowed to be a source of divergence. This is the same
discipline the plan's risk register (R6) imposes on Spark, applied to ordinary shared-memory
threads.

## Where it is safe, and why

The **record transforms** are the clean case. Their work is a map over the records: each record is
transformed independently of the others, reading only the shared, immutable header. For these:

- `par_iter_mut().for_each(...)` mutates the records in place and never reorders them, so the
  resulting record list is identical to a serial loop, element for element.
- `par_iter().filter(...).map(...).collect()` builds a `String` or `Vec`, and rayon guarantees the
  collected result is in the **same order as the sequential iterator** even through `filter` and
  `map`. So the concatenated FASTQ, or the converted record list, is byte-identical.

Neither depends on how rayon splits the work across cores; the output is a pure function of the
input. So `CleanSam`, `AddOrReplaceReadGroups`, `SamToFastq` (unpaired), and `FastqToSam`'s
conversion now run their per-record step in parallel, and their existing oracle conformance tests,
unchanged, still match Picard's goldens byte-for-byte. A determinism test additionally cleans
several thousand reads both ways and asserts the written SAM is identical.

## Where it is not (yet) safe

Three things stay serial on purpose:

- **Stateful sequential passes.** `SamToFastq`'s paired path walks a first-seen map to match mates;
  the pairing depends on encounter order, so it is not a pure per-record map.
- **Stable sorts.** Coordinate and queryname sorting must be *stable* for byte-identity (htsjdk-rs
  decision 0021), and a naive parallel sort is not. A parallel stable sort is possible but is its
  own change with its own proof, so sorts remain serial for now.
- **Floating-point folds.** The metrics collectors sum in floating point, which is not associative
  (decision 0005), so splitting a sum across cores can change the last bit. Those stay serial unless
  a commutative, exact fold is established.

The measured speedup against single-threaded Java is tracked separately in the throughput
benchmark; this decision is only about the correctness precondition that lets the speedup be taken
at all.
