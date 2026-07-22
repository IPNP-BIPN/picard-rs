# 0007. Multicore helps the transform, not the I/O, so a light transform hides the win

**Status:** accepted
**Date:** 2026-07-22

## Why this exists

Decision 0006 established that the per-record transforms may run on all cores because the parallel
output is provably the same bytes as the serial output. That is the correctness precondition. This
records what the parallelism is actually *worth* end to end, measured rather than assumed, because a
multicore edge that does not show up in wall-clock time is not the edge the program claims.

## What was measured

The throughput job runs `AddOATag` over the 2,000,000-read benchmark BAM on a real x86-64 runner,
three ways (Picard, Rust with `RAYON_NUM_THREADS=1`, Rust on all cores), and additionally dumps the
Rust phase timings at one thread and at all cores. The output is byte-identical to Picard in every
case (asserted in the same job).

The phase breakdown is the point:

| phase | time (all cores) |
|---|---|
| decode (BGZF already inflated) | ~0.82 s |
| **add_oa (the parallel transform)** | **~0.19 s** |
| encode_sam | ~1.45 s |

The `add_oa` step is the *only* thing rayon parallelizes, and it is a small fraction of the run.
Decode, SAM encoding, and writing a ~300 MB text file are serial and dominate. So the whole-program
numbers barely move between one thread (~7.4 s) and all cores (~7.1 s), even though the transform
step itself scales close to linearly with cores.

## The finding

**Multicore accelerates the per-record work, not the I/O around it.** For a *light* transform like
`AddOATag` (a little string formatting per read), the work is a few percent of an I/O- and
serialization-bound pipeline, so the end-to-end win is in the noise. Rust lands on par with Picard
here and byte-identical, which is the honest headline: not "faster because parallel", but "the same
bytes, at parity, with the transform step already parallel and idle capacity to spare".

The multicore edge becomes an *end-to-end* win only when the per-record work is heavy enough to rival
the fixed I/O cost: the CPU-bound inner loops (PairHMM and genotyping in the variant callers, the
heavier metrics folds) are where it will pay, not the thin record transforms. The `RAYON_NUM_THREADS`
knob and the phase timings are kept in the benchmark so that this ratio is visible for each new
transform rather than assumed, and so a future heavy transform can show the win it actually has.

## Consequence

The benchmark reports the transform-phase speedup (1 thread vs N) alongside the whole-program time, so
neither over-claims. Byte-identity is gated in the same job: a throughput number is never reported for
output that does not match. This keeps the multicore claim honest and measured, per the risk register's
R6 discipline applied to ordinary shared-memory threads.
