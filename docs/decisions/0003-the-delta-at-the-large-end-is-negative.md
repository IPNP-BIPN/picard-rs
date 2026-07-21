# 0003. The within-stratum delta at the large end: measured, and negative

**Status:** accepted; second data point, and it agrees with the first
**Date:** 2026-07-21
**Follows:** [0001](0001-the-metrics-archetype-is-not-homogeneous.md),
[0002](0002-the-first-within-stratum-delta.md)

## What was measured

Decision 0002 measured a pair at the **small** end of the size distribution, both tools under
400 Java lines, and found the second cost *more* per line than the first. The obvious objection
was that the pair was unrepresentative. So this measures the largest pair available in a single
stratum: `CollectInsertSizeMetrics` and `CollectAlignmentSummaryMetrics`, which share the
`histogram` + `multi_level` + `single_pass` + `needs_r` signals and differ only by
`needs_reference`.

| | Java ported | Rust written (non-test) | Rust per Java line |
|---|---:|---:|---:|
| `CollectInsertSizeMetrics` (first member) | 526 | 399 | **0.76** |
| `CollectAlignmentSummaryMetrics` (second member) | 1217 | 1030 | **0.85** |

The second member's 1030 is `alignment_summary.rs` at 876 plus `adapter.rs` at 154. Its 1217
Java lines are `CollectAlignmentSummaryMetrics`, `AlignmentSummaryMetrics`,
`AlignmentSummaryMetricsCollector`, `AdapterUtility` and `ChimeraUtil`.

Counted the other way, including the new shared machinery the second member forced into
htsjdk-rs — `alignment_block.rs`, `sequence.rs`, `fasta.rs`, and the histogram key-union fix,
about 190 further non-test lines — the second member cost **1220 Rust lines for 1217 Java
lines, a ratio of 1.00**.

## The result

**The second member of the stratum cost more per Java line than the first, at both ends of the
size distribution.** Two independent pairs, one small and one large, and the delta is negative
in both.

The conformance suites are the same size too: 198 lines for the first, 193 for the second. The
harness was not reused; it was copied and adapted, and the copy is as long as the original.

## Why, concretely

What the two tools actually share is 20 lines: `PairOrientation` and `pair_orientation`. Nothing
else. `CollectInsertSizeMetrics` bins insert sizes by orientation into three histograms.
`CollectAlignmentSummaryMetrics` fills 34 metric fields across four pairing categories from a
reference comparison, an adapter matcher and a chimera test. They have the same *shape* — accept
records, accumulate, emit a `MetricsFile` — and disjoint bodies.

The stratum signals turn out to describe the **plumbing**, not the work. `single_pass` says how
records arrive. `histogram` says the output has a histogram section. `multi_level` says the bean
extends `MultilevelMetrics`. Every one of those was paid for before either tool was written, and
none of them is where a tool's lines go.

## What this does to the plan's sizing model

The plan's model was: "port the shape once, then pay a small delta per member." Decision 0002
already replaced it with "pay the infrastructure once, then roughly full price per member." This
confirms that at four times the size, and sharpens it in one direction:

**The per-member cost tracks the member's own Java line count, not its stratum.** For sizing
purposes the useful predictor is the footprint measured by `tools/stratify/stratify.py`, at a
ratio near 0.8 to 1.0 Rust lines per Java line, and the stratum tells you almost nothing beyond
that.

That is a worse answer for the program than the archetype story, and it is the one the
measurements give. Two pairs is still two pairs, and both are in the metrics stratum, which is
the most homogeneous of the twenty-four; the variant callers and the Spark tools have had no
measurement at all. But the direction has now been checked at both ends of the one stratum where
the archetype hypothesis was most likely to hold, and it did not hold at either.

## The honest caveat, restated

Lines are a poor proxy for effort and this is the second decision to say so. The finding that
cost the most time in this port changed almost no lines: `BAD_CYCLES` is binned by the offset
within an alignment block rather than by the read cycle, which took a probe in the oracle to
establish and one comment plus one variable name to reproduce. A sizing model built on line
counts will systematically under-price exactly the work that makes the port bit-identical rather
than merely correct.

## Verification

`crates/picard-analysis/tests/alignment_summary_conformance.rs`: 22 cases, each carrying its
input BAM, its reference FASTA and Picard's own metrics file, all byte-equal under the two
declared canonicalization rules.

Sabotage-checked, and the third result changed the port:

| deliberate break | cases failed |
|---|---|
| index the bad cycle by the read position, as the parameter name asks | 3 |
| drop `MathUtil.divide`'s magnitude threshold for a plain zero guard | 1 |
| remove the supplementary-record guard in `collectReadData` | **0** |

The third failing to fail is the point. It exposed that
`AlignmentSummaryMetricsCollector.acceptRecord` filters `isSecondaryOrSupplementary()` before
any per-unit collector runs, so the guard inside `collectReadData` is unreachable through this
tool — and with it, the comment above it:

```java
// NB: for read count metrics, do not include supplementary records, but for base count
// metrics, do include supplementary records.
```

is false of its own tool. The oracle settles it: the `supplementary` case is one ordinary
20-base read plus one supplementary 20-base read, and Picard reports `PF_ALIGNED_BASES = 20`.
Were the comment true it would be 40. The port reproduces both guards, and does not repeat the
comment as though it described behaviour.


## Prior art: `fulcrumgenomics/riker`

`fulcrumgenomics/riker` is an independent Rust reimplementation of these same Picard QC tools,
MIT-licensed, from the maintainers of Picard and htsjdk. It covers `alignment`, `isize`, `basic`,
`gcbias`, `rna`, `wgs`, `hybcap` and `error` — the same ground this repository is porting.

**It cannot be a source to port from, and the reason is its own stated goal.** Riker says it "is
not intended to be a drop-in replacement for Picard": lowercase `snake_case` headers, no metadata
lines, `frac_` in place of `pct_`, and "bug fixes that yield slightly different outputs in
specific edge cases". Functional equivalence is its target. Byte equivalence is this project's,
and those are not the same target in the one place that matters. Copying riker would import its
deliberate deviations, and the licence being compatible does not make it correct.

**It is, however, the best available map of divergence candidates.** Riker's `ERRATA.md` is a
curated list of exactly the places a careful reimplementer chooses to differ from Picard, written
by people who know this codebase better than anyone. Every entry in it is a place this port must
*not* differ. So each is turned into a corpus case rather than left as prose. Two are pinned:

| riker's claim | Picard, measured | pinned as |
|---|---|---|
| "Picard computes `mean_aligned_read_length` over all PF reads, including unmapped reads which contribute zero to the sum" | one mapped and one unmapped 20-base read gives `MEAN_ALIGNED_READ_LENGTH = 10`, not 20 | `riker_mean_aligned_dilution` |
| "Picard counts all mapped, paired, non-proper reads as improperly paired, including reads whose mate is unmapped" | a mapped non-proper read with an unmapped mate gives `PF_READS_IMPROPER_PAIRS = 1`, not 0 | `riker_improper_pair_unmapped_mate` |

Both are reproduced by this port byte-identically, and neither needed any code written for it,
because the port follows the source. Riker's reading of Picard and this port's agree in both
cases, which is a genuine independent cross-check of the reading rather than of the bytes.

**One asymmetry is worth recording.** Riker's errata does not mention `BAD_CYCLES`, cycle
counting, or alignment-block offsets. The divergence measured above — that the bad-cycle
histogram is binned by the offset within an alignment block rather than by the read cycle — is
not on their list. That is not a criticism of riker, whose goal does not require finding it. It
is evidence about method: byte comparison against the reference surfaces behaviours that reading
the reference and reimplementing it carefully does not, even when the reimplementers are the
reference's own maintainers.
