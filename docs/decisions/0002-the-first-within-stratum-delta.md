# 0002. The first within-stratum delta: measured, and smaller than the archetype story assumed

**Status:** accepted; one data point, explicitly not a program-wide answer
**Date:** 2026-07-21
**Follows:** [0001](0001-the-metrics-archetype-is-not-homogeneous.md)

## What was measured

Two stratum-mates from the `single_pass` stratum, ported over **identical inputs** in one
module so any shared machinery is paid once and visible:

| | Java ported | Rust written (non-test) | Rust per Java line |
|---|---:|---:|---:|
| shared: cycle convention, record filter, `OQ` helper | — | **67** | — |
| `MeanQualityByCycle` (+ `HistogramGenerator`) | 385 | 151 | **0.39** |
| `CollectBaseDistributionByCycle` (+ its bean) | 291 | 169 | **0.58** |

Plus 144 lines of test shared by both.

Decision 0001 fixed the methodology to report **Rust lines per Java line ported**, because the
raw counts are dominated by size differences. On that measure the second tool cost **more** per
line than the first, not less.

## Why the delta did not appear, and what that means

The honest reading is that these two tools share almost nothing beyond a 67-line convention.
`MeanQualityByCycle` accumulates quality sums into two parallel arrays and emits a histogram;
`CollectBaseDistributionByCycle` accumulates base counts into a 5-by-cycle matrix and emits
metric rows. Same *shape*, same traversal, same cycle indexing — and essentially disjoint
bodies.

The amortisation that did occur is real but small and lives elsewhere:

- **The infrastructure**, paid before either: `MetricsFile`, `FormatUtil`, `Histogram`, the BAM
  reader, the oracle image, the conformance-harness pattern. That is what makes a 150-line port
  possible at all.
- **The 67 shared lines**, which is the archetype delta proper: about 18% of the second tool's
  cost avoided.

So the archetype hypothesis is **not refuted**, but its mechanism is different from the one the
plan assumed. The plan's model was "port the shape once, then pay a small delta per member".
What the measurement shows is "pay the *infrastructure* once, then pay roughly full price per
member, minus a small shared-convention discount".

## What this does and does not license

It is **one pair, in one stratum, at the small end of the size distribution**. Both tools are
under 400 Java lines; the stratum also contains `ConvertSequencingArtifactToOxoG` at 1438. A
9.34x spread means this pair is not representative of its own stratum, let alone of the 44.

It also measures line counts, which are a proxy for effort and a poor one. Two of the findings
in this session — the `EnumMap` iteration order, the FPU's NaN sign — cost hours and changed
almost no lines.

What it does establish, concretely:

1. The infrastructure investment paid off. A 291-line Java tool became 169 lines of Rust that is
   byte-identical to Picard, which would have been impossible before `MetricsFile` and the
   harness existed.
2. **The per-member cost does not approach zero.** Any sizing that assumed members 2..n of an
   archetype are nearly free is wrong. On this evidence the right planning assumption is roughly
   **0.4 to 0.6 Rust lines per Java line, per tool, with an ~18% discount for stratum-mates**.
3. The corpus, not the port, is where the risk sits. See below.

## The corpus gap that sabotage found

A sabotage making `baseToInt` case-sensitive produced **zero** divergences. The reason is not a
weak corpus in the ordinary sense: BAM stores bases as nibbles and decoding always yields upper
case (htsjdk-rs decision 0008), so **a lower-case base cannot survive the file format**. That
branch of `baseToInt` is unreachable through BAM input and is exercisable only through SAM text,
which is not yet ported.

Widening the corpus to IUPAC ambiguity codes did close a real gap: sabotaging *their* mapping
diverges all 7 cases. Both facts are recorded in the code, because "the test found nothing" and
"there is nothing to find" look identical from the outside and are not the same.

## Next

The `histogram, multi_level, single_pass` stratum still has `CollectRnaSeqMetrics` (879) and
`CollectAlignmentSummaryMetrics` (978) unported, and those are where a delta at the *large* end
would be measured. Both drag in substantial new input plumbing — refFlat annotation parsing and
reference FASTA reading respectively — which is itself the finding that the stratification's
"environment cost" column was pointing at.
