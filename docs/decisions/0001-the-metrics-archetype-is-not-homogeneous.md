# 0001. The metrics archetype is not homogeneous, so a delta from one pair would be an artefact

**Status:** accepted; calibration gate redesigned, delta still unmeasured
**Date:** 2026-07-21

## Context

The program's feasibility rests on one number. The plan commits to all 311 GATK and Picard
tools, and the honest effort range for that is 40 to 100 person-years unless the marginal cost
of tool number 300 falls far below the cost of tool number 3. The mechanism for that is
**archetypes**: port a tool shape once, then pay only the delta for each further member.

The metrics collector is the largest archetype, **57 of the 311 tools**. The plan's Phase 1
calibration gate says to port three representatives and measure the marginal cost of the second
and third before committing to fan-out.

## What was actually measured

The first member is done. `CollectQualityYieldMetrics` matches Picard 3.4.0 on all 10 goldens.

| | |
|---|---|
| Java ported | 314 lines, one file |
| Rust written | 200 non-test lines, 169 test lines |
| shared infrastructure it needed | picard-rs oracle, BAM reader, `MetricsFile`, `FormatUtil`, the conformance harness pattern |

All of that shared infrastructure now exists, so the second member should pay only its own
logic. That is the hypothesis the gate was going to test.

## The finding that changes the gate

Before porting a second member, the candidates were measured. They are not the same size, and
not by a small margin:

| tool | Java lines, including its collector and bean classes |
|---|---:|
| `CollectQualityYieldMetrics` | **314** (1 file) |
| `CollectInsertSizeMetrics` | **526** (3 files) |
| `CollectAlignmentSummaryMetrics` | **978** (3 files) |

A **3.1×** spread across three members of one archetype. Worse, the spread is not only in
volume. `CollectQualityYieldMetrics` accumulates counters over a single pass. The other two pull
in machinery it never touches: `MultiLevelCollector`, which accumulates separately per sample,
per library and per read group, and `Histogram`'s percentile logic (median, median absolute
deviation) with its own numerical behaviour.

Several members also require an R installation for their chart output
(`QualityScoreDistribution`, `MeanQualityByCycle`, `CollectBaseDistributionByCycle`), which is
an oracle-environment cost, not a porting cost, and it falls on some members and not others.

## The stratification, once actually run

`tools/stratify/stratify.py` classifies every Picard tool that writes a `MetricsFile`, using
only symbols present in the pinned source. Three results, in order of how much they change the
picture.

**There are 44 such tools, not 57.** The plan's group table lists 57 under "Metrics". That is
the same conflation the inventory correction already caught once: 57 is the count of *metric
definition classes*, which are output-file schemas with no CLI. The archetype is 44 tools.

**Scanning one file per tool misclassifies nearly half of them.** A tool's machinery usually
lives in a class it delegates to: `CollectInsertSizeMetrics` reads as a plain single-pass tool,
and its `InsertSizeMetricsCollector` extends `MultiLevelCollector`. Following referenced
collector classes one level, **20 of 44 tools (45%) gain machinery invisible in their own
file**. A stratification without that step would have put half the archetype in the wrong
stratum and flattered any delta measured inside it.

**Environment requirements must be separated from porting machinery.** Needing R, a reference
sequence or an interval list is oracle-and-input cost, paid once for the repository. Needing
`MultiLevelCollector`, `Histogram` or `MergeableMetricBase` is code that must be ported.
Stratifying on everything gives 24 strata for 44 tools, which is close to no amortisation at
all; stratifying on porting machinery alone gives **11**:

| tools | line spread | machinery |
|---:|---|---|
| 12 | 114-1074 | (none of the four) |
| 6 | 236-1047 | `histogram`, `mergeable_base` |
| 5 | 158-613 | `histogram` |
| 4 | 77-381 | `mergeable_base`, `single_pass` |
| 4 | 136-497 | `mergeable_base` |
| 3 | 191-245 | `histogram`, `multi_level`, `single_pass` |
| 3 | 154-251 | `single_pass` |
| 3 | 105-549 | `histogram`, `mergeable_base`, `single_pass` |
| 2 | 176-237 | `histogram`, `single_pass` |
| 1 | 219 | `multi_level`, `single_pass` |
| 1 | 231 | `histogram`, `multi_level` |

Stratification reduces the fragmentation. It does **not** homogenise size: the six-tool stratum
still spans 236 to 1047 lines, a 4.4x range.

## Which triple the gate should use

One stratum is tight enough for the delta to mean something:

| tool | lines |
|---|---:|
| `CollectInsertSizeMetrics` | 191 |
| `CollectRnaSeqMetrics` | 215 |
| `CollectAlignmentSummaryMetrics` | 245 |

Three members, `histogram` + `multi_level` + `single_pass`, spanning 191 to 245 lines. That is
the calibration triple the plan asked for, now **chosen by measurement rather than by guess**.
The first pays for `MultiLevelCollector` and `Histogram`; the second and third pay the delta,
and because their sizes are within 28% of each other the delta is attributable to the
amortisation rather than to the sample.

The already-ported `CollectQualityYieldMetrics` is in a different stratum
(`mergeable_base` + `single_pass`), so it is not a member of this triple and its cost does not
belong in the delta.

## Decision

**Do not report a delta from a single pair.** Porting `CollectQualityYieldMetrics` and then one
neighbour would produce a number, and that number would say more about which neighbour was
picked than about the archetype. Reporting it would be exactly the kind of precision-by-
invention the plan warns against, and it would then be used to size Phases 2 through 5.

The gate is redesigned to sample the archetype rather than sample three arbitrary points:

1. **Stratify first**, mechanically, from symbols in the pinned source, following referenced
   collector classes so delegated machinery is not missed. Done; see above.
2. **Port the tight triple** `CollectInsertSizeMetrics`, `CollectRnaSeqMetrics`,
   `CollectAlignmentSummaryMetrics`, whose sizes are within 28% of each other.
3. **Measure the delta of the second and third**, which is the only place the amortisation
   claim is actually made.

The number that decides the program is the **within-stratum delta**, and the count of tools in
each stratum weights it. On the numbers above, the largest stratum holds 6 of 44 tools, so even
a good delta amortises over a smaller group than "57 collectors, one shape" suggested.

## Why this is recorded as a decision rather than fixed quietly

The plan's own risk register lists R3, "archetype delta does not fall as hoped", as critical,
and says the gate is what measures it. Discovering that the gate as specified would have
produced a misleading number is a result of running it, and it is worth more than the number
would have been. Recording it here means the redesign is visible, and the eventual delta comes
with the stratification that makes it meaningful.

The delta remains **unmeasured**. Nothing in this repository should be read as evidence about
the cost of tools 2 through 311.
