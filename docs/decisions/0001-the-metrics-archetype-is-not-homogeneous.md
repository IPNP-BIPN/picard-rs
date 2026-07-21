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

## Decision

**Do not report a delta from a single pair.** Porting `CollectQualityYieldMetrics` and then one
neighbour would produce a number, and that number would say more about which neighbour was
picked than about the archetype. Reporting it would be exactly the kind of precision-by-
invention the plan warns against, and it would then be used to size Phases 2 through 5.

The gate is redesigned to sample the archetype rather than sample three arbitrary points:

1. **Stratify first.** Classify all 57 collectors by the machinery they need: single-pass
   counters, `MultiLevelCollector`, `Histogram` percentiles, R chart output, reference sequence.
   This is mechanical and comes from the inventory generator, not from judgement.
2. **Port one member per stratum at full price**, so each distinct piece of shared machinery is
   paid once and its cost is attributed to the stratum rather than smeared across the archetype.
3. **Then measure the delta within a stratum**, which is the only place the amortisation claim
   is actually made.

The number that decides the program is the **within-stratum delta**, and the count of tools in
each stratum weights it.

## Why this is recorded as a decision rather than fixed quietly

The plan's own risk register lists R3, "archetype delta does not fall as hoped", as critical,
and says the gate is what measures it. Discovering that the gate as specified would have
produced a misleading number is a result of running it, and it is worth more than the number
would have been. Recording it here means the redesign is visible, and the eventual delta comes
with the stratification that makes it meaningful.

The delta remains **unmeasured**. Nothing in this repository should be read as evidence about
the cost of tools 2 through 311.
