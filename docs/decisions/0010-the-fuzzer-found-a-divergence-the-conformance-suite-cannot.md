# 0010. The fuzzer found a divergence the conformance suite cannot

**Status:** accepted
**Date:** 2026-07-29

## What the fuzzer is

`tools/fuzz/FuzzDriver.java` runs the reference in one warm JVM with a JaCoCo agent attached,
seeded from the covering array's rows, and keeps a mutant when it reaches a probe nothing has
reached or produces an outcome nothing has produced. `tools/fuzz/run_fuzz.py` drives it, replays
every kept case against the port, and minimizes any divergence.

Measured, first sessions:

| tool | seeds | iterations | probes from seeds | probes after fuzzing | kept |
|---|---:|---:|---:|---:|---:|
| `CollectQualityYieldMetrics` | 10 | 60 | 286 | **314** | 8 |
| `CollectAlignmentSummaryMetrics` | 16 | 60 | 272 | **~300** | 16 |

Sixty mutations add about 10% more reference probes than the whole pairwise array reaches. That is
the number the plan wanted and could not state: it is what combinatorial coverage leaves on the
table, measured rather than argued.

The signal is JaCoCo *probe* coverage, not branch coverage in the strict sense. A probe sits at a
branch boundary, so it is a proxy: good enough to steer a search and to serve as a stopping
condition, not a figure to publish as "branch coverage". Recorded here so it is not quoted as one.

## The finding

Minimized to two arguments, the port and the reference disagree:

```
CollectAlignmentSummaryMetrics --INPUT=<fixtures>/small.bam \
                              --REFERENCE_SEQUENCE=<fixtures>/ref.fasta
```

```
-FIRST_OF_PAIR  ... 35  0.754288  0.753273  ... 21 ...
+FIRST_OF_PAIR  ... 35  0.755496  0.754411  ... 21 ...
-SECOND_OF_PAIR ... 37  0.751375  0.750845  ...  0 ...
+SECOND_OF_PAIR ... 37  0.754875  0.753380  ...  1 ...
```

Three columns: `PF_MISMATCH_RATE`, `PF_HQ_ERROR_RATE`, `BAD_CYCLES`. The port counts slightly more
mismatches, and one more bad cycle.

Without `--REFERENCE_SEQUENCE` the two agree exactly, so this is in the reference-comparison path.
The fixture reference is not the smooth one the existing corpus uses: it carries a twenty-base run
of `N` and a GC-rich stretch, deliberately (`MakeFixtures.java`). The hypothesis that fits the
shape of the difference is **ambiguous-base handling**: htsjdk's `SequenceUtil` does not count a
read base against an `N` reference base as a mismatch the way a plain byte comparison does, and a
mismatch that is counted also lands in the cycle histogram that produces `BAD_CYCLES`. That is a
hypothesis with a measurement behind it, not a conclusion; the slice that fixes it owns the proof.

## Why the conformance suite could not have found this

It is not a weakness of the suite, it is what the suite is for. The alignment-summary corpus was
built to pin the tool's behaviour on the inputs the port was written against, and those inputs
have no `N` in the reference. A conformance corpus tests what its author thought of. A fuzzer
seeded from a covering array and steered by the reference's own coverage tests what the author did
not, which is the only mechanism in the programme that can.

This is the first divergence found by something other than a person reading Java, and it was found
on the second tool it was pointed at, within sixty mutations.

## What happens next

The fix belongs to a tool slice, not here: this is the harness recording what it caught. The
minimized case becomes a conformance case in that slice, and the fix is byte-identity against the
oracle on it.

The fuzz job runs a short session on every CI run and publishes its findings as an artefact.
Findings are not committed, for the reason decision 0008 exists: a finding produced anywhere other
than real x86-64 CI cannot be quoted. The value of running it every time is the regression signal,
a divergence class that disappears or a new one that appears, rather than the depth of any one
session.

## What is still missing

The mutator edits arguments only. The other half of the search is mutating the **input**: dropping
records, corrupting flags, truncating a BAM. That reaches parser and validation branches no
argument can, and it is the natural next slice. `MakeFixtures.java` already produces the corpus it
would mutate.
