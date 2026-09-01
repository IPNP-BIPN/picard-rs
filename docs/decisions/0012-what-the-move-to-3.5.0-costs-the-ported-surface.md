# 0012. What the move to Picard 3.5.0 costs the ported surface

**Status:** accepted
**Date:** 2026-09-01

## The question this answers, and the one it does not

The pin stays at Picard `3.4.0` until the tool ports are done. That decision is recorded once, in
[IPNP-BIPN/gatk-rs#810](https://github.com/IPNP-BIPN/gatk-rs/issues/810), and nothing here changes
it: three ports are coherent only because they move together, and htsjdk 4.2.0 to 5.0.0 is a major
bump that has to land in `htsjdk-rs` first.

What was missing was the *size* of the bill. [#114](https://github.com/IPNP-BIPN/picard-rs/issues/114)
lists the upstream release notes, tool by tool, which is a list of what Broad said changed. This
document is the measurement: what changed in the source the ports actually read, and which of it
can move a byte this repository claims.

## The measurement

`git diff 3.4.0..3.5.0 -- src/main/java` in the pinned clone is **29 files**. Twenty of them are
named by a port module in `crates/picard-analysis/src`, matched by the `picard.<package>.<Class>`
references those modules carry in their headers, which is the same text that says which symbol each
one ports.

Every changed file falls into one of four groups, and only the last one can move the bytes.

### 1. The chart argument became optional (10 files, no output change)

`CollectAlignmentSummaryMetrics`, `CollectBaseDistributionByCycle`, `CollectGcBiasMetrics`,
`CollectInsertSizeMetrics`, `CollectMultipleMetrics`, `CollectRnaSeqMetrics`, `CollectRrbsMetrics`,
`CollectWgsMetricsWithNonZeroCoverage`, `MeanQualityByCycle` and `QualityScoreDistribution` all
take the same shape:

```java
-    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.")
+    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.", optional=true)
```

plus a `customCommandLineValidation` that refuses the chart when `RExecutor.runningInGatkLiteDocker()`
says R is absent, and a null guard around the R invocation. The metrics file is written before any
of that and is not touched.

So the **conformance suites do not move**: their goldens are metrics files on the default path. The
**covering arrays do**: an argument that was required becomes optional, which is a row that used to
be a refusal and is now a run. The cost is re-generating those arrays, not re-deriving those
goldens.

### 2. Documentation (1 file)

`AlignmentSummaryMetrics`'s `PCT_CHIMERAS` prose is harmonised with the code
([#2012](https://github.com/broadinstitute/picard/pull/2012)): chimeras are described as a property
of a pair at MAPQ ≥ 20, and the `SA` tag as belonging to either read. The metric is computed the
same way. Nothing to re-measure; the port's own comments should be updated to the clearer wording
when the move happens, because they were written from the old text.

### 3. Compensation for htsjdk 5.0.0 (2 files)

`ViewSam` is the one to keep, because it is a fact about the I/O layer rather than about the tool:

```java
+ // htsjdk 5.0 dropped the trailing newline from SAMRecord.getSAMString(); add one
+ // ourselves when missing so this code works with both pre- and post-5.0 htsjdk.
```

The output is therefore unchanged **because** the tool now compensates. Whatever `htsjdk-rs` does
for `getSAMString` has to make the same choice explicit, or ViewSam's port will silently write one
newline too few. `CommandLineProgram`'s two added lines are whitespace.

### 4. Behaviour, and what it re-opens (16 files)

| what | where | what it costs here |
|---|---|---|
| Physical location widens from `short` to `int` | `ReadEnds`, `PhysicalLocationForMateCigar`, `EstimateLibraryComplexity` (+ its codec's `writeShort`→`writeInt`) | optical-duplicate counts change wherever a flowcell coordinate exceeded 32767 and wrapped. Every golden holding duplicate metrics is re-opened, and the port's own `ReadEnds` equivalent changes type |
| Flow-based duplicate marking forces codec serialization | `MarkDuplicates`, `ReadEndsForMarkDuplicatesCodec` | only under `FLOW_USE_END_IN_UNPAIRED_READS`; the flow path is not ported yet, so this lands as a requirement on the port rather than as a re-measurement |
| `RESTORE_HARDCLIPS` guards on `hasAttribute` before reading `XB`/`XQ` | `RevertSam` | a record carrying one tag and not the other used to throw; it now proceeds. The rejection corpus is where this shows |
| Input-is-a-pipe is no longer re-sorted on error | `AbstractAlignmentMerger`, `SamAlignmentMerger` | a new refusal (`Input is not regular file, cannot sort and restart`) on a path the port does not reach: `MergeBamAlignment`'s inputs here are files |
| `CREATE_INDEX` is plumbed through | `FilterVcf` | a run that produced no index now produces one. `FilterVcf`'s suite lives in gatk-rs |
| Inconclusive results no longer count as errors; the non-zero-LOD check is accumulated rather than re-derived from the metrics | `CrosscheckFingerprints` | changes which comparisons are written under `OUTPUT_ERRORS_ONLY`, so the crosscheck goldens are re-opened |
| SNP matching subsets the VC to the SNP instead of taking the first record at the position | `FingerprintChecker` | reachable by `CheckFingerprint`, `CalculateFingerprintMetrics` and `CrosscheckFingerprints`; a multi-allelic site at a fingerprint position now resolves differently |
| A pair whose mate is on another contig no longer contributes an insert size | `InsertSizeMetricsCollector` | `CollectInsertSizeMetrics`, `CollectMultipleMetrics` and the targeted collectors that embed it. Every insert-size histogram over a corpus with inter-contig pairs moves, which is the one entry here that touches a metrics golden on the default path |
| Proper orientation now requires both mates on the same contig | `RnaSeqMetricsCollector` | `CollectRnaSeqMetrics` ([#2027](https://github.com/broadinstitute/picard/pull/2027)): the assert that used to abort is gone, and `CORRECT_STRAND_READS` moves on any corpus with inter-contig pairs |
| The R-absence helper itself | `RExecutor` | new `runningInGatkLiteDocker()`, read by the ten tools in group 1; no output of its own |

## What this says about the schedule

The move is not 29 files of work. It is:

* **nothing** for most metrics goldens: ten tools changed only around the chart;
* **the covering arrays regenerated**, because ten tools' argument surfaces changed shape;
* **six families re-opened** — duplicate metrics, fingerprinting, `RevertSam`'s rejections,
  `FilterVcf`'s index, and the two that a reading of the release notes alone would have missed:
  insert size and RNA-seq metrics both stop counting pairs whose mates are on different contigs.
  Each needs its golden re-derived on the new oracle and the difference explained rather than
  absorbed;
* **one fact to carry into htsjdk-rs**: `SAMRecord.getSAMString()` loses its trailing newline, and
  a port that keeps it will be wrong in exactly the way Picard 3.5.0 had to work around.

That is a bounded bill, and it is bounded *because* it was measured against the ported symbols
rather than against the release notes. The two collectors are the argument for doing it that way:
the release notes describe #2027 as protecting `CollectRnaSeqMetrics` from an assert, and what the
diff says is that both collectors now require a pair's mates to share a contig, which is a metrics
change on any corpus that has inter-contig pairs. The pin still does not move until Milestone P
closes.

## How this was measured

```sh
git -C picard diff --numstat 6c3f23bc2e0d229d75e9f9b04200396bcd067526..3.5.0 -- src/main/java
```

against the `picard.<package>.<Class>` references carried in the headers of
`crates/picard-analysis/src/*.rs`. The four groups above are 10 + 1 + 2 + 16, which is the 29.

The nine changed files that no port names by fully qualified class are `CollectGcBiasMetrics`
(ported: `gc.rs` names `picard.analysis.GcBiasUtils`, the collector it reproduces, and mentions the
tool class only unqualified), `CommandLineProgram`, `RExecutor`, `FilterVcf` (whose suite is in
gatk-rs), the four `markduplicates/util` classes and `MarkDuplicatesForFlowHelper`. Each is
accounted for above, so the name match is a way of finding what to read, not the measurement
itself: the two collectors that matter most here were found by reading all 29 diffs, and one of
them is named by a module the match did find.
