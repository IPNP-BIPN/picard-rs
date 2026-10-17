# Validation

How picard-rs establishes that its output is byte-identical to Picard 3.5.0, and on which paths.
The [README](../README.md) is the user-facing summary; this page is the evidence.

## Method

Every change ports a **named symbol** of the pinned reference source: a specific class or method,
read from the pinned clone and translated. Behaviour is never reconstructed from documentation,
papers, or memory, and each change records the source file and symbol it ports. A tool counts as
ported when its output is byte-identical to the reference on every case and every row measured.


## Status

All 121 Picard tools this repository ports (120 binaries) are measured on two grounds:

- the **conformance suites** compare the port with the pinned reference, byte for byte, on each
  suite's cases (see "Conformance suites");
- a **t=2 covering array** for each tool, regenerated and re-measured by CI on real x86-64, runs
  every row against both sides. Every one of the 121 arrays matches on every row: share **1.000**
  in [`tools/coverage/measured.json`](../tools/coverage/measured.json).

Every argument specific to a tool takes at least two values in its array, a refusal included, and
`JAVASCRIPT_FILE` runs on an embedded JavaScript engine as Picard runs it on Nashorn. What still
does not vary is common to the tools, each with its reason written in
[`fixtures.json`](../tools/coverage/fixtures.json) or the array: the ten arguments the repository
excludes outright (`COMPRESSION_LEVEL`, `CREATE_MD5_FILE`, `MAX_RECORDS_IN_RAM`,
`USE_JDK_DEFLATER`, `USE_JDK_INFLATER`, `VERBOSITY`, `QUIET`, `help`, `showHidden`, `version`),
`TMP_DIR` and `arguments_file`, the output path, and `REFERENCE_SEQUENCE`, `CREATE_INDEX` or
`VALIDATION_STRINGENCY` where a tool's corpus offers one value.

"Byte-identical" means byte-identical on the paths measured: the suites' cases and the arrays'
rows. A row is a pair of argument values, not every combination of all of them.

| Tool | Conformance suite | Covering array (t=2) | Share |
|---|---|---|---|
| `AccumulateQualityYieldMetrics` | oracle-backed | 6 rows, 6 accepted, 2 distinct | **1.000** |
| `AccumulateVariantCallingMetrics` | oracle-backed | 15 rows, 12 accepted, 4 distinct | **1.000** |
| `AddCommentsToBam` | covering array only | 13 rows, 10 accepted, 6 distinct | **1.000** |
| `AddOATag` | oracle-backed | 9 rows, 9 accepted, 4 distinct | **1.000** |
| `AddOrReplaceReadGroups` | oracle-backed | 21 rows, 13 accepted, 13 distinct | **1.000** |
| `BaitDesigner` | oracle-backed | 18 rows, 14 accepted, 14 distinct | **1.000** |
| `BamIndexStats` | oracle-backed | 9 rows, 3 accepted, 1 distinct | **1.000** |
| `BamToBfq` | oracle-backed | 17 rows, 7 accepted, 6 distinct | **1.000** |
| `BedToIntervalList` | oracle-backed | 10 rows, 10 accepted, 4 distinct | **1.000** |
| `BpmToNormalizationManifestCsv` | oracle-backed | 10 rows, 4 accepted, 2 distinct | **1.000** |
| `BuildBamIndex` | covering array only | 18 rows, 9 accepted, 3 distinct | **1.000** |
| `CalculateFingerprintMetrics` | oracle-backed | 19 rows, 19 accepted, 12 distinct | **1.000** |
| `CalculateReadGroupChecksum` | oracle-backed | 6 rows, 6 accepted, 2 distinct | **1.000** |
| `CheckDuplicateMarking` | oracle-backed | 13 rows, 4 accepted, 1 distinct | **1.000** |
| `CheckFingerprint` | oracle-backed | 19 rows, 8 accepted, 8 distinct | **1.000** |
| `CheckIlluminaDirectory` | oracle-backed | 17 rows, 6 accepted, 1 distinct | **1.000** |
| `CheckTerminatorBlock` | covering array only | 12 rows, 6 accepted, 2 distinct | **1.000** |
| `CleanSam` | oracle-backed | 9 rows, 9 accepted, 2 distinct | **1.000** |
| `ClusterCrosscheckMetrics` | oracle-backed | 10 rows, 10 accepted, 5 distinct | **1.000** |
| `CollectAlignmentSummaryMetrics` | oracle-backed | 18 rows, 18 accepted, 16 distinct | **1.000** |
| `CollectArraysVariantCallingMetrics` | oracle-backed | 11 rows, 6 accepted, 3 distinct | **1.000** |
| `CollectBaseDistributionByCycle` | oracle-backed | 12 rows, 9 accepted, 3 distinct | **1.000** |
| `CollectDuplicateMetrics` | oracle-backed | 9 rows, 6 accepted, 4 distinct | **1.000** |
| `CollectGcBiasMetrics` | oracle-backed | 23 rows, 18 accepted, 13 distinct | **1.000** |
| `CollectHiSeqXPfFailMetrics` | oracle-backed | 12 rows, 8 accepted, 4 distinct | **1.000** |
| `CollectHsMetrics` | oracle-backed | 17 rows, 10 accepted, 9 distinct | **1.000** |
| `CollectIlluminaBasecallingMetrics` | oracle-backed | 14 rows, 10 accepted, 8 distinct | **1.000** |
| `CollectIlluminaLaneMetrics` | oracle-backed | 14 rows, 8 accepted, 5 distinct | **1.000** |
| `CollectIndependentReplicateMetrics` | oracle-backed | 14 rows, 13 accepted, 11 distinct | **1.000** |
| `CollectInsertSizeMetrics` | oracle-backed | 23 rows, 13 accepted, 9 distinct | **1.000** |
| `CollectJumpingLibraryMetrics` | oracle-backed | 9 rows, 6 accepted, 4 distinct | **1.000** |
| `CollectMultipleMetrics` | oracle-backed | 19 rows, 13 accepted, 10 distinct | **1.000** |
| `CollectOxoGMetrics` | oracle-backed | 18 rows, 9 accepted, 8 distinct | **1.000** |
| `CollectQualityYieldMetrics` | oracle-backed | 10 rows, 10 accepted, 5 distinct | **1.000** |
| `CollectQualityYieldMetricsFlow` | oracle-backed | 12 rows, 8 accepted, 8 distinct | **1.000** |
| `CollectQualityYieldMetricsSNVQ` | oracle-backed | 12 rows, 8 accepted, 7 distinct | **1.000** |
| `CollectRawWgsMetrics` | oracle-backed | 22 rows, 16 accepted, 14 distinct | **1.000** |
| `CollectRnaSeqMetrics` | oracle-backed | 26 rows, 20 accepted, 16 distinct | **1.000** |
| `CollectRrbsMetrics` | oracle-backed | 21 rows, 11 accepted, 9 distinct | **1.000** |
| `CollectSamErrorMetrics` | oracle-backed | 28 rows, 12 accepted, 11 distinct | **1.000** |
| `CollectSequencingArtifactMetrics` | oracle-backed | 19 rows, 10 accepted, 9 distinct | **1.000** |
| `CollectTargetedPcrMetrics` | oracle-backed | 18 rows, 14 accepted, 10 distinct | **1.000** |
| `CollectUmiPrevalenceMetrics` | oracle-backed | 10 rows, 10 accepted, 6 distinct | **1.000** |
| `CollectVariantCallingMetrics` | oracle-backed | 10 rows, 2 accepted, 2 distinct | **1.000** |
| `CollectWgsMetrics` | oracle-backed | 22 rows, 16 accepted, 14 distinct | **1.000** |
| `CollectWgsMetricsWithNonZeroCoverage` | oracle-backed | 22 rows, 18 accepted, 16 distinct | **1.000** |
| `CombineGenotypingArrayVcfs` | oracle-backed | 15 rows, 9 accepted, 3 distinct | **1.000** |
| `CompareGtcFiles` | oracle-backed | 10 rows, 4 accepted, 1 distinct | **1.000** |
| `CompareMetrics` | oracle-backed | 32 rows, 6 accepted, 3 distinct | **1.000** |
| `CompareSAMs` | oracle-backed | 11 rows, 7 accepted, 3 distinct | **1.000** |
| `ConvertHaplotypeDatabaseToVcf` | oracle-backed | 12 rows, 9 accepted, 3 distinct | **1.000** |
| `ConvertSequencingArtifactToOxoG` | oracle-backed | 13 rows, 5 accepted, 3 distinct | **1.000** |
| `CreateBafRegressMetricsFile` | oracle-backed | 9 rows, 6 accepted, 2 distinct | **1.000** |
| `CreateExtendedIlluminaManifest` | oracle-backed | 13 rows, 4 accepted, 2 distinct | **1.000** |
| `CreateSequenceDictionary` | oracle-backed | 10 rows, 8 accepted, 6 distinct | **1.000** |
| `CreateVerifyIDIntensityContaminationMetricsFile` | oracle-backed | 9 rows, 6 accepted, 2 distinct | **1.000** |
| `CrosscheckFingerprints` | oracle-backed | 24 rows, 8 accepted, 7 distinct | **1.000** |
| `CrosscheckReadGroupFingerprints` | oracle-backed | 20 rows, 7 accepted, 6 distinct | **1.000** |
| `DownsampleSam` | oracle-backed | 13 rows, 13 accepted, 8 distinct | **1.000** |
| `EstimateLibraryComplexity` | oracle-backed | 14 rows, 14 accepted, 9 distinct | **1.000** |
| `ExtractFingerprint` | oracle-backed | 15 rows, 15 accepted, 12 distinct | **1.000** |
| `ExtractIlluminaBarcodes` | oracle-backed | 30 rows, 17 accepted, 16 distinct | **1.000** |
| `ExtractSequences` | oracle-backed | 7 rows, 7 accepted, 4 distinct | **1.000** |
| `FastqToSam` | oracle-backed | 21 rows, 7 accepted, 7 distinct | **1.000** |
| `FifoBuffer` | oracle-backed | 10 rows, 10 accepted, 4 distinct | **1.000** |
| `FilterSamReads` | oracle-backed | 26 rows, 5 accepted, 4 distinct | **1.000** |
| `FilterVcf` | covering array only | 18 rows, 11 accepted, 10 distinct | **1.000** |
| `FindMendelianViolations` | oracle-backed | 12 rows, 5 accepted, 2 distinct | **1.000** |
| `FixMateInformation` | oracle-backed | 12 rows, 7 accepted, 4 distinct | **1.000** |
| `FixVcfHeader` | covering array only | 20 rows, 8 accepted, 8 distinct | **1.000** |
| `GatherBamFiles` | covering array only | 9 rows, 9 accepted, 2 distinct | **1.000** |
| `GatherVcfs` | covering array only | 31 rows, 14 accepted, 4 distinct | **1.000** |
| `GenotypeConcordance` | oracle-backed | 12 rows, 8 accepted, 6 distinct | **1.000** |
| `GtcToVcf` | oracle-backed | 18 rows, 13 accepted, 9 distinct | **1.000** |
| `IdentifyContaminant` | oracle-backed | 13 rows, 13 accepted, 11 distinct | **1.000** |
| `IlluminaBasecallsToFastq` | oracle-backed | 30 rows, 13 accepted, 12 distinct | **1.000** |
| `IlluminaBasecallsToSam` | oracle-backed | 31 rows, 14 accepted, 12 distinct | **1.000** |
| `IntervalListToBed` | oracle-backed | 8 rows, 8 accepted, 5 distinct | **1.000** |
| `IntervalListTools` | oracle-backed | 32 rows, 15 accepted, 5 distinct | **1.000** |
| `LiftOverHaplotypeMap` | oracle-backed | 8 rows, 5 accepted, 3 distinct | **1.000** |
| `LiftOverIntervalList` | oracle-backed | 10 rows, 6 accepted, 5 distinct | **1.000** |
| `LiftoverVcf` | oracle-backed | 16 rows, 16 accepted, 14 distinct | **1.000** |
| `MakeSitesOnlyVcf` | covering array only | 10 rows, 8 accepted, 5 distinct | **1.000** |
| `MakeVcfSampleNameMap` | oracle-backed | 6 rows, 3 accepted, 1 distinct | **1.000** |
| `MarkDuplicates` | oracle-backed | 31 rows, 16 accepted, 14 distinct | **1.000** |
| `MarkDuplicatesWithMateCigar` | oracle-backed | 27 rows, 5 accepted, 3 distinct | **1.000** |
| `MarkIlluminaAdapters` | oracle-backed | 28 rows, 16 accepted, 4 distinct | **1.000** |
| `MeanQualityByCycle` | oracle-backed | 12 rows, 9 accepted, 3 distinct | **1.000** |
| `MergeBamAlignment` | oracle-backed | 29 rows, 16 accepted, 16 distinct | **1.000** |
| `MergePedIntoVcf` | oracle-backed | 19 rows, 19 accepted, 8 distinct | **1.000** |
| `MergeSamFiles` | oracle-backed | 18 rows, 1 accepted, 1 distinct | **1.000** |
| `MergeVcfs` | covering array only | 27 rows, 16 accepted, 11 distinct | **1.000** |
| `NonNFastaSize` | oracle-backed | 7 rows, 7 accepted, 2 distinct | **1.000** |
| `NormalizeFasta` | oracle-backed | 8 rows, 8 accepted, 5 distinct | **1.000** |
| `PositionBasedDownsampleSam` | oracle-backed | 12 rows, 12 accepted, 4 distinct | **1.000** |
| `QualityScoreDistribution` | oracle-backed | 15 rows, 11 accepted, 4 distinct | **1.000** |
| `RenameSampleInVcf` | covering array only | 13 rows, 8 accepted, 5 distinct | **1.000** |
| `ReorderSam` | oracle-backed | 10 rows, 10 accepted, 2 distinct | **1.000** |
| `ReplaceSamHeader` | oracle-backed | 10 rows, 2 accepted, 1 distinct | **1.000** |
| `RevertOriginalBaseQualitiesAndAddMateCigar` | oracle-backed | 15 rows, 13 accepted, 5 distinct | **1.000** |
| `RevertSam` | oracle-backed | 21 rows, 16 accepted, 14 distinct | **1.000** |
| `SamFormatConverter` | oracle-backed | 9 rows, 7 accepted, 2 distinct | **1.000** |
| `SamToFastq` | oracle-backed | 19 rows, 3 accepted, 3 distinct | **1.000** |
| `SamToFastqWithTags` | oracle-backed | 18 rows, 2 accepted, 2 distinct | **1.000** |
| `ScatterIntervalsByNs` | oracle-backed | 10 rows, 10 accepted, 6 distinct | **1.000** |
| `SetNmAndUqTags` | covering array only | 11 rows, 8 accepted, 4 distinct | **1.000** |
| `SetNmMdAndUqTags` | oracle-backed | 11 rows, 8 accepted, 4 distinct | **1.000** |
| `SimpleMarkDuplicatesWithMateCigar` | oracle-backed | 30 rows, 11 accepted, 7 distinct | **1.000** |
| `SortGff` | oracle-backed | 14 rows, 5 accepted, 4 distinct | **1.000** |
| `SortSam` | oracle-backed | 10 rows, 7 accepted, 2 distinct | **1.000** |
| `SortVcf` | covering array only | 10 rows, 7 accepted, 4 distinct | **1.000** |
| `SplitSamByLibrary` | oracle-backed | 9 rows, 9 accepted, 3 distinct | **1.000** |
| `SplitSamByNumberOfReads` | oracle-backed | 13 rows, 12 accepted, 9 distinct | **1.000** |
| `SplitVcfs` | covering array only | 22 rows, 12 accepted, 11 distinct | **1.000** |
| `UmiAwareMarkDuplicatesWithMateCigar` | oracle-backed | 26 rows, 16 accepted, 13 distinct | **1.000** |
| `UpdateVcfSequenceDictionary` | covering array only | 17 rows, 14 accepted, 4 distinct | **1.000** |
| `ValidateSamFile` | oracle-backed | 15 rows, 8 accepted, 1 distinct | **1.000** |
| `VcfFormatConverter` | covering array only | 16 rows, 7 accepted, 5 distinct | **1.000** |
| `VcfToAdpc` | oracle-backed | 15 rows, 12 accepted, 7 distinct | **1.000** |
| `VcfToIntervalList` | covering array only | 11 rows, 8 accepted, 6 distinct | **1.000** |
| `ViewSam` | oracle-backed | 14 rows, 5 accepted, 4 distinct | **1.000** |

## Bit-identity contract

Goldens are produced by the pinned reference in a digest-pinned `linux/amd64` container on
JDK 17, on real x86-64 CI, invoked through `gatk <Tool>` so the argument parser is fixed to
Barclay semantics.

Two settings are pinned in the oracle contract and are part of the claim:

- `USE_JDK_DEFLATER=true` / `USE_JDK_INFLATER=true`. Picard defaults to the Intel GKL
  (ISA-L/igzip) deflater, which emits different bytes than zlib for identical input. Pinning
  the JDK deflater makes the claim reproducible.
- The oracle asserts `IntelDeflaterFactory.usingIntelDeflater()` and fails the run on silent
  native-library degradation, rather than producing a golden that matches no real machine.

Fields legitimately allowed to vary (timestamps, version strings, command lines) are
canonicalized under explicitly declared rules. Adding a canonicalization rule is a reviewable,
CI-gated event, because canonicalization is how a bit-identity claim quietly weakens.

Values that cannot be matched exactly are quarantined and reported with their measured
divergence rate, and the output is described as **bio-identical** rather than
**bit-identical**, following the vocabulary of
[broadinstitute/gatk#9384](https://github.com/broadinstitute/gatk/pull/9384).

## Conformance suites

Every suite is declared once in [`tools/conformance/manifest.json`](../tools/conformance/manifest.json):
the harness that regenerates it, the goldens it produces, and every rule that canonicalizes a field
away, each with the reason it is there. The oracle jobs of `.github/workflows/ci.yml` are generated
from that manifest, and the same run happens locally with

```sh
python3 tools/conformance/run_suite.py --list
python3 tools/conformance/run_suite.py --suites metrics
```

A suite's `status` is part of the claim it supports: **oracle-backed** means CI regenerates the
golden in the pinned container and compares it on every run; **unchecked** means the golden is
committed and read by the Rust tests but has never been re-derived. Decision 0008 records how
sixteen goldens came to be in the second category, and decision 0011 records the answer it left
open: those sixteen came back green, so every suite is oracle-backed today.

## Argument coverage

The conformance suites cover each tool's default path. Coverage of its **argument surface** is
measured separately, with t-wise covering arrays generated in gatk-rs from the pinned inventory and
run here against the oracle:

```sh
python3 tools/coverage/run_array.py --tool CollectAlignmentSummaryMetrics \
    --port target/release/CollectAlignmentSummaryMetrics
```

`tools/coverage/MakeFixtures.java` builds the deterministic corpus the rows run against (a covering
array cannot invent a file path), `fixtures.json` declares which values this repository is willing
to pass and which arguments it refuses to vary, and the run writes a corpus in the same dump format
the conformance suites use. Decision 0009 records what the first two runs measured.

Which tools the CI job runs is declared once, in the `coverage` block of
[`tools/conformance/manifest.json`](../tools/conformance/manifest.json), and the job is generated from
it alongside the oracle jobs. A tool's `status` there carries the same distinction the suites' does:
**measured** means the number in `tools/coverage/measured.json` was produced by that job on real
x86-64 and is re-derived on every run, so a port that gains or loses argument surface fails the run;
**pending** means the array runs and the measurement is published as an artefact, but nothing is
committed, because a coverage number produced on a developer machine is the same mistake as a golden
produced there. Committing the published file and flipping the status is what turns the run into a
claim.

Beyond the array, `tools/fuzz/` mutates from its rows and keeps whatever reaches a probe of the
reference that nothing has reached, with JaCoCo measuring inside the pinned container:

```sh
python3 tools/fuzz/run_fuzz.py --tool CollectAlignmentSummaryMetrics --iterations 200 \
    --port target/release/CollectAlignmentSummaryMetrics
```

Sixty mutations reach about 10% more of the reference than the whole pairwise array does, and the
first two sessions found a numeric divergence the conformance corpus cannot contain: decision 0010.

## Decisions

| # | Title |
|---|---|
| [0001](decisions/0001-the-metrics-archetype-is-not-homogeneous.md) | The metrics archetype is not homogeneous |
| [0002](decisions/0002-the-first-within-stratum-delta.md) | The first within-stratum delta, smaller than the archetype story assumed |
| [0003](decisions/0003-the-delta-at-the-large-end-is-negative.md) | The within-stratum delta at the large end is negative |
| [0004](decisions/0004-the-speedup-is-mostly-jvm-startup-until-it-is-not.md) | The speedup is mostly JVM startup, until the input is big enough that it is not |
| [0005](decisions/0005-the-rnaseq-coverage-fold-order-is-unobservable.md) | The RnaSeq coverage fold order is unobservable |
| [0006](decisions/0006-per-record-transforms-parallelize-without-touching-the-bytes.md) | Per-record transforms parallelize without touching the bytes |
| [0007](decisions/0007-multicore-helps-the-transform-not-the-io.md) | Multicore helps the transform, not the I/O |
| [0008](decisions/0008-sixteen-goldens-had-never-been-produced-by-the-oracle.md) | Sixteen goldens had never been produced by the oracle |
| [0009](decisions/0009-the-first-covering-array-run-measures-the-argument-surface.md) | The first covering-array run measures the argument surface, and it is 0% |
| [0010](decisions/0010-the-fuzzer-found-a-divergence-the-conformance-suite-cannot.md) | The fuzzer found a divergence the conformance suite cannot |
| [0011](decisions/0011-the-sixteen-unchecked-goldens-came-back-green.md) | The sixteen unchecked goldens came back green |
| [0012](decisions/0012-what-the-move-to-3.5.0-costs-the-ported-surface.md) | What the move to Picard 3.5.0 costs the ported surface |
| [0013](decisions/0013-picard-rs-moves-to-3.5.0-on-its-own.md) | picard-rs moves to Picard 3.5.0 on its own |

## Relationship to `fulcrumgenomics/riker`

[riker](https://github.com/fulcrumgenomics/riker) is an independent, MIT-licensed Rust
reimplementation of these same Picard QC tools, from the maintainers of Picard and htsjdk. It is
the closest existing work to this repository, so the distinction matters and is worth stating
plainly.

**riker targets functional equivalence; this repository targets byte equivalence.** riker's own
README says it "is not intended to be a drop-in replacement for Picard": lowercase `snake_case`
headers, no metadata lines, `frac_` for `pct_`, and "bug fixes that yield slightly different
outputs". So riker is the better tool for new work. picard-rs is a byte-for-byte reproduction
of Picard itself, bugs included, so that existing GATK and Picard pipelines can be reproduced
exactly.

That makes riker **a source of divergence candidates, never a source to port from**: copying it
would import its deliberate deviations, and the licence being compatible does not make it correct.
Its `ERRATA.md` is a curated list of exactly where a careful reimplementer differs from Picard,
and every entry is a place this port must *not* differ. Two are already pinned as conformance
cases (`riker_mean_aligned_dilution`, `riker_improper_pair_unmapped_mate` in the alignment-summary
suite), measured against the reference rather than trusted. Where riker's reading of Picard and
this port's agree, that is an independent cross-check of the reading; where riker's errata is
silent on something byte comparison catches (Picard's alignment-block cycle binning, recorded in
`docs/decisions/0003`), that is evidence for the method rather than against riker.

## Commit attribution

Commits are co-authored with the model that wrote them. On 2026-07-21 the history of all three
repositories was rewritten to add that trailer uniformly, at the maintainer's request, changing
every commit SHA. The **current** htsjdk-rs pin was moved to the rewritten commit, but pins in
this repository's *historical* commits name pre-rewrite htsjdk-rs SHAs that no longer exist, so a
checkout of an old picard-rs commit can no longer fetch its exact dependency. Historical builds
before that date are therefore no longer bit-reproducible; current and future ones are. This was
a deliberate trade.
