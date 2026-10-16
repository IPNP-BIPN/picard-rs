# picard-rs

**Picard 3.5.0, in Rust, with byte-identical output.** All 121 command-line programs of the
Broad Institute's [Picard](https://github.com/broadinstitute/picard), reimplemented as native
binaries. Each one writes the same bytes Picard writes on the same input, so an existing pipeline
can switch to it without changing a single result, and without a JVM.

> **This is not the official Picard.** It is an independent reimplementation, not affiliated
> with or endorsed by the Broad Institute.

[![CI](https://github.com/IPNP-BIPN/picard-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/IPNP-BIPN/picard-rs/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/IPNP-BIPN/picard-rs)](https://github.com/IPNP-BIPN/picard-rs/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Install

### Prebuilt binaries

Each [release](https://github.com/IPNP-BIPN/picard-rs/releases) carries one archive per platform
(Linux x86-64, Linux arm64, macOS arm64) with every tool's binary and a SHA-256 checksum:

```bash
curl -LO https://github.com/IPNP-BIPN/picard-rs/releases/download/v0.1.0/picard-rs-v0.1.0-x86_64-unknown-linux-gnu.tar.gz
tar -xzf picard-rs-v0.1.0-x86_64-unknown-linux-gnu.tar.gz
export PATH="$PWD/picard-rs-v0.1.0-x86_64-unknown-linux-gnu/bin:$PATH"
```

### From source

Rust 1.97.1, pinned in `rust-toolchain.toml` (rustup installs it on first build):

```bash
git clone https://github.com/IPNP-BIPN/picard-rs
cd picard-rs
cargo build --release --locked
```

The binaries land in `target/release/`. A clean build takes a few minutes.

## Usage

Each Picard tool is its own binary, named in kebab-case, and takes Picard's own `NAME=value`
arguments, long or short names alike. Where you ran

```bash
java -jar picard.jar CollectAlignmentSummaryMetrics I=in.bam O=metrics.txt R=ref.fasta
```

you run

```bash
collect-alignment-summary-metrics I=in.bam O=metrics.txt R=ref.fasta
```

and get the same `metrics.txt`, byte for byte. The arguments, their defaults and their refusals
are Picard's: the [Picard documentation](https://broadinstitute.github.io/picard/) describes them.

Three names are not the mechanical kebab-case: `CompareSAMs` is `compare-sams`,
`ConvertSequencingArtifactToOxoG` is `convert-sequencing-artifact-to-oxog`, and the deprecated
`SetNmAndUqTags` runs as `set-nm-md-and-uq-tags`, which is what Picard's subclass does.

## Tools

All 121 command-line programs of Picard 3.5.0, by Picard's own program groups:

| Family | Count | Tools |
|---|---|---|
| Diagnostics and QC metrics | 47 | `AccumulateQualityYieldMetrics`, `AccumulateVariantCallingMetrics`, `BamIndexStats`, `CalculateFingerprintMetrics`, `CalculateReadGroupChecksum`, `CheckDuplicateMarking`, `CheckFingerprint`, `CheckTerminatorBlock`, `ClusterCrosscheckMetrics`, `CollectAlignmentSummaryMetrics`, `CollectArraysVariantCallingMetrics`, `CollectBaseDistributionByCycle`, `CollectGcBiasMetrics`, `CollectHiSeqXPfFailMetrics`, `CollectHsMetrics`, `CollectIndependentReplicateMetrics`, `CollectInsertSizeMetrics`, `CollectJumpingLibraryMetrics`, `CollectMultipleMetrics`, `CollectOxoGMetrics`, `CollectQualityYieldMetrics`, `CollectQualityYieldMetricsFlow`, `CollectQualityYieldMetricsSNVQ`, `CollectRawWgsMetrics`, `CollectRnaSeqMetrics`, `CollectRrbsMetrics`, `CollectSamErrorMetrics`, `CollectSequencingArtifactMetrics`, `CollectTargetedPcrMetrics`, `CollectUmiPrevalenceMetrics`, `CollectVariantCallingMetrics`, `CollectWgsMetrics`, `CollectWgsMetricsWithNonZeroCoverage`, `CompareMetrics`, `CompareSAMs`, `ConvertHaplotypeDatabaseToVcf`, `ConvertSequencingArtifactToOxoG`, `CrosscheckFingerprints`, `CrosscheckReadGroupFingerprints` (deprecated), `EstimateLibraryComplexity`, `ExtractFingerprint`, `IdentifyContaminant`, `LiftOverHaplotypeMap`, `MeanQualityByCycle`, `QualityScoreDistribution`, `ValidateSamFile`, `ViewSam` |
| Read data (SAM/BAM) | 31 | `AddCommentsToBam`, `AddOATag`, `AddOrReplaceReadGroups`, `BamToBfq`, `BuildBamIndex`, `CleanSam`, `CollectDuplicateMetrics`, `DownsampleSam`, `FastqToSam`, `FilterSamReads`, `FixMateInformation`, `GatherBamFiles`, `MarkDuplicates`, `MarkDuplicatesWithMateCigar`, `MergeBamAlignment`, `MergeSamFiles`, `PositionBasedDownsampleSam`, `ReorderSam`, `ReplaceSamHeader`, `RevertOriginalBaseQualitiesAndAddMateCigar`, `RevertSam`, `SamFormatConverter`, `SamToFastq`, `SamToFastqWithTags`, `SetNmAndUqTags` (deprecated), `SetNmMdAndUqTags`, `SimpleMarkDuplicatesWithMateCigar`, `SortSam`, `SplitSamByLibrary`, `SplitSamByNumberOfReads`, `UmiAwareMarkDuplicatesWithMateCigar` |
| Variant manipulation (VCF) | 12 | `FixVcfHeader`, `GatherVcfs`, `LiftoverVcf`, `MakeSitesOnlyVcf`, `MakeVcfSampleNameMap`, `MergeVcfs`, `RenameSampleInVcf`, `SortVcf`, `SplitVcfs`, `UpdateVcfSequenceDictionary`, `VcfFormatConverter`, `VcfToIntervalList` |
| Variant evaluation | 2 | `FindMendelianViolations`, `GenotypeConcordance` |
| Variant filtering | 1 | `FilterVcf` |
| Intervals | 4 | `BedToIntervalList`, `IntervalListToBed`, `IntervalListTools`, `LiftOverIntervalList` |
| Reference | 6 | `BaitDesigner`, `CreateSequenceDictionary`, `ExtractSequences`, `NonNFastaSize`, `NormalizeFasta`, `ScatterIntervalsByNs` |
| Base calling (Illumina) | 7 | `CheckIlluminaDirectory`, `CollectIlluminaBasecallingMetrics`, `CollectIlluminaLaneMetrics`, `ExtractIlluminaBarcodes`, `IlluminaBasecallsToFastq`, `IlluminaBasecallsToSam`, `MarkIlluminaAdapters` |
| Genotyping arrays | 9 | `BpmToNormalizationManifestCsv`, `CombineGenotypingArrayVcfs`, `CompareGtcFiles`, `CreateBafRegressMetricsFile`, `CreateExtendedIlluminaManifest`, `CreateVerifyIDIntensityContaminationMetricsFile`, `GtcToVcf`, `MergePedIntoVcf`, `VcfToAdpc` |
| Other | 2 | `FifoBuffer`, `SortGff` |

## How identical is identical

Byte-identical on every path measured, and CI measures them on every change, on real x86-64,
against Picard 3.5.0 itself running in a digest-pinned container:

- **107 conformance suites** compare the port's output with Picard's, byte for byte;
- **a covering array per tool** runs every pair of argument values on both sides. All 121 arrays
  match on every row (share 1.000, [`tools/coverage/measured.json`](tools/coverage/measured.json)).

`USE_JDK_DEFLATER=true` and `USE_JDK_INFLATER=true` are part of the claim: Picard's default Intel
deflater writes different bytes for the same data, so the comparison pins the JDK one. Fields that legitimately vary between runs (timestamps, version strings, command lines)
are canonicalized under declared, reviewed rules.

[`docs/VALIDATION.md`](docs/VALIDATION.md) has the per-tool table, the method, the arguments that
are not varied and why, and the decisions that record each measurement.

## Throughput

Measured, not claimed. `tools/benchmark/run.sh` times the port against Picard on the same input,
inside the same container, and asserts byte-identity of the two outputs in the same run.

At 2,000,000 reads on real x86-64: **1.46x** wall clock, **1.36x** net of JVM startup, output
byte-identical. No optimisation work has been done yet, so 1.46x is a floor;
[decision 0004](docs/decisions/0004-the-speedup-is-mostly-jvm-startup-until-it-is-not.md) records
how the number was measured.

## Reference version

Picard `3.5.0`, commit `a53d7ebd975606b79ca0eaddb7bf93f555ec4db4`, whose release jar bundles
htsjdk `5.0.0`. The three htsjdk 5.0.0 behaviours Picard reaches are ported here, beside the
Picard code that reaches them
([decision 0013](docs/decisions/0013-picard-rs-moves-to-3.5.0-on-its-own.md)).

## Part of a three-repository program

| Repo | Ports | Depends on |
|---|---|---|
| [`htsjdk-rs`](https://github.com/IPNP-BIPN/htsjdk-rs) | htsjdk 4.2.0 | (none) |
| `picard-rs` | Picard 3.5.0 (htsjdk 5.0.0) | `htsjdk-rs` |
| [`gatk-rs`](https://github.com/IPNP-BIPN/gatk-rs) | GATK 4.6.2.0 (in progress) | `picard-rs`, `htsjdk-rs` |

## Who

picard-rs is developed by [BIPN](https://github.com/IPNP-BIPN) (Bioinformatics for Psychiatry &
Neuroscience) at the Institute of Psychiatry and Neuroscience of Paris (IPNP), Inserm U1266,
Université Paris Cité. Commits are co-authored with the AI model that wrote them, and held to the
same CI as any other change.

Issues are welcome, above all a file that does not come out identical to Picard's: that is
exactly what this project exists to find.

## License

MIT, matching Picard. See `LICENSE`.
