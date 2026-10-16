# picard-rs

Picard 3.5.0's 121 tools as native binaries, with byte-identical output. If you run
[Picard](https://github.com/broadinstitute/picard) today, you can run picard-rs with the same
arguments on the same files and get the same output files, without Java.

> **This is not the official Picard.** It is an independent reimplementation, not affiliated
> with or endorsed by the Broad Institute.

[![CI](https://github.com/IPNP-BIPN/picard-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/IPNP-BIPN/picard-rs/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/IPNP-BIPN/picard-rs)](https://github.com/IPNP-BIPN/picard-rs/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Install

### Prebuilt binaries

Download the archive for your platform from the
[latest release](https://github.com/IPNP-BIPN/picard-rs/releases/latest). There is one per
platform (Linux x86-64, Linux arm64, macOS arm64), each with a `.sha256` file beside it.

```bash
VERSION=v0.1.0
TARGET=x86_64-unknown-linux-gnu   # or aarch64-unknown-linux-gnu, aarch64-apple-darwin
curl -LO https://github.com/IPNP-BIPN/picard-rs/releases/download/$VERSION/picard-rs-$VERSION-$TARGET.tar.gz
tar -xzf picard-rs-$VERSION-$TARGET.tar.gz
export PATH="$PWD/picard-rs-$VERSION-$TARGET/bin:$PATH"
```

### From source

You need [rustup](https://rustup.rs). The toolchain (Rust 1.97.1) is pinned in
`rust-toolchain.toml` and installed on the first build.

```bash
git clone https://github.com/IPNP-BIPN/picard-rs
cd picard-rs
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"
```

## Running a tool

Every Picard tool is its own binary, named in kebab-case: `SortSam` is `sort-sam`,
`CollectAlignmentSummaryMetrics` is `collect-alignment-summary-metrics`. Arguments are Picard's,
written the way Picard takes them, `NAME=value`, with the long or the short name:

```bash
# Picard
java -jar picard.jar MarkDuplicates I=sorted.bam O=marked.bam M=dup_metrics.txt

# picard-rs
mark-duplicates I=sorted.bam O=marked.bam M=dup_metrics.txt
```

The records, the metrics and the rest of the header in `marked.bam` and `dup_metrics.txt` come
out identical to Picard's; [Getting the same bytes as Picard](#getting-the-same-bytes-as-picard)
lists the setting that matters and the header lines that do not match yet. The arguments, their defaults and what they do are Picard's, so the
[Picard documentation](https://broadinstitute.github.io/picard/command-line-overview.html) is the
reference for every tool.

### Keeping your scripts as they are

If your scripts call `picard <Tool> ...` (the Bioconda wrapper's form), this shell function gives
them the same command with picard-rs underneath. Put it in your `~/.bashrc` or `~/.zshrc`:

```bash
picard() {
  local tool="$1"; shift
  case "$tool" in
    CompareSAMs) tool=compare-sams ;;
    ConvertSequencingArtifactToOxoG) tool=convert-sequencing-artifact-to-oxog ;;
    SetNmAndUqTags) tool=set-nm-md-and-uq-tags ;;
    *) tool=$(printf '%s' "$tool" | sed -E 's/([a-z0-9])([A-Z])/\1-\2/g; s/([A-Z])([A-Z][a-z])/\1-\2/g' | tr '[:upper:]' '[:lower:]') ;;
  esac
  "$tool" "$@"
}
```

```bash
picard SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate
picard CollectInsertSizeMetrics I=sorted.bam O=insert_size.txt H=insert_size.pdf
```

Three names are not the mechanical kebab-case, and the function handles them: `CompareSAMs` is
`compare-sams`, `ConvertSequencingArtifactToOxoG` is `convert-sequencing-artifact-to-oxog`, and
the deprecated `SetNmAndUqTags` runs as `set-nm-md-and-uq-tags`, which is what Picard's subclass
does.

### A typical run

```bash
add-or-replace-read-groups I=aligned.bam O=rg.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=sample1
sort-sam I=rg.bam O=sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
mark-duplicates I=sorted.bam O=marked.bam M=dup_metrics.txt CREATE_INDEX=true
collect-alignment-summary-metrics I=marked.bam O=alignment_summary.txt R=reference.fasta
collect-insert-size-metrics I=marked.bam O=insert_size.txt H=insert_size.pdf
```

### What differs from running Picard

The output files are the same. The command line around them has a few differences, listed here so
none of them surprises you:

- **Only `NAME=value` is read.** Picard's newer `--INPUT x` syntax is not: `-h`, `--help` and
  `--version` work, and an argument the tool does not declare is refused, but a declared one written
  that way is ignored.
- **`COMPRESSION_LEVEL`, `MAX_RECORDS_IN_RAM`, `TMP_DIR`, `USE_JDK_DEFLATER`, `USE_JDK_INFLATER`,
  `VERBOSITY` and `QUIET` are accepted and have no effect.** Output is always written the way Picard
  writes it at its default compression level with the JDK deflater, and memory is not bounded by a
  record count.
- **Headers**: metrics files carry placeholder command-line and date lines, and no tool adds its
  own `@PG` record (see [Getting the same bytes as Picard](#getting-the-same-bytes-as-picard)).
- An unknown or misspelled argument is refused with Picard's own message, usage and exit status.
  Any other error prints `Error: ...` on stderr and exits with status 1.

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

## Getting the same bytes as Picard

picard-rs reproduces Picard run with `USE_JDK_DEFLATER=true USE_JDK_INFLATER=true`. By default
Picard compresses with the Intel deflater, which writes different bytes for the same data, so a
BAM or a compressed VCF from a default Picard run decompresses to the same content but is not the
same file. To compare file for file, run Picard with those two arguments:

```bash
java -jar picard.jar SortSam I=in.bam O=picard.bam SORT_ORDER=coordinate USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
sort-sam I=in.bam O=picard-rs.bam SORT_ORDER=coordinate
cmp picard.bam picard-rs.bam && echo identical
```

Three kinds of header line do not match Picard's yet, and the comparisons in CI remove them
before comparing:

- a metrics file's `# <Tool> <arguments>` and `# Started on: <date>` lines, which picard-rs writes
  as the placeholders `<command line>` and `<timestamp>`;
- the `@PG` record a tool adds to the SAM/BAM header in Picard (`MarkDuplicates`, for one, adds
  `@PG ID:MarkDuplicates` with its version and command line). picard-rs does not add it, although
  the reads still carry `PG:Z:MarkDuplicates`;
- `UR:` paths in `@SQ` lines and sequence dictionaries, which are where the files were on disk.

Everything else, records and metrics included, is compared byte for byte. So a SortSam BAM compares with
`cmp`; for a metrics file or a MarkDuplicates BAM, compare everything below
those lines.

How this is measured: on every change, CI runs Picard 3.5.0 in a pinned container on real x86-64
and compares byte for byte, on 107 conformance suites and on a covering array per tool that runs
every pair of argument values on both sides. All 121 arrays match on every row.
[`docs/VALIDATION.md`](docs/VALIDATION.md) has the per-tool results and the method.

## Performance

At 2,000,000 reads on real x86-64, picard-rs runs **1.46x** faster than Picard in wall-clock time
(1.36x net of JVM startup), with byte-identical output. On small inputs the gap is larger, because
there is no JVM to start. No optimisation work has been done yet.
[Decision 0004](docs/decisions/0004-the-speedup-is-mostly-jvm-startup-until-it-is-not.md) has the
measurements.

## Reporting a difference

If a file from picard-rs differs from Picard's (run with the two arguments above), that is a bug,
and the most useful issue this project can receive. Please
[open an issue](https://github.com/IPNP-BIPN/picard-rs/issues/new) with the tool, the full
command line, and if possible a small input that reproduces it.

## Who

picard-rs is developed by [BIPN](https://github.com/IPNP-BIPN) (Bioinformatics for Psychiatry &
Neuroscience) at the Institute of Psychiatry and Neuroscience of Paris (IPNP), Inserm U1266,
Université Paris Cité. It is part of a wider effort with
[htsjdk-rs](https://github.com/IPNP-BIPN/htsjdk-rs) and [gatk-rs](https://github.com/IPNP-BIPN/gatk-rs)
(GATK, in progress). Commits are co-authored with the AI model that wrote them, and go through the
same CI as any other change.

## License

MIT, matching Picard. See `LICENSE`.
