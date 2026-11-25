# picard-rs against Picard 3.5.0: CPU, peak RSS, wall

Every picard-rs tool (121) timed against Picard 3.5.0 at chr20 scale by `tools/speed`. The targets are CPU at most 0.5x Picard's and peak RSS at most 0.25x Picard's, with output bytes unchanged.

## How it is measured

- Inputs: the NA12878 chr20 30x BAM reduced to chr20's dictionary, the chr20 FASTA, dbSNP and gnomAD chr20, and everything `tools/speed/make_scale_inputs.py` derives from them (`rows.txt` says which row reads what).
- Both sides in one container (`picard-rs-speed:3.5.0`: the oracle's temurin digest, the same `picard.jar`), 4 CPUs, 6 GB, `RAYON_NUM_THREADS=4`, Picard with its default heap. `wait4` gives each run's CPU (user + system) and peak RSS.
- Port: origin/main c7921ed, `--profile dist`, rust 1.98.0. Host: arm64, so Picard runs without GKL (JDK Inflater and Deflater).
- This table is one cold run per side (`--reps 0`) on a host shared with other jobs, so small ratios are noisy. Each perf PR carries its own medians of warm runs.

## Where it stands

132 rows (121 tools, plus a plate of 384 barcodes and the real dbSNP and gnomAD VCFs):

| Verdict | Rows |
|---|---:|
| at target | 38 |
| CPU and RSS over | 53 |
| CPU over | 16 |
| RSS over | 11 |
| refused by one side or both | 14 |

Open PRs improve eight rows, and four of them (BamIndexStats, CalculateReadGroupChecksum, CollectJumpingLibraryMetrics, AddCommentsToBam) reach the target: #491 (CollectHsMetrics, CollectTargetedPcrMetrics), #492 (BamIndexStats, CalculateReadGroupChecksum), #493 (CollectSequencingArtifactMetrics), #494 (BuildBamIndex), #495 (CollectJumpingLibraryMetrics) and #496 (AddCommentsToBam). The last column gives each one's branch CPU and RSS, and its two ratios to Picard.

What holds the rest:

- **CPU of the tools that write a BAM**: deflate at level 5. Picard's exact bytes fix the level and the block cut, so this is the floor unless the deflater itself gets faster.
- **CPU of the tools that read one**: about 45% of their CPU is htsjdk-rs's `BamRecord::decode` (tags and bases), measured on CollectDuplicateMetrics and CollectSequencingArtifactMetrics. That code belongs to htsjdk-rs.
- **RSS of the streaming tools (100-220 MiB on 4 cores)**: the decode-ahead batches and the blocks in flight. The window has not been tuned against the RSS target yet.
- **MarkDuplicates and its kin (2.7-6 GB)**: these move to `crates/picard-markduplicates` in #484 and are worked on after it merges. MarkDuplicatesWithMateCigar is killed at the 6 GB limit, where Picard uses 375 MiB.
- **Rows marked "row to fix"**: both sides refuse the input, so the row needs another input, not a fix to the port.
- **Divergence found**: CollectSamErrorMetrics refuses a `.vcf.gz` VCF (`stream did not contain valid UTF-8`) that Picard reads.

## Every row, worst RSS ratio first

CPU s is user + system time of every thread, RSS MiB the peak resident set, Picard against the port (origin/main) on the same input; for a row an open PR changes, the last column gives the port after it, CPU and RSS, and the two ratios to Picard.

| Row | Picard CPU s | port CPU s | CPU ratio | Picard RSS MiB | port RSS MiB | RSS ratio | Verdict | Why, or what is open |
|---|---:|---:|---:|---:|---:|---:|---|---|
| MarkDuplicatesWithMateCigar | 126.6 | 181.2 | 1.43 | 375 | 6120 | 16.33 | port failed | port killed at the 6 GB limit (Picard 375 MiB); MarkDuplicates* move to picard-markduplicates in #484, optimised after it |
| DownsampleSam | 60.0 | 72.2 | 1.20 | 262 | 1735 | 6.61 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| PositionBasedDownsampleSam | 80.8 | 95.5 | 1.18 | 383 | 1715 | 4.47 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| FastqToSam | 144.1 | 99.9 | 0.69 | 1339 | 5929 | 4.43 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| UmiAwareMarkDuplicatesWithMateCigar | 125.7 | 27.5 | 0.22 | 1702 | 5850 | 3.44 | port failed | both refuse: the reads carry no RX UMI; no UMI data at scale |
| EstimateLibraryComplexity | 58.6 | 40.1 | 0.68 | 1285 | 4399 | 3.42 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| SimpleMarkDuplicatesWithMateCigar | 295.0 | 100.0 | 0.34 | 1739 | 5447 | 3.13 | RSS over | RSS: not yet diagnosed; after #484 moves it to picard-markduplicates |
| CheckDuplicateMarking | 69.3 | 25.2 | 0.36 | 1647 | 3317 | 2.01 | port failed | both exit non-zero on this input; row to fix |
| CollectJumpingLibraryMetrics | 50.0 | 27.3 | 0.54 | 506 | 856 | 1.69 | CPU / RSS over | #495: 20.4 s / 105 MiB (0.41 / 0.21) |
| BuildBamIndex | 19.6 | 25.3 | 1.29 | 518 | 867 | 1.67 | CPU / RSS over | #494: 23.6 s / 21 MiB (1.20 / 0.04) |
| MarkDuplicates | 161.4 | 110.2 | 0.68 | 1742 | 2762 | 1.59 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed; after #484 moves it to picard-markduplicates |
| VcfToIntervalList@gnomad | 6.4 | 6.9 | 1.08 | 256 | 337 | 1.32 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| CollectIndependentReplicateMetrics | 28.8 | 21.6 | 0.75 | 1690 | 1875 | 1.11 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| BaitDesigner | 3.2 | 0.3 | 0.08 | 234 | 214 | 0.92 | RSS over | RSS: not yet diagnosed |
| CreateSequenceDictionary | 1.0 | 0.1 | 0.14 | 219 | 196 | 0.90 | RSS over | RSS: not yet diagnosed |
| NormalizeFasta | 1.3 | 0.1 | 0.09 | 213 | 182 | 0.86 | RSS over | RSS: not yet diagnosed |
| NonNFastaSize | 1.0 | 0.1 | 0.06 | 219 | 182 | 0.83 | RSS over | RSS: not yet diagnosed |
| ScatterIntervalsByNs | 2.0 | 0.1 | 0.04 | 219 | 180 | 0.82 | RSS over | RSS: not yet diagnosed |
| CollectGcBiasMetrics | 32.5 | 30.2 | 0.93 | 452 | 355 | 0.78 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| MakeVcfSampleNameMap | 0.6 | 0.0 | 0.03 | 157 | 110 | 0.70 | RSS over | RSS: not yet diagnosed |
| ConvertHaplotypeDatabaseToVcf | 2.3 | 0.1 | 0.03 | 405 | 278 | 0.69 | RSS over | RSS: not yet diagnosed |
| SetNmAndUqTags | 90.7 | 75.5 | 0.83 | 334 | 225 | 0.67 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| SetNmMdAndUqTags | 89.7 | 75.6 | 0.84 | 338 | 223 | 0.66 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| MergeSamFiles | 136.2 | 128.5 | 0.94 | 259 | 165 | 0.64 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| CollectMultipleMetrics | 106.6 | 76.2 | 0.71 | 370 | 228 | 0.62 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectSequencingArtifactMetrics | 69.6 | 113.4 | 1.63 | 461 | 281 | 0.61 | CPU / RSS over | #493: 44.7 s / 276 MiB (0.64 / 0.60) |
| CrosscheckFingerprints | 472.0 | 47.6 | 0.10 | 1666 | 1014 | 0.61 | RSS over | RSS: not yet diagnosed |
| CollectRrbsMetrics | 71.3 | 46.6 | 0.65 | 741 | 417 | 0.56 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectWgsMetricsWithNonZeroCoverage | 121.7 | 65.7 | 0.54 | 380 | 213 | 0.56 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectVariantCallingMetrics@dbsnp | 0.8 | 12.0 | 14.61 | 178 | 99 | 0.56 | Picard failed | Picard refuses the sites-only VCF; row to fix |
| CollectVariantCallingMetrics@gnomad | 1.1 | 11.4 | 10.75 | 182 | 99 | 0.54 | Picard failed | Picard refuses the sites-only VCF; row to fix |
| AddOATag | 157.7 | 116.9 | 0.74 | 303 | 150 | 0.50 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| RevertSam | 184.2 | 108.6 | 0.59 | 1399 | 681 | 0.49 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| SplitSamByLibrary | 72.0 | 69.1 | 0.96 | 233 | 113 | 0.49 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CleanSam | 104.8 | 85.6 | 0.82 | 285 | 138 | 0.48 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| ViewSam | 28.3 | 58.2 | 2.05 | 387 | 182 | 0.47 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| MergeVcfs | 4.9 | 2.8 | 0.57 | 364 | 163 | 0.45 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| VcfToIntervalList | 3.0 | 1.6 | 0.55 | 399 | 176 | 0.44 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| FixMateInformation | 164.0 | 91.6 | 0.56 | 229 | 100 | 0.44 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| MarkIlluminaAdapters | 149.7 | 92.2 | 0.62 | 216 | 93 | 0.43 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| BamToBfq | 85.5 | 79.0 | 0.92 | 254 | 109 | 0.43 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| FindMendelianViolations | 0.7 | 1.6 | 2.38 | 168 | 71 | 0.42 | Picard failed | Picard refuses the row; row to fix |
| IdentifyContaminant | 141.8 | 27.0 | 0.19 | 1207 | 508 | 0.42 | RSS over | RSS: not yet diagnosed |
| CheckFingerprint | 33.1 | 21.4 | 0.65 | 1246 | 522 | 0.42 | port failed | both answer with no result for this sample; row to fix |
| SplitVcfs | 0.8 | 1.5 | 1.80 | 169 | 71 | 0.42 | port failed | both refuse: dbSNP has MNPs, which SplitVcfs refuses; row to fix |
| CalculateFingerprintMetrics | 74.7 | 62.8 | 0.84 | 1076 | 450 | 0.42 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CrosscheckReadGroupFingerprints | 67.2 | 35.6 | 0.53 | 1099 | 459 | 0.42 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| MergeBamAlignment | 216.7 | 139.1 | 0.64 | 1690 | 705 | 0.42 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| CollectOxoGMetrics | 110.4 | 53.0 | 0.48 | 531 | 215 | 0.41 | RSS over | RSS: not yet diagnosed |
| CollectSamErrorMetrics | 3406.9 | 0.1 | 0.00 | 696 | 267 | 0.38 | port failed | Picard over the 1 h timeout with the default aggregations; the port refuses a .vcf.gz VCF (divergence to fix) |
| GenotypeConcordance | 9.9 | 6.0 | 0.61 | 343 | 131 | 0.38 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| AddOrReplaceReadGroups | 119.0 | 91.4 | 0.77 | 301 | 113 | 0.37 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| CombineGenotypingArrayVcfs | 5.8 | 4.6 | 0.79 | 353 | 131 | 0.37 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| CollectQualityYieldMetrics | 31.3 | 23.5 | 0.75 | 287 | 104 | 0.36 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| ExtractFingerprint | 25.7 | 20.3 | 0.79 | 1413 | 501 | 0.35 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| SortSam | 152.1 | 85.0 | 0.56 | 1730 | 612 | 0.35 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| SplitSamByNumberOfReads | 85.3 | 83.5 | 0.98 | 329 | 116 | 0.35 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| FilterSamReads | 0.8 | 0.1 | 0.13 | 162 | 55 | 0.34 | port failed | both refuse FILTER=includeAligned on this input; row to fix |
| CollectAlignmentSummaryMetrics | 73.9 | 44.7 | 0.61 | 636 | 212 | 0.33 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectHsMetrics | 52.8 | 231.9 | 4.39 | 616 | 203 | 0.33 | CPU / RSS over | #491: 35.1 s / 212 MiB (0.66 / 0.34) |
| SamToFastqWithTags | 0.7 | 0.4 | 0.62 | 163 | 54 | 0.33 | port failed | both refuse: the reads already carry OA; row to fix |
| VcfFormatConverter@gnomad | 14.3 | 10.4 | 0.73 | 246 | 81 | 0.33 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| MakeSitesOnlyVcf@gnomad | 16.4 | 10.5 | 0.64 | 242 | 79 | 0.33 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| CollectRnaSeqMetrics | 40.3 | 35.2 | 0.87 | 344 | 111 | 0.32 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| ValidateSamFile | 26.7 | 21.4 | 0.80 | 362 | 111 | 0.31 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| QualityScoreDistribution | 23.8 | 31.6 | 1.33 | 345 | 105 | 0.30 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| SamFormatConverter | 28.3 | 40.2 | 1.42 | 349 | 106 | 0.30 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| SplitVcfs@gnomad | 13.9 | 10.2 | 0.73 | 276 | 83 | 0.30 | CPU / RSS over | CPU: not yet diagnosed; RSS: not yet diagnosed |
| CollectBaseDistributionByCycle | 23.7 | 28.2 | 1.19 | 369 | 110 | 0.30 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectDuplicateMetrics | 16.6 | 22.6 | 1.36 | 369 | 107 | 0.29 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectRawWgsMetrics | 125.8 | 80.7 | 0.64 | 730 | 211 | 0.29 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| SamToFastq | 36.8 | 49.1 | 1.34 | 284 | 82 | 0.29 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| CollectTargetedPcrMetrics | 91.5 | 227.0 | 2.48 | 695 | 201 | 0.29 | CPU / RSS over | #491: 32.8 s / 198 MiB (0.36 / 0.28) |
| CollectQualityYieldMetricsFlow | 0.8 | 0.2 | 0.23 | 160 | 45 | 0.28 | port failed | both refuse: no flow-based reads at scale |
| CollectQualityYieldMetricsSNVQ | 0.8 | 0.3 | 0.35 | 162 | 44 | 0.27 | port failed | both refuse: no SNVQ tags at scale |
| CollectWgsMetrics | 123.0 | 74.4 | 0.60 | 749 | 199 | 0.26 | CPU / RSS over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics); RSS: not yet diagnosed |
| RevertOriginalBaseQualitiesAndAddMateCigar | 0.8 | 0.5 | 0.69 | 164 | 42 | 0.25 | CPU / RSS over | CPU: deflate at level 5 for the same bytes (measured on AddOATag); RSS: not yet diagnosed |
| MeanQualityByCycle | 21.2 | 27.9 | 1.31 | 422 | 102 | 0.24 | CPU over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics) |
| CollectVariantCallingMetrics | 4.2 | 2.1 | 0.50 | 375 | 91 | 0.24 | at target |  |
| ExtractIlluminaBarcodes@384 | 1.7 | 0.3 | 0.16 | 205 | 48 | 0.24 | at target |  |
| CollectHiSeqXPfFailMetrics | 1.9 | 0.2 | 0.09 | 222 | 50 | 0.23 | at target |  |
| CompareSAMs | 35.3 | 54.7 | 1.55 | 441 | 98 | 0.22 | CPU over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics) |
| SortVcf@dbsnp | 34.6 | 11.5 | 0.33 | 1762 | 391 | 0.22 | at target |  |
| UpdateVcfSequenceDictionary@gnomad | 14.0 | 10.0 | 0.72 | 283 | 63 | 0.22 | CPU over | CPU: not yet diagnosed |
| CollectArraysVariantCallingMetrics | 3.1 | 1.4 | 0.44 | 433 | 95 | 0.22 | at target |  |
| MakeSitesOnlyVcf | 4.0 | 1.8 | 0.46 | 380 | 83 | 0.22 | at target |  |
| VcfToAdpc | 3.5 | 2.7 | 0.78 | 437 | 94 | 0.22 | CPU over | CPU: not yet diagnosed |
| MergePedIntoVcf | 5.2 | 3.5 | 0.66 | 595 | 123 | 0.21 | CPU over | CPU: not yet diagnosed |
| FilterVcf | 6.1 | 2.8 | 0.47 | 389 | 79 | 0.20 | at target |  |
| VcfFormatConverter | 4.9 | 2.6 | 0.54 | 395 | 77 | 0.19 | CPU over | CPU: not yet diagnosed |
| RenameSampleInVcf | 2.5 | 1.6 | 0.62 | 376 | 71 | 0.19 | CPU over | CPU: not yet diagnosed |
| CollectUmiPrevalenceMetrics | 19.9 | 31.6 | 1.59 | 377 | 70 | 0.19 | CPU over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics) |
| GatherVcfs | 5.2 | 2.7 | 0.52 | 388 | 71 | 0.18 | CPU over | CPU: not yet diagnosed |
| FixVcfHeader | 6.8 | 4.0 | 0.59 | 401 | 71 | 0.18 | CPU over | CPU: not yet diagnosed |
| ExtractIlluminaBarcodes | 2.1 | 0.5 | 0.22 | 229 | 40 | 0.18 | at target |  |
| LiftOverHaplotypeMap | 2.0 | 0.0 | 0.02 | 343 | 60 | 0.18 | at target |  |
| UpdateVcfSequenceDictionary | 4.8 | 2.6 | 0.55 | 404 | 71 | 0.18 | CPU over | CPU: not yet diagnosed |
| SortVcf | 13.8 | 3.8 | 0.28 | 1772 | 309 | 0.17 | at target |  |
| CollectIlluminaBasecallingMetrics | 1.7 | 0.1 | 0.05 | 224 | 36 | 0.16 | at target |  |
| ReplaceSamHeader | 1.2 | 0.3 | 0.27 | 162 | 24 | 0.15 | at target |  |
| CollectIlluminaLaneMetrics | 1.5 | 0.0 | 0.01 | 195 | 28 | 0.15 | at target |  |
| CollectInsertSizeMetrics | 36.5 | 44.0 | 1.21 | 642 | 90 | 0.14 | CPU over | CPU: htsjdk-rs record decode (measured on CollectDuplicateMetrics) |
| CheckIlluminaDirectory | 1.0 | 0.0 | 0.01 | 161 | 22 | 0.14 | at target |  |
| CheckTerminatorBlock | 0.8 | 0.0 | 0.01 | 155 | 21 | 0.14 | at target |  |
| GatherBamFiles | 1.1 | 0.3 | 0.26 | 154 | 20 | 0.13 | at target |  |
| IlluminaBasecallsToSam | 8.7 | 2.6 | 0.30 | 613 | 81 | 0.13 | at target |  |
| AddCommentsToBam | 1.2 | 6.0 | 4.99 | 161 | 21 | 0.13 | CPU over | #496: 0.3 s / 25 MiB (0.25 / 0.16) |
| CalculateReadGroupChecksum | 0.8 | 4.9 | 5.91 | 156 | 20 | 0.13 | CPU over | #492: 0.0 s / 22 MiB (0.00 / 0.14) |
| CompareMetrics | 0.8 | 0.0 | 0.01 | 161 | 21 | 0.13 | at target |  |
| AccumulateVariantCallingMetrics | 0.7 | 0.0 | 0.01 | 162 | 21 | 0.13 | at target |  |
| ClusterCrosscheckMetrics | 0.9 | 0.0 | 0.02 | 171 | 22 | 0.13 | at target |  |
| BamIndexStats | 0.9 | 5.6 | 5.99 | 162 | 20 | 0.13 | CPU over | #492: 0.0 s / 20 MiB (0.00 / 0.12) |
| LiftoverVcf@gnomad | 56.3 | 17.2 | 0.31 | 1759 | 221 | 0.13 | at target |  |
| ConvertSequencingArtifactToOxoG | 1.0 | 0.0 | 0.01 | 163 | 20 | 0.12 | at target |  |
| LiftoverVcf | 14.0 | 4.1 | 0.29 | 1705 | 207 | 0.12 | at target |  |
| CreateVerifyIDIntensityContaminationMetricsFile | 0.8 | 0.0 | 0.02 | 187 | 22 | 0.12 | at target |  |
| CreateBafRegressMetricsFile | 1.0 | 0.0 | 0.02 | 208 | 23 | 0.11 | at target |  |
| AccumulateQualityYieldMetrics | 0.8 | 0.0 | 0.01 | 193 | 20 | 0.10 | at target |  |
| CreateExtendedIlluminaManifest | 9.8 | 0.4 | 0.04 | 1022 | 90 | 0.09 | at target |  |
| IlluminaBasecallsToFastq | 10.3 | 1.6 | 0.16 | 696 | 58 | 0.08 | at target |  |
| BedToIntervalList | 1.6 | 0.1 | 0.04 | 311 | 24 | 0.08 | at target |  |
| IlluminaBasecallsToFastq@384 | 15.5 | 0.7 | 0.05 | 804 | 61 | 0.08 | at target |  |
| BpmToNormalizationManifestCsv | 3.3 | 0.3 | 0.08 | 710 | 50 | 0.07 | at target |  |
| ReorderSam | 141.7 | 75.2 | 0.53 | 1746 | 113 | 0.06 | CPU over | CPU: deflate at level 5 for the same bytes (measured on AddOATag) |
| LiftOverIntervalList | 4.2 | 0.2 | 0.04 | 399 | 24 | 0.06 | at target |  |
| IntervalListTools | 5.4 | 0.1 | 0.02 | 439 | 24 | 0.06 | at target |  |
| SortGff | 5.0 | 1.3 | 0.25 | 876 | 46 | 0.05 | at target |  |
| CompareGtcFiles | 2.0 | 0.5 | 0.25 | 1118 | 56 | 0.05 | port failed | exit 1 is the answer (the two GTCs differ), on both sides |
| IntervalListToBed | 3.4 | 0.1 | 0.03 | 535 | 24 | 0.05 | at target |  |
| GtcToVcf | 10.2 | 1.6 | 0.16 | 1688 | 64 | 0.04 | at target |  |
| ExtractSequences | 4.3 | 0.6 | 0.14 | 587 | 21 | 0.04 | at target |  |
| FifoBuffer | 0.7 | 0.2 | 0.23 | 665 | 21 | 0.03 | at target |  |
