# Changelog

All notable changes to this repository, newest first. The workspace carries one version
(`[workspace.package]` in `Cargo.toml`). The format is
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); 0.x is what an unreleased port is.

## [Unreleased]

The reference stays Picard 3.5.0, on htsjdk 5.0.0 through htsjdk-rs `2116a33`.

### Added

- Every measured tool is fuzzed against Picard in CI, in eight shards, frozen into
  `tools/fuzz/measured.json` (#481).
- SortSam timed against Picard on real x86-64 in CI.
- Release archives carry `NOTICE`, `THIRD-PARTY-NOTICES` and the licence of every linked crate
  (`THIRD-PARTY-LICENSES.txt`); the release workflow drafts the release itself, and its binaries are
  run on QEMU's baseline `qemu64` CPU before they are packaged.

### Changed

- SortSam reads and sorts a clean BAM on every core; a streamed BAM's records are decoded on every
  core, and each BAM is read once.
- Picard's legacy parser's refusals are refused here too; the covering arrays stay at 100%, and
  MergeBamAlignment and IlluminaBasecallsToFastq are frozen at 1.000.

### Fixed

- Damaged Illumina, BPM and tile-metrics files, and ignored write errors, fail as in Picard.
- `.clocs` and `_pos.txt` cluster positions are read.
- Audit findings: the write-time BAI, MarkDuplicates barcodes and UMIs, SortVcf inputs, interval
  names, BED parsing.

## [0.2.0] - 2026-10-21

Every tool streams and is faster than Picard 3.5.0 on all 121. See the
[release notes](https://github.com/IPNP-BIPN/picard-rs/releases/tag/v0.2.0).

## [0.1.0] - 2026-10-16

The first release. See the [release notes](https://github.com/IPNP-BIPN/picard-rs/releases/tag/v0.1.0).
