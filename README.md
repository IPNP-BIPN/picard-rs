# picard-rs

Native Rust reimplementation of [Picard](https://github.com/broadinstitute/picard), targeting
**byte-identical** output against a pinned reference build. Work in progress.

> **This is not the official Picard.** It is an independent reimplementation, not affiliated
> with or endorsed by the Broad Institute.

## Reference version

Ported from Picard `3.4.0`, commit `6c3f23bc2e0d229d75e9f9b04200396bcd067526`. This is the
exact version pinned by GATK 4.6.2.0 (`picardVersion = '3.4.0'` in its `build.gradle`), so the
two ports are coherent by construction rather than by coincidence.

## Method

Every feature branch ports a **named symbol** of the pinned reference source: a specific class
or method, read from the pinned clone and translated. Behavior is never reconstructed from
documentation, papers, or memory. Each branch records the source file and symbol it ports.

Work is organized three levels deep:

```
main                                    shared infrastructure only
└── tool/picard-<toolname>               one per tool
    └── feat/picard-<toolname>-<feature> one per ported symbol
```

A feature merges into its tool branch when it is byte-identical for its scope. A tool merges
into `main` when the whole tool is byte-identical across its corpus.

## Status

| Tool | Status |
|---|---|
| `CollectQualityYieldMetrics` | **byte-identical** |
| `MeanQualityByCycle` | **byte-identical** |
| `CollectBaseDistributionByCycle` | **byte-identical** |
| `CollectInsertSizeMetrics` | **byte-identical**, first member of the calibration pair |
| `CollectAlignmentSummaryMetrics` | **byte-identical**, 22 cases, decision 0003 |

## Decisions

| # | Title |
|---|---|
| [0001](docs/decisions/0001-the-metrics-archetype-is-not-homogeneous.md) | The metrics archetype is not homogeneous |
| [0002](docs/decisions/0002-the-first-within-stratum-delta.md) | The first within-stratum delta, smaller than the archetype story assumed |
| [0003](docs/decisions/0003-the-delta-at-the-large-end-is-negative.md) | The within-stratum delta at the large end is negative |

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

## Part of a three-repository program

| Repo | Ports | Depends on |
|---|---|---|
| `htsjdk-rs` | htsjdk 4.2.0 | (none) |
| `picard-rs` | Picard 3.4.0 | `htsjdk-rs` |
| `gatk-rs` | GATK 4.6.2.0 | `picard-rs`, `htsjdk-rs` |

## License

MIT, matching Picard. See `LICENSE`.
