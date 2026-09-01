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

"Byte-identical" here means byte-identical **on the paths its conformance suite covers**, which is
the default path plus the cases listed. Over the tool's whole argument surface it is measured
separately and is much lower: see "Argument coverage" below and decision 0009.

## Throughput

Measured, not claimed. `tools/benchmark/run.sh` times the port against Picard on the same input,
inside the same container, and asserts byte-identity of the two outputs in the same run.

At 2,000,000 reads on real x86-64: **1.46x** wall clock, **1.36x** net of JVM startup, output
byte-identical.

The same code measured 5.2x on a 200,000-read input and 2.0x under x86-64 emulation. Each time
the measurement got more careful the number got smaller, which is the direction benchmarks move
when they start optimistic;
[decision 0004](docs/decisions/0004-the-speedup-is-mostly-jvm-startup-until-it-is-not.md) records
all three and why. No optimisation work has been done, so 1.46x is a floor.

## Decisions

| # | Title |
|---|---|
| [0001](docs/decisions/0001-the-metrics-archetype-is-not-homogeneous.md) | The metrics archetype is not homogeneous |
| [0002](docs/decisions/0002-the-first-within-stratum-delta.md) | The first within-stratum delta, smaller than the archetype story assumed |
| [0003](docs/decisions/0003-the-delta-at-the-large-end-is-negative.md) | The within-stratum delta at the large end is negative |
| [0004](docs/decisions/0004-the-speedup-is-mostly-jvm-startup-until-it-is-not.md) | The speedup is mostly JVM startup, until the input is big enough that it is not |
| [0005](docs/decisions/0005-the-rnaseq-coverage-fold-order-is-unobservable.md) | The RnaSeq coverage fold order is unobservable |
| [0006](docs/decisions/0006-per-record-transforms-parallelize-without-touching-the-bytes.md) | Per-record transforms parallelize without touching the bytes |
| [0007](docs/decisions/0007-multicore-helps-the-transform-not-the-io.md) | Multicore helps the transform, not the I/O |
| [0008](docs/decisions/0008-sixteen-goldens-had-never-been-produced-by-the-oracle.md) | Sixteen goldens had never been produced by the oracle |
| [0009](docs/decisions/0009-the-first-covering-array-run-measures-the-argument-surface.md) | The first covering-array run measures the argument surface, and it is 0% |
| [0010](docs/decisions/0010-the-fuzzer-found-a-divergence-the-conformance-suite-cannot.md) | The fuzzer found a divergence the conformance suite cannot |

## Conformance suites

Every suite is declared once in [`tools/conformance/manifest.json`](tools/conformance/manifest.json):
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
    --port target/release/collect-alignment-summary-metrics
```

`tools/coverage/MakeFixtures.java` builds the deterministic corpus the rows run against (a covering
array cannot invent a file path), `fixtures.json` declares which values this repository is willing
to pass and which arguments it refuses to vary, and the run writes a corpus in the same dump format
the conformance suites use. Decision 0009 records what the first two runs measured.

Which tools the CI job runs is declared once, in the `coverage` block of
[`tools/conformance/manifest.json`](tools/conformance/manifest.json), and the job is generated from
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
    --port target/release/collect-alignment-summary-metrics
```

Sixty mutations reach about 10% more of the reference than the whole pairwise array does, and the
first two sessions found a numeric divergence the conformance corpus cannot contain: decision 0010.

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

## Relationship to `fulcrumgenomics/riker`

[riker](https://github.com/fulcrumgenomics/riker) is an independent, MIT-licensed Rust
reimplementation of these same Picard QC tools, from the maintainers of Picard and htsjdk. It is
the closest existing work to this repository, so the distinction matters and is worth stating
plainly.

**riker targets functional equivalence; this repository targets byte equivalence.** riker's own
README says it "is not intended to be a drop-in replacement for Picard": lowercase `snake_case`
headers, no metadata lines, `frac_` for `pct_`, and "bug fixes that yield slightly different
outputs". So riker is the better tool to use; this is a byte-for-byte reproduction of the
existing one, bugs included, so that GATK and Picard pipelines can be reproduced exactly.

That makes riker **a source of divergence candidates, never a source to port from** — copying it
would import its deliberate deviations, and the licence being compatible does not make it correct.
Its `ERRATA.md` is a curated list of exactly where a careful reimplementer differs from Picard,
and every entry is a place this port must *not* differ. Two are already pinned as conformance
cases (`riker_mean_aligned_dilution`, `riker_improper_pair_unmapped_mate` in the alignment-summary
suite), measured against the reference rather than trusted. Where riker's reading of Picard and
this port's agree, that is an independent cross-check of the reading; where riker's errata is
silent on something byte comparison catches — Picard's alignment-block cycle binning, recorded in
`docs/decisions/0003` — that is evidence for the method rather than against riker.

## Commit attribution

Commits are co-authored with the model that wrote them. On 2026-07-21 the history of all three
repositories was rewritten to add that trailer uniformly, at the maintainer's request, changing
every commit SHA. The **current** htsjdk-rs pin was moved to the rewritten commit, but pins in
this repository's *historical* commits name pre-rewrite htsjdk-rs SHAs that no longer exist, so a
checkout of an old picard-rs commit can no longer fetch its exact dependency. Historical builds
before that date are therefore no longer bit-reproducible; current and future ones are. This was
a deliberate trade.

## License

MIT, matching Picard. See `LICENSE`.
