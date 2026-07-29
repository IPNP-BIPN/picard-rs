# Fanning the port out across parallel agents

The programme is 311 tools. The porting loop is now mechanical enough to run many at once, and
this is the protocol that keeps that from turning into a merge queue full of conflicts. It applies
to `htsjdk-rs`, `picard-rs` and `gatk-rs` alike.

## The loop, per tool

1. Read the symbol in the pinned clone. Never reconstruct behaviour from documentation, papers, or
   memory; the branch records the file and symbol it ports.
2. Write the Rust module and the `*Dump.java` harness.
3. Add the suite to `tools/conformance/manifest.json`: goldens, compare mode, and every rule that
   canonicalizes a field away with the reason it is there.
4. Regenerate the workflow (`tools/conformance/generate_ci.py`) and commit it: the `guard` job
   fails if the two disagree.
5. Generate the covering array in `gatk-rs`, run it (`tools/coverage/run_array.py`), and declare
   any constraint the reference's own rejections reveal.
6. One PR per slice. Merge only on green.

## Branches and worktrees

```
main                                    shared infrastructure only
└── tool/picard-<toolname>               one per tool
    └── feat/picard-<toolname>-<feature> one per ported symbol
```

One worktree per tool, so N agents do not share an index:

```sh
git worktree add ../picard-rs-sortsam -b tool/picard-sortsam
```

## The rule that makes parallelism work

**Shared infrastructure changes only on `main`, one at a time.**

Shared means: the crates every tool depends on (`htsjdk-*`), `tools/conformance/*`,
`tools/coverage/*`, the oracle image, and the generated workflow. An agent that needs a new shared
symbol opens a separate infrastructure PR first, and its tool PR depends on that landing.

Without this, N agents each add a function to `htsjdk-bam`, each regenerate `ci.yml`, and every
one of the N pull requests conflicts with the other N-1 in a file none of them is really about.
The tool modules themselves are disjoint by construction, so once the shared layer is serialized
the fan-out is genuinely parallel.

The generated workflow makes this cheap to check: a tool PR that touches only its own module, its
harness, and one manifest entry produces a `ci.yml` diff of exactly one matrix row.

## What an agent may decide alone, and what it may not

**Alone:** how to structure the Rust module, which symbols to port first, what the fixture corpus
needs, whether a divergence is a fixture defect or a real one.

**Not alone**, because each one weakens or redefines a claim the whole programme rests on:

- **Adding a canonicalization rule.** It is how a bit-identity claim quietly narrows. It is
  declared in the manifest with a reason, and reviewed as such.
- **Committing a golden that CI did not produce.** Decision 0008 is what that costs: sixteen
  goldens turned out to have been produced on a laptop, and nine still carry the macOS temp paths
  that prove it. A suite with no golden yet is declared `golden-pending`; CI publishes the
  candidate and committing that artefact is a deliberate act.
- **Declaring a constraint on a covering array.** It removes rows from the coverage claim.
  Decision 0009 shows the failure mode: the first shape of a constraint was wrong twice, and only
  running it against the oracle said so.
- **Quarantining a field.** It downgrades the tool from bit-identical to bio-identical, and the
  measured divergence rate has to be attached.

## Budget

Every conformance job restores the ~1 GB oracle image once. Suites are grouped
(`ci.suites_per_job` in the manifest) so the restore is amortized rather than paid per suite. Raise
the parallelism only after measuring what the matrix costs on real x86-64 runners: the constraint
is runner minutes, not agents.
