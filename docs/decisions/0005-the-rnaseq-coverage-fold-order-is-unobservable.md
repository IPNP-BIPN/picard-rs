# 0005. The RnaSeq coverage-histogram fold order is unobservable at printed precision

**Status:** accepted; measured, not assumed
**Date:** 2026-07-22

## Why this exists

`CollectRnaSeqMetrics` is the first ported collector whose output depends, in principle, on the
iteration order of a Java `HashMap`. `RnaSeqMetricsCollector.computeCoverageMetrics`
(`RnaSeqMetricsCollector.java:373`, `:411`) accumulates `normalized / transcriptCount` into the
`normalized_coverage` histogram **in floating point**, iterating the picked transcripts in the
order a `HashMap<Transcript, int[]>` yields them. Floating-point addition is not associative, so a
different iteration order can change the last bit of a bin, which `FormatUtil` might then print
differently. This is the same class of hazard as the SAM header's `LinkedHashMap` (htsjdk-rs
decision 0009) and the VCF header's `TreeSet` (decision 0016): an ordering rule leaking into bytes.

Unlike those, this one is a *floating-point* fold, so whether it leaks is a question about rounding,
not about the order itself. It was measured rather than assumed.

## What was measured

The `coverage` conformance case is built to make the fold order matter: three genes, each a single
600bp coding transcript, with deliberately different coverage depth (45, 30 and 22 fragments at
different spacings), so the three transcripts' normalized-coverage arrays differ and the per-bin sum
of three distinct doubles is genuinely order-sensitive.

The picked transcripts were folded in three different orders against the same oracle golden:

| Fold order | `.rna_metrics` bytes |
|---|---|
| deterministic content order (the four coordinates then the name) | identical |
| that order reversed | identical |
| Java `HashMap` bucket order (`spread(hash) & (cap-1)`, content-based `Transcript.hashCode`) | identical |

All three produce the same file. The last-ULP differences the reorderings introduce vanish under
`FormatUtil`'s formatting. So for this corpus the fold order is **unobservable at printed
precision**, the RnaSeq analogue of htsjdk-rs decision 0020 (where the overlap set's order was shown
not to escape).

## What the port does, and what it does not claim

The collector folds the transcripts in a single **deterministic content order** (sorted by
`transcription_start, transcription_end, coding_start, coding_end, name`). This is not a claim to
reproduce Java's `HashMap` order: that order was measured not to matter, and claiming to match an
unobservable thing would be unfalsifiable. The deterministic order exists for a concrete reason of
its own: `OverlapDetector::get_all` iterates a Rust `HashMap`, whose order is randomized per run, so
without an explicit sort the output would vary run to run even though it stays byte-identical to the
oracle. The sort removes that nondeterminism.

An earlier draft did reproduce the `HashMap` bucket order exactly, complete with a no-collision
assertion. It was removed once the measurement showed it made no difference to the bytes: carrying
an elaborate "matches Java" mechanism that no test can confirm is the kind of unfalsifiable
complexity this project's method exists to avoid.

## Residual risk

The measurement covers three transcripts. A corpus with many high-coverage transcripts (the tool
keeps up to the 1000 most highly expressed) sums many more doubles per bin, where the rounding may
no longer absorb the reordering. If such a corpus ever shows the fold order reaching the bytes, the
exact `HashMap` order would have to be rebuilt and, crucially, **verified against that corpus**
rather than trusted. Until then, the deterministic order stands, with this decision recording that
it is a measured convenience and not a reproduction of Java's map.

The four `MEDIAN_*` metrics are unaffected either way: they come from `Histogram.getMedian`, which
sorts, so they are order-independent by construction.
