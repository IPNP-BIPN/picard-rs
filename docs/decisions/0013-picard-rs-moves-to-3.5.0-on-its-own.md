# 0013. picard-rs moves to Picard 3.5.0 on its own

**Status:** accepted
**Date:** 2026-10-15

## The decision

The reference is Picard `3.5.0` (commit `a53d7ebd975606b79ca0eaddb7bf93f555ec4db4`), whose
release jar bundles htsjdk `5.0.0`. gatk-rs stays at GATK 4.6.2.0, which pins Picard 3.4.0, so
the three repositories no longer name one set of pins. That reverses the order
[IPNP-BIPN/gatk-rs#810](https://github.com/IPNP-BIPN/gatk-rs/issues/810) set (htsjdk-rs first,
then this repository, then gatk-rs), and the reason is that the bill was measured first:
[decision 0012](0012-what-the-move-to-3.5.0-costs-the-ported-surface.md) bounded what 3.5.0 costs
the ported surface, and nothing in it needs htsjdk-rs to move. What picard-rs needs from htsjdk
5.0.0 is three behaviours, and each is ported here, beside the Picard code that reaches it:

* `SamPairUtil.getPairOrientation` reads the reverse mate's five-prime position from `MC` when it
  is set, and calls a tie `FR` (`insert_size::pair_orientation`, and the three copies the
  alignment mergers carry);
* `SAMRecord.toString()` is the whole SAM line, so an exception that names a record names all
  eleven fields and every tag (`view_sam::sam_string`);
* `getSAMString()` no longer ends in a newline, which only the dumps and `ViewSam` see.

htsjdk-rs's own move ([IPNP-BIPN/htsjdk-rs#182](https://github.com/IPNP-BIPN/htsjdk-rs/issues/182))
is still the large piece, CRAM 3.1 writing and the libdeflate default, and is unchanged by this.

## How it was measured

The oracle image was moved first, alone, so the first CI run compared Picard 3.5.0's output with
every 3.4.0 golden before any port code changed. An oracle-backed case that no longer matches now
publishes its regenerated dump at the golden's own path, so a re-frozen golden is the x86-64 job's
artefact, compressed the way the others are, and never a local run (decision 0008).

That run failed 27 of 107 suites. 25 of them were one change: the dumps write records through
`getSAMString()`, which lost its trailing newline. The dumps now add it back exactly as Picard
3.5.0's `ViewSam` does, which is the reference's own answer to the same change, and those 25
goldens came back identical.

## What moved, and why

Seven goldens were re-frozen, each with its difference explained:

| golden | the difference | the cause |
|---|---|---|
| `merge_bam_alignment_hits` | the `EarliestFragment` refusal names the record by its SAM line | htsjdk 5.0.0 `SAMRecord.toString` |
| `mate_cigar_duplicates` | `Mate CIGAR (Tag MC) not found:` names the record by its SAM line | the same |
| `position_based_downsample_sam` | the previous `@PG` says `VN=3.5.0` | the version |
| `crosscheck_fingerprints` | under `OUTPUT_ERRORS_ONLY`, a run with only expected results exits 0 with an empty table instead of refusing itself | Picard #2018: the non-zero-LOD check covers every comparison, not the rows written |
| `collect_wgs_metrics_with_non_zero_coverage` | a run without `CHART_OUTPUT` writes its metrics instead of being refused | the chart argument became optional |
| `pf_fail_metrics` | the unmatched data types are listed `QualityScores, BaseCalls, PF` | see below |
| `collect_multiple_metrics` | two of five `chart-stability` lines say `same` | see below |

Two of those are not behaviour, and are recorded as what they are:

* **`pf_fail_metrics`.** `findUnmatchedTypes` returns a `HashSet<IlluminaDataType>`, and an enum's
  hash is its identity hash, so the order is whatever this JVM's allocation sequence gives. A
  different jar moved it. It is stable for one image, which is all a golden claims.
* **`chart-stability`.** The line says whether two runs of one fixture write identical PDFs, and
  a PDF carries its creation time to the second. Under 3.4.0 the two runs never shared a second;
  under 3.5.0 some do. The line measures the clock, so `multiplemetrics` now records it without
  comparing it.

Behaviour the port took on without a golden moving, because the corpora do not reach it yet:

* `CollectInsertSizeMetrics` skips a pair whose mates are on different contigs, and
  `CollectRnaSeqMetrics` requires one contig before calling a pair `FR`;
* `FingerprintChecker` takes, at each SNP of an indexed VCF, the first record
  `subsetVCToMatchSnp` accepts;
* `RevertSam` restores hard clips only when both `XB` and `XQ` are present, which the port already
  did;
* the ten chart arguments are optional and refused only in the GATK Lite image, and
  `CollectMultipleMetrics` and `CollectRrbsMetrics` take `DO_NOT_CREATE_PLOTS`.

## The argument surface

The arrays are generated in gatk-rs from an inventory of GATK 4.6.2.0's documentation, which
describes Picard 3.4.0. Picard 3.5.0's surface differs from it in exactly the twelve places
decision 0012 lists (ten chart arguments optional, two `DO_NOT_CREATE_PLOTS` added), so the ten
affected arrays were regenerated from that inventory with those twelve edits, and the chart
fixtures gained the absent value so the arrays run both with and without a chart. The other
arrays are unchanged, and `SortSam`'s was already out of step with its generator before this move.

## What this leaves

* Module headers that say "Ported from Picard 3.4.0" stay: the symbols they name did not change
  between the tags (the 29 files that did are in decision 0012), so the claim is still true.
* gatk-rs reaches Picard through this repository. When it moves to GATK 4.7.0.0 the pins meet
  again; until then the README states them separately.
