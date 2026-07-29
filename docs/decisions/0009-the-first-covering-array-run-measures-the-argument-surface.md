# 0009. The first covering-array run measures the argument surface, and it is 0%

**Status:** accepted
**Date:** 2026-07-29

## What was run

`tools/coverage/run_array.py` turns each row of a t=2 covering array (generated in gatk-rs from the
pinned inventory) into a command line, runs it against the oracle in the pinned container, and
optionally against the port. The first two tools it was pointed at:

| tool | rows | accepted by the reference | rejected | distinct outputs | rejection classes |
|---|---:|---:|---:|---:|---:|
| `CollectQualityYieldMetrics` | 11 | 3 | 8 | 3 | 2 |
| `CollectAlignmentSummaryMetrics` | 16 | 12 | 4 | 9 | 2 |

Nine distinct outputs from twelve accepted rows means the arguments are **observable on this
corpus**: a covering array over arguments that all produce the same file would be green and
worthless, and that had to be measured rather than assumed.

## The port matches the reference on 0 of 16 rows

That number is the point of the exercise, and it is not a regression. Every divergence attributes
to an argument the binary does not implement:

| divergence | rows | cause |
|---|---:|---|
| category column differs (fields 0-2), extra rows | 8 | `--METRIC_ACCUMULATION_LEVEL` other than `ALL_READS`: the reference emits one row per sample, library or read group; the port emits `ALL_READS` only |
| fields 11-13 and 26 differ | 6 | `--IS_BISULFITE_SEQUENCED=true`, which the port ignores |
| trailing empty columns | some rows | the reference writes the empty `SAMPLE`/`LIBRARY`/`READ_GROUP` columns; the port stops at the last populated one |

None of them is arithmetic. `crates/picard-analysis/src/bin/collect-alignment-summary-metrics.rs`
says so in its own header: it takes three arguments, exists for the benchmark, and is not the
Barclay command line the programme commits to. The array simply makes the consequence countable.

So the README's "byte-identical" for this tool means **byte-identical on the default path**, which
is what its conformance suite covers. Pairwise over the argument surface, it is 0%. Both numbers
are true and they are not the same claim; the dashboard now carries the weaker one as
"not measured", and this is the first tool where it *is* measured.

The trailing-empty-columns difference is the one finding here that is not explained by a missing
argument, and it is worth its own slice: the conformance suite does not catch it because that
suite compares the metrics body of a default run, where the columns happen to line up.

## Rejected rows are outputs, not failures

The reference rejects rows, and the message is behaviour the port has to reproduce:

* `FLOW_MODE is obsolete. Flow support now provided by CollectQualityYieldMetricsFlow`
* `File ... should be coordinate sorted but the header says the sort order is queryname. If you
  believe the file to be coordinate sorted you may pass ASSUME_SORTED=true`
* `Requesting earlier reference sequence: 0 < 1`

The runner records the exit code and the message as that row's output, so a port that accepts what
the reference refuses is a divergence rather than a pass.

## Two things this exposed about the method

**A fixture defect looks exactly like a divergence.** The first run had nine of eleven rows
failing, one class being `Supplementary alignment flag should not be set for unaligned read`: the
generated corpus set the supplementary flag on a read that another rule had made unmapped, so the
rows were testing htsjdk's validator. Fixed in `MakeFixtures.java`, with the reason recorded there.

**A held-at default must not be passed back.** The first run of `CollectAlignmentSummaryMetrics`
rejected all sixteen rows: `--ADAPTER_SEQUENCE`'s default is a `List[String]`, and re-serializing
it produced one token with brackets and commas that the parser read as positional arguments. The
array now records, per excluded argument, whether its value came from a fixture (must be passed:
the tool has no default) or from the tool's own default (must be omitted: that is what a default
is).

## Constraints closed the rejection problem

Eight of `CollectQualityYieldMetrics`'s eleven rows were rejected before the tool ran, most of them
for two reasons that are properties of the tool rather than of the array. Declared constraints
(`forbid`, with a reason, scoped to the tools they apply to) are now honoured during generation:
`covering.py` never emits a forbidden tuple, and `--verify` counts forbidden tuples separately, so
a coverage figure is computed over combinations that can exist.

Re-measured against the oracle, same fixtures:

| tool | before | after |
|---|---|---|
| `CollectQualityYieldMetrics` | 11 rows, 3 accepted, 3 distinct outputs | 10 rows, **10 accepted**, 4 distinct outputs |
| `CollectAlignmentSummaryMetrics` | 16 rows, 12 accepted, 9 distinct outputs | 16 rows, **16 accepted**, 11 distinct outputs |

Every row now runs, and both tools produce more distinct outputs than before: rows spent on
rejections were rows not spent on the tool.

The shape of a constraint had to be measured too. The queryname clause was first declared as two
argument-pair constraints (`with ASSUME_SORTED=false`, and `with CREATE_INDEX=true`), and the rerun
showed rows still failing with `ASSUME_SORTED=true, CREATE_INDEX=false`. The truth was simpler and
different: a queryname input is invalid for these collectors under *any* combination, because
`ASSUME_SORTED=false` is refused up front and `true` lets it through until the reference walk fails
mid-stream. The clause is scoped to the two collectors, since the same file is an ordinary input
for `SortSam` or `SamToFastq`, and a global clause would have deleted a dimension those tools need.

## The rejections became their own suite

They are behaviour the port owes and the constrained arrays no longer reach them, so
`tools/rejection-conformance/RejectionDump.java` runs the four refused combinations through
`PicardCommandLine` and records the exception class and message. The fixture path is replaced with
a token: the message is the behaviour, the path is not.

The suite is declared `golden-pending`, a third status beside `oracle-backed` and `unchecked`. It
has no golden and none may be committed from here: the bit-identity contract accepts goldens
produced on real x86-64 CI, and this harness has so far only run in an emulated container, which is
exactly how the sixteen goldens of decision 0008 came to exist. The generated job runs the suite and
uploads its dump as an artefact; committing that artefact is what turns the suite `oracle-backed`.

A pending suite still has to assert something, and the first attempt shows why the something must
be named. With no fixtures mounted, the harness produced four rows of
`Cannot read non-existent file` and the declared row count of four accepted them. The suite now
declares the phrases it exists for (`FLOW_MODE is obsolete`,
`should be coordinate sorted but the header says the sort order is queryname`,
`Requesting earlier reference sequence`), and that caught a second defect immediately: the fixture
mount is read-only, so writing `--OUTPUT` into it failed before the tool reached the reference walk.
The harness now writes outside the mount.

The port's 0/16 stands until the tools take more than three arguments. The number to watch is not
whether it is zero but whether it moves when a slice lands.
