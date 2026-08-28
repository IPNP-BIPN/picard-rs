# 0011. The sixteen unchecked goldens came back green

**Status:** accepted
**Date:** 2026-08-28

## What decision 0008 left open

Decision 0008 found sixteen goldens that no CI step had ever regenerated against the oracle, gave
them suites of their own, and marked them `unchecked` so a green run could not be read as the
stronger claim. It closed on a question:

> The first CI run of the sixteen may well be red, and a red result there is the point of the
> exercise: it distinguishes "the golden was fine, only its paths were machine-specific" from "the
> golden was never right".

That run has happened, and many since. Every one of the fourteen suites the sixteen goldens sit in
regenerates its dump in the pinned container on real x86-64 CI and compares it, and every one is
green:

| suite | compared |
|---|---:|
| `intervallisttools` (three corpora) | 51 |
| `mergebamalignment` | 23 |
| `bedtointervallist` | 21 |
| `liftoverintervallist` | 18 |
| `scatterintervalsbyns` | 18 |
| `gcbias` | 27 |
| `normalizefasta` | 15 |
| `mergebamalignment-full` | 15 |
| `qualityscoredistribution` | 11 |
| `extractsequences` | 10 |
| `intervallisttobed` | 10 |
| `samtofastqwithtags` | 10 |
| `nonnfastasize` | 6 |
| `accumulatequalityyield` | 4 |

The answer is the first of the two: the goldens were fine, and only their paths were
machine-specific.

## The criterion is decision 0008's own

Decision 0008 already settled what the status means, in the section that admitted four
`oracle-backed` goldens carrying macOS paths:

> The distinction that matters is not where a file was authored but whether anything re-derives it.

Fourteen suites now re-derive their goldens on every run. By that criterion they are
`oracle-backed`, and this change says so.

## What the promotion rests on

Nine of the suites canonicalize something away, and the promotion is only as strong as those
rules, so they are named here rather than left to be read out of the manifest:

- `strip_ur` on `bedtointervallist`, `intervallisttobed`, `intervallisttools`,
  `liftoverintervallist`, `extractsequences`, `scatterintervalsbyns`, `samtofastqwithtags` and
  `mergebamalignment-full`. The `@SQ UR:` field carries the fixture's own path, which
  `Files.createTempDirectory` makes new on every run: it is not something a byte comparison can be
  asked to reproduce;
- `strip_pg` on `intervallisttools` as well, for the same reason applied to the command line;
- `strip_line_prefixes` on `qualityscoredistribution` and `gcbias`, for the `# <Tool>` header and
  the `# Started on:` line, which are the command line and the wall clock.

Five suites canonicalize nothing at all: `accumulatequalityyield`, `nonnfastasize`,
`normalizefasta` and both `mergebamalignment` corpora bar the `UR:` of the full one. Those are
compared byte for byte.

## What is deliberately left as it was

The `unchecked` status stays in the manifest and in `generate_ci.py`'s job grouping. It is what a
new golden gets when its suite has not yet run against the oracle, and the mechanism decision 0008
built is the reason this promotion could be made on evidence rather than on assertion. Deleting
the weaker label would take that away.
