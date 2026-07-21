#!/usr/bin/env bash
# Probe: is BAD_CYCLES computed from the read index or from the alignment-block offset?
#
# AlignmentSummaryMetricsCollector.collectQualityData walks the read's alignment blocks and, on a
# mismatch, does:
#
#     badCycleHistogram.increment(CoordMath.getCycle(negativeStrand, readBases.length, i));
#
# where `i` is the offset *within the current alignment block*, while getCycle's third parameter
# is declared `final int readBaseIndex`. Two bases in different blocks that share a block offset
# therefore land in the same cycle bin, and a read is charged one bad cycle instead of two.
#
# The two inputs differ only in where the reference mismatches sit:
#
#   collide   mismatches at read indices 3 and 13, block offsets 3 and 3
#   distinct  mismatches at read indices 3 and 15, block offsets 3 and 5
#
# Both reads are 20 bases with CIGAR 10M5D10M, so both have exactly two alignment blocks and
# exactly two mismatches. If the cycle were the read index, both cases would give BAD_CYCLES=2.
#
# Expected, and what the pinned oracle prints:
#
#   collide:  BAD_CYCLES=1  PF_HQ_MEDIAN_MISMATCHES=2
#   distinct: BAD_CYCLES=2  PF_HQ_MEDIAN_MISMATCHES=2
#
# PF_HQ_MEDIAN_MISMATCHES=2 in both is the control: both mismatches are seen in both runs, so the
# difference in BAD_CYCLES is the binning and not a missed mismatch.
set -euo pipefail
cd "$(dirname "$0")"

for case in collide distinct; do
  cp "$case.fasta" ref.fasta
  cp "$case.fasta.fai" ref.fasta.fai
  cp "$case.dict" ref.dict
  docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work picard-rs-oracle:3.4.0 \
    "java -jar \$PICARD_JAR CollectAlignmentSummaryMetrics I=probe.sam O=out_$case.txt R=ref.fasta >/dev/null 2>&1"
  python3 - "$case" <<'PY'
import sys
case = sys.argv[1]
lines = open(f'out_{case}.txt').read().split('\n')
h = next(i for i, l in enumerate(lines) if l.startswith('CATEGORY'))
d = dict(zip(lines[h].split('\t'), lines[h + 1].split('\t')))
print(f"{case}: BAD_CYCLES={d['BAD_CYCLES']}  PF_HQ_MEDIAN_MISMATCHES={d['PF_HQ_MEDIAN_MISMATCHES']}")
PY
done
rm -f ref.fasta ref.fasta.fai ref.dict
