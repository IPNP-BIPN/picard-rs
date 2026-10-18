#!/usr/bin/env bash
# Wall clock and peak RSS of the release binaries on the benchmark BAM, natively, best of N.
# Usage: tools/benchmark/local.sh <bin dir> [runs]   (writes outputs under target/benchmark/out-<tag>)
set -uo pipefail
cd "$(dirname "$0")/../.."
BIN="$1"; RUNS="${2:-3}"; B=$PWD/target/benchmark; O=$B/out-$(basename "$(dirname "$BIN")")-$(basename "$BIN"); mkdir -p "$O"
# Binaries are named as Picard names the tool, or in kebab case before that rename.
bin() { if [ -x "$BIN/$1" ]; then echo "$BIN/$1"; else echo "$BIN/$(echo "$1" | sed -E 's/([a-z0-9])([A-Z])/\1-\2/g; s/([A-Z])([A-Z][a-z])/\1-\2/g' | tr 'A-Z' 'a-z')"; fi; }
run() { # name, args...
  local name=$1; shift 2; set -- "$(bin "$name")" "$@"; local best=999 rss=0
  for _ in $(seq "$RUNS"); do
    /usr/bin/time -l "$@" >/dev/null 2>"$O/$name.err"
    local s=$(awk '/real/{print $1}' "$O/$name.err"); local r=$(awk '/maximum resident/{print $1}' "$O/$name.err")
    if awk "BEGIN{exit !($s < $best)}"; then best=$s; rss=$r; fi
  done
  printf "%-34s %7.2f s %8.0f MB\n" "$name" "$best" "$((rss/1048576))"
}
run CollectAlignmentSummaryMetrics "$BIN/CollectAlignmentSummaryMetrics" I=$B/bench.bam O=$O/casm.txt R=$B/bench.fasta
run AddOATag "$BIN/AddOATag" I=$B/bench.bam O=$O/oa.bam
run CollectInsertSizeMetrics "$BIN/CollectInsertSizeMetrics" I=$B/bench.bam O=$O/ism.txt H=$O/ism.pdf
run CollectQualityYieldMetrics "$BIN/CollectQualityYieldMetrics" I=$B/bench.bam O=$O/qy.txt
run MeanQualityByCycle "$BIN/MeanQualityByCycle" I=$B/bench.bam O=$O/mqc.txt CHART=$O/mqc.pdf
run SortSam "$BIN/SortSam" I=$B/bench.bam O=$O/sorted.bam SORT_ORDER=queryname
# MarkDuplicates is quadratic somewhere: more than 15 minutes at 2M reads. Its own investigation.
# run MarkDuplicates "$BIN/MarkDuplicates" I=$B/bench.bam O=$O/md.bam M=$O/md.txt
run SamToFastq "$BIN/SamToFastq" I=$B/bench.bam F=$O/r1.fq F2=$O/r2.fq
run ValidateSamFile "$BIN/ValidateSamFile" I=$B/bench.bam MODE=SUMMARY
run CollectGcBiasMetrics "$BIN/CollectGcBiasMetrics" I=$B/bench.bam O=$O/gc.txt CHART=$O/gc.pdf S=$O/gcs.txt R=$B/bench.fasta
