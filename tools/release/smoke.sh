#!/usr/bin/env bash
# The release binaries read a BAM, write one through the BGZF deflater, and read it back.
# `RUN` prefixes every call (e.g. "qemu-x86_64 -cpu qemu64"), so the same check proves the binaries
# need nothing a baseline x86-64 lacks: an unsupported instruction is SIGILL, not a wrong count.
set -euo pipefail
bin=$1
here=$(cd "$(dirname "$0")" && pwd)
out=$(mktemp -d)
run=${RUN:-}

$run "$bin/SortSam" I="$here/smoke.bam" O="$out/out.bam" SORT_ORDER=queryname 2>/dev/null
test "$($run "$bin/ViewSam" I="$out/out.bam" ALIGNMENT_STATUS=All PF_STATUS=All 2>/dev/null | grep -vc '^@')" = 2
echo "smoke: ok${run:+ under $run}"
