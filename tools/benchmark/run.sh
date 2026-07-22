#!/usr/bin/env bash
# Times picard-rs against Picard on the same input, on the same platform, and asserts that the
# two outputs are byte-identical in the same run.
#
# Both binaries run inside the pinned oracle container. Timing a native binary against an
# emulated JVM would produce a large number and no information.
#
# Usage: tools/benchmark/run.sh [read count]   (default 2000000)
set -euo pipefail
cd "$(dirname "$0")/../.."
READS="${1:-2000000}"
OUT="$PWD/target/benchmark"
mkdir -p "$OUT"

echo "== generating $READS reads"
docker run --rm --platform linux/amd64 \
  -v "$PWD/tools/benchmark":/harness:ro -v "$OUT":/out -w /work picard-rs-oracle:3.4.0 \
  "cp /harness/MakeBenchmarkBam.java . && javac -cp \$PICARD_JAR -d . MakeBenchmarkBam.java \
   && java -cp .:\$PICARD_JAR MakeBenchmarkBam $READS /out/bench"

echo "== cross-building the port for linux/amd64"
docker run --rm --platform linux/amd64 -v "$PWD":/src -v "$OUT":/out -w /src rust:1.97-slim \
  bash -c 'cargo build --release --bin collect-alignment-summary-metrics --bin add-oa-tag --target-dir /out/amd64 2>&1 | tail -1'

cat > "$OUT/bench.sh" <<'INNER'
cd /out
TIMEFORMAT='%3R'
echo "jvm_startup_floor_s"
for i in 1 2; do { time java -jar $PICARD_JAR CollectAlignmentSummaryMetrics --version >/dev/null 2>&1 ; } 2>&1; done
echo "picard_gkl_default_s"
for i in 1 2 3; do { time java -jar $PICARD_JAR CollectAlignmentSummaryMetrics I=bench.bam O=/out/picard.txt R=bench.fasta >/dev/null 2>&1 ; } 2>&1; done
echo "picard_jdk_inflater_s"
for i in 1 2 3; do { time java -jar $PICARD_JAR CollectAlignmentSummaryMetrics I=bench.bam O=/out/picard_jdk.txt R=bench.fasta USE_JDK_INFLATER=true USE_JDK_DEFLATER=true >/dev/null 2>&1 ; } 2>&1; done
echo "rust_s"
for i in 1 2 3; do { time /out/amd64/release/collect-alignment-summary-metrics I=bench.bam O=/out/rust.txt R=bench.fasta >/dev/null 2>&1 ; } 2>&1; done
echo "rust_phases"
PICARD_RS_TIMING=1 /out/amd64/release/collect-alignment-summary-metrics I=bench.bam O=/out/rust.txt R=bench.fasta
echo "== AddOATag: a per-record transform that parallelizes =="
echo "picard_addoatag_s"
for i in 1 2 3; do { time java -jar $PICARD_JAR AddOATag I=bench.bam O=/out/picard_oa.sam VALIDATION_STRINGENCY=SILENT >/dev/null 2>&1 ; } 2>&1; done
echo "rust_addoatag_1t_s"
for i in 1 2 3; do { time env RAYON_NUM_THREADS=1 /out/amd64/release/add-oa-tag I=bench.bam O=/out/rust_oa.sam >/dev/null 2>&1 ; } 2>&1; done
echo "rust_addoatag_nt_s"
for i in 1 2 3; do { time /out/amd64/release/add-oa-tag I=bench.bam O=/out/rust_oa.sam >/dev/null 2>&1 ; } 2>&1; done
echo "rust_addoatag_1t_phases"
PICARD_RS_TIMING=1 RAYON_NUM_THREADS=1 /out/amd64/release/add-oa-tag I=bench.bam O=/out/rust_oa.sam
echo "rust_addoatag_nt_phases"
PICARD_RS_TIMING=1 /out/amd64/release/add-oa-tag I=bench.bam O=/out/rust_oa.sam
INNER

echo "== timing"
docker run --rm --platform linux/amd64 -v "$OUT":/out -w /out picard-rs-oracle:3.4.0 'bash /out/bench.sh'

echo "== byte-identity"
python3 - "$OUT" <<'PY'
import sys
out = sys.argv[1]
def body(path):
    return ''.join(l for l in open(path)
                   if not l.startswith('# CollectAlignmentSummaryMetrics')
                   and not l.startswith('# Started on:'))
picard, rust = body(f'{out}/picard.txt'), body(f'{out}/rust.txt')
ok = picard == rust
print("picard vs picard-rs byte-identical:", ok)
print("GKL vs JDK inflater identical:", picard == body(f'{out}/picard_jdk.txt'))
oa_ok = open(f'{out}/picard_oa.sam','rb').read() == open(f'{out}/rust_oa.sam','rb').read()
print("AddOATag picard vs picard-rs byte-identical:", oa_ok)
if not ok or not oa_ok:
    if not ok:
        for i, (a, b) in enumerate(zip(picard.split('\n'), rust.split('\n'))):
            if a != b:
                print(f"line {i}\n  picard: {a[:200]}\n  ours  : {b[:200]}")
                break
    sys.exit(1)
PY
