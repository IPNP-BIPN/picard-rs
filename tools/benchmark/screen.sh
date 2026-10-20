#!/usr/bin/env bash
# Every tool that reads a BAM, once, on the benchmark input, with a timeout: a tool that does not
# finish is a quadratic path to look at, and the RSS column is where the memory goes.
# Usage: tools/benchmark/screen.sh <bin dir> [timeout seconds]   (inputs from tools/benchmark/local.sh's dir)
set -uo pipefail
BIN="$1"; LIMIT="${2:-120}"
B=$(cd "$(dirname "$0")/../.." && pwd)/target/benchmark
[ -f "$B/bench.bam" ] || B=$HOME/GATK/picard-rs/target/benchmark
O=$B/screen-$(basename "$(dirname "$BIN")"); rm -rf "$O"; mkdir -p "$O"
I=$B/bench.bam; R=$B/bench.fasta; Q=$B/qname.bam; T=$B/targets.interval_list
# The binaries are named as Picard names the tool; the list below keeps the kebab names.
exe() { local want; want=$(echo "$1" | tr -d - | tr 'A-Z' 'a-z'); for f in "$BIN"/*; do [ "$(basename "$f" | tr 'A-Z' 'a-z')" = "$want" ] && { echo "$f"; return; }; done; echo "$BIN/$1"; }
run() {
  local name=$1; shift
  /usr/bin/time -l timeout "$LIMIT" "$(exe "$name")" "$@" >"$O/$name.out" 2>"$O/$name.err"; local code=$?
  # CPU seconds (user + sys) as well as wall clock: other work on the machine inflates the latter.
  local s=$(awk '/ real /{print $1}' "$O/$name.err"); local c=$(awk '/ real /{print $3+$5}' "$O/$name.err"); local r=$(awk '/maximum resident/{printf "%d", $1/1048576}' "$O/$name.err")
  local why=$(grep -v -E '^ +[0-9]|real|^$' "$O/$name.err" | tail -1 | cut -c1-70)
  printf "%-46s %3s %7s s %7s cpu %6s MB  %s\n" "$name" "$code" "$s" "$c" "$r" "$why"
}
run add-comments-to-bam I=$I O=$O/a.bam C=hello
run add-oa-tag I=$I O=$O/b.bam
run add-or-replace-read-groups I=$I O=$O/c.bam RGLB=l RGPL=ILLUMINA RGPU=u RGSM=s
run build-bam-index I=$I O=$O/bench.bai
run bam-index-stats I=$B/coord.bam  # wants coord.bai beside it
run calculate-read-group-checksum I=$I O=$O/d.md5
run check-terminator-block I=$I
run clean-sam I=$I O=$O/e.bam
run collect-alignment-summary-metrics I=$I O=$O/f.txt R=$R
run collect-base-distribution-by-cycle I=$I O=$O/g.txt CHART=$O/g.pdf
run collect-duplicate-metrics I=$I M=$O/h.txt
run collect-gc-bias-metrics I=$I O=$O/i.txt CHART=$O/i.pdf S=$O/i2.txt R=$R
run collect-hs-metrics I=$I O=$O/j.txt R=$R BAIT_INTERVALS=$T TARGET_INTERVALS=$T
run collect-insert-size-metrics I=$I O=$O/k.txt H=$O/k.pdf
run collect-jumping-library-metrics I=$I O=$O/l.txt
run collect-multiple-metrics I=$I O=$O/m R=$R
run collect-oxo-g-metrics I=$I O=$O/n.txt R=$R
run collect-quality-yield-metrics I=$I O=$O/o.txt
run collect-quality-yield-metrics-flow I=$I O=$O/o2.txt
run collect-quality-yield-metrics-snvq I=$I O=$O/o3.txt
run collect-raw-wgs-metrics I=$I O=$O/p.txt R=$R
run collect-rna-seq-metrics I=$I O=$O/q.txt REF_FLAT=$B/refflat.txt STRAND_SPECIFICITY=NONE
run collect-rrbs-metrics I=$I M=$O/r R=$R
run collect-sam-error-metrics I=$I O=$O/s R=$R ERROR_METRICS=ERROR ERROR_METRICS=ERROR:READ_ORDINALITY ERROR_METRICS=ERROR:BASE_QUALITY ERROR_METRICS=ERROR:CYCLE
run collect-sequencing-artifact-metrics I=$I O=$O/t R=$R
run collect-targeted-pcr-metrics I=$I O=$O/u.txt R=$R AMPLICON_INTERVALS=$T TARGET_INTERVALS=$T
run collect-umi-prevalence-metrics I=$I O=$O/v.txt
run collect-wgs-metrics I=$I O=$O/w.txt R=$R
run collect-wgs-metrics-with-non-zero-coverage I=$I O=$O/w2.txt CHART=$O/w2.pdf R=$R
run compare-sams $I $B/coord.bam O=$O/x.txt
run downsample-sam I=$I O=$O/y.bam P=0.5
run estimate-library-complexity I=$I O=$O/z.txt
run filter-sam-reads I=$I O=$O/fa.bam FILTER=includeAligned
run fix-mate-information I=$Q O=$O/fb.bam
run gather-bam-files I=$I O=$O/fc.bam
run mark-duplicates I=$I O=$O/fd.bam M=$O/fd.txt
run mark-duplicates-with-mate-cigar I=$I O=$O/fe.bam M=$O/fe.txt
run simple-mark-duplicates-with-mate-cigar I=$I O=$O/ff.bam M=$O/ff.txt
run umi-aware-mark-duplicates-with-mate-cigar I=$I O=$O/fg.bam M=$O/fg.txt UMI_METRICS=$O/fg2.txt
run mark-illumina-adapters I=$Q O=$O/fh.bam M=$O/fh.txt
run mean-quality-by-cycle I=$I O=$O/fi.txt CHART=$O/fi.pdf
run merge-sam-files I=$I I=$B/coord.bam O=$O/fj.bam
run position-based-downsample-sam I=$I O=$O/fk.bam FRACTION=0.5
run quality-score-distribution I=$I O=$O/fl.txt CHART=$O/fl.pdf
run reorder-sam I=$I O=$O/fm.bam SEQUENCE_DICTIONARY=$B/bench.dict
run replace-sam-header I=$I HEADER=$B/coord.bam O=$O/fn.bam
run revert-original-base-qualities-and-add-mate-cigar I=$Q O=$O/fo.bam
run revert-sam I=$I O=$O/fp.bam
# The benchmark BAM names every record apart, so a mate, a secondary or a supplementary record has
# no unmapped read to merge with, and MergeBamAlignment refuses it as Picard does. Its unpaired
# primary reads that keep every base, in queryname order, are the alignment; RevertSam's copy of
# them is the unmapped input. Prepared here, not measured.
"$(exe view-sam)" I=$I | awk '/^@/ || ($2 % 2 == 0 && $2 % 512 < 256 && $2 < 2048 && $6 !~ /H/)' >$O/fp1.sam
"$(exe sort-sam)" I=$O/fp1.sam O=$O/fp1.bam SORT_ORDER=queryname 2>/dev/null; rm -f $O/fp1.sam
"$(exe revert-sam)" I=$O/fp1.bam O=$O/fp1u.bam 2>/dev/null
run merge-bam-alignment UNMAPPED_BAM=$O/fp1u.bam ALIGNED_BAM=$O/fp1.bam R=$R O=$O/fp2.bam
run sam-format-converter I=$I O=$O/fq.sam
run sam-to-fastq I=$Q F=$O/fr1.fq F2=$O/fr2.fq VALIDATION_STRINGENCY=SILENT
run sam-to-fastq-with-tags I=$Q F=$O/fs1.fq F2=$O/fs2.fq QUALITY_TAG_GROUP=OQ SEQUENCE_TAG_GROUP=OA
run set-nm-md-and-uq-tags I=$I O=$O/ft.bam R=$R
run set-nm-and-uq-tags I=$I O=$O/ft2.bam R=$R
run sort-sam I=$I O=$O/fu.bam SORT_ORDER=queryname
run split-sam-by-library I=$I O=$O/fv
run split-sam-by-number-of-reads I=$I OUTPUT=$O/fw SPLIT_TO_N_FILES=4
run validate-sam-file I=$I MODE=SUMMARY
run view-sam I=$I
run check-duplicate-marking I=$O/fd.bam
