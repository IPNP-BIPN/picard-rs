#!/usr/bin/env bash
# Every tool that does not read a BAM (VCF, arrays, basecalls, fingerprints, interval lists, FASTA,
# FASTQ, GFF, metrics files), once, on inputs large enough to measure the work, with a timeout: a
# tool that does not finish is a quadratic path to look at, and the RSS column is where the memory
# goes. The few that read a BAM besides (the fingerprint tools, BamToBfq) read bench.bam.
# Inputs: tools/benchmark/make_nonbam_inputs.py, after MakeBigArrays.java in the oracle container.
# Usage: tools/benchmark/screen_nonbam.sh <bin dir | oracle> [timeout seconds] [tool ...]
# `oracle` runs Picard 3.5.0 in the oracle container (linux/amd64, emulated on this laptop, so its
# times say nothing) on the same command lines, into screen-nonbam-oracle, for a diff of outputs.
set -uo pipefail
BIN="$1"; LIMIT="${2:-300}"; shift; shift || true; ONLY=" $* "
B=$HOME/GATK/picard-rs/target/benchmark
N=$B/nonbam; G=$B/big; A=$N/arrays
if [ "$BIN" = oracle ]; then O=$B/screen-nonbam-oracle; else O=$B/screen-nonbam-$(basename "$(dirname "$BIN")"); fi; [ "$ONLY" = "  " ] && rm -rf "$O"; mkdir -p "$O"
mkdir -p $N/ill_out $O/bfq $O/baits $O/bc
R=$B/bench.fasta; D=$B/bench.dict; I=$B/bench.bam; RUN=$N/ill_big/Data/Intensities/BaseCalls
# A row is `run <tool>[@<label>] <arguments>`: the label tells two rows of one tool apart.
# The binaries are named as Picard names the tool; the list below keeps the kebab names.
exe() { local want; want=$(echo "$1" | tr -d - | tr 'A-Z' 'a-z'); for f in "$BIN"/*; do [ "$(basename "$f" | tr 'A-Z' 'a-z')" = "$want" ] && { echo "$f"; return; }; done; echo "$BIN/$1"; }
run() {
  local name=$1; local bin=${1%@*}; shift
  [ "$ONLY" != "  " ] && [[ "$ONLY" != *" $name "* ]] && return
  if [ "$BIN" = oracle ]; then
    # Picard's class name: the kebab-case parts capitalised, but for the few with a capital inside.
    local tool=$(echo "$bin" | awk -F- '{for(i=1;i<=NF;i++) printf "%s", toupper(substr($i,1,1)) substr($i,2)}' \
      | sed -e 's/VerifyIdIntensity/VerifyIDIntensity/' -e 's/ToOxog$/ToOxoG/')
    /usr/bin/time -l timeout "$LIMIT" docker run -i --rm --platform linux/amd64 -v "$B:$B" -w "$O" picard-rs-oracle:3.5.0 \
      "java -Xmx8g -Dsamjdk.use_libdeflate=false -jar \$PICARD_JAR $tool $*" >"$O/$name.out" 2>"$O/$name.err"
  else
    /usr/bin/time -l timeout "$LIMIT" "$(exe "$bin")" "$@" >"$O/$name.out" 2>"$O/$name.err"
  fi
  local code=$?
  # CPU seconds (user + sys) as well as wall clock: other work on the machine inflates the latter.
  local s=$(awk '/ real /{print $1}' "$O/$name.err"); local c=$(awk '/ real /{print $3+$5}' "$O/$name.err"); local r=$(awk '/maximum resident/{printf "%d", $1/1048576}' "$O/$name.err")
  local why=$(grep -v -E '^ +[0-9]|real|^$' "$O/$name.err" | tail -1 | cut -c1-70)
  printf "%-48s %3s %7s s %7s cpu %6s MB  %s\n" "$name" "$code" "$s" "$c" "$r" "$why"
}
# VCF
run collect-variant-calling-metrics I=$G/big.vcf DBSNP=$N/dbsnp.vcf O=$O/cvcm SD=$D
run accumulate-variant-calling-metrics I=$O/cvcm I=$O/cvcm O=$O/avcm
run filter-vcf I=$G/big.vcf O=$O/filter.vcf MIN_DP=10 MIN_GQ=20 MIN_AB=0.2
run find-mendelian-violations I=$G/big.vcf TRIOS=$G/trio.ped O=$O/mendel.txt
run fix-vcf-header I=$N/big_missing.vcf O=$O/fixed.vcf
run gather-vcfs I=$N/part1.vcf I=$N/part2.vcf I=$N/part3.vcf I=$N/part4.vcf O=$O/gathered.vcf
run genotype-concordance TRUTH_VCF=$G/big.vcf CALL_VCF=$N/big_call.vcf TRUTH_SAMPLE=child CALL_SAMPLE=child O=$O/gc OUTPUT_VCF=true
run liftover-vcf I=$G/big.vcf O=$O/lifted.vcf REJECT=$O/rejected.vcf CHAIN=$G/identity.chain R=$R
run make-sites-only-vcf I=$G/big.vcf O=$O/sites.vcf
run make-vcf-sample-name-map I=$N/big_one.vcf O=$O/names.txt
run merge-vcfs I=$N/part1.vcf I=$N/part2.vcf I=$N/part3.vcf I=$N/part4.vcf O=$O/merged.vcf
run rename-sample-in-vcf I=$N/big_one.vcf O=$O/renamed.vcf NEW_SAMPLE_NAME=renamed
run sort-vcf I=$G/shuffled.vcf O=$O/sorted.vcf
run split-vcfs I=$G/big.vcf SNP_OUTPUT=$O/snp.vcf INDEL_OUTPUT=$O/indel.vcf
run update-vcf-sequence-dictionary I=$G/big.vcf O=$O/updated.vcf SD=$D
run vcf-format-converter I=$G/big.vcf O=$O/converted.vcf REQUIRE_INDEX=false
run vcf-to-interval-list I=$G/big.vcf O=$O/vcf.interval_list
run collect-independent-replicate-metrics I=$B/coord.bam VCF=$G/big.vcf SAMPLE=child O=$O/cirm.txt MATRIX_OUTPUT=$O/cirm_matrix.txt
# Genotyping arrays
run gtc-to-vcf I=$A/arr_s1.gtc O=$O/arr.vcf EXTENDED_ILLUMINA_MANIFEST=$A/arr_extended.csv CLUSTER_FILE=$A/arr.egt ILLUMINA_BEAD_POOL_MANIFEST_FILE=$A/arr.bpm R=$A/arr_ref.fasta EXPECTED_GENDER=Female SAMPLE_ALIAS=sample1 PIPELINE_VERSION=1.0 ANALYSIS_VERSION_NUMBER=1
run bpm-to-normalization-manifest-csv I=$A/arr.bpm CLUSTER_FILE=$A/arr.egt O=$O/norm.csv
run compare-gtc-files I=$A/arr_s1.gtc I=$A/arr_s2.gtc ILLUMINA_BEAD_POOL_MANIFEST_FILE=$A/arr.bpm
run create-extended-illumina-manifest I=$A/arr_manifest.csv O=$O/ext.csv REPORT_FILE=$O/ext_report.txt BAD_ASSAYS_FILE=$O/ext_bad.txt CLUSTER_FILE=$A/arr.egt DBSNP_FILE=$A/arr_dbsnp.vcf TARGET_BUILD=37 R=$A/arr_ref.fasta
run collect-arrays-variant-calling-metrics I=$A/arr_s1.vcf DBSNP=$A/arr_dbsnp.vcf O=$O/avcm_arr
run combine-genotyping-array-vcfs I=$A/arr_s1.vcf I=$A/arr_s2.vcf O=$O/combined.vcf
run vcf-to-adpc VCF=$A/arr_s1.vcf VCF=$A/arr_s2.vcf O=$O/adpc.bin SAMPLES_FILE=$O/adpc_samples.txt NUM_MARKERS_FILE=$O/adpc_markers.txt
run merge-ped-into-vcf ORIGINAL_VCF=$A/arr_s1.vcf PED_FILE=$N/zcall.ped MAP_FILE=$N/zcall.map ZCALL_THRESHOLDS_FILE=$N/zcall_thresholds.txt ZCALL_VERSION=1.0.0.0 O=$O/zcall.vcf
run create-baf-regress-metrics-file I=$N/bafregress.txt O=$O/baf
run create-verify-id-intensity-contamination-metrics-file I=$N/verifyidintensity.txt O=$O/verifyid
# Fingerprints
run convert-haplotype-database-to-vcf I=$N/hap.txt O=$O/hap.vcf R=$R
run lift-over-haplotype-map I=$N/hap.txt O=$O/hap_lifted.txt CHAIN=$G/identity.chain SD=$D
run check-fingerprint I=$I GENOTYPES=$G/big.vcf HAPLOTYPE_MAP=$N/hap.txt O=$O/cfp EXPECTED_SAMPLE_ALIAS=child R=$R
run extract-fingerprint I=$I O=$O/efp.vcf HAPLOTYPE_MAP=$N/hap.txt R=$R
run identify-contaminant I=$I O=$O/contam.vcf HAPLOTYPE_MAP=$N/hap.txt R=$R
run calculate-fingerprint-metrics I=$I O=$O/cfm.txt HAPLOTYPE_MAP=$N/hap.txt R=$R
run crosscheck-fingerprints I=$N/fp_many.vcf CROSSCHECK_BY=SAMPLE HAPLOTYPE_MAP=$N/hap.txt O=$O/xc.txt LOD_THRESHOLD=-5 MATRIX_OUTPUT=$O/xc_matrix.txt EXIT_CODE_WHEN_MISMATCH=0
run crosscheck-read-group-fingerprints I=$I HAPLOTYPE_MAP=$N/hap.txt O=$O/xcrg.txt R=$R
run cluster-crosscheck-metrics I=$O/xc.txt O=$O/cluster.txt LOD_THRESHOLD=5
# Illumina basecalls
run check-illumina-directory BASECALLS_DIR=$RUN READ_STRUCTURE=25T8B25T LANES=1
run extract-illumina-barcodes BASECALLS_DIR=$RUN READ_STRUCTURE=25T8B25T BARCODE_FILE=$N/ill_barcodes.txt LANE=1 OUTPUT_DIR=$O/bc METRICS_FILE=$O/bc_metrics.txt
run collect-illumina-basecalling-metrics BASECALLS_DIR=$RUN BARCODES_DIR=$O/bc LANE=1 READ_STRUCTURE=25T8B25T INPUT=$N/ill_barcodes.txt O=$O/ibm.txt
run collect-illumina-lane-metrics RUN_DIRECTORY=$N/ill_lanes_many OUTPUT_DIRECTORY=$O OUTPUT_PREFIX=lanes READ_STRUCTURE=25T8B25T
run collect-hi-seq-x-pf-fail-metrics BASECALLS_DIR=$RUN O=$O/pffail LANE=1
run illumina-basecalls-to-fastq BASECALLS_DIR=$RUN BARCODES_DIR=$O/bc LANE=1 READ_STRUCTURE=25T8B25T MULTIPLEX_PARAMS=$N/ill_multiplex.txt RUN_BARCODE=R1 MACHINE_NAME=M1 FLOWCELL_BARCODE=FC1
run illumina-basecalls-to-sam BASECALLS_DIR=$RUN BARCODES_DIR=$O/bc LANE=1 READ_STRUCTURE=25T8B25T LIBRARY_PARAMS=$N/ill_library.txt RUN_BARCODE=R1 SEQUENCING_CENTER=BI TMP_DIR=$O
# A plate of 384 barcodes on one tile: what grows with the barcode count rather than the clusters.
R384=$N/ill_384/Data/Intensities/BaseCalls; mkdir -p $O/bc384 $N/ill_out_384
run extract-illumina-barcodes@384 BASECALLS_DIR=$R384 READ_STRUCTURE=25T8B25T BARCODE_FILE=$N/ill_barcodes_384.txt LANE=1 OUTPUT_DIR=$O/bc384 METRICS_FILE=$O/bc384_metrics.txt
run illumina-basecalls-to-fastq@384 BASECALLS_DIR=$R384 BARCODES_DIR=$O/bc384 LANE=1 READ_STRUCTURE=25T8B25T MULTIPLEX_PARAMS=$N/ill_multiplex_384.txt RUN_BARCODE=R1 MACHINE_NAME=M1 FLOWCELL_BARCODE=FC1
# Interval lists, FASTA, GFF, FASTQ
run bed-to-interval-list I=$G/many.bed O=$O/bed.interval_list SD=$D
run interval-list-to-bed I=$G/many.interval_list O=$O/il.bed
run interval-list-tools I=$G/many.interval_list O=$O/uniq.interval_list UNIQUE=true SORT=true
run lift-over-interval-list I=$G/many.interval_list O=$O/lifted.interval_list CHAIN=$G/identity.chain SD=$D REJECT=$O/rejected.interval_list
run bait-designer TARGETS=$N/bait_targets.interval_list R=$R DESIGN_NAME=design OUTPUT_DIRECTORY=$O/baits
run scatter-intervals-by-ns R=$N/genome.fasta O=$O/ns.interval_list
run non-n-fasta-size I=$N/genome.fasta O=$O/nonn.txt
run normalize-fasta I=$N/genome.fasta O=$O/normalized.fasta LINE_LENGTH=100
run create-sequence-dictionary R=$N/genome.fasta O=$O/genome.dict
run extract-sequences INTERVAL_LIST=$G/many.interval_list R=$R O=$O/extracted.fasta
run sort-gff I=$G/big.gff3 O=$O/sorted.gff3
run fastq-to-sam FASTQ=$N/reads.fq OUTPUT=$O/fastq.bam SAMPLE_NAME=s1 READ_GROUP_NAME=rg1
run bam-to-bfq I=$I ANALYSIS_DIR=$O/bfq FLOWCELL_BARCODE=FC1 LANE=1 READ_NAME_PREFIX=x PAIRED_RUN=false
run fifo-buffer < $N/genome.fasta
# Metrics files
run accumulate-quality-yield-metrics $(for f in $N/qy/q1[0-9][0-9].metrics; do printf 'I=%s ' $f; done) O=$O/aqy.txt
run compare-metrics I=$B/screen-screenbin/k.txt I=$B/screen-screenbin/k.txt O=$O/cmp.txt
run convert-sequencing-artifact-to-oxog INPUT_BASE=$B/screen-screenbin/t OUTPUT_BASE=$O/oxog
