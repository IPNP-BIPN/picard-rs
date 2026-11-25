#!/usr/bin/env python3
"""Derive every input tools/speed/rows.txt names from the chr20 scale data, reproducibly.

    python3 tools/speed/make_scale_inputs.py <scale-data dir> <derived dir>

The scale data (gatk-rs tools/speed/scale.json pins it: NA12878 chr20 30x BAM, chr20 FASTA,
dbSNP 138 and gnomAD chr20 VCFs) is mounted read-only; everything else a row needs is built here
from it with fixed seeds, laid out as tools/benchmark/make_nonbam_inputs.py expects (`bench.*` and
`big/`), and that script is then run on it for the arrays, Illumina, fingerprint and metrics inputs.

The queryname BAM and the duplicate-marked BAM are made by Picard itself (in the speed image),
so no input is shaped by the port. Real VCFs are sites-only, so the multi-sample callset is
dbSNP's biallelic sites with a seeded trio of genotypes.
"""

import gzip
import os
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCALE = Path(sys.argv[1]).resolve()
D = Path(sys.argv[2]).resolve()
BIG = D / "big"
rng = random.Random(20)
CONTIG, LENGTH = "chr20", 64444167


def link(name, target):
    """A symlink resolved inside the container (scale data at /data), or a copy of a small file."""
    path = D / name
    if path.exists() or path.is_symlink():
        return
    if target.stat().st_size < 1 << 20:
        path.write_bytes(target.read_bytes())
    else:
        path.symlink_to(Path("/data") / target.name)


def write(path, lines):
    if not path.exists():
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "w") as f:
            f.writelines(lines)
        tmp.rename(path)


def picard(command):
    """One Picard run in the speed image, scale data at /data, derived inputs at /derived."""
    subprocess.run(
        ["docker", "run", "--rm", "-m", "6g", "-v", f"{SCALE}:/data:ro", "-v", f"{D}:/derived",
         "picard-rs-speed:3.5.0", f"java -jar $PICARD_JAR {command}"],
        check=True,
    )


BIG.mkdir(parents=True, exist_ok=True)
for name in ("bench.bam", "bench.bam.bai", "bench.fasta", "bench.fasta.fai", "bench.dict"):
    source = {"bench.bam": "NA12878.chr20.bam", "bench.bam.bai": "NA12878.chr20.bam.bai",
              "bench.fasta": "chr20.fasta", "bench.fasta.fai": "chr20.fasta.fai",
              "bench.dict": "chr20.dict"}[name]
    link(name, SCALE / source)
link("bench.bai", SCALE / "NA12878.chr20.bam.bai")
link("coord.bam", SCALE / "NA12878.chr20.bam")
link("coord.bai", SCALE / "NA12878.chr20.bam.bai")
link("big/genome.fasta", SCALE / "chr20.fasta")
link("big/genome.fasta.fai", SCALE / "chr20.fasta.fai")
link("big/genome.dict", SCALE / "chr20.dict")
dict_text = (SCALE / "chr20.dict").read_text()
header_il = "@HD\tVN:1.6\n" + "".join(l + "\n" for l in dict_text.splitlines() if l.startswith("@SQ"))

# The trio callset: every biallelic dbSNP site, GT:AD:DP:GQ per sample, AF in INFO.
if not (BIG / "big.vcf").exists():
    head = ["##fileformat=VCFv4.2\n",
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n',
            f"##contig=<ID={CONTIG},length={LENGTH}>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmother\tfather\tchild\n"]
    body = []
    gts = ("0/0", "0/1", "1/1")
    with gzip.open(SCALE / "dbsnp138.chr20.vcf.gz", "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c, pos, ident, ref, alt = line.split("\t", 5)[:5]
            if "," in alt:
                continue
            row = [c, pos, ident, ref, alt, str(rng.randrange(10, 900)), "PASS",
                   f"AF={rng.random():.3f}", "GT:AD:DP:GQ"]
            for _ in range(3):
                g = rng.randrange(3)
                dp = rng.randrange(0, 60)
                alt_reads = (0, dp // 2, dp)[g]
                row.append(f"{gts[g]}:{dp - alt_reads},{alt_reads}:{dp}:{rng.randrange(0, 99)}")
            body.append("\t".join(row) + "\n")
    write(BIG / "big.vcf", head + body)
    shuffled = body[:]
    rng.shuffle(shuffled)
    write(BIG / "shuffled.vcf", head + shuffled)
write(BIG / "trio.ped", ["fam1\tchild\tfather\tmother\t1\t0\n", "fam1\tfather\t0\t0\t1\t0\n",
                         "fam1\tmother\t0\t0\t2\t0\n"])

# 200,000 intervals of 50-1000 bases, unsorted and overlapping, as a list and as a BED.
spans = []
for i in range(200_000):
    start = rng.randrange(1, LENGTH - 1000)
    spans.append((start, start + rng.randrange(50, 1000), "+-"[i % 2], f"i{i}"))
write(BIG / "many.interval_list", [header_il] + [f"{CONTIG}\t{s}\t{e}\t{st}\t{n}\n" for s, e, st, n in spans])
write(BIG / "many.bed", [f"{CONTIG}\t{s - 1}\t{e}\t{n}\t0\t{st}\n" for s, e, st, n in spans])

# An identity chain (one block over the whole contig): LiftOver's work without a moved base.
write(BIG / "identity.chain", [f"chain 1000 {CONTIG} {LENGTH} + 0 {LENGTH} {CONTIG} {LENGTH} + 0 {LENGTH} 1\n",
                               f"{LENGTH}\n", "\n"])

# 100,000 GFF3 genes with one mRNA and two exons each, shuffled, for SortGff.
gff = []
for i in range(100_000):
    start = rng.randrange(1, LENGTH - 5000)
    end = start + rng.randrange(1000, 5000)
    gff.append(f"{CONTIG}\tsim\tgene\t{start}\t{end}\t.\t+\t.\tID=g{i}\n")
    gff.append(f"{CONTIG}\tsim\tmRNA\t{start}\t{end}\t.\t+\t.\tID=m{i};Parent=g{i}\n")
    gff.append(f"{CONTIG}\tsim\texon\t{start}\t{start + 200}\t.\t+\t.\tParent=m{i}\n")
    gff.append(f"{CONTIG}\tsim\texon\t{end - 200}\t{end}\t.\t+\t.\tParent=m{i}\n")
rng.shuffle(gff)
write(BIG / "big.gff3", ["##gff-version 3\n"] + gff)

# An exome-like target set (200 bases every 5 kb, 12,888 targets) and genes made of them.
targets = [(s, s + 199) for s in range(60_000, LENGTH - 60_000, 5000)]
write(D / "targets.interval_list", [header_il] + [f"{CONTIG}\t{s}\t{e}\t+\tt{i}\n" for i, (s, e) in enumerate(targets)])
refflat = []
for g in range(0, len(targets) - 4, 4):
    exons = targets[g:g + 4]
    refflat.append(f"gene{g}\ttx{g}\t{CONTIG}\t{'+-'[g % 8 == 0]}\t{exons[0][0] - 1}\t{exons[-1][1]}\t"
                   f"{exons[0][0] - 1}\t{exons[-1][1]}\t4\t{','.join(str(s - 1) for s, _ in exons)},\t"
                   f"{','.join(str(e) for _, e in exons)},\n")
write(D / "refflat.txt", refflat)

# Picard's own queryname sort and duplicate marking of the BAM, once.
if not (D / "qname.bam").exists():
    picard("SortSam I=/data/NA12878.chr20.bam O=/derived/qname.tmp.bam SORT_ORDER=queryname")
    os.rename(D / "qname.tmp.bam", D / "qname.bam")
if not (D / "md.bam").exists():
    picard("MarkDuplicates I=/data/NA12878.chr20.bam O=/derived/md.tmp.bam M=/derived/md.metrics")
    os.rename(D / "md.tmp.bam", D / "md.bam")

# The genotyping-array chip (500,000 loci) from the corpus's own writers, in the speed image.
arrays = D / "nonbam" / "arrays"
if not (arrays / "arr_s1.vcf").exists():
    arrays.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["docker", "run", "--rm", "-m", "6g", "-v", f"{HERE.parent}:/tools:ro", "-v", f"{D}:/derived",
         "picard-rs-speed:3.5.0",
         "mkdir -p /tmp/mba && cd /tmp/mba && cp /tools/coverage/*.java /tools/benchmark/MakeBigArrays.java . "
         "&& javac -cp $PICARD_JAR -d . *.java && java -cp .:$PICARD_JAR MakeBigArrays "
         "/derived/nonbam/arrays 500000"],
        check=True,
    )

subprocess.run([sys.executable, str(HERE.parent / "benchmark" / "make_nonbam_inputs.py"), str(D)], check=True)

# make_nonbam_inputs.py names its Illumina output prefixes by host path; a row writes in /work/out.
for path in (D / "nonbam").glob("ill_*.txt"):
    text = path.read_text()
    for suffix in ("_384", ""):
        text = text.replace(f"{D / 'nonbam'}/ill_out{suffix}/", "/work/out/")
    path.write_text(text)
