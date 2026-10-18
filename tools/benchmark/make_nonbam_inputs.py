#!/usr/bin/env python3
"""The inputs tools/benchmark/screen_nonbam.sh runs the tools that do not read a BAM on.

Scaled-up versions of the conformance corpus (tools/coverage/MakeFixtures.java and friends), large
enough that a run measures the work and not the startup. Reads the benchmark inputs already in
<bench>/ (bench.bam, qname.bam, bench.fasta) and <bench>/big/ (big.vcf, genome.fasta, many.*,
big.gff3, identity.chain), and writes <bench>/nonbam/. The genotyping-array chip comes from
MakeBigArrays.java, run in the oracle container (it reuses the corpus's binary writers and runs the
reference's GtcToVcf); this script then derives the zCall files from its VCF when it is there.

Usage: tools/benchmark/make_nonbam_inputs.py [bench dir]
"""
import os
import random
import shutil
import struct
import subprocess
import sys
from pathlib import Path

B = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home() / "GATK/picard-rs/target/benchmark")
BIG = B / "big"
N = B / "nonbam"
N.mkdir(exist_ok=True)
rng = random.Random(42)


def vcf_split(path):
    header, body = [], []
    with open(path) as f:
        for line in f:
            (header if line.startswith("#") else body).append(line)
    return header, body


header, body = vcf_split(BIG / "big.vcf")
print(f"big.vcf: {len(body)} records")


def write(name, lines):
    with open(N / name, "w") as f:
        f.writelines(lines)


# dbSNP: every other site of the callset, sites only, so known and novel are both non-empty.
meta = [h for h in header if not h.startswith("##FORMAT") and not h.startswith("#CHROM")]
write("dbsnp.vcf", meta + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
      + ["\t".join(r.split("\t")[:5]) + "\t.\t.\t.\n" for r in body[::2]])

# A header that lost two of the lines its records use, for FixVcfHeader.
write("big_missing.vcf", [h for h in header if "ID=AF," not in h and "ID=AD," not in h] + body)

# Four shards for GatherVcfs and MergeVcfs.
quarter = (len(body) + 3) // 4
for i in range(4):
    write(f"part{i + 1}.vcf", header + body[i * quarter:(i + 1) * quarter])

# One sample (the child) for RenameSampleInVcf.
one = header[:-1] + ["\t".join(header[-1].rstrip("\n").split("\t")[:10]) + "\n"]
write("big_one.vcf", one + ["\t".join(r.rstrip("\n").split("\t")[:10]) + "\n" for r in body])

# A second callset for GenotypeConcordance: every seventh child genotype flipped.
calls = []
for i, r in enumerate(body):
    if i % 7 == 3:
        f = r.split("\t")
        g = f[9].split(":")
        g[0] = {"0/0": "0/1", "0/1": "1/1", "1/1": "0/0"}.get(g[0], g[0])
        f[9] = ":".join(g)
        r = "\t".join(f)
    calls.append(r)
write("big_call.vcf", header + calls)

# A haplotype map over the callset's biallelic SNPs: one in twenty, in blocks of up to three
# (the first site of a block is the anchor of the others).
dict_lines = [l for l in open(B / "bench.dict") if l.startswith("@SQ")]
hap = ["@HD\tVN:1.6\tSO:coordinate\n"] + dict_lines
hap.append("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS\n")
sites = []
anchor = None
last = -10**9
for i, r in enumerate(body):
    if i % 20:
        continue
    f = r.split("\t", 6)
    if len(f[3]) != 1 or len(f[4]) != 1 or int(f[1]) - last < 2:
        continue
    pos = int(f[1])
    name = f"snp{pos}"
    anchor = None if (len(sites) % 3 == 0) else anchor
    maf = 0.05 + (pos % 41) / 100
    hap.append(f"{f[0]}\t{pos}\t{name}\t{f[3]}\t{f[4]}\t{maf:.2f}\t{anchor or ''}\t\n")
    sites.append((f[0], pos, name, f[3], f[4]))
    anchor = anchor or name
    last = pos
write("hap.txt", hap)
print(f"hap.txt: {len(sites)} sites")

# A fingerprint VCF of a hundred samples over the map's sites, GT and PL, for the crosscheck
# tools and, through their output, ClusterCrosscheckMetrics. Samples come in pairs that agree.
samples = [f"s{i:03d}" for i in range(100)]
fp = ["##fileformat=VCFv4.2\n", '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n',
      '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n']
fp += [l.replace("@SQ\tSN:", "##contig=<ID=").replace("\tLN:", ",length=").rstrip("\n") + ">\n"
       for l in dict_lines]
fp.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
gts = ["0/0", "0/1", "1/1"]
pls = ["0,30,300", "30,0,30", "300,30,0"]
for c, pos, name, ref, alt in sites:
    row = [c, str(pos), name, ref, alt, ".", "PASS", ".", "GT:PL"]
    for s in range(0, 100, 2):
        g = rng.randrange(3)
        row += [f"{gts[g]}:{pls[g]}"] * 2
    fp.append("\t".join(row) + "\n")
write("fp_many.vcf", fp)

# BaitDesigner targets: five thousand 300-base targets, one every 190 bases.
with open(N / "bait_targets.interval_list", "w") as f:
    f.write("@HD\tVN:1.6\n" + "".join(dict_lines))
    for i in range(5000):
        start = 1000 + i * 190
        f.write(f"chr1\t{start}\t{start + 299}\t+\tt{i}\n")

# Metrics files for the tools that read other tools' metrics.
for name in ("quality_yield_one.metrics",):
    src = Path(__file__).resolve().parents[2] / "target/corpus/fixtures" / name
    if src.exists():
        qd = N / "qy"
        qd.mkdir(exist_ok=True)
        for i in range(2000):
            shutil.copy(src, qd / f"q{i}.metrics")
write("bafregress.txt", ["sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom\n"] + [
    f"sample{i}\t{rng.uniform(-0.01, 0.02):.6f}\t0.0001\t{rng.uniform(-5, 5):.4f}\t1e-06\t0.99\t{1000 + i}\n"
    for i in range(20000)])
write("verifyidintensity.txt", ["ID\t%Mix\t\tLLK\t\tLLK0\t\n", "-" * 65 + "\n"] + [
    f"{i}\t{rng.uniform(0, 0.3):.6f}\t{rng.uniform(1e4, 2e5):.1f}\t{rng.uniform(1e4, 2e5):.1f}\n"
    for i in range(20000)])

# FASTQ from the benchmark's queryname BAM.
if not (N / "reads.fq").exists() or (N / "reads.fq").stat().st_size == 0:
    subprocess.run(f"samtools fastq -@4 -F 0x900 -1 {N}/r1.fq -2 {N}/r2.fq -0 {N}/reads.fq "
                   f"-s /dev/null {B}/qname.bam", shell=True, check=True, stderr=subprocess.DEVNULL)

# An Illumina run: one lane, two tiles of 200,000 clusters, 25T8B25T, on the corpus's layout
# (tools/coverage/MakeIlluminaCoverageFixtures.java): BCLs per cycle, filter and locs per tile,
# tile metrics in InterOp. Twenty-four barcodes; one cluster in ten has a mismatch in its barcode
# and one in fifty a no-call, and one in eight fails the filter.
CLUSTERS, STRUCTURE, CYCLES = 200_000, "25T8B25T", 58


def barcode_set(count, draw):
    """`count` eight-base barcodes, each at least three mismatches from every other."""
    found = []
    while len(found) < count:
        b = "".join(draw.choice("ACGT") for _ in range(8))
        if all(sum(x != y for x, y in zip(b, o)) >= 3 for o in found):
            found.append(b)
    return found


def illumina_run(run, tiles, barcodes):
    if (run / "Data/Intensities/BaseCalls/L001").exists():
        return
    root = run / "Data/Intensities/BaseCalls"
    lane = root / "L001"
    (run / "InterOp").mkdir(parents=True, exist_ok=True)
    # Phasing and prephasing (codes 200 + 2r, 201 + 2r) for the template reads, reads 1 and 3.
    tm = bytearray(b"\x02\x0a")
    for t in tiles:
        for code, value in ((100, CLUSTERS), (101, CLUSTERS * 7 // 8), (102, 1.5e6), (103, 1.3e6),
                            (200, 0.1), (201, 0.2), (204, 0.11), (205, 0.21)):
            tm += struct.pack("<HHHf", 1, t, code, value)
    (run / "InterOp/TileMetricsOut.bin").write_bytes(bytes(tm))
    (run / "RunInfo.xml").write_text(
        '<?xml version="1.0"?>\n<RunInfo><Run><Reads><Read Number="1" NumCycles="25" IsIndexedRead="N"/>'
        '<Read Number="2" NumCycles="8" IsIndexedRead="Y"/><Read Number="3" NumCycles="25" '
        'IsIndexedRead="N"/></Reads></Run></RunInfo>\n')
    locs = struct.pack("<ifi", 1, 1.0, CLUSTERS) + b"".join(
        struct.pack("<ff", 10.0 + (c % 2000) * 1.5, 10.0 + (c // 2000) * 1.5) for c in range(CLUSTERS))
    (root.parent / "s.locs").parent.mkdir(parents=True, exist_ok=True)
    (root.parent / "s.locs").write_bytes(locs)
    for t in tiles:
        trng = random.Random(t)
        reads = []
        for c in range(CLUSTERS):
            bc = list(barcodes[trng.randrange(len(barcodes))])
            if c % 10 == 3:
                bc[trng.randrange(8)] = trng.choice("ACGT")
            if c % 50 == 7:
                bc[trng.randrange(8)] = "N"
            r = [trng.choice("ACGT") for _ in range(25)] + bc + [trng.choice("ACGT") for _ in range(25)]
            reads.append(r)
        for cycle in range(CYCLES):
            calls = bytearray(4 + CLUSTERS)
            struct.pack_into("<i", calls, 0, CLUSTERS)
            for c in range(CLUSTERS):
                base = reads[c][cycle]
                if base != "N":
                    q = 2 + (c * 7 + cycle * 3) % 39
                    calls[4 + c] = "ACGT".index(base) | (q << 2)
            d = lane / f"C{cycle + 1}.1"
            d.mkdir(parents=True, exist_ok=True)
            (d / f"s_1_{t}.bcl").write_bytes(bytes(calls))
        filt = struct.pack("<iii", 0, 3, CLUSTERS) + bytes(0 if c % 8 == 5 else 1 for c in range(CLUSTERS))
        (lane / f"s_1_{t}.filter").write_bytes(filt)
        (lane / f"s_1_{t}.locs").write_bytes(locs)
    (root / "config.xml").write_text("<BaseCallAnalysis/>\n")


def barcode_files(suffix, barcodes):
    out = N / f"ill_out{suffix}"
    write(f"ill_barcodes{suffix}.txt", ["barcode_sequence_1\tbarcode_name\tlibrary_name\n"]
          + [f"{b}\tbc{i}\tlib{i}\n" for i, b in enumerate(barcodes)])
    write(f"ill_multiplex{suffix}.txt", ["OUTPUT_PREFIX\tBARCODE_1\n"]
          + [f"{out}/fq{i}\t{b}\n" for i, b in enumerate(barcodes)] + [f"{out}/fqrest\tN\n"])
    write(f"ill_library{suffix}.txt", ["OUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tBARCODE_1\n"]
          + [f"{out}/s{i}.bam\tsample{i}\tlib{i}\t{b}\n" for i, b in enumerate(barcodes)]
          + [f"{out}/rest.bam\trest\trest\tN\n"])


barcodes = barcode_set(24, rng)
run = N / "ill_big"
illumina_run(run, (1101, 1102), barcodes)
barcode_files("", barcodes)
# The same layout with a plate of 384 barcodes on one tile: the per-cluster work that grows with
# the barcode count is what the twenty-four of ill_big would not show.
plate = barcode_set(384, random.Random(384))
illumina_run(N / "ill_384", (1101,), plate)
barcode_files("_384", plate)

# A run directory with nothing but tile metrics, for CollectIlluminaLaneMetrics: four lanes of
# 2,500 tiles with eight codes each (80,000 records), the size of a NovaSeq flow cell's InterOp.
many = N / "ill_lanes_many"
(many / "InterOp").mkdir(parents=True, exist_ok=True)
tm = bytearray(b"\x02\x0a")
for lane in range(1, 5):
    for t in range(2500):
        tile = 1101 + (t // 50) * 100 + t % 50
        for code, value in ((100, 1000 + t), (101, 800 + t), (102, 1e4 + t), (103, 8e3 + t),
                            (200, 0.1 + t / 1e5), (201, 0.2), (204, 0.11), (205, 0.21 + t / 1e5)):
            tm += struct.pack("<HHHf", lane, tile, code, value)
(many / "InterOp/TileMetricsOut.bin").write_bytes(bytes(tm))
(many / "RunInfo.xml").write_text((run / "RunInfo.xml").read_text())

# zCall files over the array chip's first sample, when MakeBigArrays has written it.
arr = N / "arrays/arr_s1.vcf"
if arr.exists():
    _, arows = vcf_split(arr)
    mp, ped, thr = [], ["FAM1\tsample1\t0\t0\t2\t-9"], ["Name\tThr_X\tThr_Y\n"]
    for i, r in enumerate(arows):
        f = r.split("\t")
        mp.append(f"{f[0]}\t{f[2]}\t0\t{f[1]}\n")
        ped.append(["A\tA", "A\tB", "B\tB", "0\t0"][i % 4])
        thr.append(f"{f[2]}\t{'NA' if i % 9 == 0 else 0.2 + (i % 5) / 10}\t{'NA' if i % 9 == 0 else 0.3 + (i % 3) / 10}\n")
    write("zcall.map", mp)
    write("zcall.ped", ["\t".join(ped) + "\n"])
    write("zcall_thresholds.txt", thr)
print("done")
