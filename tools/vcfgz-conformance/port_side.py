#!/usr/bin/env python3
"""The port's side of VcfGzDump: the same inputs, the same cases, the same rows.

    python3 tools/vcfgz-conformance/port_side.py target/release > port.txt
    diff <(grep -v '^#' tools/conformance/pending/vcfgz.VcfGzDump.txt) port.txt

The inputs are written by the formula VcfGzDump.java uses, and each one's MD5 is a row, so a
generator that drifted shows as an `input` row before any output is compared. Each case runs the
binary named after the tool, under the arguments the dump hands `instanceMain`.
"""

import hashlib
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

CONTIGS = ["chr1", "chr2", "chr3"]
LENGTH = 60000
PER_CONTIG = 4000
STEP = 13
BASES = "ACGT"
ROWS = []


def emit(kind, name, payload):
    payload = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")
    ROWS.append(f"{kind}\t{name}\t{payload}")


def base(contig, position):
    return BASES[(position * 31 + position // 7 + contig * 5) % 4]


def nxt(b):
    return BASES[(BASES.index(b) + 1) % 4]


def record(i, contig, position, samples):
    indel = i % 5 == 0
    first = nxt(base(contig, position)) if i % 17 == 0 else base(contig, position)
    ref = first + base(contig, position + 1) if indel else first
    alt = first if indel else nxt(first)
    line = (
        f"{CONTIGS[contig]}\t{position}\t{'rs' + str(i) if i % 3 == 0 else '.'}\t{ref}\t{alt}\t"
        f"{i % 997}.{i % 10}\t{'LowQual' if i % 11 == 0 else 'PASS'}\t"
        f"DP={i * 7 % 200};AC=1;AF=0.167;AN=6\tGT:AD:DP:GQ:PL"
    )
    gts = ["0/0", "0/1", "1/1", "./."]
    for s in range(samples):
        g = (i + s) % 4
        line += "\t" + gts[g]
        if g != 3:
            ref_depth, alt_depth = 10 + i % 30, i % 20
            line += (
                f":{ref_depth},{alt_depth}:{ref_depth + alt_depth}:{i * 13 % 99}"
                f":0,{i % 50},{i % 90 + 10}"
            )
    return line + "\n"


def header(contig_lines, samples):
    h = "##fileformat=VCFv4.2\n"
    for c in contig_lines:
        h += f"##contig=<ID={c},length={LENGTH}>\n"
    h += (
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    )
    for name in ["child", "father", "mother"][:samples]:
        h += "\t" + name
    return h + "\n"


def body(contigs, samples, parity=-1, swap=-1):
    lines = []
    i = 0
    for c in contigs:
        for k in range(PER_CONTIG):
            if parity < 0 or i % 2 == parity:
                lines.append(record(i, c, 11 + k * STEP, samples))
            i += 1
    if swap >= 0:
        lines[swap], lines[swap + 1] = lines[swap + 1], lines[swap]
    return "".join(lines)


def md5(data):
    return hashlib.md5(data).hexdigest()


def inputs(d):
    def put(name, text):
        data = text.encode()
        (d / name).write_bytes(data)
        emit("input", name, f"{md5(data)}\t{len(data)}")

    order = [0, 1, 2]
    put("in.vcf", header(CONTIGS, 3) + body(order, 3))
    put("unsorted.vcf", header(CONTIGS, 3) + body([2, 0, 1], 3))
    put("swapped.vcf", header(CONTIGS, 3) + body(order, 3, swap=100))
    put("nodict.vcf", header([], 3) + body(order, 3))
    put("missing.vcf", header(["chr1", "chr2"], 3) + body(order, 3))
    put("single.vcf", header(CONTIGS, 1) + body(order, 1))
    put("single-nodict.vcf", header([], 1) + body(order, 1))
    put("m1.vcf", header(CONTIGS, 3) + body(order, 3, parity=0))
    put("m2.vcf", header(CONTIGS, 3) + body(order, 3, parity=1))
    for c in range(3):
        put(f"p{c + 1}.vcf", header(CONTIGS, 3) + body([c], 3))
    fasta, dictionary, other, chain = "", "@HD\tVN:1.6\n", "@HD\tVN:1.6\n", ""
    for c, name in enumerate(CONTIGS):
        fasta += f">{name}\n"
        seq = "".join(base(c, p) for p in range(1, LENGTH + 1))
        fasta += "".join(seq[i : i + 60] + "\n" for i in range(0, LENGTH, 60))
        dictionary += f"@SQ\tSN:{name}\tLN:{LENGTH}\n"
        other += f"@SQ\tSN:{name}\tLN:{2 * LENGTH}\n"
        chain += f"chain 1000 {name} {LENGTH} + 0 {LENGTH} {name} {LENGTH} + 0 {LENGTH} {c + 1}\n"
        chain += f"{LENGTH}\n\n"
    put("ref.fasta", fasta)
    put("ref.dict", dictionary)
    put("other.dict", other)
    put("id.chain", chain)


def inflate_bgzf(data):
    """A BGZF file's contents: each member inflated in turn."""
    out = b""
    while data:
        d = zlib.decompressobj(31)
        out += d.decompress(data)
        data = d.unused_data
    return out


def java_exception(stderr):
    """The uncaught exception's `Class: message`, or the logged error's first line."""
    for line in stderr.splitlines():
        if line.startswith('Exception in thread "main" '):
            return line[len('Exception in thread "main" ') :]
    return None


def run(d, bins, name, tool, *args):
    out = d / f"case-{name}"
    out.mkdir()
    argv = [a.replace("@", str(d)).replace("%", str(out)) for a in args]
    done = subprocess.run([str(bins / tool), *argv], capture_output=True, text=True, cwd=out)
    exception = java_exception(done.stderr)
    emit("exit", name, exception if exception is not None else str(done.returncode))
    for f in sorted(out.iterdir()):
        key = f"{name}/{f.name}"
        if f.name.endswith(".idx"):
            emit("file", key, "present")
            continue
        data = f.read_bytes()
        emit("file", key, f"{md5(data)}\t{len(data)}")
        if f.name.endswith(".tbi"):
            emit("tbi", key, inflate_bgzf(data).hex())


def main():
    bins = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="vcfgz") as tmp:
        d = Path(tmp)
        inputs(d)
        for name, extra in [
            ("sort", []),
            ("sort-no-index", ["CREATE_INDEX=false"]),
            ("sort-level-3", ["COMPRESSION_LEVEL=3"]),
            ("sort-level-9", ["COMPRESSION_LEVEL=9"]),
            ("sort-jdk-deflater", ["USE_JDK_DEFLATER=true"]),
            ("sort-jdk-deflater-level-1", ["USE_JDK_DEFLATER=true", "COMPRESSION_LEVEL=1"]),
        ]:
            run(d, bins, name, "SortVcf", "I=@/unsorted.vcf", "O=%/out.vcf.gz", *extra)
        run(d, bins, "sort-bgz", "SortVcf", "I=@/unsorted.vcf", "O=%/out.vcf.bgz")
        run(d, bins, "sort-gz-input", "SortVcf", "I=@/case-sort/out.vcf.gz", "O=%/out.vcf.gz")

        run(d, bins, "sites-only", "MakeSitesOnlyVcf", "I=@/in.vcf", "O=%/out.vcf.gz")
        run(d, bins, "sites-only-out-of-order", "MakeSitesOnlyVcf", "I=@/swapped.vcf",
            "O=%/out.vcf.gz")
        run(d, bins, "sites-only-missing-contig", "MakeSitesOnlyVcf", "I=@/missing.vcf",
            "O=%/out.vcf.gz")

        run(d, bins, "rename", "RenameSampleInVcf", "I=@/single.vcf", "O=%/out.vcf.gz",
            "NEW_SAMPLE_NAME=renamed")
        run(d, bins, "rename-no-dictionary", "RenameSampleInVcf", "I=@/single-nodict.vcf",
            "O=%/out.vcf.gz", "NEW_SAMPLE_NAME=renamed", "CREATE_INDEX=true")

        run(d, bins, "update-dictionary", "UpdateVcfSequenceDictionary", "I=@/in.vcf",
            "O=%/out.vcf.gz", "SEQUENCE_DICTIONARY=@/other.dict", "CREATE_INDEX=true")

        run(d, bins, "filter", "FilterVcf", "I=@/in.vcf", "O=%/out.vcf.gz", "MIN_DP=30")
        run(d, bins, "filter-no-dictionary", "FilterVcf", "I=@/nodict.vcf", "O=%/out.vcf.gz",
            "MIN_DP=30")

        run(d, bins, "fix-header", "FixVcfHeader", "I=@/in.vcf", "O=%/out.vcf.gz")

        run(d, bins, "merge", "MergeVcfs", "I=@/m1.vcf", "I=@/m2.vcf", "O=%/out.vcf.gz")

        run(d, bins, "split", "SplitVcfs", "I=@/in.vcf", "SNP_OUTPUT=%/snp.vcf.gz",
            "INDEL_OUTPUT=%/indel.vcf.gz", "STRICT=false")

        run(d, bins, "liftover", "LiftoverVcf", "I=@/in.vcf", "O=%/out.vcf.gz",
            "REJECT=%/reject.vcf.gz", "CHAIN=@/id.chain", "R=@/ref.fasta")

        run(d, bins, "convert", "VcfFormatConverter", "I=@/in.vcf", "O=%/out.vcf.gz",
            "REQUIRE_INDEX=false")
        run(d, bins, "convert-gz-to-gz", "VcfFormatConverter", "I=@/case-sort/out.vcf.gz",
            "O=%/out.vcf.gz")
        run(d, bins, "convert-gz-to-vcf", "VcfFormatConverter", "I=@/case-sort/out.vcf.gz",
            "O=%/out.vcf")

        for c in range(1, 4):
            run(d, bins, f"part{c}", "SortVcf", f"I=@/p{c}.vcf", f"O=%/p{c}.vcf.gz")
        run(d, bins, "gather-block-copy", "GatherVcfs", "I=@/case-part1/p1.vcf.gz",
            "I=@/case-part2/p2.vcf.gz", "I=@/case-part3/p3.vcf.gz", "O=%/out.vcf.gz")
        run(d, bins, "gather-conventional", "GatherVcfs", "I=@/p1.vcf", "I=@/p2.vcf",
            "I=@/p3.vcf", "O=%/out.vcf.gz")
    print("\n".join(ROWS))


if __name__ == "__main__":
    main()
