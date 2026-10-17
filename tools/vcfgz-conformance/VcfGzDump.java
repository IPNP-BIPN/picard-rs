/*
 * Every Picard tool that writes a VCF, asked for a `.vcf.gz`, taken from the reference.
 *
 * `VariantContextWriterBuilder` decides the container from the name: `.vcf.gz` and `.vcf.bgz` are a
 * BLOCK_COMPRESSED_VCF, a `VCFWriter` over a `BlockCompressedOutputStream` at the level and deflater
 * `instanceMain` installed, and the builder always gives it a `TabixIndexCreator`, which indexes
 * when `INDEX_ON_THE_FLY` is set and writes a `.tbi` at `close()`. What is measured is the bytes of
 * both files, so the deflater, the level, the block boundaries and the index's chunks all count.
 *
 * Behaviours this is built to catch.
 *
 *   - THE `.vcf.gz` AND ITS `.tbi` ARE BYTE-IDENTICAL for SortVcf, MakeSitesOnlyVcf,
 *     RenameSampleInVcf, UpdateVcfSequenceDictionary, FilterVcf, FixVcfHeader, MergeVcfs,
 *     SplitVcfs, LiftoverVcf (its REJECT a `.vcf.gz` too) and VcfFormatConverter;
 *   - `CREATE_INDEX=false` WRITES NO `.tbi`, and `COMPRESSION_LEVEL` AND `USE_JDK_DEFLATER` MOVE THE
 *     BYTES OF BOTH FILES, the index being written through the same static deflater;
 *   - A `.vcf.bgz` IS THE SAME CONTAINER under another name;
 *   - THE TABIX CREATOR NEEDS NO DICTIONARY where the `.idx` one does: FilterVcf sets one only for a
 *     `.vcf` name, and RenameSampleInVcf's builder refuses an index without one only for `VCF`, so
 *     a file with no contig lines is still indexed;
 *   - IT CHECKS THE ORDER, which the `.idx` creator does not: MakeSitesOnlyVcf over a record that
 *     goes back on its contig is an IllegalArgumentException naming both features by their virtual
 *     offsets, and the file holds what the 128 KiB buffer had passed on;
 *   - A CONTIG THE DICTIONARY LACKS IS A NullPointerException when its first record opens it;
 *   - GatherVcfs OVER `.vcf.gz` INPUTS COPIES BLOCKS: the first file as it is, each later file from
 *     the block its header ends in, deflated again, and one terminator, with no `.tbi` whatever
 *     CREATE_INDEX says; over plain inputs it gathers conventionally and indexes.
 *
 * The inputs are written by a formula, not a random generator, so the port's side can write the
 * same bytes (tools/vcfgz-conformance/inputs.py); each input's MD5 is a row, so a drift in either
 * generator shows before any output is compared.
 *
 * Output:
 *
 *     input\t<file>\t<md5>\t<size>
 *     exit\t<case>\t<status, or the exception class and message>
 *     file\t<case>/<file>\t<md5>\t<size>       (a `.idx` is `present`: it records the file's mtime)
 *     tbi\t<case>/<file>\t<the index, inflated, in hex>
 *
 * Usage: VcfGzDump
 */

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public class VcfGzDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final String[] CONTIGS = {"chr1", "chr2", "chr3"};
    static final int LENGTH = 60000;
    static final int PER_CONTIG = 4000;
    static final int STEP = 13;
    static final String BASES = "ACGT";

    static char base(final int contig, final int position) {
        return BASES.charAt((position * 31 + position / 7 + contig * 5) % 4);
    }

    static char next(final char b) {
        return BASES.charAt((BASES.indexOf(b) + 1) % 4);
    }

    /** One record line. `i` runs over every record of the file; `contig` indexes CONTIGS. */
    static String record(final int i, final int contig, final int position, final int samples) {
        final boolean indel = i % 5 == 0;
        // Every seventeenth reference base is wrong, so LiftoverVcf has rejects to write.
        final char first = i % 17 == 0 ? next(base(contig, position)) : base(contig, position);
        final String ref = indel ? "" + first + base(contig, position + 1) : "" + first;
        final String alt = indel ? "" + first : "" + next(first);
        final StringBuilder line = new StringBuilder();
        line.append(CONTIGS[contig]).append('\t').append(position).append('\t')
                .append(i % 3 == 0 ? "rs" + i : ".").append('\t').append(ref).append('\t')
                .append(alt).append('\t').append(i % 997).append('.').append(i % 10).append('\t')
                .append(i % 11 == 0 ? "LowQual" : "PASS").append('\t')
                .append("DP=").append(i * 7 % 200).append(";AC=1;AF=0.167;AN=6")
                .append("\tGT:AD:DP:GQ:PL");
        final String[] gts = {"0/0", "0/1", "1/1", "./."};
        for (int s = 0; s < samples; s++) {
            final int g = (i + s) % 4;
            line.append('\t').append(gts[g]);
            if (g != 3) {
                final int ref_depth = 10 + i % 30;
                final int alt_depth = i % 20;
                line.append(':').append(ref_depth).append(',').append(alt_depth)
                        .append(':').append(ref_depth + alt_depth)
                        .append(':').append(i * 13 % 99)
                        .append(":0,").append(i % 50).append(',').append(i % 90 + 10);
            }
        }
        return line.append('\n').toString();
    }

    static String header(final List<String> contigLines, final int samples) {
        final StringBuilder h = new StringBuilder("##fileformat=VCFv4.2\n");
        for (final String c : contigLines) {
            h.append("##contig=<ID=").append(c).append(",length=").append(LENGTH).append(">\n");
        }
        h.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n")
                .append("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\">\n")
                .append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n")
                .append("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"AN\">\n")
                .append("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
                .append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
                .append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n")
                .append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n")
                .append("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n")
                .append("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n")
                .append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        final String[] names = {"child", "father", "mother"};
        for (int s = 0; s < samples; s++) {
            h.append('\t').append(names[s]);
        }
        return h.append('\n').toString();
    }

    /**
     * The records of the given contigs in the given order, `PER_CONTIG` each; `keep` selects by the
     * running index (`i % 2 == parity` for MergeVcfs' halves, -1 for all), and `swap` exchanges
     * the record at that running index with the next one.
     */
    static String body(final int[] contigs, final int samples, final int parity, final int swap) {
        final List<String> lines = new ArrayList<>();
        int i = 0;
        for (final int c : contigs) {
            for (int k = 0; k < PER_CONTIG; k++, i++) {
                if (parity >= 0 && i % 2 != parity) {
                    continue;
                }
                lines.add(record(i, c, 11 + k * STEP, samples));
            }
        }
        if (swap >= 0) {
            final String a = lines.get(swap);
            lines.set(swap, lines.get(swap + 1));
            lines.set(swap + 1, a);
        }
        return String.join("", lines);
    }

    static final List<String> ALL = Arrays.asList(CONTIGS);
    static final int[] ORDER = {0, 1, 2};

    static Path dir;

    static void input(final String name, final String text) throws Exception {
        final Path file = dir.resolve(name);
        Files.writeString(file, text, StandardCharsets.UTF_8);
        emit("input", name, md5(Files.readAllBytes(file)) + "\t" + Files.size(file));
    }

    static void inputs() throws Exception {
        input("in.vcf", header(ALL, 3) + body(ORDER, 3, -1, -1));
        input("unsorted.vcf", header(ALL, 3) + body(new int[]{2, 0, 1}, 3, -1, -1));
        input("swapped.vcf", header(ALL, 3) + body(ORDER, 3, -1, 100));
        input("nodict.vcf", header(List.of(), 3) + body(ORDER, 3, -1, -1));
        input("missing.vcf", header(List.of("chr1", "chr2"), 3) + body(ORDER, 3, -1, -1));
        input("single.vcf", header(ALL, 1) + body(ORDER, 1, -1, -1));
        input("single-nodict.vcf", header(List.of(), 1) + body(ORDER, 1, -1, -1));
        input("m1.vcf", header(ALL, 3) + body(ORDER, 3, 0, -1));
        input("m2.vcf", header(ALL, 3) + body(ORDER, 3, 1, -1));
        for (int c = 0; c < 3; c++) {
            input("p" + (c + 1) + ".vcf", header(ALL, 3) + body(new int[]{c}, 3, -1, -1));
        }
        final StringBuilder fasta = new StringBuilder();
        final StringBuilder dict = new StringBuilder("@HD\tVN:1.6\n");
        final StringBuilder other = new StringBuilder("@HD\tVN:1.6\n");
        final StringBuilder chain = new StringBuilder();
        for (int c = 0; c < 3; c++) {
            fasta.append('>').append(CONTIGS[c]).append('\n');
            for (int p = 1; p <= LENGTH; p++) {
                fasta.append(base(c, p));
                if (p % 60 == 0 || p == LENGTH) {
                    fasta.append('\n');
                }
            }
            dict.append("@SQ\tSN:").append(CONTIGS[c]).append("\tLN:").append(LENGTH).append('\n');
            other.append("@SQ\tSN:").append(CONTIGS[c]).append("\tLN:").append(2 * LENGTH)
                    .append('\n');
            chain.append("chain 1000 ").append(CONTIGS[c]).append(' ').append(LENGTH)
                    .append(" + 0 ").append(LENGTH).append(' ').append(CONTIGS[c]).append(' ')
                    .append(LENGTH).append(" + 0 ").append(LENGTH).append(' ').append(c + 1)
                    .append('\n').append(LENGTH).append("\n\n");
        }
        input("ref.fasta", fasta.toString());
        input("ref.dict", dict.toString());
        input("other.dict", other.toString());
        input("id.chain", chain.toString());
    }

    static String md5(final byte[] bytes) throws Exception {
        final StringBuilder hex = new StringBuilder();
        for (final byte b : MessageDigest.getInstance("MD5").digest(bytes)) {
            hex.append(String.format("%02x", b));
        }
        return hex.toString();
    }

    static String inflated(final Path file) throws Exception {
        try (InputStream in = new BlockCompressedInputStream(new FileInputStream(file.toFile()))) {
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            in.transferTo(out);
            final StringBuilder hex = new StringBuilder();
            for (final byte b : out.toByteArray()) {
                hex.append(String.format("%02x", b));
            }
            return hex.toString();
        }
    }

    interface Tool {
        int run(String[] argv);
    }

    /**
     * One run in its own directory: `@` in an argument is the inputs' directory, `%` this case's.
     * The deflater factory is put back first, because `instanceMain` installs GKL only when
     * USE_JDK_DEFLATER is false and never puts the default back.
     */
    static void run(final String name, final Tool tool, final String... args) throws Exception {
        final Path out = Files.createDirectory(dir.resolve("case-" + name));
        final String[] argv = Arrays.stream(args)
                .map(a -> a.replace("@", dir.toString()).replace("%", out.toString()))
                .toArray(String[]::new);
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        String status;
        try {
            status = String.valueOf(tool.run(argv));
        } catch (final Throwable e) {
            status = e.getClass().getName() + ": " + e.getMessage();
        }
        emit("exit", name, status);
        try (Stream<Path> files = Files.list(out)) {
            for (final Path f : files.sorted().toList()) {
                final String key = name + "/" + f.getFileName();
                if (f.toString().endsWith(".idx")) {
                    emit("file", key, "present");
                    continue;
                }
                emit("file", key, md5(Files.readAllBytes(f)) + "\t" + Files.size(f));
                if (f.toString().endsWith(".tbi")) {
                    emit("tbi", key, inflated(f));
                }
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        dir = Files.createTempDirectory("vcfgz");
        inputs();

        final Tool sort = a -> new picard.vcf.SortVcf().instanceMain(a);
        run("sort", sort, "I=@/unsorted.vcf", "O=%/out.vcf.gz");
        run("sort-no-index", sort, "I=@/unsorted.vcf", "O=%/out.vcf.gz", "CREATE_INDEX=false");
        // GKL routes levels 1 and 2 through ISA-L, which only an x86-64 port build reaches; level 3
        // is GKL's own deflate, and level 1 is measured through libdeflate instead.
        run("sort-level-3", sort, "I=@/unsorted.vcf", "O=%/out.vcf.gz", "COMPRESSION_LEVEL=3");
        run("sort-level-9", sort, "I=@/unsorted.vcf", "O=%/out.vcf.gz", "COMPRESSION_LEVEL=9");
        run("sort-jdk-deflater", sort, "I=@/unsorted.vcf", "O=%/out.vcf.gz",
                "USE_JDK_DEFLATER=true");
        run("sort-jdk-deflater-level-1", sort, "I=@/unsorted.vcf", "O=%/out.vcf.gz",
                "USE_JDK_DEFLATER=true", "COMPRESSION_LEVEL=1");
        run("sort-bgz", sort, "I=@/unsorted.vcf", "O=%/out.vcf.bgz");
        run("sort-gz-input", sort, "I=@/case-sort/out.vcf.gz", "O=%/out.vcf.gz");

        final Tool sitesOnly = a -> new picard.vcf.MakeSitesOnlyVcf().instanceMain(a);
        run("sites-only", sitesOnly, "I=@/in.vcf", "O=%/out.vcf.gz");
        run("sites-only-out-of-order", sitesOnly, "I=@/swapped.vcf", "O=%/out.vcf.gz");
        run("sites-only-missing-contig", sitesOnly, "I=@/missing.vcf", "O=%/out.vcf.gz");

        final Tool rename = a -> new picard.vcf.RenameSampleInVcf().instanceMain(a);
        run("rename", rename, "I=@/single.vcf", "O=%/out.vcf.gz", "NEW_SAMPLE_NAME=renamed");
        run("rename-no-dictionary", rename, "I=@/single-nodict.vcf", "O=%/out.vcf.gz",
                "NEW_SAMPLE_NAME=renamed", "CREATE_INDEX=true");

        run("update-dictionary", a -> new picard.vcf.UpdateVcfSequenceDictionary().instanceMain(a),
                "I=@/in.vcf", "O=%/out.vcf.gz", "SEQUENCE_DICTIONARY=@/other.dict",
                "CREATE_INDEX=true");

        final Tool filter = a -> new picard.vcf.filter.FilterVcf().instanceMain(a);
        run("filter", filter, "I=@/in.vcf", "O=%/out.vcf.gz", "MIN_DP=30");
        run("filter-no-dictionary", filter, "I=@/nodict.vcf", "O=%/out.vcf.gz", "MIN_DP=30");

        run("fix-header", a -> new picard.vcf.FixVcfHeader().instanceMain(a),
                "I=@/in.vcf", "O=%/out.vcf.gz");

        run("merge", a -> new picard.vcf.MergeVcfs().instanceMain(a),
                "I=@/m1.vcf", "I=@/m2.vcf", "O=%/out.vcf.gz");

        run("split", a -> new picard.vcf.SplitVcfs().instanceMain(a), "I=@/in.vcf",
                "SNP_OUTPUT=%/snp.vcf.gz", "INDEL_OUTPUT=%/indel.vcf.gz", "STRICT=false");

        run("liftover", a -> new picard.vcf.LiftoverVcf().instanceMain(a), "I=@/in.vcf",
                "O=%/out.vcf.gz", "REJECT=%/reject.vcf.gz", "CHAIN=@/id.chain",
                "R=@/ref.fasta");

        final Tool converter = a -> new picard.vcf.VcfFormatConverter().instanceMain(a);
        run("convert", converter, "I=@/in.vcf", "O=%/out.vcf.gz", "REQUIRE_INDEX=false");
        run("convert-gz-to-gz", converter, "I=@/case-sort/out.vcf.gz", "O=%/out.vcf.gz");
        run("convert-gz-to-vcf", converter, "I=@/case-sort/out.vcf.gz", "O=%/out.vcf");

        // The gathered parts are this run's own SortVcf outputs, so the port's side makes them with
        // its own SortVcf, which the `sort` case checks.
        for (int c = 1; c <= 3; c++) {
            run("part" + c, sort, "I=@/p" + c + ".vcf", "O=%/p" + c + ".vcf.gz");
        }
        final Tool gather = a -> new picard.vcf.GatherVcfs().instanceMain(a);
        run("gather-block-copy", gather, "I=@/case-part1/p1.vcf.gz", "I=@/case-part2/p2.vcf.gz",
                "I=@/case-part3/p3.vcf.gz", "O=%/out.vcf.gz");
        run("gather-conventional", gather, "I=@/p1.vcf", "I=@/p2.vcf", "I=@/p3.vcf",
                "O=%/out.vcf.gz");

        System.out.print(buf);
    }
}
