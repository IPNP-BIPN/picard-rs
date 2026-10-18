/*
 * Oracle dump for records htsjdk did not write: a BAM built byte by byte the way samtools writes
 * one (tags in insertion order, small positive integers as unsigned types, an H tag), run through
 * the tools that pass records through, some of which set something on them and some of which do
 * not. htsjdk writes a record nothing was set on from the bytes it read, so its tags keep the
 * input's order and types; a record that was touched is re-encoded, sorted and canonical.
 *
 * Each case emits the tool, its arguments and the records of the output, hex, without the header
 * (whose @PG carries paths). The input BAM itself is the first row.
 *
 *   java -cp picard.jar:. RawRecordsDump
 */
import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.GZIPInputStream;

public class RawRecordsDump {
    static final String HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n"
        + "@RG\tID:rg1\tSM:s1\tLB:lib1\tPL:ILLUMINA\n";

    public static void main(final String[] args) throws Exception {
        final File input = File.createTempFile("raw-in-", ".bam");
        input.deleteOnExit();
        writeInput(input);
        emit("input_bam", "input", hex(Files.readAllBytes(input.toPath())));
        final File dict = File.createTempFile("raw-dict-", ".dict");
        dict.deleteOnExit();
        Files.write(dict.toPath(), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n".getBytes(StandardCharsets.US_ASCII));

        run("sortsam_queryname", "SortSam", input, "SORT_ORDER=queryname");
        run("sortsam_coordinate", "SortSam", input, "SORT_ORDER=coordinate");
        run("markduplicates_no_pg_tag", "markduplicates.MarkDuplicates", input, "METRICS_FILE=" + temp(".txt"), "ADD_PG_TAG_TO_READS=false");
        run("markduplicates_pg_tag", "markduplicates.MarkDuplicates", input, "METRICS_FILE=" + temp(".txt"));
        run("mergesamfiles", "MergeSamFiles", input);
        run("samformatconverter", "SamFormatConverter", input);
        run("cleansam", "CleanSam", input);
        run("reordersam", "ReorderSam", input, "SEQUENCE_DICTIONARY=" + dict.getAbsolutePath());
        run("downsamplesam", "DownsampleSam", input, "P=1");
        run("fixmateinformation", "FixMateInformation", input);
        runSplit("splitsambylibrary", "SplitSamByLibrary", input);
        runSplit("splitsambynumberofreads", "SplitSamByNumberOfReads", input, "SPLIT_TO_N_FILES=3");
        run("addorreplacereadgroups", "AddOrReplaceReadGroups", input, "RGID=rg2", "RGLB=lib2", "RGPL=ILLUMINA", "RGPU=u", "RGSM=s2");
    }

    static String temp(final String ext) throws Exception {
        final File f = File.createTempFile("raw-side-", ext);
        f.deleteOnExit();
        return f.getAbsolutePath();
    }

    static void run(final String kase, final String tool, final File input, final String... extra) throws Exception {
        final File out = File.createTempFile("raw-out-", ".bam");
        out.delete();
        out.deleteOnExit();
        final String[] argv = new String[2 + extra.length];
        argv[0] = "INPUT=" + input.getAbsolutePath();
        argv[1] = "OUTPUT=" + out.getAbsolutePath();
        System.arraycopy(extra, 0, argv, 2, extra.length);
        emit("tool", kase, tool.substring(tool.lastIndexOf('.') + 1));
        final StringBuilder shown = new StringBuilder();
        for (final String a : extra) {
            shown.append(shown.length() == 0 ? "" : " ").append(a.contains("raw-") ? a.substring(0, a.indexOf('=') + 1) + "<side>" : a);
        }
        emit("args", kase, shown.toString());
        final Object program = Class.forName("picard.sam." + tool).getDeclaredConstructor().newInstance();
        final int rc = ((picard.cmdline.CommandLineProgram) program).instanceMain(argv);
        emit("exit", kase, Integer.toString(rc));
        emit("records", kase, out.exists() ? hex(records(out)) : "");
    }

    /** A tool that writes a directory of BAMs: each file's records are a row. */
    static void runSplit(final String kase, final String tool, final File input, final String... extra) throws Exception {
        final File dir = Files.createTempDirectory("raw-split-").toFile();
        final String[] argv = new String[2 + extra.length];
        argv[0] = "INPUT=" + input.getAbsolutePath();
        argv[1] = "OUTPUT=" + dir.getAbsolutePath();
        System.arraycopy(extra, 0, argv, 2, extra.length);
        emit("tool", kase, tool);
        emit("args", kase, String.join(" ", extra));
        final Object program = Class.forName("picard.sam." + tool).getDeclaredConstructor().newInstance();
        emit("exit", kase, Integer.toString(((picard.cmdline.CommandLineProgram) program).instanceMain(argv)));
        final File[] files = dir.listFiles();
        java.util.Arrays.sort(files);
        for (final File f : files) {
            emit("file:" + f.getName(), kase, hex(records(f)));
            f.delete();
        }
        dir.delete();
    }

    /** The decompressed BAM after its header and reference list. */
    static byte[] records(final File bam) throws Exception {
        final byte[] all;
        try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(bam))) { all = in.readAllBytes(); }
        final ByteBuffer b = ByteBuffer.wrap(all).order(ByteOrder.LITTLE_ENDIAN);
        int p = 4;
        final int lText = b.getInt(p); p += 4 + lText;
        final int nRef = b.getInt(p); p += 4;
        for (int i = 0; i < nRef; i++) { final int l = b.getInt(p); p += 4 + l + 4; }
        final byte[] out = new byte[all.length - p];
        System.arraycopy(all, p, out, 0, out.length);
        return out;
    }

    /** A coordinate-sorted BAM of mapped pairs, every record's tags written as samtools would. */
    static void writeInput(final File file) throws Exception {
        final ByteArrayOutputStream body = new ByteArrayOutputStream();
        final ByteBuffer h = le(12 + HEADER.length() + 4 + 4 + 5 + 4);
        h.put("BAM\1".getBytes(StandardCharsets.US_ASCII)).putInt(HEADER.length()).put(HEADER.getBytes(StandardCharsets.US_ASCII));
        h.putInt(1).putInt(5).put("chr1\0".getBytes(StandardCharsets.US_ASCII)).putInt(10000);
        body.write(h.array(), 0, h.position());
        long seed = 11;
        final java.util.List<Object[]> records = new java.util.ArrayList<>();
        for (int i = 0; i < 40; i++) {
            final int pair = i / 2;
            final boolean first = i % 2 == 0;
            final int start = 100 + pair * 37 + (first ? 0 : 180);
            final int mateStart = 100 + pair * 37 + (first ? 180 : 0);
            final int flag = 1 | 2 | (first ? 64 | 32 : 128 | 16);
            final String name = "M01:7:FC1:1:1101:" + (1000 + pair * 13 % 7) + ":" + (2000 + pair);
            final StringBuilder seq = new StringBuilder();
            final byte[] qual = new byte[30];
            for (int k = 0; k < 30; k++) {
                seed = seed * 6364136223846793005L + 1442695040888963407L;
                seq.append("ACGT".charAt((int) ((seed >>> 33) & 3)));
                qual[k] = (byte) (20 + (int) ((seed >>> 40) % 20));
            }
            // samtools order: whatever the aligner appended, unsorted, with unsigned small types.
            final ByteArrayOutputStream tags = new ByteArrayOutputStream();
            tag(tags, "NM", 'C', 1 + pair % 3);
            tag(tags, "MQ", 'C', 58);
            if (pair % 3 == 0) tag(tags, "XS", 'S', 300 + pair);
            tags.write("MCZ30M\0".getBytes(StandardCharsets.US_ASCII));
            if (pair % 5 == 0) tags.write("XHH0A1F\0".getBytes(StandardCharsets.US_ASCII));
            tags.write("RGZrg1\0".getBytes(StandardCharsets.US_ASCII));
            final byte[] t = tags.toByteArray();
            final byte[] nameBytes = (name + "\0").getBytes(StandardCharsets.US_ASCII);
            final int lSeq = seq.length();
            final int blockSize = 32 + nameBytes.length + 4 + (lSeq + 1) / 2 + lSeq + t.length;
            final ByteBuffer r = le(4 + blockSize);
            r.putInt(blockSize).putInt(0).putInt(start - 1).put((byte) nameBytes.length).put((byte) 60)
                .putShort((short) htsjdk.samtools.GenomicIndexUtil.regionToBin(start - 1, start - 1 + 30))
                .putShort((short) 1).putShort((short) flag).putInt(lSeq).putInt(0).putInt(mateStart - 1)
                .putInt(first ? 210 : -210);
            r.put(nameBytes).putInt((30 << 4) | 0);
            for (int k = 0; k < lSeq; k += 2) {
                final int hi = "=ACMGRSVTWYHKDBN".indexOf(seq.charAt(k));
                final int lo = k + 1 < lSeq ? "=ACMGRSVTWYHKDBN".indexOf(seq.charAt(k + 1)) : 0;
                r.put((byte) ((hi << 4) | lo));
            }
            r.put(qual).put(t);
            records.add(new Object[]{start, i, java.util.Arrays.copyOf(r.array(), r.position())});
        }
        // Coordinate order, the file's declared one; ties keep their creation order.
        records.sort((a, b) -> (Integer) a[0] != (int) (Integer) b[0] ? Integer.compare((Integer) a[0], (Integer) b[0]) : Integer.compare((Integer) a[1], (Integer) b[1]));
        for (final Object[] record : records) {
            final byte[] bytes = (byte[]) record[2];
            body.write(bytes, 0, bytes.length);
        }
        try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(new FileOutputStream(file), (File) null)) {
            out.write(body.toByteArray());
        }
    }

    static void tag(final ByteArrayOutputStream out, final String name, final char type, final int value) {
        out.write(name.charAt(0));
        out.write(name.charAt(1));
        out.write(type);
        if (type == 'C' || type == 'c') out.write(value);
        else { out.write(value & 0xff); out.write((value >> 8) & 0xff); }
    }

    static ByteBuffer le(final int n) { return ByteBuffer.allocate(n).order(ByteOrder.LITTLE_ENDIAN); }

    static String hex(final byte[] bytes) {
        final StringBuilder s = new StringBuilder(bytes.length * 2);
        for (final byte b : bytes) s.append(String.format("%02x", b));
        return s.toString();
    }

    static void emit(final String kind, final String kase, final String payload) {
        System.out.println(kind + "\t" + kase + "\t" + payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }
}
