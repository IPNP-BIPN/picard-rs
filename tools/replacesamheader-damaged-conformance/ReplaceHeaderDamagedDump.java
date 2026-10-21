/*
 * Oracle dump for ReplaceSamHeader's block-copy path on damaged BAMs, against Picard 3.5.0.
 *
 * HEADER is read by SamReaderFactory.getFileHeader, which inflates no further than the header: a
 * HEADER whose records or terminator are damaged is accepted, and only a header that does not
 * inflate throws. INPUT goes to BamFileIoUtils.reheaderBamFile, which refuses a file whose last
 * block is not a valid one, reads the header to find the first record, and copies the record
 * blocks without inflating them: damaged records reach the output as they are.
 *
 * Each case emits the two files it ran on (`header_bam`/`header_sam`, `input_bam`, as hex), then
 * `exit` or `thrown` (Throwable.toString(), the work directory written as <dir>), then the output
 * as hex or `output_absent`.
 *
 *   java -cp picard.jar:. ReplaceHeaderDamagedDump
 */
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;

public class ReplaceHeaderDamagedDump {
    static File W;

    static final String SAM =
            "@HD\tVN:1.6\tSO:coordinate\n"
            + "@SQ\tSN:chr1\tLN:1000\n"
            + "@RG\tID:rg1\tSM:s\n"
            + "r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            + "r2\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    static final String HEADER =
            "@HD\tVN:1.6\tSO:coordinate\n"
            + "@SQ\tSN:chr1\tLN:1000\n"
            + "@RG\tID:rg1\tSM:replaced\tLB:lib1\n"
            + "@CO\tedited by hand\n";
    /** Not a BGZF block: htsjdk's reader would refuse it, the block copy does not look. */
    static final byte[] JUNK = "\u001f\u008bJUNKJUNKJUNKJUNKJUNKJUNKJUNK".getBytes(java.nio.charset.StandardCharsets.ISO_8859_1);
    static final int EOF = 28;

    public static void main(final String[] args) throws Exception {
        W = Files.createTempDirectory("rshdamaged").toFile();
        final byte[] bam = bam(SAM);
        final int first = blockSize(bam, 0);
        final byte[] headerOnly = Arrays.copyOf(bam, first);
        final byte[] records = Arrays.copyOfRange(bam, first, bam.length - EOF);
        final byte[] eof = Arrays.copyOfRange(bam, bam.length - EOF, bam.length);
        final byte[] header = HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8);

        // HEADER damaged, INPUT clean.
        run("header_clean", bam, bam);
        run("header_truncated_in_header", Arrays.copyOf(bam, first / 2), bam);
        run("header_without_terminator", headerOnly, bam);
        run("header_records_cut", Arrays.copyOf(bam, first + records.length / 2), bam);
        run("header_records_junk", concat(headerOnly, JUNK, eof), bam);
        // INPUT damaged, HEADER clean.
        runSam("input_clean", header, bam);
        runSam("input_without_terminator", header, Arrays.copyOf(bam, bam.length - EOF));
        runSam("input_records_cut", header, Arrays.copyOf(bam, first + records.length / 2));
        runSam("input_records_junk", header, concat(headerOnly, JUNK, eof));
        runSam("input_records_junk_unterminated", header, concat(headerOnly, records, JUNK));
    }

    static byte[] bam(final String sam) throws Exception {
        final File in = new File(W, "build.sam");
        Files.write(in.toPath(), sam.getBytes(java.nio.charset.StandardCharsets.UTF_8));
        final File out = new File(W, "build.bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(in);
                SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, out)) {
            for (final SAMRecord rec : reader) writer.addAlignment(rec);
        }
        final byte[] bytes = Files.readAllBytes(out.toPath());
        in.delete();
        out.delete();
        return bytes;
    }

    /** BSIZE + 1, from the BC subfield of the block at `at`. */
    static int blockSize(final byte[] b, final int at) {
        return ((b[at + 16] & 0xff) | (b[at + 17] & 0xff) << 8) + 1;
    }

    static byte[] concat(final byte[]... parts) {
        int n = 0;
        for (final byte[] p : parts) n += p.length;
        final byte[] out = new byte[n];
        int at = 0;
        for (final byte[] p : parts) {
            System.arraycopy(p, 0, out, at, p.length);
            at += p.length;
        }
        return out;
    }

    static void run(final String kase, final byte[] header, final byte[] input) throws Exception {
        emit("header_bam", kase, hex(header));
        go(kase, new File(W, kase + ".header.bam"), header, input);
    }

    static void runSam(final String kase, final byte[] header, final byte[] input) throws Exception {
        emit("header_sam", kase, hex(header));
        go(kase, new File(W, kase + ".header.sam"), header, input);
    }

    static void go(final String kase, final File headerFile, final byte[] header, final byte[] input) throws Exception {
        emit("input_bam", kase, hex(input));
        Files.write(headerFile.toPath(), header);
        final File in = new File(W, kase + ".in.bam");
        Files.write(in.toPath(), input);
        final File out = new File(W, kase + ".out.bam");
        try {
            final int rc = new picard.sam.ReplaceSamHeader().instanceMain(new String[] {
                    "INPUT=" + in.getAbsolutePath(),
                    "HEADER=" + headerFile.getAbsolutePath(),
                    "OUTPUT=" + out.getAbsolutePath()});
            emit("exit", kase, Integer.toString(rc));
        } catch (final Throwable t) {
            emit("thrown", kase, t.toString().replace(W.getAbsolutePath(), "<dir>"));
        }
        if (out.exists()) emit("output", kase, hex(Files.readAllBytes(out.toPath())));
        else emit("output_absent", kase, "");
    }

    static String hex(final byte[] b) {
        final StringBuilder s = new StringBuilder(b.length * 2);
        for (final byte x : b) s.append(String.format("%02x", x & 0xff));
        return s.toString();
    }

    static void emit(final String kind, final String kase, final String payload) {
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
