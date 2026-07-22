/*
 * Oracle dump harness for SortSam (SAM -> BAM) conformance in picard-rs.
 *
 * Writes an unsorted SAM, sorts it to a BAM with SortSam SORT_ORDER=coordinate and
 * USE_JDK_DEFLATER=true (so the BGZF blocks come from zlib, which the port's writer matches), and
 * emits a `sam` row (the input) and a `bam` row (the coordinate-sorted BAM, hex). SortSam adds no @PG,
 * so the whole BAM is compared raw.
 *
 *   java -cp picard-fat.jar:. SortSamBamDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class SortSamBamDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:unsorted\n");
        sam.append("@SQ\tSN:chr1\tLN:100000\n");
        sam.append("@SQ\tSN:chr2\tLN:100000\n");
        sam.append("@RG\tID:rg1\tSM:s\n");
        // Deliberately out of coordinate order across two contigs.
        sam.append("r3\t0\tchr2\t500\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        sam.append("r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        sam.append("r4\t16\tchr1\t100\t30\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        sam.append("r2\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        sam.append("u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");

        final File in = File.createTempFile("ss-", ".sam");
        in.deleteOnExit();
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam.toString()); }

        final File bam = File.createTempFile("ss-", ".bam");
        bam.deleteOnExit();
        final int rc = new picard.sam.SortSam().instanceMain(new String[] {
                "INPUT=" + in.getAbsolutePath(),
                "OUTPUT=" + bam.getAbsolutePath(),
                "SORT_ORDER=coordinate",
                "USE_JDK_DEFLATER=true",
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("SortSam exited " + rc); System.exit(rc); }

        emit("sam", "case", esc(sam.toString()));
        emit("bam", "case", hex(Files.readAllBytes(bam.toPath())));
    }

    static String hex(final byte[] b) {
        final StringBuilder sb = new StringBuilder(b.length * 2);
        for (final byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static String esc(final String p) {
        return p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static void emit(final String kind, final String kase, final String payload) {
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
