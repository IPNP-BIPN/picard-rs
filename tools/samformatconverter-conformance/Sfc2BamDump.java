/*
 * Oracle dump harness for SamFormatConverter (SAM -> BAM) conformance in picard-rs.
 *
 * Writes a small coordinate-sorted SAM (mixed flags/cigars/tags plus an unmapped read), converts it
 * to BAM with SamFormatConverter and USE_JDK_DEFLATER=true (so the BGZF blocks come from java.util.zip
 * / zlib, which the port's writer matches; Picard's default GKL/igzip would not), and emits a `sam`
 * row (the input) and a `bam` row (the BAM bytes, hex). SamFormatConverter adds no @PG, so the whole
 * BAM is compared raw.
 *
 *   java -cp picard-fat.jar:. Sfc2BamDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class Sfc2BamDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:100000\n");
        sam.append("@RG\tID:rg1\tSM:s\n");
        sam.append("r0\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        sam.append("r1\t147\tchr1\t200\t60\t2M2D2M\t=\t100\t-104\tACGT\tIIII\tRG:Z:rg1\tMQ:i:42\n");
        sam.append("r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");

        final File in = File.createTempFile("s2b-", ".sam");
        in.deleteOnExit();
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam.toString()); }

        final File bam = File.createTempFile("s2b-", ".bam");
        bam.deleteOnExit();
        final int rc = new picard.sam.SamFormatConverter().instanceMain(new String[] {
                "INPUT=" + in.getAbsolutePath(),
                "OUTPUT=" + bam.getAbsolutePath(),
                "USE_JDK_DEFLATER=true",
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("SamFormatConverter exited " + rc); System.exit(rc); }

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
