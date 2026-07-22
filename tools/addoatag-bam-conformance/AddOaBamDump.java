/*
 * Oracle dump harness for AddOATag (SAM -> BAM) conformance in picard-rs.
 *
 * Writes a coordinate-sorted SAM (mapped, reverse, and unmapped reads, one with an NM tag and one an
 * existing OA), runs AddOATag to a BAM with USE_JDK_DEFLATER=true (so the BGZF blocks come from zlib,
 * which the port's writer matches), and emits a `sam` row (the input) and a `bam` row (the tagged BAM,
 * hex). AddOATag adds no @PG, so the whole BAM is compared raw.
 *
 *   java -cp picard-fat.jar:. AddOaBamDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class AddOaBamDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("@SQ\tSN:chr2\tLN:1000\n");
        sam.append("m1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:2\n");
        sam.append("m3\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\tOA:Z:old,1,+,4M,10,;\n");
        sam.append("m2\t16\tchr2\t200\t30\t2M2S\t*\t0\t0\tACGT\tIIII\n");
        sam.append("u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");

        final File in = File.createTempFile("oab-", ".sam");
        in.deleteOnExit();
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam.toString()); }

        final File bam = File.createTempFile("oab-", ".bam");
        bam.deleteOnExit();
        final int rc = new picard.sam.AddOATag().instanceMain(new String[] {
                "INPUT=" + in.getAbsolutePath(),
                "OUTPUT=" + bam.getAbsolutePath(),
                "USE_JDK_DEFLATER=true",
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("AddOATag exited " + rc); System.exit(rc); }

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
