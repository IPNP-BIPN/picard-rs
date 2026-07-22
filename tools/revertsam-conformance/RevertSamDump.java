/*
 * Oracle dump harness for RevertSam (SAM I/O, default option path) conformance in picard-rs.
 *
 * Writes a coordinate-sorted SAM whose queryname order differs from its coordinate order (zeb@100,
 * amy@200, mid@300), covering a duplicate read with an OQ to restore and NM/MD/AS to clear, a
 * negative-strand read (reverse-complement), and a proper pair with an MC tag (mate info cleared),
 * then runs RevertSam with default options and emits `input` and `output` rows. RevertSam writes a
 * bare header (@HD + @RG only) with no @PG and no timestamp, so both SAMs are compared raw.
 *
 *   java -cp picard-fat.jar:. RevertSamDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class RevertSamDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("@RG\tID:rg1\tSM:s\n");
        // duplicate(0x400), OQ to restore, NM/MD/AS to clear.
        sam.append("zeb\t1024\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tOQ:Z:5555\tNM:i:1\tMD:Z:4\tAS:i:40\tRG:Z:rg1\n");
        // negative strand(0x10): reverse-complement of bases + reversal of quals.
        sam.append("amy\t16\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tABCD\tRG:Z:rg1\n");
        // proper pair (0x1|0x2|0x20|0x40 = 99), mate mapped, MC tag to clear.
        sam.append("mid\t99\tchr1\t300\t60\t4M\t=\t350\t54\tACGT\tIIII\tMC:Z:4M\tRG:Z:rg1\n");

        final File input = File.createTempFile("revert-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("revert-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.RevertSam().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("RevertSam exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
