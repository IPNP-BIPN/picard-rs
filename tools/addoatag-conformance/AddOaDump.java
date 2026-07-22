/*
 * Oracle dump harness for AddOATag (SAM I/O, default path) conformance in picard-rs.
 *
 * Writes a raw SAM (controlled bytes) with mapped/reverse/unmapped reads, some carrying an NM tag
 * and one an existing OA tag, runs AddOATag, and emits an `input` and an `output` row. AddOATag
 * adds no @PG and no timestamp, so both SAMs are compared raw.
 *
 *   java -cp picard-fat.jar:. AddOaDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class AddOaDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("@SQ\tSN:chr2\tLN:1000\n");
        sam.append("m1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:2\n");
        sam.append("m3\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\tOA:Z:old,1,+,4M,10,;\n"); // existing OA
        sam.append("m2\t16\tchr2\t200\t30\t2M2S\t*\t0\t0\tACGT\tIIII\n");       // reverse, no NM
        sam.append("u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");                  // unmapped

        final File input = File.createTempFile("addoa-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("addoa-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.AddOATag().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("AddOATag exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
