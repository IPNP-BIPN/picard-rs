/*
 * Oracle dump harness for RevertOriginalBaseQualitiesAndAddMateCigar (SAM I/O, default path)
 * conformance in picard-rs.
 *
 * Writes a coordinate-sorted proper FR pair, both mapped on-reference, each carrying an OQ to restore
 * and no MC, runs the tool with default options, and emits `input` and `output` rows. The tool clones
 * the header, sets only its sort order, and adds no @PG and no timestamp, so both SAMs are compared
 * raw.
 *
 *   java -cp picard-fat.jar:. RevertOrigDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class RevertOrigDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        // proper pair, forward first (0x1|0x2|0x20|0x40 = 99), each with an OQ, no MC yet.
        sam.append("p1\t99\tchr1\t100\t60\t4M\t=\t300\t204\tACGT\tIIII\tOQ:Z:5555\n");
        // reverse second (0x1|0x2|0x10|0x80 = 147).
        sam.append("p1\t147\tchr1\t300\t50\t4M\t=\t100\t-204\tACGT\tJJJJ\tOQ:Z:AAAA\n");

        final File input = File.createTempFile("revorig-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("revorig-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("RevertOriginal... exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
