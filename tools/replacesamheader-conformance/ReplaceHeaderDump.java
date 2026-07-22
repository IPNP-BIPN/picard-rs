/*
 * Oracle dump harness for ReplaceSamHeader (SAM I/O, standardReheader path) conformance in picard-rs.
 *
 * Writes a coordinate-sorted INPUT SAM (with an @RG the records reference) and a HEADER stub SAM that
 * keeps the same @SQ block and sort order but carries a different @RG (new SM/LB) and an added @CO,
 * runs ReplaceSamHeader, and emits `input`, `header`, and `output` rows. standardReheader writes the
 * replacement header unchanged (no @PG, no timestamp) and the records presorted (input order), so both
 * SAMs are compared raw.
 *
 *   java -cp picard-fat.jar:. ReplaceHeaderDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class ReplaceHeaderDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder in = new StringBuilder();
        in.append("@HD\tVN:1.6\tSO:coordinate\n");
        in.append("@SQ\tSN:chr1\tLN:1000\n");
        in.append("@SQ\tSN:chr2\tLN:1000\n");
        in.append("@RG\tID:old\tSM:sampleA\n");
        in.append("r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n");
        in.append("r2\t0\tchr2\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n");

        final StringBuilder hdr = new StringBuilder();
        hdr.append("@HD\tVN:1.6\tSO:coordinate\n");
        hdr.append("@SQ\tSN:chr1\tLN:1000\n");
        hdr.append("@SQ\tSN:chr2\tLN:1000\n");
        hdr.append("@RG\tID:old\tSM:sampleB\tLB:lib1\n");
        hdr.append("@CO\tedited by hand\n");

        final File input = File.createTempFile("reheader-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(in.toString()); }

        final File header = File.createTempFile("reheader-hdr-", ".sam");
        header.deleteOnExit();
        try (PrintStream ps = new PrintStream(header)) { ps.print(hdr.toString()); }

        final File out = File.createTempFile("reheader-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.ReplaceSamHeader().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "HEADER=" + header.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("ReplaceSamHeader exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("header", "case", hdr.toString());
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
