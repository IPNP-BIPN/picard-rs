/*
 * Oracle dump harness for FixMateInformation (single SAM input, default options) conformance in
 * picard-rs.
 *
 * Writes a raw coordinate-sorted SAM whose pairs have absent mate info (RNEXT *, PNEXT 0, TLEN 0,
 * no MC/MQ), runs FixMateInformation, and emits an `input` and an `output` row. FixMateInformation
 * adds no @PG and no timestamp, so both SAMs are compared raw. Read leniently, since the input is by
 * definition inconsistent.
 *
 *   java -cp picard-fat.jar:. FixMateDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class FixMateDump {
    static final String S = "ACGTACGTAC";
    static final String Q = "IIIIIIIIII";

    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("@SQ\tSN:chr2\tLN:1000\n");
        // Two FR pairs on chr1, and a pair split across chr1/chr2, all with absent mate info.
        sam.append(rec("pA", 65,  "chr1", 100, "10M"));
        sam.append(rec("pB", 65,  "chr1", 150, "10M"));
        sam.append(rec("pA", 145, "chr1", 300, "10M"));
        sam.append(rec("pB", 145, "chr1", 350, "10M"));
        sam.append(rec("pC", 65,  "chr1", 500, "10M"));
        sam.append(rec("pC", 129, "chr2", 600, "10M")); // mate on another contig, forward

        final File input = File.createTempFile("fixmate-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("fixmate-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.FixMateInformation().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("FixMateInformation exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static String rec(final String name, final int flag, final String rname, final int pos, final String cigar) {
        return name + "\t" + flag + "\t" + rname + "\t" + pos + "\t60\t" + cigar + "\t*\t0\t0\t" + S + "\t" + Q + "\n";
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
