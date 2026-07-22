/*
 * Oracle dump harness for CleanSam (SAM I/O) conformance in picard-rs.
 *
 * Writes a raw input SAM (so the exact bytes are controlled) over a short 150bp contig, with reads
 * that hang off the end of the reference by various amounts, a read that already carries a trailing
 * soft clip, a reverse-strand overhang, and an unmapped read with a nonzero MAPQ. Runs CleanSam and
 * emits an `input` row and an `output` row. CleanSam adds no @PG and no timestamp, so both SAMs are
 * compared raw.
 *
 *   java -cp picard-fat.jar:. CleanSamDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class CleanSamDump {
    static final String SEQ = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT";  // 36
    static final String QUAL = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"; // 36

    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:150\n");
        sam.append(rec("inside",   0,  "chr1", 10,  "36M"));    // ends 45, fine
        sam.append(rec("withclip", 0,  "chr1", 130, "30M6S"));  // 30M ends 159, overhang 9, trailing 6S
        sam.append(rec("offend",   0,  "chr1", 130, "36M"));    // ends 165, overhang 15
        sam.append(rec("neg",      16, "chr1", 135, "36M"));    // reverse, ends 170, overhang 20
        sam.append("unmapped\t4\t*\t0\t30\t*\t*\t0\t0\tACGT\tIIII\n"); // MAPQ 30 on unmapped

        final File input = File.createTempFile("cleansam-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("cleansam-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.CleanSam().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
        });
        if (rc != 0) { System.err.println("CleanSam exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static String rec(final String name, final int flag, final String rname, final int pos, final String cigar) {
        return name + "\t" + flag + "\t" + rname + "\t" + pos + "\t60\t" + cigar + "\t*\t0\t0\t" + SEQ + "\t" + QUAL + "\n";
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
