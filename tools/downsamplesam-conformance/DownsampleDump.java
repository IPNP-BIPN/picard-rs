/*
 * Oracle dump harness for DownsampleSam (SAM I/O, ConstantMemory strategy) conformance in picard-rs.
 *
 * Writes a coordinate-sorted SAM of 12 distinctly-named reads, runs DownsampleSam with PROBABILITY=0.5
 * and RANDOM_SEED=1 (the default seed) to SAM output, and emits an `input` row and an `output` row.
 *
 * DownsampleSam adds a @PG provenance record whose CL: is the command line (temp paths, every option),
 * which is canonicalized away in comparison; the conformance strips @PG lines from both sides and
 * compares the surviving records plus the rest of the header raw.
 *
 *   java -cp picard-fat.jar:. DownsampleDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class DownsampleDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:100000\n");
        sam.append("@RG\tID:rg1\tSM:s\n");
        int pos = 100;
        for (int i = 0; i < 12; i++) {
            sam.append("read").append(i).append("\t0\tchr1\t").append(pos)
               .append("\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
            pos += 100;
        }

        final File input = File.createTempFile("ds-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("ds-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.DownsampleSam().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "PROBABILITY=0.5",
                "RANDOM_SEED=1",
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("DownsampleSam exited " + rc); System.exit(rc); }

        emit("input", "case", sam.toString());
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
