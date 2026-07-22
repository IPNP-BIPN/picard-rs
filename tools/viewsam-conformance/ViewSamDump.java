import java.io.File;
import java.io.PrintStream;
import java.io.ByteArrayOutputStream;
import java.nio.file.Files;

public class ViewSamDump {
    static String input;

    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("m1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n");           // mapped, PF
        sam.append("qc\t512\tchr1\t150\t60\t4M\t*\t0\t0\tACGT\tIIII\n");         // mapped, vendor-fail (0x200)
        sam.append("m2\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\n");           // mapped, PF
        sam.append("u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");                  // unmapped, PF
        input = sam.toString();

        emit("input", "base", input);
        emit("default", "case", run("ALIGNMENT_STATUS=All", "PF_STATUS=All"));
        emit("aligned", "case", run("ALIGNMENT_STATUS=Aligned"));
        emit("unaligned", "case", run("ALIGNMENT_STATUS=Unaligned"));
        emit("pf", "case", run("PF_STATUS=PF"));
        emit("nonpf", "case", run("PF_STATUS=NonPF"));
    }

    static String run(final String... opts) throws Exception {
        final File in = File.createTempFile("viewsam-in-", ".sam");
        in.deleteOnExit();
        try (PrintStream ps = new PrintStream(in)) { ps.print(input); }

        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        final PrintStream captured = new PrintStream(bytes);
        final PrintStream saved = System.out;
        final String[] argv = new String[opts.length + 2];
        argv[0] = "INPUT=" + in.getAbsolutePath();
        argv[1] = "VALIDATION_STRINGENCY=SILENT";
        System.arraycopy(opts, 0, argv, 2, opts.length);
        int rc;
        try {
            System.setOut(captured);
            rc = new picard.sam.ViewSam().instanceMain(argv);
        } finally {
            System.setOut(saved);
        }
        if (rc != 0) { System.err.println("ViewSam exited " + rc); System.exit(rc); }
        return bytes.toString();
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
