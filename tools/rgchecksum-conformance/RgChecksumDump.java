import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class RgChecksumDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("@RG\tID:rg1\tSM:sampleA\tLB:lib1\tPL:illumina\n");
        sam.append("@RG\tID:rg2\tSM:sampleB\tLB:lib2\tPL:illumina\n");
        sam.append("r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");

        final File input = File.createTempFile("rgck-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File out = File.createTempFile("rgck-out-", ".read_group_md5");
        out.deleteOnExit();
        final int rc = new picard.sam.CalculateReadGroupChecksum().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("exited " + rc); System.exit(rc); }

        emit("input", "case", sam.toString());
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
