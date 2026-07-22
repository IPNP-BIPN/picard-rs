/*
 * Oracle dump harness for CreateSequenceDictionary conformance in picard-rs.
 *
 * Writes a small reference FASTA (with a soft-masked lower-case run and an N run), runs
 * CreateSequenceDictionary, and emits a `fasta` row (the reference bytes) and a `dict` row (the .dict).
 * The .dict's UR field is the reference's file: URI, which is path-dependent and is stripped in
 * comparison; @HD and each SN/LN/M5 are compared raw.
 *
 *   java -cp picard-fat.jar:. DictDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class DictDump {
    public static void main(final String[] args) throws Exception {
        final String fasta = ">chr1 some description\nACGTacgtNNNN\nACGTACGT\n>chr2\nTTTTGGGGCCCCAAAA\n";

        final File ref = File.createTempFile("dict-ref-", ".fasta");
        ref.deleteOnExit();
        try (PrintStream ps = new PrintStream(ref)) { ps.print(fasta); }

        final File out = new File(ref.getAbsolutePath().replace(".fasta", ".dict"));
        out.deleteOnExit();
        if (out.exists()) out.delete();
        final int rc = new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "REFERENCE=" + ref.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
        });
        if (rc != 0) { System.err.println("CreateSequenceDictionary exited " + rc); System.exit(rc); }

        emit("fasta", "case", fasta);
        emit("dict", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
