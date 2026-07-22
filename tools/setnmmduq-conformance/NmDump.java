/*
 * Oracle dump harness for SetNmMdAndUqTags (SAM I/O, default path) conformance in picard-rs.
 *
 * Writes a 100bp reference FASTA (+ .fai + .dict), a coordinate-sorted SAM with a perfect read, a
 * read with one mismatch (low quality, so UQ is small), and a read with a 2bp deletion, runs
 * SetNmMdAndUqTags, and emits `fasta`, `input`, and `output` rows. The tool recomputes MD/NM/UQ and
 * adds no @PG and no timestamp, so the whole SAM is compared raw.
 *
 *   java -cp picard-fat.jar:. NmDump
 */
import htsjdk.samtools.*;
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.Random;

public class NmDump {
    public static void main(final String[] args) throws Exception {
        final Random rng = new Random(5);
        final StringBuilder rb = new StringBuilder();
        for (int i = 0; i < 100; i++) rb.append("ACGT".charAt(rng.nextInt(4)));
        final String ref = rb.toString();

        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < ref.length(); i += 60) fasta.append(ref, i, Math.min(i + 60, ref.length())).append('\n');

        final File dir = Files.createTempDirectory("nmref").toFile();
        dir.deleteOnExit();
        final File fa = new File(dir, "ref.fasta");
        try (PrintStream ps = new PrintStream(fa)) { ps.print(fasta.toString()); }
        try (PrintStream ps = new PrintStream(new File(dir, "ref.fasta.fai"))) { ps.print("chr1\t100\t6\t60\t61\n"); }
        try (PrintStream ps = new PrintStream(new File(dir, "ref.dict"))) {
            ps.print("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n");
        }

        final char[] mm = ref.substring(19, 29).toCharArray();
        mm[3] = mm[3] == 'A' ? 'C' : 'A';
        final StringBuilder sam = new StringBuilder("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n");
        sam.append("perfect\t0\tchr1\t1\t60\t10M\t*\t0\t0\t").append(ref, 0, 10).append("\tIIIIIIIIII\n");
        sam.append("mm\t0\tchr1\t20\t60\t10M\t*\t0\t0\t").append(new String(mm)).append("\t##########\n");
        sam.append("del\t0\tchr1\t40\t60\t4M2D6M\t*\t0\t0\t")
           .append(ref, 39, 43).append(ref, 45, 51).append("\tIIIIIIIIII\n");

        final File in = new File(dir, "in.sam");
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam.toString()); }
        final File out = new File(dir, "out.sam");
        final int rc = new picard.sam.SetNmMdAndUqTags().instanceMain(new String[] {
                "INPUT=" + in.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + fa.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("SetNmMdAndUqTags exited " + rc); System.exit(rc); }

        emit("fasta", "case", fasta.toString());
        emit("input", "case", sam.toString());
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
