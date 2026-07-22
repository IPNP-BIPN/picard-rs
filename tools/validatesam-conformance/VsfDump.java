/*
 * Oracle corpus generator for ValidateSamFile (VERBOSE mode, SAM input, no reference) conformance.
 *
 * Emits, for a set of hand-built SAM inputs exercising the header and per-record checks that need
 * no reference and no cross-record (mate/sort) state, the exact stdout of `ValidateSamFile
 * MODE=VERBOSE`. The verbose output is raw error lines (SAMValidationError.toString) with no
 * timestamp and no banner, plus "No errors found" when clean, so the whole payload is compared raw.
 *
 *   java -cp picard-fat.jar:. VsfDump > corpus.txt
 */
import java.io.*;
import java.nio.file.*;

public class VsfDump {
    static StringBuilder out = new StringBuilder();

    static void emit(String kind, String kase, String payload) {
        String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        out.append(kind).append('\t').append(kase).append('\t').append(esc).append('\n');
    }

    static void run(String kase, String sam) throws Exception {
        File dir = Files.createTempDirectory("vsf").toFile();
        dir.deleteOnExit();
        File in = new File(dir, "in.sam");
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam); }
        File res = new File(dir, "out.txt");
        new picard.sam.ValidateSamFile().instanceMain(new String[]{
            "INPUT=" + in.getAbsolutePath(),
            "OUTPUT=" + res.getAbsolutePath(),
            "MODE=VERBOSE",
        });
        emit("input", kase, sam);
        emit("output", kase, new String(Files.readAllBytes(res.toPath())));
    }

    public static void main(String[] a) throws Exception {
        String hd = "@HD\tVN:1.6\tSO:coordinate\n";
        String sq = "@SQ\tSN:chr1\tLN:100\n";
        String rg = "@RG\tID:rg1\tSM:s\tPL:illumina\n";

        // 1. Clean file: unpaired, coordinate-sorted, RG present, NM present.
        run("clean", hd + sq + rg
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n"
            + "b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // 2. Mapped read missing NM tag (no reference -> presence check only).
        run("missing_nm", hd + sq + rg
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");

        // 3. Empty sequence dictionary (no @SQ) with only an unmapped read: no error is emitted,
        //    because the "sequence dictionary is empty" warning fires only on the first mapped read.
        run("dict_empty_unmapped", "@HD\tVN:1.6\tSO:coordinate\n" + rg
            + "a\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");

        // 4. Header with no read groups + a record with no RG tag.
        run("no_rg", hd + sq
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n");

        // 5. QUAL field is "*" (unspecified quality scores).
        run("qual_not_stored", hd + sq + rg
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\t*\tRG:Z:rg1\tNM:i:0\n");

        // 6. Read group missing its PL (platform) value.
        run("missing_pl", hd + sq + "@RG\tID:rg1\tSM:s\n"
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // 7. Read group with an invalid PL value.
        run("invalid_pl", hd + sq + "@RG\tID:rg1\tSM:s\tPL:bogusplatform\n"
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // 8. Invalid version number in header.
        run("bad_version", "@HD\tVN:9.9\tSO:coordinate\n" + sq + rg
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // 9. Several problems in one file, to pin interleaving/order.
        run("mixed", hd + sq
            + "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            + "b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\n");

        System.out.print(out);
    }
}
