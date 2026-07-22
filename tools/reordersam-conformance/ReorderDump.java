/*
 * Oracle dump harness for ReorderSam (SAM I/O, unindexed path) conformance in picard-rs.
 *
 * Writes a coordinate-sorted SAM over a read dictionary chr1, chr2, chr3, plus a SEQUENCE_DICTIONARY
 * (.dict) that swaps chr1/chr2 and drops chr3, then runs ReorderSam with
 * ALLOW_INCOMPLETE_DICT_CONCORDANCE so the chr3 read and the chr3-mate become unmapped. It emits an
 * `input` row (the SAM), a `dict` row (the .dict), and an `output` row (the reordered SAM). ReorderSam
 * clones the header and adds no @PG and no timestamp, and the presorted=false writer re-sorts by the
 * coordinate order, so both SAMs are compared raw.
 *
 *   java -cp picard-fat.jar:. ReorderDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class ReorderDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:1000\n");
        sam.append("@SQ\tSN:chr2\tLN:1000\n");
        sam.append("@SQ\tSN:chr3\tLN:1000\n");
        // r1: paired, mate on chr3 (to be dropped) with an MC tag that must be removed.
        sam.append("r1\t1\tchr1\t100\t60\t4M\tchr3\t500\t0\tACGT\tIIII\tMC:Z:4M\n");
        sam.append("r2\t0\tchr2\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        // r3: on chr3, which the new reference drops => becomes unmapped.
        sam.append("r3\t0\tchr3\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        sam.append("u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");

        final String dict = "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n";

        final File input = File.createTempFile("reorder-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }

        final File dictFile = File.createTempFile("reorder-dict-", ".dict");
        dictFile.deleteOnExit();
        try (PrintStream ps = new PrintStream(dictFile)) { ps.print(dict); }

        final File out = File.createTempFile("reorder-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.ReorderSam().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "SEQUENCE_DICTIONARY=" + dictFile.getAbsolutePath(),
                "ALLOW_INCOMPLETE_DICT_CONCORDANCE=true",
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("ReorderSam exited " + rc); System.exit(rc); }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("dict", "case", dict);
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
