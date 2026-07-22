/*
 * Oracle dump harness for FilterSamReads (SAM I/O, read-list filters) conformance in picard-rs.
 *
 * Writes a coordinate-sorted SAM and a READ_LIST_FILE with a subset of the read names, runs
 * FilterSamReads with FILTER=includeReadList and again with FILTER=excludeReadList, and emits an
 * `input` row, a `list` row, an `include` row, and an `exclude` row. FilterSamReads keeps the input
 * sort order (no re-sort) and adds no @PG and no timestamp, so every SAM is compared raw.
 *
 *   java -cp picard-fat.jar:. FilterDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class FilterDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n");
        final String[] names = {"keepA", "dropB", "keepC", "dropD", "keepE"};
        int pos = 100;
        for (final String n : names) {
            sam.append(n).append("\t0\tchr1\t").append(pos).append("\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
            pos += 100;
        }
        final String list = "keepA\nkeepC\nkeepE\n";

        final File input = File.createTempFile("fs-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam.toString()); }
        final File listFile = File.createTempFile("fs-list-", ".txt");
        listFile.deleteOnExit();
        try (PrintStream ps = new PrintStream(listFile)) { ps.print(list); }

        emit("input", "case", sam.toString());
        emit("list", "case", list);
        emit("include", "case", run(input, listFile, "includeReadList"));
        emit("exclude", "case", run(input, listFile, "excludeReadList"));
    }

    static String run(final File input, final File listFile, final String filter) throws Exception {
        final File out = File.createTempFile("fs-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.FilterSamReads().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "READ_LIST_FILE=" + listFile.getAbsolutePath(),
                "FILTER=" + filter,
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("FilterSamReads " + filter + " exited " + rc); System.exit(rc); }
        return new String(Files.readAllBytes(out.toPath()));
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
