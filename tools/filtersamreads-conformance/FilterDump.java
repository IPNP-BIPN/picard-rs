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

        // A second input with an RG tag on some reads (and one without), for the tag-value filters.
        final StringBuilder tagged = new StringBuilder();
        tagged.append("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:rg1\tSM:s\n@RG\tID:rg2\tSM:t\n");
        tagged.append("a\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        tagged.append("b\t0\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n");
        tagged.append("c\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
        tagged.append("d\t0\tchr1\t400\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        final File tagInput = File.createTempFile("fs-tag-", ".sam");
        tagInput.deleteOnExit();
        try (PrintStream ps = new PrintStream(tagInput)) { ps.print(tagged.toString()); }

        emit("tag_input", "case", tagged.toString());
        emit("include_tag", "case", runTag(tagInput, "includeTagValues"));
        emit("exclude_tag", "case", runTag(tagInput, "excludeTagValues"));

        // A queryname-sorted input for the aligned filters: a both-aligned pair, a one-unmapped pair,
        // a both-unmapped pair, an aligned singleton, and an unmapped singleton.
        final StringBuilder pairs = new StringBuilder();
        pairs.append("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100000\n");
        pairs.append("pairAA\t99\tchr1\t100\t60\t4M\t=\t200\t104\tACGT\tIIII\n");
        pairs.append("pairAA\t147\tchr1\t200\t60\t4M\t=\t100\t-104\tACGT\tIIII\n");
        pairs.append("pairAU\t97\tchr1\t300\t60\t4M\t=\t300\t0\tACGT\tIIII\n");
        pairs.append("pairAU\t141\t*\t0\t0\t*\t=\t300\t0\tACGT\tIIII\n");
        pairs.append("pairUU\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
        pairs.append("pairUU\t141\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
        pairs.append("singleA\t0\tchr1\t500\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        pairs.append("singleU\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
        final File pairsInput = File.createTempFile("fs-pairs-", ".sam");
        pairsInput.deleteOnExit();
        try (PrintStream ps = new PrintStream(pairsInput)) { ps.print(pairs.toString()); }

        emit("pairs_input", "case", pairs.toString());
        emit("include_aligned", "case", runPlain(pairsInput, "includeAligned"));
        emit("exclude_aligned", "case", runPlain(pairsInput, "excludeAligned"));
    }

    static String runPlain(final File input, final String filter) throws Exception {
        final File out = File.createTempFile("fs-aln-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.FilterSamReads().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "FILTER=" + filter,
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("FilterSamReads " + filter + " exited " + rc); System.exit(rc); }
        return new String(Files.readAllBytes(out.toPath()));
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

    static String runTag(final File input, final String filter) throws Exception {
        final File out = File.createTempFile("fs-tag-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.FilterSamReads().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "FILTER=" + filter,
                "TAG=RG",
                "TAG_VALUE=rg1",
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
