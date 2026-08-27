/*
 * CompareMetrics' verdict, taken from the reference.
 *
 * Two metrics files in, one verdict out, plus a report of the differences and a table of them.
 * What is measured is what counts as a difference, what each of the four tolerance arguments
 * forgives, and what the tool writes when it finds one.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE EXIT CODE IS THE VERDICT: zero when the two files agree and one when they do not, and
 *     the report says "equal" or "NOT equal" in the same words;
 *   - --METRICS_TO_IGNORE DROPS A COLUMN FROM THE COMPARISON entirely, so a file that differs
 *     only there compares equal;
 *   - --METRICS_NOT_REQUIRED DOES THE SAME THING TO A COLUMN THAT IS PRESENT IN BOTH: it is not
 *     only about a column one file lacks, and a run that names a differing column under it
 *     compares equal exactly as --METRICS_TO_IGNORE would;
 *   - --METRIC_ALLOWABLE_RELATIVE_CHANGE TAKES A COLON-SEPARATED PAIR and forgives a numeric
 *     difference within that relative change;
 *   - THE CHANGE IS RELATIVE TO THE FIRST FILE'S VALUE: 0.1 against 0.11 is a change of 0.1 and
 *     0.11 against 0.1 is a change of about 0.0909, so a tolerance of 0.095 forgives the second
 *     ordering and not the first;
 *   - --IGNORE_HISTOGRAM_DIFFERENCES FORGIVES THE HISTOGRAM AND NOT THE TABLE;
 *   - --KEY MATCHES ROWS BY A COLUMN rather than by position, so two files whose rows are in
 *     different orders compare equal under it and not without it;
 *   - A ROW ONE FILE HAS AND THE OTHER DOES NOT IS A DIFFERENCE, keyed or not;
 *   - TWO FILES OF DIFFERENT METRIC CLASSES ARE AN EXCEPTION rather than a difference: the
 *     reader is asked for a field the second class does not have and throws;
 *   - A FILE THAT IS NOT A METRICS FILE AT ALL THROWS TOO, from the header parser;
 *   - AND THE REPORT AND THE TABLE ARE BOTH WRITTEN WHATEVER THE VERDICT, the table being empty
 *     when the two files agree.
 *
 * Output:
 *
 *     input\t<name>\t<that metrics file, escaped>
 *     verdict\t<case>\t<exit code>
 *     report\t<case>\t<the OUTPUT file, escaped>
 *     table\t<case>\t<the OUTPUT_TABLE file without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CompareMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /**
     * One metrics file, written by hand.
     *
     * The header comments carry a timestamp in a real file; here they are fixed, which is what
     * makes the fixture reproducible.
     */
    static String metricsFile(final String metricClass, final List<String> columns,
                              final List<List<String>> rows, final List<String> histogram) {
        final List<String> lines = new ArrayList<>();
        lines.add("## htsjdk.samtools.metrics.StringHeader");
        lines.add("# fixture");
        lines.add("");
        lines.add("## METRICS CLASS\t" + metricClass);
        lines.add(String.join("\t", columns));
        for (final List<String> row : rows) {
            lines.add(String.join("\t", row));
        }
        lines.add("");
        if (histogram != null) {
            lines.add("## HISTOGRAM\tjava.lang.Integer");
            lines.addAll(histogram);
            lines.add("");
        }
        return String.join("\n", lines) + "\n";
    }

    static final String CLASS = "picard.sam.DuplicationMetrics";
    static final List<String> COLUMNS = List.of("LIBRARY", "UNPAIRED_READS_EXAMINED",
            "READ_PAIRS_EXAMINED", "PERCENT_DUPLICATION");

    static String base() {
        return metricsFile(CLASS, COLUMNS,
                List.of(List.of("libA", "100", "200", "0.1"),
                        List.of("libB", "300", "400", "0.2")),
                List.of("bin\tvalue", "1\t10", "2\t20"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final String name, final String left, final String right,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("cmp");
        final Path a = write(dir, "a.metrics", left);
        final Path b = write(dir, "b.metrics", right);
        final File report = new File(dir.toFile(), "report.txt");
        final File table = new File(dir.toFile(), "table.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "INPUT=" + a, "INPUT=" + b,
                "OUTPUT=" + report.getAbsolutePath(),
                "OUTPUT_TABLE=" + table.getAbsolutePath()));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.CompareMetrics()
                    .instanceMain(argv.toArray(new String[0]));
            emit("verdict", name, Integer.toString(code));
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":" + cause.getMessage());
            return;
        }
        if (report.exists()) {
            // The report names the two input files by their absolute paths, which differ from
            // run to run, so the temporary directory is masked out.
            emit("report", name,
                    Files.readString(report.toPath()).replace(dir.toString(), "<dir>"));
        } else {
            emit("report", name, "");
        }
        if (table.exists()) {
            final StringBuilder body = new StringBuilder();
            for (final String line : Files.readString(table.toPath()).split("\n", -1)) {
                if (!line.startsWith("#") && !line.isEmpty()) {
                    body.append(line).append('\n');
                }
            }
            emit("table", name, body.toString());
        } else {
            emit("table", name, "");
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final String base = base();
        emit("input", "base", base);

        // The same file twice.
        run("identical", base, base);

        // One numeric column differing by a tenth.
        final String changed = metricsFile(CLASS, COLUMNS,
                List.of(List.of("libA", "100", "200", "0.11"),
                        List.of("libB", "300", "400", "0.2")),
                List.of("bin\tvalue", "1\t10", "2\t20"));
        emit("input", "changed", changed);
        run("one-column-differs", base, changed);
        run("ignored-column", base, changed, "METRICS_TO_IGNORE=PERCENT_DUPLICATION");
        // A tenth of 0.1 is 0.1 relative; a tenth of 0.11 is about 0.0909.
        run("relative-change-generous", base, changed,
                "METRIC_ALLOWABLE_RELATIVE_CHANGE=PERCENT_DUPLICATION:0.2");
        run("relative-change-tight", base, changed,
                "METRIC_ALLOWABLE_RELATIVE_CHANGE=PERCENT_DUPLICATION:0.05");
        // The two files the other way round, which answers which value the change is relative to.
        run("relative-change-reversed", changed, base,
                "METRIC_ALLOWABLE_RELATIVE_CHANGE=PERCENT_DUPLICATION:0.095");
        run("relative-change-forward", base, changed,
                "METRIC_ALLOWABLE_RELATIVE_CHANGE=PERCENT_DUPLICATION:0.095");

        // A column one file does not have.
        final String fewerColumns = metricsFile(CLASS,
                List.of("LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED"),
                List.of(List.of("libA", "100", "200"), List.of("libB", "300", "400")),
                List.of("bin\tvalue", "1\t10", "2\t20"));
        emit("input", "fewer-columns", fewerColumns);
        run("missing-column", base, fewerColumns);
        run("missing-column-not-required", base, fewerColumns,
                "METRICS_NOT_REQUIRED=PERCENT_DUPLICATION");
        // Not-required turns out to be the SAME as ignored for a column present in both: naming
        // a differing column under it compares equal.
        run("not-required-but-present", base, changed,
                "METRICS_NOT_REQUIRED=PERCENT_DUPLICATION");

        // The histogram alone differing.
        final String otherHistogram = metricsFile(CLASS, COLUMNS,
                List.of(List.of("libA", "100", "200", "0.1"),
                        List.of("libB", "300", "400", "0.2")),
                List.of("bin\tvalue", "1\t10", "2\t99"));
        emit("input", "other-histogram", otherHistogram);
        run("histogram-differs", base, otherHistogram);
        run("histogram-ignored", base, otherHistogram, "IGNORE_HISTOGRAM_DIFFERENCES=true");

        // The rows in a different order.
        final String reordered = metricsFile(CLASS, COLUMNS,
                List.of(List.of("libB", "300", "400", "0.2"),
                        List.of("libA", "100", "200", "0.1")),
                List.of("bin\tvalue", "1\t10", "2\t20"));
        emit("input", "reordered", reordered);
        run("rows-reordered", base, reordered);
        run("rows-reordered-keyed", base, reordered, "KEY=LIBRARY");

        // A row one file does not have.
        final String oneRow = metricsFile(CLASS, COLUMNS,
                List.of(List.of("libA", "100", "200", "0.1")),
                List.of("bin\tvalue", "1\t10", "2\t20"));
        emit("input", "one-row", oneRow);
        run("missing-row", base, oneRow);
        run("missing-row-keyed", base, oneRow, "KEY=LIBRARY");

        // Two different metric classes.
        final String otherClass = metricsFile("picard.analysis.AlignmentSummaryMetrics", COLUMNS,
                List.of(List.of("libA", "100", "200", "0.1"),
                        List.of("libB", "300", "400", "0.2")),
                List.of("bin\tvalue", "1\t10", "2\t20"));
        emit("input", "other-class", otherClass);
        run("different-classes", base, otherClass);

        // A file that is not a metrics file at all.
        run("not-a-metrics-file", base, "this is not a metrics file\n");

        System.out.print(buf);
    }
}
