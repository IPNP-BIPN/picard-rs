/*
 * ClusterCrosscheckMetrics' clusters, taken from the reference.
 *
 * The tool reads CrosscheckFingerprints' metrics, builds a graph whose edges are the comparisons
 * above a LOD threshold, and writes the metrics back grouped by connected component. There is no
 * sequence data in it: it is a graph over one table.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - AN EDGE IS A LOD STRICTLY ABOVE THE THRESHOLD, so a comparison exactly at it makes none;
 *   - THE CLUSTERS ARE CONNECTED COMPONENTS, so A related to B and B to C puts all three
 *     together even where A and C were never related;
 *   - THE OUTPUT CARRIES EVERY COMPARISON WHOSE BOTH SIDES ARE IN ONE CLUSTER, whatever its own
 *     LOD, so the low-scoring A-to-C row comes back inside the cluster its neighbours built;
 *   - A GROUP IN NO EDGE AT ALL REACHES NO CLUSTER, so its comparisons vanish from the output
 *     entirely rather than forming a cluster of one;
 *   - CLUSTER_SIZE IS THE NUMBER OF GROUPS AND NOT THE NUMBER OF ROWS;
 *   - THE CLUSTER IDENTIFIER IS A NODE INDEX AND NOT A COUNTER, so two clusters of two are
 *     numbered 0 and 2 rather than 0 and 1: the ids are not contiguous;
 *   - RAISING THE THRESHOLD DROPS AN EDGE AND WITH IT A GROUP: the chain A-B-C under a higher
 *     threshold leaves C in no edge at all, so both of C's rows vanish and the cluster of three
 *     becomes a cluster of two;
 *   - A GROUP COMPARED WITH ITSELF KEEPS ITS ROW and adds no group to the cluster, so a file of
 *     `A-A` and `A-B` reports both rows in a cluster of size two;
 *   - THE ROWS ARE COLLECTED INTO A SET, so a duplicated comparison appears once;
 *   - THE ORDER OF THE OUTPUT IS A HASH ORDER, which is why this dump sorts before emitting;
 *   - AND A FILE WHOSE EVERY COMPARISON IS UNDER THE THRESHOLD WRITES A TABLE OF NO ROWS.
 *
 * Output:
 *
 *     in\t<name>\t<that metrics file, escaped>
 *     metrics\t<case>\t<the output table, its rows sorted, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ClusterCrosscheckMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final String HEADER = String.join("\t",
            "LEFT_GROUP_VALUE", "RIGHT_GROUP_VALUE", "RESULT", "DATA_TYPE", "LOD_SCORE",
            "LOD_SCORE_TUMOR_NORMAL", "LOD_SCORE_NORMAL_TUMOR", "LEFT_RUN_BARCODE", "LEFT_LANE",
            "LEFT_MOLECULAR_BARCODE_SEQUENCE", "LEFT_LIBRARY", "LEFT_SAMPLE", "LEFT_FILE",
            "RIGHT_RUN_BARCODE", "RIGHT_LANE", "RIGHT_MOLECULAR_BARCODE_SEQUENCE",
            "RIGHT_LIBRARY", "RIGHT_SAMPLE", "RIGHT_FILE");

    static String row(final String left, final String right, final double lod) {
        return String.join("\t", left, right,
                lod > 0 ? "EXPECTED_MATCH" : "UNEXPECTED_MISMATCH", "SAMPLE",
                Double.toString(lod), "0", "0",
                "", "", "", "", left, "", "", "", "", "", right, "");
    }

    static String metricsFile(final List<String> rows) {
        final List<String> lines = new ArrayList<>();
        lines.add("## htsjdk.samtools.metrics.StringHeader");
        lines.add("# a fixture");
        lines.add("");
        lines.add("## METRICS CLASS\tpicard.fingerprint.CrosscheckMetric");
        lines.add(HEADER);
        lines.addAll(rows);
        lines.add("");
        return String.join("\n", lines);
    }

    /** The output table without its comments, its data rows SORTED. */
    static String table(final String text) {
        final List<String> head = new ArrayList<>();
        final List<String> rows = new ArrayList<>();
        boolean seenHeader = false;
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            if (!seenHeader) {
                head.add(line);
                seenHeader = true;
            } else {
                rows.add(line);
            }
        }
        Collections.sort(rows);
        head.addAll(rows);
        return String.join("\n", head);
    }

    static void run(final String name, final List<String> rows, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("clustercrosscheck");
        final Path in = dir.resolve("in.crosscheck_metrics");
        Files.writeString(in, metricsFile(rows), StandardCharsets.UTF_8);
        final Path out = dir.resolve("out.txt");
        final List<String> argv = new ArrayList<>(List.of("I=" + in, "O=" + out));
        argv.addAll(List.of(extra));
        try {
            final int code = new picard.fingerprint.ClusterCrosscheckMetrics()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                emit("error", name, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        }
        emit("metrics", name, table(Files.readString(out, StandardCharsets.UTF_8)));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // A is related to B and B to C, but A and C are not: one cluster of three all the same,
        // and the low-scoring A-to-C row comes back inside it.
        final List<String> chain = List.of(
                row("A", "B", 10.0),
                row("B", "C", 10.0),
                row("A", "C", -5.0));
        emit("in", "chain", metricsFile(chain));
        run("transitive-chain", chain, "LOD=3");

        // Raising the threshold past one of the two edges splits the chain.
        run("threshold-splits", List.of(
                row("A", "B", 10.0),
                row("B", "C", 4.0),
                row("A", "C", -5.0)), "LOD=5");

        // A comparison exactly at the threshold makes no edge.
        run("exactly-at-the-threshold", List.of(row("A", "B", 3.0)), "LOD=3");
        run("just-above-the-threshold", List.of(row("A", "B", 3.001)), "LOD=3");

        // Two clusters that never meet.
        run("two-clusters", List.of(
                row("A", "B", 10.0),
                row("C", "D", 10.0),
                row("A", "C", -5.0)), "LOD=3");

        // A group in no edge at all, whose comparisons vanish.
        run("orphan-group", List.of(
                row("A", "B", 10.0),
                row("Z", "A", -5.0)), "LOD=3");

        // A group compared with itself.
        run("self-comparison", List.of(
                row("A", "A", 10.0),
                row("A", "B", 10.0)), "LOD=3");

        // The same comparison twice over.
        run("duplicated-row", List.of(
                row("A", "B", 10.0),
                row("A", "B", 10.0)), "LOD=3");

        // The default threshold, which is nought.
        run("default-threshold", List.of(
                row("A", "B", 0.5),
                row("C", "D", -0.5)));

        // Every comparison under the threshold.
        run("nothing-above-the-threshold", List.of(
                row("A", "B", 1.0),
                row("C", "D", 2.0)), "LOD=3");

        // A file with no comparisons at all.
        run("empty", List.of(), "LOD=3");

        System.out.print(buf);
    }
}
