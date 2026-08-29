/*
 * `CollectIlluminaBasecallingMetrics`, taken from the reference.
 *
 * The tool walks a lane's clusters and counts them per barcode: how many there were, how many
 * passed the filter, and how many bases of each were called at all. It is the summary a run gets
 * after `ExtractIlluminaBarcodes` has decided which barcode each cluster is.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE COUNTS ARE PER BARCODE, with a row for the clusters that matched none;
 *   - `--INPUT` IS THE BARCODE FILE and is optional: without it the whole lane is one row, and the
 *     barcode column is empty rather than absent;
 *   - THE PF COUNTS COME FROM THE FILTER FILE, so the fixture's one failing cluster moves the
 *     PF columns and not the totals;
 *   - THE MEAN CLUSTERS PER TILE IS A MEAN OVER TILES, so one tile makes it the count itself;
 *   - THE READ STRUCTURE DECIDES WHICH CYCLES ARE BASES, so `2T2B` counts two bases a cluster and
 *     `4T` counts four;
 *   - AND A LANE THE RUN DOES NOT HAVE IS A REFUSAL rather than an empty table.
 *
 * Output:
 *
 *     metrics\t<case>\t<the table without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BasecallingMetricsDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BasecallingMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static String table(final Path file) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final boolean withBarcodes, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("basecallingmetrics");
        final Path basecalls = MakeIlluminaFixtures.write(dir.resolve("run"));
        final Path out = dir.resolve("metrics.txt");

        final List<String> argv = new ArrayList<>(List.of(
                "BASECALLS_DIR=" + basecalls, "LANE=1", "OUTPUT=" + out));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("READ_STRUCTURE="))) {
            tail.add("READ_STRUCTURE=4T");
        }
        if (withBarcodes) {
            // The barcode each cluster was assigned is `ExtractIlluminaBarcodes`' answer, so the
            // pipeline is the reference's own: extract, then count.
            final Path barcodes = dir.resolve("barcodes.tsv");
            Files.writeString(barcodes,
                    "barcode_sequence_1\tbarcode_name\tlibrary_name\nAG\tfirst\tlibraryA\n"
                            + "CT\tsecond\tlibraryB\n", StandardCharsets.UTF_8);
            new picard.illumina.ExtractIlluminaBarcodes().instanceMain(new String[]{
                    "BASECALLS_DIR=" + basecalls, "LANE=1", "READ_STRUCTURE=2T2B",
                    "BARCODE_FILE=" + barcodes,
                    "METRICS_FILE=" + dir.resolve("barcode-metrics.txt"),
                    "OUTPUT_DIR=" + basecalls});
            tail.add("BARCODES_DIR=" + basecalls);
            tail.add("INPUT=" + barcodes);
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = new picard.illumina.CollectIlluminaBasecallingMetrics()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                System.setOut(original);
                System.setErr(originalError);
                emit("error", name, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            System.setOut(original);
            System.setErr(originalError);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        } finally {
            System.setOut(original);
            System.setErr(originalError);
        }
        emit("metrics", name, table(out));
    }

    public static void main(final String[] args) throws Exception {
        // The lane as one row, then split by the barcode the extractor assigned.
        run("the-whole-lane", false);
        run("by-barcode", true, "READ_STRUCTURE=2T2B");
        // The read structure decides how many bases a cluster contributes.
        run("a-barcode-segment-without-barcodes", false, "READ_STRUCTURE=2T2B");
        run("a-skipped-segment", false, "READ_STRUCTURE=2T2S");
        // A lane the run does not have.
        run("a-lane-that-is-not-there", false, "LANE=2");

        System.out.print(buf);
    }
}
