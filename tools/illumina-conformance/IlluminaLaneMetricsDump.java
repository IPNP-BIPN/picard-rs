/*
 * `CollectIlluminaLaneMetrics`, taken from the reference.
 *
 * The tool reads one binary file, `InterOp/TileMetricsOut.bin`, and reports what it says about a
 * lane: how many clusters it had, how many of them passed the filter, and the density of both.
 * Nothing else in the run directory is read, which is what makes it the second cheapest of the
 * Illumina tools and the one that says whether a fixture's tile metrics are well formed at all.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE METRICS COME OFF THE METRIC CODES, not off the basecalls: code 100 is the cluster count
 *     and 101 is the count that passed the filter, and a file carrying only one of the two reports
 *     the other as zero rather than refusing;
 *   - THE DENSITY IS A CODE OF ITS OWN, 102 and 103, so a file with counts and no densities
 *     reports counts and no densities;
 *   - A LANE IS THE SUM OVER ITS TILES, so two tiles of four clusters are one lane of eight;
 *   - TWO LANES ARE TWO ROWS, and each is summed on its own;
 *   - `--OUTPUT_PREFIX` NAMES THE FILE and the suffix is the tool's, so the name is a claim about
 *     the shape of the output rather than a convenience;
 *   - AND THE READ STRUCTURE IS REQUIRED where no `RunInfo.xml` is there to supply it, which is a
 *     refusal rather than a default.
 *
 * Output:
 *
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     metrics\t<case>.<file>\t<the table without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: IlluminaLaneMetricsDump
 */

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class IlluminaLaneMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One record of a tile metrics file: a lane, a tile, a metric code and a value. */
    record Metric(int lane, int tile, int code, float value) {}

    /** A tile metrics file of exactly the records given, in the order given. */
    static void writeTileMetrics(final Path file, final List<Metric> metrics) throws Exception {
        final ByteBuffer buffer = ByteBuffer.allocate(2 + metrics.size() * 10)
                .order(ByteOrder.LITTLE_ENDIAN);
        buffer.put((byte) 2);
        buffer.put((byte) 10);
        for (final Metric metric : metrics) {
            buffer.putShort((short) metric.lane());
            buffer.putShort((short) metric.tile());
            buffer.putShort((short) metric.code());
            buffer.putFloat(metric.value());
        }
        Files.createDirectories(file.getParent());
        Files.write(file, buffer.array());
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

    static void run(final String name, final List<Metric> metrics, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("lanemetrics");
        final Path run = dir.resolve("run");
        writeTileMetrics(run.resolve("InterOp").resolve("TileMetricsOut.bin"), metrics);
        // The tool reads the InterOp file alone, but it asks the run directory to exist as one.
        Files.createDirectories(run.resolve("Data").resolve("Intensities").resolve("BaseCalls"));
        final Path out = Files.createDirectories(dir.resolve("out"));

        final List<String> argv = new ArrayList<>(List.of(
                "RUN_DIRECTORY=" + run, "OUTPUT_DIRECTORY=" + out));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("O="))) {
            tail.add("O=metrics");
        }
        if (tail.stream().noneMatch(argument -> argument.startsWith("READ_STRUCTURE="))) {
            tail.add("READ_STRUCTURE=4T");
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = new picard.illumina.CollectIlluminaLaneMetrics()
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

        final List<String> written = new ArrayList<>();
        for (final File file : out.toFile().listFiles()) {
            written.add(file.getName());
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            emit("metrics", name + "." + file, table(out.resolve(file)));
        }
    }

    /**
     * The codes one tile reports: the four counts, and a phasing pair per template read.
     *
     * The phasing codes are `200 + 2i` and `201 + 2i` for the i-th read DESCRIPTOR, and both of a
     * pair have to be there: a file with one of the two is refused by name rather than reported
     * with a gap. `templates` is how many template reads the read structure declares, which is
     * what decides how many pairs are asked for.
     */
    static List<Metric> tile(final int lane, final int tile, final float clusters,
                             final float passing, final int templates) {
        final List<Metric> metrics = new ArrayList<>(List.of(
                new Metric(lane, tile, 100, clusters),
                new Metric(lane, tile, 101, passing),
                new Metric(lane, tile, 102, clusters * 10),
                new Metric(lane, tile, 103, passing * 10)));
        // The index in the code is the read DESCRIPTOR's, not the template read's: in `4T8B4T` the
        // second template read is the third descriptor, so its codes are 204 and 205 rather than
        // 202 and 203. `descriptors` says which descriptor indices carry a template read.
        for (final int index : descriptors(templates)) {
            metrics.add(new Metric(lane, tile, 200 + index * 2, 0.1f * (index + 1)));
            metrics.add(new Metric(lane, tile, 201 + index * 2, 0.2f * (index + 1)));
        }
        return metrics;
    }

    /** Which descriptor indices the template reads sit at: `4T` is {0} and `4T8B4T` is {0, 2}. */
    static List<Integer> descriptors(final int templates) {
        return templates == 1 ? List.of(0) : List.of(0, 2);
    }

    /** The same, for a run whose read structure has one template read. */
    static List<Metric> tile(final int lane, final int tile, final float clusters,
                             final float passing) {
        return tile(lane, tile, clusters, passing, 1);
    }

    public static void main(final String[] args) throws Exception {
        // One tile, then two, then two lanes: a lane is the sum over its tiles.
        run("one-tile", tile(1, 1101, 1000f, 800f));
        final List<Metric> twoTiles = new ArrayList<>(tile(1, 1101, 1000f, 800f));
        twoTiles.addAll(tile(1, 1102, 500f, 400f));
        run("two-tiles", twoTiles);
        final List<Metric> twoLanes = new ArrayList<>(twoTiles);
        twoLanes.addAll(tile(2, 1101, 200f, 100f));
        run("two-lanes", twoLanes);

        // A file with the counts and no densities, and one with the densities and no counts.
        run("counts-without-densities", List.of(
                new Metric(1, 1101, 100, 1000f),
                new Metric(1, 1101, 101, 800f)));
        run("densities-without-counts", List.of(
                new Metric(1, 1101, 102, 10000f),
                new Metric(1, 1101, 103, 8000f)));
        // And one carrying a code the tool does not read at all.
        run("a-code-nobody-reads", List.of(
                new Metric(1, 1101, 100, 1000f),
                new Metric(1, 1101, 200, 42f)));

        // The output's name, and the read structure that has no RunInfo.xml to come from.
        run("a-named-prefix", tile(1, 1101, 1000f, 800f), "O=run17");
        run("no-read-structure", tile(1, 1101, 1000f, 800f), "READ_STRUCTURE=null");
        // A read structure with two template reads asks for two phasing pairs, and a file with one
        // is refused by name.
        run("two-template-reads", tile(1, 1101, 1000f, 800f, 2), "READ_STRUCTURE=4T8B4T");
        run("two-template-reads-with-one-phasing-pair", tile(1, 1101, 1000f, 800f, 1),
                "READ_STRUCTURE=4T8B4T");
        // And half a phasing pair, which is the refusal the codes are named in.
        final List<Metric> halfAPair = new ArrayList<>(List.of(
                new Metric(1, 1101, 100, 1000f),
                new Metric(1, 1101, 101, 800f),
                new Metric(1, 1101, 102, 10000f),
                new Metric(1, 1101, 103, 8000f),
                new Metric(1, 1101, 200, 0.1f)));
        run("half-a-phasing-pair", halfAPair);

        System.out.print(buf);
    }
}
