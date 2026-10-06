import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * The Illumina corpus: run directories written byte by byte by IlluminaRun, which is
 * tools/illumina-conformance's MakeIlluminaFixtures renamed so the two copies never meet on one
 * classpath, and the barcode and parameter files the tools read beside them.
 *
 * ill_run is four cycles of one tile of one lane with four clusters, three of which pass the
 * filter; ill_run24 is the same run carried to twenty-four cycles, which CollectHiSeqXPfFailMetrics
 * needs. Cycles three and four read AG for the first two clusters and CT for the last two, which is
 * what the barcode files declare.
 */
public class MakeIlluminaCoverageFixtures {

    /** A run directory with nothing but its tile metrics and, when given, a RunInfo.xml. */
    static void laneRun(final Path run, final float[][] metrics, final String reads) throws IOException {
        final ByteBuffer buffer = ByteBuffer.allocate(2 + metrics.length * 10).order(ByteOrder.LITTLE_ENDIAN);
        buffer.put((byte) 2);
        buffer.put((byte) 10);
        for (final float[] m : metrics) {
            buffer.putShort((short) m[0]);
            buffer.putShort((short) m[1]);
            buffer.putShort((short) m[2]);
            buffer.putFloat(m[3]);
        }
        Files.createDirectories(run.resolve("InterOp"));
        Files.write(run.resolve("InterOp").resolve("TileMetricsOut.bin"), buffer.array());
        if (reads != null) {
            Files.writeString(run.resolve("RunInfo.xml"),
                    "<?xml version=\"1.0\"?>\n<RunInfo><Run><Reads>" + reads + "</Reads></Run></RunInfo>\n",
                    StandardCharsets.UTF_8);
        }
    }

    static void write(final File into) throws Exception {
        final Path dir = into.toPath();
        // CollectIlluminaLaneMetrics reads only the tile metrics. Lane 1 of two tiles with the
        // phasing of the first template read; lanes 1 and 2 with the phasing of the two template
        // reads of 4T8B4T (descriptors 0 and 2); and counts and densities with no phasing at all.
        laneRun(dir.resolve("ill_lanes_a"), new float[][]{
                {1, 1101, 100, 1000}, {1, 1101, 101, 800}, {1, 1101, 102, 10000}, {1, 1101, 103, 8000},
                {1, 1101, 200, 0.1f}, {1, 1101, 201, 0.2f},
                {1, 1102, 100, 1200}, {1, 1102, 101, 900}, {1, 1102, 102, 12500}, {1, 1102, 103, 9000},
                {1, 1102, 200, 0.15f}, {1, 1102, 201, 0.25f}},
                "<Read Number=\"1\" NumCycles=\"4\" IsIndexedRead=\"N\"/>");
        laneRun(dir.resolve("ill_lanes_b"), new float[][]{
                {1, 1101, 100, 1000}, {1, 1101, 101, 800}, {1, 1101, 102, 10000}, {1, 1101, 103, 8000},
                {1, 1101, 200, 0.1f}, {1, 1101, 201, 0.2f}, {1, 1101, 204, 0.3f}, {1, 1101, 205, 0.4f},
                {2, 1101, 100, 2000}, {2, 1101, 101, 1500}, {2, 1101, 102, 20000}, {2, 1101, 103, 15000},
                {2, 1101, 200, 0.5f}, {2, 1101, 201, 0.6f}, {2, 1101, 204, 0.7f}, {2, 1101, 205, 0.8f}},
                "<Read Number=\"1\" NumCycles=\"4\" IsIndexedRead=\"N\"/><Read Number=\"2\" NumCycles=\"8\" IsIndexedRead=\"Y\"/><Read Number=\"3\" NumCycles=\"4\" IsIndexedRead=\"N\"/>");
        laneRun(dir.resolve("ill_lanes_bare"), new float[][]{
                {1, 1101, 100, 1000}, {1, 1101, 101, 800}, {1, 1101, 102, 10000}, {1, 1101, 103, 8000}},
                null);
        IlluminaRun.write(dir.resolve("ill_run"));
        IlluminaRun.write(dir.resolve("ill_run24"), 24);
        // What ExtractIlluminaBarcodes writes per tile, by hand: the observed barcode, whether it
        // matched, the barcode matched, and the mismatch counts. The last cluster matches nothing.
        Files.createDirectories(dir.resolve("ill_bcdir"));
        Files.writeString(dir.resolve("ill_bcdir").resolve("s_1_1101_barcode.txt"),
                "AG\tY\tAG\t0\t0\nAG\tY\tAG\t0\t0\nCT\tY\tCT\t0\t0\nCT\tN\tNN\t1\t1\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("ill_barcodes.txt"),
                "barcode_sequence_1\tbarcode_name\tlibrary_name\nAG\tfirst\tlibraryA\nCT\tsecond\tlibraryB\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("ill_barcodes_one.txt"),
                "barcode_sequence_1\tbarcode_name\tlibrary_name\nAG\tfirst\tlibraryA\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("ill_multiplex.txt"),
                "OUTPUT_PREFIX\tBARCODE_1\n/work/out/first\tAG\n/work/out/second\tCT\n/work/out/rest\tN\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("ill_library.txt"),
                "OUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tBARCODE_1\n/work/out/first.bam\tsampleA\tlibraryA\tAG\n"
                        + "/work/out/second.bam\tsampleB\tlibraryB\tCT\n/work/out/rest.bam\trest\trest\tN\n",
                StandardCharsets.UTF_8);
    }
}

final class IlluminaRun {

    /** The bases a cycle's file carries, one character per cluster. */
    static final String[] CYCLES = {"ACGT", "ACGT", "AACC", "GGTT"};
    static final int CLUSTERS = 4;
    static final int LANE = 1;
    static final int TILE = 1101;

    /** `A`, `C`, `G` and `T` are 0, 1, 2 and 3 in the low two bits of a basecall byte. */
    static int code(final char base) {
        switch (base) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: throw new IllegalArgumentException("not a base: " + base);
        }
    }

    static byte basecall(final char base, final int quality) {
        // The quality occupies the six high bits, so a byte of zero can only be a no-call.
        return (byte) (code(base) | (quality << 2));
    }

    static void writeBcl(final Path file, final String bases, final int quality) throws IOException {
        final ByteBuffer buffer = ByteBuffer.allocate(4 + bases.length())
                .order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(bases.length());
        for (final char base : bases.toCharArray()) {
            buffer.put(basecall(base, quality));
        }
        Files.createDirectories(file.getParent());
        Files.write(file, buffer.array());
    }

    static void writeFilter(final Path file, final boolean[] passed) throws IOException {
        final ByteBuffer buffer = ByteBuffer.allocate(12 + passed.length)
                .order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(0);
        buffer.putInt(3);
        buffer.putInt(passed.length);
        for (final boolean pass : passed) {
            buffer.put((byte) (pass ? 1 : 0));
        }
        Files.createDirectories(file.getParent());
        Files.write(file, buffer.array());
    }

    static void writeLocs(final Path file, final int clusters) throws IOException {
        final ByteBuffer buffer = ByteBuffer.allocate(12 + clusters * 8)
                .order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(1);
        buffer.putFloat(1.0f);
        buffer.putInt(clusters);
        for (int cluster = 0; cluster < clusters; cluster++) {
            // Coordinates a hundred apart, so no two clusters share a position.
            buffer.putFloat(100.0f * (cluster + 1));
            buffer.putFloat(200.0f * (cluster + 1));
        }
        Files.createDirectories(file.getParent());
        Files.write(file, buffer.array());
    }

    /**
     * A tile metrics file, which is what tells the tools which tiles a lane HAS.
     *
     * Version two: a byte of version, a byte of record size, and then ten bytes per record, which
     * are the lane, the tile and the metric code as unsigned shorts and the value as a float.
     */
    static void writeTileMetrics(final Path file, final int lane, final int[] tiles)
            throws IOException {
        final ByteBuffer buffer = ByteBuffer.allocate(2 + tiles.length * 10)
                .order(ByteOrder.LITTLE_ENDIAN);
        buffer.put((byte) 2);
        buffer.put((byte) 10);
        for (final int tile : tiles) {
            buffer.putShort((short) lane);
            buffer.putShort((short) tile);
            // Code 100 is the cluster count, which is the metric the tools read a tile's existence
            // off; the value itself is only used by the metrics tools.
            buffer.putShort((short) 100);
            buffer.putFloat(CLUSTERS);
        }
        Files.createDirectories(file.getParent());
        Files.write(file, buffer.array());
    }

    /**
     * The whole RUN directory, which is what the tools take rather than a basecalls directory on
     * its own: `<run>/Data/Intensities/BaseCalls` is what `--BASECALLS_DIR` names, `s.locs` sits
     * beside it in `Intensities`, and `InterOp/TileMetricsOut.bin` is two levels above that.
     *
     * Four cycles, one lane, one tile, four clusters.
     */
    static Path write(final Path run) throws IOException {
        return write(run, CYCLES.length);
    }

    /**
     * The same directory with a chosen number of CYCLES, for a tool that needs more than four.
     *
     * `CollectHiSeqXPfFailMetrics` builds its read structure from `N_CYCLES` in a field
     * initialiser, which runs before the parser assigns the argument, so it always asks for
     * twenty-four cycles whatever the command line says. A fixture for it has to have them.
     */
    static Path write(final Path run, final int cycles) throws IOException {
        final Path intensities = run.resolve("Data").resolve("Intensities");
        final Path root = intensities.resolve("BaseCalls");
        writeTileMetrics(run.resolve("InterOp").resolve("TileMetricsOut.bin"), LANE,
                new int[]{TILE});
        writeLocs(intensities.resolve("s.locs"), CLUSTERS);
        final Path lane = root.resolve(String.format("L%03d", LANE));
        for (int cycle = 1; cycle <= cycles; cycle++) {
            // Past the fourth cycle the four patterns repeat, so a longer run carries the same
            // four clusters and nothing new to reason about.
            writeBcl(lane.resolve("C" + cycle + ".1")
                            .resolve(String.format("s_%d_%d.bcl", LANE, TILE)),
                    CYCLES[(cycle - 1) % CYCLES.length], 30);
        }
        // Three of the four clusters pass the filter, which is what makes a `PF` count worth
        // reading: a tool that ignored the filter would report four.
        writeFilter(lane.resolve(String.format("s_%d_%d.filter", LANE, TILE)),
                new boolean[]{true, true, true, false});
        writeLocs(lane.resolve(String.format("s_%d_%d.locs", LANE, TILE)), CLUSTERS);
        Files.writeString(root.resolve("config.xml"), "<BaseCallAnalysis/>\n",
                StandardCharsets.UTF_8);
        return root;
    }
}
