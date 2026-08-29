/*
 * A basecalls directory, written byte by byte.
 *
 * Picard's Illumina tools read a directory rather than a file, and the reference's own test data is
 * not in the pinned clone, so the fixture is built here from the formats the readers document:
 *
 *   - a `.bcl` is an unsigned int of clusters and then one byte per cluster, whose low two bits are
 *     the base (`A`, `C`, `G`, `T`) and whose remaining six are the quality. A byte of zero is a
 *     no-call, whatever the quality bits would have said;
 *   - a `.filter` is three unsigned ints (zero, the version, which must be three, and the number of
 *     clusters) and then one byte per cluster, one for a cluster that passed;
 *   - a `.locs` is a little-endian one, a float version of 1.0, an unsigned int of clusters, and
 *     then two floats per cluster.
 *
 * Everything is little-endian, which is what `ByteBuffer.order(LITTLE_ENDIAN)` is for here.
 *
 * Usage: MakeIlluminaFixtures <directory>
 */

import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class MakeIlluminaFixtures {

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
        final Path intensities = run.resolve("Data").resolve("Intensities");
        final Path root = intensities.resolve("BaseCalls");
        writeTileMetrics(run.resolve("InterOp").resolve("TileMetricsOut.bin"), LANE,
                new int[]{TILE});
        writeLocs(intensities.resolve("s.locs"), CLUSTERS);
        final Path lane = root.resolve(String.format("L%03d", LANE));
        for (int cycle = 1; cycle <= CYCLES.length; cycle++) {
            writeBcl(lane.resolve("C" + cycle + ".1")
                            .resolve(String.format("s_%d_%d.bcl", LANE, TILE)),
                    CYCLES[cycle - 1], 30);
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

    public static void main(final String[] args) throws Exception {
        final Path root = Paths.get(args[0]);
        Files.createDirectories(root);
        write(root);

        try (final OutputStream out = System.out) {
            out.write(("wrote " + root + "\n").getBytes(StandardCharsets.UTF_8));
        }
    }
}
