/*
 * `BaitDesigner`, taken from the reference.
 *
 * The tool turns a list of targets into a list of baits: fixed-length probes laid along each
 * target, with primers on each end, ready to be ordered. What decides the answer is the strategy,
 * the bait's size and offset, and what the tool does with a target too small to tile.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE BAITS ARE LAID AT A FIXED OFFSET by default, so a target longer than one bait is covered
 *     by several that overlap;
 *   - A TARGET SHORTER THAN A BAIT STILL GETS `--MINIMUM_BAITS_PER_TARGET` of them, which means
 *     baits that run off both ends of it;
 *   - THE THREE STRATEGIES LAY THE SAME TARGET DIFFERENTLY, and `CenteredConstrained` is the one
 *     that keeps a bait centred on a short target rather than tiling past it;
 *   - THE PRIMERS ARE PREPENDED AND APPENDED to every bait's sequence in the Agilent file and to
 *     none of the interval lists;
 *   - `--MERGE_NEARBY_TARGETS` JOINS TWO TARGETS a bait could cover together, so two intervals
 *     become one design;
 *   - THE DESIGN'S SUMMARY IS A METRICS FILE, with the bait and target territory and the
 *     percentage of one covered by the other;
 *   - AND THE OUTPUT IS A DIRECTORY of files named after the design.
 *
 * Output:
 *
 *     files\t<case>\t<the files written, sorted, space separated>
 *     baits\t<case>\t<the bait interval list's lines, escaped>
 *     pool\t<case>\t<the pool file: one row per bait, its sequence between the primers>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BaitDesignerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;

public class BaitDesignerDump {

    static final StringBuilder buf = new StringBuilder();
    static final int LENGTH = 2000;

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A contig of `ACGT` repeating, so a bait's sequence is known by arithmetic. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int index = 0; index < LENGTH; index++) {
            bases.append("ACGT".charAt(index % 4));
        }
        final StringBuilder text = new StringBuilder(">chr1\n");
        for (int index = 0; index < bases.length(); index += 60) {
            text.append(bases, index, Math.min(index + 60, bases.length())).append('\n');
        }
        return text.toString();
    }

    /** An interval list over that contig, with the intervals given as `start end name`. */
    static String targets(final List<int[]> intervals) {
        final StringBuilder text = new StringBuilder("@HD\tVN:1.6\tSO:coordinate\n");
        text.append("@SQ\tSN:chr1\tLN:").append(LENGTH).append('\n');
        int index = 0;
        for (final int[] interval : intervals) {
            text.append("chr1\t").append(interval[0]).append('\t').append(interval[1])
                    .append("\t+\ttarget").append(++index).append('\n');
        }
        return text.toString();
    }

    /** A file's lines, without the SAM-style header an interval list carries. */
    static String lines(final Path file, final boolean dropHeader) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.isEmpty() || (dropHeader && (line.startsWith("@") || line.startsWith("#")))) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final List<int[]> intervals, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("baitdesigner");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path targets = dir.resolve("targets.interval_list");
        Files.writeString(targets, targets(intervals), StandardCharsets.UTF_8);
        final Path out = dir.resolve("design");

        final List<String> argv = new ArrayList<>(List.of(
                "TARGETS=" + targets, "R=" + reference, "DESIGN_NAME=design",
                "OUTPUT_DIRECTORY=" + out));
        argv.addAll(Arrays.asList(extra));

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.util.BaitDesigner()
                    .instanceMain(argv.toArray(new String[0]));
            System.setOut(original);
            System.setErr(originalError);
            if (code != 0) {
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
        try (final Stream<Path> walk = Files.walk(out)) {
            walk.filter(Files::isRegularFile)
                    .forEach(path -> written.add(path.getFileName().toString()));
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            if (file.endsWith("baits.interval_list")) {
                emit("baits", name, lines(out.resolve(file), true));
            } else if (file.endsWith("pool0.design.txt")) {
                // The pool file is the order: one row per bait, its sequence between the primers.
                // `design_parameters.txt` is skipped, because it is the tool's own usage text
                // rather than anything the design decided.
                emit("pool", name, lines(out.resolve(file), true));
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        // A target longer than a bait, tiled at the default offset.
        final List<int[]> one = List.of(new int[]{201, 500});
        run("a-target-longer-than-a-bait", one);
        // A target shorter than a bait, which still gets the minimum number of them.
        run("a-target-shorter-than-a-bait", List.of(new int[]{201, 260}));
        run("a-short-target-with-one-bait", List.of(new int[]{201, 260}),
                "MINIMUM_BAITS_PER_TARGET=1");

        // The four strategies over the same short target.
        // The three strategies the reference declares, over the same short target.
        for (final String strategy : new String[]{
                "CenteredConstrained", "FixedOffset", "Simple"}) {
            run("strategy-" + strategy.toLowerCase(), List.of(new int[]{201, 260}),
                    "DESIGN_STRATEGY=" + strategy);
        }

        // The bait's shape.
        run("a-smaller-bait", one, "BAIT_SIZE=60");
        run("a-wider-offset", one, "BAIT_OFFSET=120");
        run("with-padding", one, "PADDING=50");

        // Two targets a bait could cover together, merged and not.
        final List<int[]> two = List.of(new int[]{201, 260}, new int[]{301, 360});
        run("two-nearby-targets", two);
        run("two-nearby-targets-unmerged", two, "MERGE_NEARBY_TARGETS=false");

        System.out.print(buf);
    }
}
