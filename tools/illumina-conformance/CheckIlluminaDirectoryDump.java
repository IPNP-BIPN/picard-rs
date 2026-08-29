/*
 * `CheckIlluminaDirectory`, taken from the reference.
 *
 * The tool answers one question about a basecalls directory: are the files a run would need there,
 * for every tile and every cycle the read structure asks for. It is the cheapest of the Illumina
 * tools and the one that says what the others will refuse.
 *
 * Six behaviours this is built to catch.
 *
 *   - A COMPLETE DIRECTORY IS ZERO and says so quietly;
 *   - A MISSING CYCLE IS AN EXIT CODE, and the count of what is missing is the code itself rather
 *     than one: the tool returns the NUMBER of failures;
 *   - THE READ STRUCTURE DECIDES HOW MANY CYCLES ARE ASKED FOR, so the same directory passes under
 *     `4T` and fails under `6T`;
 *   - `--FAKE_FILES` WRITES WHAT IS MISSING rather than refusing, which turns a failing run into a
 *     passing one and leaves files behind;
 *   - `--TILE_NUMBERS` NARROWS THE QUESTION to a tile, so asking about a tile that is not there is
 *     a different failure from a cycle that is not;
 *   - AND `--DATA_TYPES` DECIDES WHICH FILES COUNT: a directory with no `.locs` passes when the
 *     positions are not asked for and fails when they are.
 *
 * Output:
 *
 *     code\t<case>\t<the exit status>
 *     files\t<case>\t<the directory's files afterwards, sorted, space separated>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CheckIlluminaDirectoryDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class CheckIlluminaDirectoryDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** Every file under the directory, relative to it and sorted. */
    static String listing(final Path root) throws Exception {
        final List<String> names = new ArrayList<>();
        try (final Stream<Path> walk = Files.walk(root)) {
            walk.filter(Files::isRegularFile)
                    .map(path -> root.relativize(path).toString())
                    .sorted(Comparator.naturalOrder())
                    .forEach(names::add);
        }
        return String.join(" ", names);
    }

    static void run(final String name, final List<String> remove, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("checkillumina");
        final Path basecalls = MakeIlluminaFixtures.write(dir.resolve("run"));
        for (final String relative : remove) {
            Files.deleteIfExists(basecalls.resolve(relative));
        }

        final List<String> argv = new ArrayList<>(List.of("B=" + basecalls, "L=1"));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("READ_STRUCTURE="))) {
            tail.add("READ_STRUCTURE=4T");
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        int code;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            code = new picard.illumina.CheckIlluminaDirectory()
                    .instanceMain(argv.toArray(new String[0]));
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
        emit("code", name, String.valueOf(code));
        emit("files", name, listing(basecalls));
    }

    public static void main(final String[] args) throws Exception {
        // A directory with everything in it.
        run("a-complete-directory", List.of());
        // A cycle that is not there, and the same directory asked for fewer cycles.
        run("a-missing-cycle", List.of("L001/C4.1/s_1_1101.bcl"));
        run("a-missing-cycle-not-asked-for", List.of("L001/C4.1/s_1_1101.bcl"),
                "READ_STRUCTURE=3T");
        // More cycles than the directory has.
        run("more-cycles-than-there-are", List.of(), "READ_STRUCTURE=6T");
        // The filter and the positions, each asked for and each not.
        run("a-missing-filter", List.of("L001/s_1_1101.filter"));
        run("a-missing-locs", List.of("L001/s_1_1101.locs"));
        run("a-missing-locs-with-basecalls-only", List.of("L001/s_1_1101.locs"),
                "DATA_TYPES=BaseCalls");
        // The positions are not in the default set at all, so the file has to be ASKED for before
        // its absence is a failure.
        // The per-tile `.locs` is not the only position file: `s.locs` sits beside the basecalls
        // directory and is what a run with one position file per LANE uses, so a case about the
        // positions has to take that one away.
        run("a-missing-s-locs-with-positions-asked-for", List.of("../s.locs"),
                "DATA_TYPES=Position");
        run("a-missing-s-locs", List.of("../s.locs"));
        run("a-complete-directory-with-positions-asked-for", List.of(), "DATA_TYPES=Position");
        // A tile the directory does not carry, and the one it does.
        run("a-tile-that-is-there", List.of(), "T=1101");
        run("a-tile-that-is-not", List.of(), "T=1102");
        // And the argument that writes what is missing rather than refusing it.
        run("faking-what-is-missing", List.of("L001/C4.1/s_1_1101.bcl"), "F=true");

        System.out.print(buf);
    }
}
