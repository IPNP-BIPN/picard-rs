/*
 * `CollectHiSeqXPfFailMetrics`, taken from the reference.
 *
 * The tool asks WHY a cluster failed the filter. A HiSeq X calls a cluster non-PF for one of a few
 * reasons, and the tool separates them by looking at the first cycles of the read: a cluster whose
 * bases are missing entirely is one kind of failure, one whose bases are there but poor is
 * another, and one that looks fine is a third.
 *
 * Five behaviours this is built to catch.
 *
 *   - THE METRICS ARE PER TILE, not per lane, so a run of one tile is one row;
 *   - THE NON-PF CLUSTERS ARE CLASSIFIED, and the classes sum to the count the filter file gives;
 *   - `--N_CYCLES` IS READ BEFORE IT IS ASSIGNED: the read structure is a FINAL field initialised
 *     from the argument's default, and a field initialiser runs before the parser sets anything,
 *     so a run always looks at twenty-four cycles whatever the command line asked for. That is the
 *     reference's behaviour and the reason this fixture has twenty-four cycles at all;
 *   - `--PROB_EXPLICIT_READS` WRITES A SECOND FILE, sampling the clusters it explains;
 *   - AND THE LANE MUST EXIST, which is a refusal like the other Illumina tools'.
 *
 * Output:
 *
 *     files\t<case>\t<the files written, sorted, space separated>
 *     metrics\t<case>.<file>\t<the table without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: PfFailMetricsDump
 */

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

public class PfFailMetricsDump {

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

    static void run(final String name, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("pffail");
        // Twenty-four cycles, because the tool builds its read structure from `N_CYCLES` in a
        // field initialiser: the assignment from the command line happens later, so the structure
        // is always the default's twenty-four whatever `--N_CYCLES` says.
        final Path basecalls = MakeIlluminaFixtures.write(dir.resolve("run"), 24);
        final Path out = Files.createDirectories(dir.resolve("out"));

        final List<String> argv = new ArrayList<>(List.of(
                "BASECALLS_DIR=" + basecalls, "OUTPUT=" + out.resolve("pf")));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("LANE="))) {
            tail.add("LANE=1");
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.illumina.quality.CollectHiSeqXPfFailMetrics()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                System.setOut(original);
                System.setErr(originalError);
                final List<String> reasons = new ArrayList<>();
                for (final String line : said.toString(StandardCharsets.UTF_8).split("\n", -1)) {
                    if (line.startsWith("ERROR:")) {
                        reasons.add(line);
                    }
                }
                emit("error", name, "exit " + code
                        + (reasons.isEmpty() ? "" : " " + String.join(" ", reasons)));
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
            emit("metrics", name + "." + file, table(out.resolve(file)));
        }
    }

    public static void main(final String[] args) throws Exception {
        run("the-whole-lane");
        // `--N_CYCLES` is READ BEFORE IT IS ASSIGNED, so neither of these changes the answer: the
        // read structure was built from the default in a field initialiser.
        run("two-cycles", "N_CYCLES=2");
        run("forty-cycles", "N_CYCLES=40");
        // The sampled file, which only exists when the probability is not zero.
        run("with-explicit-reads", "PROB_EXPLICIT_READS=1");
        // A lane the run does not have.
        run("a-lane-that-is-not-there", "LANE=2");

        System.out.print(buf);
    }
}
