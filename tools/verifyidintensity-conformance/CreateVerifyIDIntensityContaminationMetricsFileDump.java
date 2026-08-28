/*
 * CreateVerifyIDIntensityContaminationMetricsFile's metrics, taken from the reference.
 *
 * The tool turns VerifyIDIntensity's own stdout into a Picard metrics file. It parses, and what
 * is measured is what its three regexes accept and what they do with everything else.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE OUTPUT ARGUMENT IS A BASENAME AND NOT A FILE, the tool appending
 *     `.verifyidintensity_metrics` to whatever it is given;
 *   - THE FIRST TWO LINES ARE FIXED: a header naming the four columns and then a run of dashes,
 *     each matched by its own pattern;
 *   - THE DASHES ARE NOT COUNTED, so one dash is as good as forty;
 *   - THE COLUMNS ARE SPLIT ON RUNS OF WHITESPACE, so tabs and single spaces both parse;
 *   - THE FRACTION MAY OPEN ON A DOT, `.5` parsing as five tenths;
 *   - THE LIKELIHOODS MAY BE NEGATIVE and the fraction may NOT: a negative `%Mix` is refused
 *     while a negative `LLK` is the ordinary case;
 *   - THE ID MUST BE AN UNSIGNED INTEGER, so a negative one is refused;
 *   - AN UNRECOGNISED LINE IS REFUSED BY A MESSAGE QUOTING IT and naming the input's path;
 *   - A FILE THAT ENDS EARLY IS A NullPointerException AND NOT A PicardException: the reader
 *     answers null and the matcher is handed it without a guard;
 *   - THE NUMBERS ARE REFORMATTED ON THE WAY OUT, `-1.0` coming back as `-1`;
 *   - AND A FILE OF ONLY A HEADER AND DASHES WRITES A METRICS FILE WITH NO TABLE AT ALL, not even
 *     the column line: the writer emits its comments and stops.
 *
 * Output:
 *
 *     in\t<case>\t<the input file, escaped>
 *     metrics\t<case>\t<the metrics file without its comments, escaped>
 *     name\t<case>\t<the output file's name>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CreateVerifyIDIntensityContaminationMetricsFileDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final String HEADER = "ID\t%Mix\tLLK\tLLK0";
    static final String DASHES = "----------------------------";

    static String file(final List<String> lines) {
        return String.join("\n", lines) + "\n";
    }

    /** The metrics table without its comment lines. */
    static String table(final String text) {
        final List<String> kept = new ArrayList<>();
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final String content) throws Exception {
        final Path dir = Files.createTempDirectory("verifyidintensity");
        final Path in = dir.resolve("in.txt");
        Files.writeString(in, content, StandardCharsets.UTF_8);
        emit("in", name, content);
        final Path base = dir.resolve("out");
        try {
            final int code = new picard.arrays.CreateVerifyIDIntensityContaminationMetricsFile()
                    .instanceMain(new String[]{"I=" + in, "O=" + base});
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
        final Path written = dir.resolve("out.verifyidintensity_metrics");
        if (!Files.exists(written)) {
            emit("error", name, "no output at " + written.getFileName());
            return;
        }
        emit("name", name, written.getFileName().toString());
        emit("metrics", name, table(Files.readString(written, StandardCharsets.UTF_8)));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        run("one-row", file(List.of(HEADER, DASHES, "0\t0.05\t-1234.5\t-2345.6")));
        run("three-rows", file(List.of(HEADER, DASHES,
                "0\t0.05\t-1234.5\t-2345.6",
                "1\t0.5\t-1.0\t-2.0",
                "2\t0\t0\t0")));

        // The columns split on runs of whitespace, so single spaces parse as well as tabs.
        run("spaces-not-tabs", file(List.of("ID   %Mix   LLK   LLK0", DASHES,
                "0   0.05   -1234.5   -2345.6")));
        // And so does a single dash.
        run("one-dash", file(List.of(HEADER, "-", "0\t0.05\t-1234.5\t-2345.6")));

        // A fraction that opens on a dot.
        run("leading-dot", file(List.of(HEADER, DASHES, "0\t.5\t-1.0\t-2.0")));
        // A positive likelihood, which the pattern also allows.
        run("positive-likelihood", file(List.of(HEADER, DASHES, "0\t0.05\t1234.5\t2345.6")));

        // A negative fraction, which it does not.
        run("negative-fraction", file(List.of(HEADER, DASHES, "0\t-0.05\t-1.0\t-2.0")));
        // A negative id, likewise.
        run("negative-id", file(List.of(HEADER, DASHES, "-1\t0.05\t-1.0\t-2.0")));
        // A row with a column too few.
        run("short-row", file(List.of(HEADER, DASHES, "0\t0.05\t-1.0")));
        // A header that is not the header.
        run("wrong-header", file(List.of("SAMPLE\t%Mix\tLLK\tLLK0", DASHES, "0\t0.05\t-1.0\t-2.0")));

        // A file of only a header and dashes.
        run("no-rows", file(List.of(HEADER, DASHES)));
        // A file of only a header.
        run("header-only", file(List.of(HEADER)));
        // A file with nothing in it at all.
        run("empty", "");

        System.out.print(buf);
    }
}
