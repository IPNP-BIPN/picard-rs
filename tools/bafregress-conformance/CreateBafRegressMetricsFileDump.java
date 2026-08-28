/*
 * CreateBafRegressMetricsFile's metrics, taken from the reference.
 *
 * The tool turns bafRegress's own stdout into a Picard metrics file. It parses, and it derives one
 * column the input does not carry.
 *
 * Thirteen behaviours this is built to catch.
 *
 *   - THE OUTPUT ARGUMENT IS A BASENAME, the tool appending `.bafregress_metrics` to it;
 *   - THE HEADER IS COMPARED AS A WHOLE STRING and not by a pattern, so it must be exactly the
 *     seven names separated by single tabs: spaces where the tabs should be are refused;
 *   - THE HEADER'S REFUSAL QUOTES IT IN SINGLE QUOTES, which the rows' refusal does not;
 *   - THE ROWS SPLIT ON RUNS OF WHITESPACE all the same, so a row of spaces parses under a header
 *     of tabs;
 *   - LOG10_PVAL IS DERIVED AND NOT READ, being the base-ten logarithm of the p-value column;
 *   - A P-VALUE OF ZERO GIVES A LOG10_PVAL OF -Infinity, which the writer renders as `-?`;
 *   - A NEGATIVE ESTIMATE IS ACCEPTED, the columns being parsed as doubles rather than matched;
 *   - SO IS EXPONENT NOTATION, `1e-5` parsing where a regex over digits and dots would refuse it;
 *   - A ROW WITH THE WRONG NUMBER OF COLUMNS IS AN IOException wrapped in a PicardException, so
 *     it is a different class from the header's refusal, and its message counts the columns it
 *     found and quotes the row unquoted;
 *   - A FRACTIONAL NHOM IS A NumberFormatException, which is neither of those two: it escapes the
 *     IOException the parser wraps and reaches the caller by itself;
 *   - EXPONENT NOTATION IS REWRITTEN ON THE WAY OUT, `1e-5` coming back as `0.00001`;
 *   - AND A FILE WITH NOTHING IN IT AT ALL IS A NullPointerException, the header comparison being
 *     called on the null the reader answered.
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

public class CreateBafRegressMetricsFileDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final String HEADER = "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom";

    static String file(final List<String> lines) {
        return String.join("\n", lines) + "\n";
    }

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
        final Path dir = Files.createTempDirectory("bafregress");
        final Path in = dir.resolve("in.txt");
        Files.writeString(in, content, StandardCharsets.UTF_8);
        emit("in", name, content);
        final Path base = dir.resolve("out");
        try {
            final int code = new picard.arrays.CreateBafRegressMetricsFile()
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
        final Path written = dir.resolve("out.bafregress_metrics");
        if (!Files.exists(written)) {
            emit("error", name, "no output at " + written.getFileName());
            return;
        }
        emit("name", name, written.getFileName().toString());
        emit("metrics", name, table(Files.readString(written, StandardCharsets.UTF_8)));
    }

    static String row(final String sample, final String estimate, final String stderr,
                      final String tval, final String pval, final String callrate,
                      final String nhom) {
        return String.join("\t", sample, estimate, stderr, tval, pval, callrate, nhom);
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        run("one-row", file(List.of(HEADER,
                row("SAMPLE_A", "0.05", "0.01", "5.0", "0.001", "0.99", "1000"))));
        run("three-rows", file(List.of(HEADER,
                row("SAMPLE_A", "0.05", "0.01", "5.0", "0.001", "0.99", "1000"),
                row("SAMPLE_B", "0.10", "0.02", "5.0", "0.5", "0.98", "2000"),
                row("SAMPLE_C", "0.00", "0.00", "0.0", "1.0", "1.0", "3000"))));

        // The rows split on runs of whitespace even under a header of tabs.
        run("spaces-in-the-row", file(List.of(HEADER,
                "SAMPLE_A   0.05   0.01   5.0   0.001   0.99   1000")));
        // The header is compared as a whole string, so spaces in it are refused.
        run("spaces-in-the-header", file(List.of(
                "sample estimate stderr tval pval callrate Nhom",
                row("SAMPLE_A", "0.05", "0.01", "5.0", "0.001", "0.99", "1000"))));
        // And so is a header naming a column differently.
        run("wrong-header", file(List.of(
                "SAMPLE\testimate\tstderr\ttval\tpval\tcallrate\tNhom",
                row("SAMPLE_A", "0.05", "0.01", "5.0", "0.001", "0.99", "1000"))));

        // A p-value of zero, whose logarithm is not finite.
        run("zero-pvalue", file(List.of(HEADER,
                row("SAMPLE_A", "0.05", "0.01", "5.0", "0", "0.99", "1000"))));
        // A negative estimate and exponent notation, which a double parser takes.
        run("negative-estimate", file(List.of(HEADER,
                row("SAMPLE_A", "-0.05", "0.01", "-5.0", "0.001", "0.99", "1000"))));
        run("exponent-notation", file(List.of(HEADER,
                row("SAMPLE_A", "5e-2", "1e-2", "5.0", "1e-5", "0.99", "1000"))));

        // A fractional Nhom, which the integer parser refuses.
        run("fractional-nhom", file(List.of(HEADER,
                row("SAMPLE_A", "0.05", "0.01", "5.0", "0.001", "0.99", "1000.5"))));
        // A row with a column too few and one too many.
        run("short-row", file(List.of(HEADER, "SAMPLE_A\t0.05\t0.01\t5.0\t0.001\t0.99")));
        run("long-row", file(List.of(HEADER,
                "SAMPLE_A\t0.05\t0.01\t5.0\t0.001\t0.99\t1000\t1")));

        // A file of only a header.
        run("header-only", file(List.of(HEADER)));
        // A file with nothing in it at all.
        run("empty", "");

        System.out.print(buf);
    }
}
