/*
 * ConvertSequencingArtifactToOxoG's metrics, taken from the reference.
 *
 * The tool reads CollectSequencingArtifactMetrics' two detail files and rewrites them as the
 * older OxoG file. There is no sequence data involved at all: it is arithmetic over two tables.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE OUTPUT REPORTS THE `C` CONTEXTS ONLY, so the contexts of the output are those whose
 *     middle base is C in the PRE-ADAPTER input;
 *   - THE PRE-ADAPTER FIGURES COME FROM THE REVERSE COMPLEMENT CONTEXT: a row for `ACA` is built
 *     from the input's `TGT`, because OxoG reverse-complements its contexts;
 *   - THE BAIT-BIAS FIGURES COME FROM THE CONTEXT ITSELF, so one output row draws on two
 *     different input rows;
 *   - ONLY THE C>A AND G>T TRANSITIONS ARE READ, every other row of both inputs being ignored;
 *   - TOTAL_SITES IS ALWAYS ZERO, the input not carrying it;
 *   - THE OXIDATION ERROR RATE HAS A FLOOR OF ONE BASE AND NOT OF A SMALL NUMBER, so a context
 *     with fewer oxidised alternates than unoxidised ones reports `1 / TOTAL_BASES` rather than a
 *     negative rate or a tiny one;
 *   - THE TWO BAIT-BIAS RATES HAVE A FLOOR OF 1e-10 INSTEAD, and they are opposite differences of
 *     the same two numbers, so at most one of them is above that floor;
 *   - THE Q SCORES ARE THE PHRED OF THOSE FLOORED RATES, so the floor of 1e-10 shows up as 100;
 *   - --OUTPUT_BASE DEFAULTS TO --INPUT_BASE, and each of the three file names is derived from a
 *     basename by a fixed extension;
 *   - NAMING NEITHER A BASENAME NOR THE FILE IT WOULD DERIVE IS REFUSED by a message naming both
 *     arguments;
 *   - AND THE ROWS ARE WRITTEN IN A HASH ORDER, both the libraries and the contexts coming out of
 *     HashSets, which is why this dump SORTS them before emitting.
 *
 * Output:
 *
 *     in\t<name>\t<that input file, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, its rows sorted, escaped>
 *     name\t<case>\t<the output file's name>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ConvertSequencingArtifactToOxoGDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final String PRE_ADAPTER_HEADER =
            "SAMPLE_ALIAS\tLIBRARY\tREF_BASE\tALT_BASE\tCONTEXT\tPRO_REF_BASES\tPRO_ALT_BASES\t"
                    + "CON_REF_BASES\tCON_ALT_BASES\tERROR_RATE\tQSCORE";
    static final String BAIT_BIAS_HEADER =
            "SAMPLE_ALIAS\tLIBRARY\tREF_BASE\tALT_BASE\tCONTEXT\tFWD_CXT_REF_BASES\t"
                    + "FWD_CXT_ALT_BASES\tREV_CXT_REF_BASES\tREV_CXT_ALT_BASES\tFWD_ERROR_RATE\t"
                    + "REV_ERROR_RATE\tERROR_RATE\tQSCORE";

    /** A Picard metrics file: the banner the reader wants, then the class line and the table. */
    static String metricsFile(final String beanClass, final String header, final List<String> rows) {
        final List<String> lines = new ArrayList<>();
        lines.add("## htsjdk.samtools.metrics.StringHeader");
        lines.add("# a fixture");
        lines.add("");
        lines.add("## METRICS CLASS\t" + beanClass);
        lines.add(header);
        lines.addAll(rows);
        lines.add("");
        return String.join("\n", lines);
    }

    static String preAdapterRow(final String library, final char ref, final char alt,
                                final String context, final long proRef, final long proAlt,
                                final long conRef, final long conAlt) {
        return String.join("\t", "sample1", library, String.valueOf(ref), String.valueOf(alt),
                context, Long.toString(proRef), Long.toString(proAlt), Long.toString(conRef),
                Long.toString(conAlt), "0", "0");
    }

    static String baitBiasRow(final String library, final char ref, final char alt,
                              final String context, final long fwdRef, final long fwdAlt,
                              final long revRef, final long revAlt) {
        return String.join("\t", "sample1", library, String.valueOf(ref), String.valueOf(alt),
                context, Long.toString(fwdRef), Long.toString(fwdAlt), Long.toString(revRef),
                Long.toString(revAlt), "0", "0", "0", "0");
    }

    /** The metrics table without its comments, its data rows SORTED. */
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

    static void run(final String name, final List<String> preAdapter, final List<String> baitBias,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("artifacttooxog");
        final Path base = dir.resolve("in");
        Files.writeString(dir.resolve("in.pre_adapter_detail_metrics"),
                metricsFile("picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterDetailMetrics",
                        PRE_ADAPTER_HEADER, preAdapter), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("in.bait_bias_detail_metrics"),
                metricsFile("picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasDetailMetrics",
                        BAIT_BIAS_HEADER, baitBias), StandardCharsets.UTF_8);
        final List<String> argv = new ArrayList<>();
        if (extra.length == 0 || !Arrays.asList(extra).contains("NO_INPUT_BASE")) {
            argv.add("INPUT_BASE=" + base);
        }
        for (final String value : extra) {
            if (!value.equals("NO_INPUT_BASE")) {
                argv.add(value.replace("<dir>", dir.toString()));
            }
        }
        try {
            final int code = new picard.analysis.artifacts.ConvertSequencingArtifactToOxoG()
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
        Path written = dir.resolve("in.oxog_metrics");
        for (final String value : extra) {
            if (value.startsWith("OXOG_OUT=")) {
                written = Path.of(value.substring("OXOG_OUT=".length()).replace("<dir>", dir.toString()));
            } else if (value.startsWith("OUTPUT_BASE=")) {
                written = Path.of(value.substring("OUTPUT_BASE=".length())
                        .replace("<dir>", dir.toString()) + ".oxog_metrics");
            }
        }
        if (!Files.exists(written)) {
            emit("error", name, "no output at " + written.getFileName());
            return;
        }
        emit("name", name, written.getFileName().toString());
        emit("metrics", name, table(Files.readString(written, StandardCharsets.UTF_8)));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // One C context, `ACA`, whose pre-adapter data sits under its reverse complement `TGT`.
        // The C>A row of ACA is what names the context; the G>T row of TGT is what is read.
        final List<String> preAdapter = List.of(
                preAdapterRow("lib1", 'C', 'A', "ACA", 100, 5, 200, 3),
                preAdapterRow("lib1", 'G', 'T', "TGT", 1000, 40, 2000, 10),
                // A transition the tool ignores.
                preAdapterRow("lib1", 'A', 'G', "AAA", 7, 7, 7, 7));
        final List<String> baitBias = List.of(
                baitBiasRow("lib1", 'C', 'A', "ACA", 900, 30, 800, 10),
                baitBiasRow("lib1", 'G', 'T', "TGT", 1, 1, 1, 1),
                baitBiasRow("lib1", 'A', 'G', "AAA", 7, 7, 7, 7));
        emit("in", "pre-adapter", metricsFile(
                "picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterDetailMetrics",
                PRE_ADAPTER_HEADER, preAdapter));
        emit("in", "bait-bias", metricsFile(
                "picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasDetailMetrics",
                BAIT_BIAS_HEADER, baitBias));
        run("one-context", preAdapter, baitBias);

        // The oxidation floor: fewer oxidised alternates than unoxidised ones.
        run("oxidation-floor",
                List.of(preAdapterRow("lib1", 'C', 'A', "ACA", 100, 5, 200, 3),
                        preAdapterRow("lib1", 'G', 'T', "TGT", 1000, 2, 2000, 50)),
                baitBias);

        // The bait-bias floor: the two rates equal, so both differences are nought.
        run("bait-bias-floor", preAdapter,
                List.of(baitBiasRow("lib1", 'C', 'A', "ACA", 900, 100, 900, 100),
                        baitBiasRow("lib1", 'G', 'T', "TGT", 1, 1, 1, 1)));
        // And the other direction, where the reverse rate is the larger.
        run("bait-bias-reversed", preAdapter,
                List.of(baitBiasRow("lib1", 'C', 'A', "ACA", 900, 10, 800, 30),
                        baitBiasRow("lib1", 'G', 'T', "TGT", 1, 1, 1, 1)));

        // Two contexts and two libraries.
        final List<String> two = new ArrayList<>();
        final List<String> twoBait = new ArrayList<>();
        for (final String library : List.of("lib1", "lib2")) {
            two.add(preAdapterRow(library, 'C', 'A', "ACA", 100, 5, 200, 3));
            two.add(preAdapterRow(library, 'G', 'T', "TGT", 1000, 40, 2000, 10));
            two.add(preAdapterRow(library, 'C', 'A', "TCT", 100, 5, 200, 3));
            two.add(preAdapterRow(library, 'G', 'T', "AGA", 500, 20, 900, 5));
            twoBait.add(baitBiasRow(library, 'C', 'A', "ACA", 900, 30, 800, 10));
            twoBait.add(baitBiasRow(library, 'C', 'A', "TCT", 700, 20, 600, 5));
            twoBait.add(baitBiasRow(library, 'G', 'T', "TGT", 1, 1, 1, 1));
            twoBait.add(baitBiasRow(library, 'G', 'T', "AGA", 1, 1, 1, 1));
        }
        run("two-libraries-two-contexts", two, twoBait);

        // The file names: a basename for the output, and the file named outright.
        run("output-base", preAdapter, baitBias, "OUTPUT_BASE=<dir>/other");
        run("oxog-out", preAdapter, baitBias, "OXOG_OUT=<dir>/named.txt");

        // Neither a basename nor the files it would derive.
        run("no-basename", preAdapter, baitBias, "NO_INPUT_BASE");

        System.out.print(buf);
    }
}
