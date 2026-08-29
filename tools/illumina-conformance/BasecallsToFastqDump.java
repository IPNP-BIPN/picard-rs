/*
 * `IlluminaBasecallsToFastq`, taken from the reference.
 *
 * This is the tool the whole Illumina family exists for: it turns a basecalls directory into
 * FASTQ. Everything upstream of it decides what a cluster IS; this decides what a read looks like
 * once written.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE READ STRUCTURE CUTS THE CYCLES: `4T` is one read of four bases and `2T2T` is two reads of
 *     two, out of the same four cycles and the same clusters;
 *   - A BARCODE SEGMENT IS NOT WRITTEN as a read, so `2T2B` writes two bases per cluster and not
 *     four;
 *   - THE READ NAME CARRIES THE RUN, THE LANE, THE TILE AND THE CLUSTER'S COORDINATES, which is
 *     where the `.locs` file's floats end up;
 *   - `--INCLUDE_NON_PF_READS` DECIDES WHETHER THE FILTERED CLUSTER IS THERE AT ALL, so the same
 *     tile is four reads or three;
 *   - `--READ_NAME_FORMAT` CHOOSES BETWEEN TWO NAMINGS of the same cluster;
 *   - `--MULTIPLEX_PARAMS` SPLITS THE OUTPUT BY BARCODE, one file per declared barcode plus one
 *     for the rest, and it is mutually exclusive with `--OUTPUT_PREFIX`;
 *   - A BARCODE THE PARAMS DO NOT DECLARE IS A REFUSAL unless `--IGNORE_UNEXPECTED_BARCODES` says
 *     otherwise;
 *   - AND THE QUALITIES ARE THE BCL'S SIX BITS, written as ASCII 33 plus the value.
 *
 * Output:
 *
 *     files\t<case>\t<the files written, sorted, space separated>
 *     fastq\t<case>.<file>\t<its text, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BasecallsToFastqDump
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

public class BasecallsToFastqDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A file's text, decompressed where the tool gzipped it. */
    static String text(final Path path) throws Exception {
        final byte[] bytes = Files.readAllBytes(path);
        if (bytes.length > 1 && (bytes[0] & 0xff) == 0x1f && (bytes[1] & 0xff) == 0x8b) {
            try (final java.util.zip.GZIPInputStream in =
                         new java.util.zip.GZIPInputStream(Files.newInputStream(path))) {
                return new String(in.readAllBytes(), StandardCharsets.UTF_8);
            }
        }
        return new String(bytes, StandardCharsets.UTF_8);
    }

    static void run(final String name, final String multiplexParams, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("basecallstofastq");
        final Path basecalls = MakeIlluminaFixtures.write(dir.resolve("run"));
        final Path out = Files.createDirectories(dir.resolve("out"));

        final List<String> argv = new ArrayList<>(List.of(
                "BASECALLS_DIR=" + basecalls, "LANE=1",
                "RUN_BARCODE=run17", "MACHINE_NAME=machine", "FLOWCELL_BARCODE=flowcell",
                "COMPRESS_OUTPUTS=false"));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("READ_STRUCTURE="))) {
            tail.add("READ_STRUCTURE=4T");
        }
        if (multiplexParams == null) {
            tail.add("OUTPUT_PREFIX=" + out.resolve("reads"));
        } else {
            // A multiplexed run reads the barcode each cluster was ASSIGNED, which is the file
            // `ExtractIlluminaBarcodes` writes rather than anything in the basecalls themselves.
            // So the pipeline is the reference's own: extract, then convert.
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
            final Path params = dir.resolve("params.tsv");
            Files.writeString(params,
                    multiplexParams.replace("<out>", out.toString()), StandardCharsets.UTF_8);
            tail.add("MULTIPLEX_PARAMS=" + params);
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = new picard.illumina.IlluminaBasecallsToFastq()
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
        try (final Stream<Path> walk = Files.walk(out)) {
            walk.filter(Files::isRegularFile)
                    .forEach(path -> written.add(path.getFileName().toString()));
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            emit("fastq", name + "." + file, text(out.resolve(file)));
        }
    }

    public static void main(final String[] args) throws Exception {
        // The four cycles as one read, then as two, then with two of them a barcode.
        run("one-read-of-four", null);
        run("two-reads-of-two", null, "READ_STRUCTURE=2T2T");
        run("a-barcode-segment", null, "READ_STRUCTURE=2T2B");
        // A skip is not written either.
        run("a-skipped-segment", null, "READ_STRUCTURE=2T2S");

        // The cluster that failed the filter, in and out.
        run("without-the-non-pf-reads", null, "INCLUDE_NON_PF_READS=false");

        // The two namings of a cluster.
        run("the-other-read-name-format", null, "READ_NAME_FORMAT=ILLUMINA");

        // Split by barcode: the fixture's clusters carry `AG` and `CT` at cycles three and four.
        run("split-by-barcode",
                "OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n<out>/second\tCT\n",
                "READ_STRUCTURE=2T2B");
        // One barcode declared, and the clusters carrying the other one.
        run("an-undeclared-barcode",
                "OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n",
                "READ_STRUCTURE=2T2B");
        run("an-undeclared-barcode-ignored",
                "OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n",
                "READ_STRUCTURE=2T2B", "IGNORE_UNEXPECTED_BARCODES=true");
        // And the params file's own row for everything that matched nothing, which is the row
        // whose barcode column is EMPTY rather than one carrying a placeholder.
        run("a-row-for-the-rest",
                "OUTPUT_PREFIX\tBARCODE_1\n<out>/first\tAG\n<out>/rest\t\n",
                "READ_STRUCTURE=2T2B");

        System.out.print(buf);
    }
}
