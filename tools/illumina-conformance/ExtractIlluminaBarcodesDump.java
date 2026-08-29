/*
 * `ExtractIlluminaBarcodes`, taken from the reference.
 *
 * The tool reads the barcode cycles of every cluster in a tile and decides which declared barcode
 * each one is, if any. What it writes is a `_barcode.txt` per tile, one line per cluster, and a
 * metrics file counting how many clusters each barcode took.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE MATCH IS BY EDIT DISTANCE, not by equality: `--MAX_MISMATCHES` decides how far a cluster
 *     may be from a declared barcode and still be it;
 *   - `--MIN_MISMATCH_DELTA` IS THE SECOND TEST: a cluster near TWO barcodes is matched to neither
 *     unless the better match wins by that margin, so two barcodes one base apart make every
 *     cluster ambiguous;
 *   - A CLUSTER THAT MATCHES NOTHING IS STILL WRITTEN, with an `N` verdict rather than a line
 *     missing;
 *   - THE METRICS COUNT PF AND NON-PF SEPARATELY, so the filter file decides two of the columns;
 *   - `--BARCODE` AND `--BARCODE_FILE` ARE TWO WAYS TO SAY THE SAME THING, and the second carries
 *     a name and a library the first cannot;
 *   - A BARCODE OF THE WRONG LENGTH IS NOT REFUSED BUT CRASHES: the read structure decides how
 *     many cycles a barcode has, and a longer one walks off the end of the cluster's bases with an
 *     `ArrayIndexOutOfBoundsException`. That is the reference's behaviour and not a tidy refusal,
 *     which is exactly why it is measured;
 *   - `--MINIMUM_BASE_QUALITY` DROPS A CLUSTER whose barcode bases are too poor, which is a
 *     different rejection from a mismatch;
 *   - AND `--OUTPUT_DIR` DECIDES WHERE THE PER-TILE FILES GO, the basecalls directory being the
 *     default and not a requirement.
 *
 * Output:
 *
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     barcodes\t<case>\t<the per-tile file, escaped>
 *     files\t<case>\t<the names written beside the basecalls, sorted, space separated>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: ExtractIlluminaBarcodesDump
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

public class ExtractIlluminaBarcodesDump {

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

    /** A barcode file: the sequence, a name and a library. */
    static String barcodeFile(final List<String[]> rows) {
        final StringBuilder text = new StringBuilder("barcode_sequence_1\tbarcode_name\tlibrary_name\n");
        for (final String[] row : rows) {
            text.append(String.join("\t", row)).append('\n');
        }
        return text.toString();
    }

    static void run(final String name, final List<String[]> barcodes, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("extractbarcodes");
        final Path basecalls = MakeIlluminaFixtures.write(dir.resolve("run"));
        final Path metrics = dir.resolve("metrics.txt");
        final Path out = Files.createDirectories(dir.resolve("out"));

        final List<String> argv = new ArrayList<>(List.of(
                "BASECALLS_DIR=" + basecalls, "LANE=1", "METRICS_FILE=" + metrics,
                "OUTPUT_DIR=" + out));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("READ_STRUCTURE="))) {
            tail.add("READ_STRUCTURE=2T2B");
        }
        if (!barcodes.isEmpty()) {
            final Path file = dir.resolve("barcodes.tsv");
            Files.writeString(file, barcodeFile(barcodes), StandardCharsets.UTF_8);
            tail.add("BARCODE_FILE=" + file);
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = new picard.illumina.ExtractIlluminaBarcodes()
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

        emit("metrics", name, table(metrics));
        final List<String> written = new ArrayList<>();
        try (final Stream<Path> walk = Files.walk(out)) {
            walk.filter(Files::isRegularFile).forEach(path -> written.add(path.getFileName().toString()));
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            // The per-tile file is gzipped where the tool chose to compress it, and plain
            // otherwise; what is compared is the text either way.
            final Path path = out.resolve(file);
            final byte[] bytes = Files.readAllBytes(path);
            final String text;
            if (bytes.length > 1 && (bytes[0] & 0xff) == 0x1f && (bytes[1] & 0xff) == 0x8b) {
                try (final java.util.zip.GZIPInputStream in =
                             new java.util.zip.GZIPInputStream(Files.newInputStream(path))) {
                    text = new String(in.readAllBytes(), StandardCharsets.UTF_8);
                }
            } else {
                text = new String(bytes, StandardCharsets.UTF_8);
            }
            emit("barcodes", name + "." + file, text);
        }
    }

    public static void main(final String[] args) throws Exception {
        // The fixture's four clusters carry `AG`, `AG`, `CT` and `CT` in cycles three and four.
        final List<String[]> two = List.<String[]>of(
                new String[]{"AG", "first", "libraryA"},
                new String[]{"CT", "second", "libraryB"});
        run("two-barcodes", two);
        // One declared barcode, so half the clusters match nothing.
        run("one-barcode", List.<String[]>of(new String[]{"AG", "first", "libraryA"}));
        // A barcode nothing carries.
        run("a-barcode-nobody-has", List.<String[]>of(new String[]{"TT", "absent", "libraryC"}));

        // The two thresholds, each side of their cut.
        run("one-mismatch-allowed", List.<String[]>of(new String[]{"AT", "near", "libraryA"}));
        run("no-mismatch-allowed", List.<String[]>of(new String[]{"AT", "near", "libraryA"}),
                "MAX_MISMATCHES=0");
        // Two barcodes EQUIDISTANT from the clusters, which is what the delta is for: `AG` is one
        // base from `AA` and one from `GG`, so the better match wins by nothing.
        run("two-equidistant-barcodes", List.<String[]>of(
                new String[]{"AA", "first", "libraryA"},
                new String[]{"GG", "second", "libraryB"}));
        run("two-equidistant-barcodes-with-no-delta", List.<String[]>of(
                new String[]{"AA", "first", "libraryA"},
                new String[]{"GG", "second", "libraryB"}), "MIN_MISMATCH_DELTA=0");

        // The quality floor, which rejects rather than mismatches.
        run("a-quality-floor-below-the-bases", two, "MINIMUM_BASE_QUALITY=10");
        run("a-quality-floor-above-the-bases", two, "MINIMUM_BASE_QUALITY=40");

        // A barcode whose length is not the read structure's.
        run("a-barcode-of-the-wrong-length", List.<String[]>of(new String[]{"AGT", "long", "libraryA"}));
        // And the argument that names barcodes without a file.
        run("barcodes-without-a-file", List.of(), "BARCODE=AG", "BARCODE=CT");

        System.out.print(buf);
    }
}
