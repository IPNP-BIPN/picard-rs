/*
 * `CompareGtcFiles`, taken from the reference.
 *
 * The tool reads two genotype call files and reports every field that differs between them, which
 * is what a caller runs after re-running a chip. What it compares is every GETTER on the file
 * object, by reflection, so the answer names Java methods rather than file fields.
 *
 * Five behaviours this is built to catch.
 *
 *   - TWO IDENTICAL FILES ARE ZERO and say nothing;
 *   - A DIFFERENT SAMPLE NAME IS NOT A DIFFERENCE: the fields expected to differ between two runs
 *     of one chip are excluded by name;
 *   - A DIFFERENT GENOTYPE IS ONE, and the message counts the elements of the array that differ
 *     rather than naming them;
 *   - A DIFFERENT ARRAY LENGTH IS A DIFFERENT MESSAGE, naming both lengths;
 *   - AND THE COMPARISON NEEDS THE MANIFEST, because normalized intensities are computed from it
 *     rather than stored.
 *
 * Output:
 *
 *     code\t<case>\t<the exit status>
 *     (the differences themselves go through log4j and cannot be captured in-process; the status
 *     is the answer a caller reads)
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CompareGtcDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class CompareGtcDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static void run(final String name, final MakeGtc.Sample left, final MakeGtc.Sample right)
            throws Exception {
        final Path dir = Files.createTempDirectory("comparegtc");
        final Path bpm = MakeBpm.write(dir.resolve("fixture.bpm"), "fixture.bpm", MakeBpm.LOCI);
        // One transformation per UNIQUE normalization id the manifest declares, which is the id
        // plus a hundred times the assay type: the reader indexes the transformations by that.
        final java.util.Set<Integer> unique = new java.util.TreeSet<>();
        for (final MakeBpm.Locus locus : MakeBpm.LOCI) {
            unique.add(locus.normalizationId() + 100 * locus.assayType());
        }
        final int normalizations = unique.size();
        final Path one = MakeGtc.write(dir.resolve("one.gtc"), left, "fixture.egt", "fixture.bpm",
                normalizations);
        final Path two = MakeGtc.write(dir.resolve("two.gtc"), right, "fixture.egt", "fixture.bpm",
                normalizations);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.arrays.illumina.CompareGtcFiles().instanceMain(
                    new String[]{"I=" + one, "I=" + two, "BPM_FILE=" + bpm});
            System.setOut(original);
            System.setErr(originalError);
            emit("code", name, String.valueOf(code));
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
        // The differences themselves are not recorded, and the reason is worth stating: the tool
        // writes them through log4j, which holds the stream it was initialised with, so a dump
        // that redirects `System.err` in-process captures nothing. What IS the tool's answer to a
        // caller is the STATUS, which is what every case here compares.
    }

    public static void main(final String[] args) throws Exception {
        final MakeGtc.Sample base = MakeGtc.fixture("sample1");
        run("two-identical-files", base, base);
        run("a-different-sample-name", base, MakeGtc.fixture("sample2"));
        run("a-different-genotype", base, new MakeGtc.Sample(base.name(),
                List.of(1, 2, 3, 3), base.rawX(), base.rawY(), base.scores(), base.callRate()));
        run("two-different-genotypes", base, new MakeGtc.Sample(base.name(),
                List.of(2, 2, 3, 3), base.rawX(), base.rawY(), base.scores(), base.callRate()));
        run("different-intensities", base, new MakeGtc.Sample(base.name(),
                base.genotypes(), List.of(1000, 2000, 3000, 9999), base.rawY(), base.scores(),
                base.callRate()));
        run("a-different-call-rate", base, new MakeGtc.Sample(base.name(),
                base.genotypes(), base.rawX(), base.rawY(), base.scores(), 0.5f));
        run("a-shorter-file", base, new MakeGtc.Sample(base.name(),
                List.of(1, 2, 3), List.of(1000, 2000, 3000), List.of(1100, 2100, 3100),
                List.of(0.7f, 0.8f, 0.9f), base.callRate()));

        System.out.print(buf);
    }
}
