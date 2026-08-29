/*
 * `BpmToNormalizationManifestCsv`, taken from the reference.
 *
 * The tool reads Illumina's binary bead pool manifest and writes the part of it a normalization
 * needs: one row per locus, with the address, the normalization id and the assay type. It is the
 * smallest of the genotyping-array tools and the one that says whether a manifest parses at all.
 *
 * Five behaviours this is built to catch.
 *
 *   - THE ROWS ARE THE LOCI, in the manifest's own order rather than sorted;
 *   - THE NORMALIZATION ID IS NOT THE ONE IN THE FILE: the parser adds a hundred times the assay
 *     type to it, so a locus with id 1 and assay type 2 reports 201;
 *   - AN ASSAY TYPE AND A B ADDRESS HAVE TO AGREE, and a manifest where they do not is refused by
 *     name;
 *   - A NORMALIZATION ID ABOVE A HUNDRED IS REFUSED, which is what the addition above would
 *     otherwise hide;
 *   - AND THE HEADER CARRIES THE MANIFEST'S OWN NAME, taken from the file rather than from the
 *     command line.
 *
 * Output:
 *
 *     csv\t<case>\t<the written file, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BpmToCsvDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BpmToCsvDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static void run(final String name, final List<MakeBpm.Locus> loci) throws Exception {
        final Path dir = Files.createTempDirectory("bpm");
        final Path bpm = MakeBpm.write(dir.resolve("fixture.bpm"), "fixture.bpm", loci);
        // The tool reads a CLUSTER file beside the manifest, which is a different format again.
        final List<String> names = new ArrayList<>();
        for (final MakeBpm.Locus locus : loci) {
            names.add(locus.name());
        }
        final Path egt = MakeEgt.write(dir.resolve("fixture.egt"), "fixture.bpm", names);
        final Path out = dir.resolve("normalization.csv");

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.arrays.illumina.BpmToNormalizationManifestCsv()
                    .instanceMain(new String[]{"I=" + bpm, "O=" + out, "CLUSTER_FILE=" + egt});
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
        emit("csv", name, Files.readString(out, StandardCharsets.UTF_8));
    }

    public static void main(final String[] args) throws Exception {
        run("four-loci", MakeBpm.LOCI);
        // One locus, of each assay type.
        run("one-locus-of-assay-type-zero", List.of(MakeBpm.LOCI.get(0)));
        run("one-locus-of-assay-type-one", List.of(MakeBpm.LOCI.get(1)));

        // An assay type that disagrees with its B address, each way round.
        run("assay-type-zero-with-a-b-address", List.of(new MakeBpm.Locus(
                "rs1", "[A/G]", "1", 1000, 11, 12, 0, 1, "TOP", "TOP", "+")));
        run("assay-type-one-without-a-b-address", List.of(new MakeBpm.Locus(
                "rs1", "[A/G]", "1", 1000, 11, 0, 1, 1, "TOP", "TOP", "+")));
        // A normalization id the parser refuses.
        run("a-normalization-id-above-a-hundred", List.of(new MakeBpm.Locus(
                "rs1", "[A/G]", "1", 1000, 11, 0, 0, 101, "TOP", "TOP", "+")));

        System.out.print(buf);
    }
}
