/*
 * A cluster file, written byte by byte.
 *
 * The `.egt` is Illumina's cluster definition: for every locus in a manifest, where the three
 * genotypes sit in intensity space and how tight each cluster is. The genotyping-array tools read
 * it beside the manifest, and the reference's own test data is not in the pinned clone, so it is
 * written here from the format `InfiniumEGTFile` documents by reading it:
 *
 *   - a header of an int version, five strings, a byte, and the manifest's name;
 *   - then an int gentrain version, the manifest's name AGAIN, and the number of loci;
 *   - then per locus: three counts, then four triples of floats, then fifteen floats nobody reads;
 *   - then per locus four more numbers of which only the second is kept;
 *   - then two strings apiece, the second of which is the locus's name.
 *
 * The names have to be UNIQUE: the parser refuses a repeat by name.
 *
 * Usage: MakeEgt <file>
 */

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class MakeEgt {

    static void writeString(final ByteArrayOutputStream out, final String text) {
        MakeBpm.writeString(out, text);
    }

    static void writeInt(final ByteArrayOutputStream out, final int value) {
        MakeBpm.writeInt(out, value);
    }

    static void writeFloat(final ByteArrayOutputStream out, final float value) {
        MakeBpm.writeFloat(out, value);
    }

    /** The whole file, over the locus names given. */
    static byte[] bytes(final String manifestName, final List<String> names) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        // The header.
        writeInt(out, 3);
        writeString(out, "gencall-1");
        writeString(out, "cluster-1");
        writeString(out, "call-1");
        writeString(out, "normalization-1");
        writeString(out, "2020-01-01");
        out.write(1);
        writeString(out, manifestName);

        // The data.
        writeInt(out, 3);
        writeString(out, manifestName);
        writeInt(out, names.size());
        for (int index = 0; index < names.size(); index++) {
            // The three genotype counts, then the deviations and means of R and theta.
            for (final int count : new int[]{10 + index, 20 + index, 30 + index}) {
                writeInt(out, count);
            }
            for (final float value : new float[]{0.1f, 0.2f, 0.3f}) {
                writeFloat(out, value);
            }
            for (final float value : new float[]{1.0f, 1.1f, 1.2f}) {
                writeFloat(out, value);
            }
            for (final float value : new float[]{0.01f, 0.02f, 0.03f}) {
                writeFloat(out, value);
            }
            for (final float value : new float[]{0.2f, 0.5f, 0.8f}) {
                writeFloat(out, value);
            }
            // Fifteen floats nobody reads.
            for (int unused = 0; unused < 15; unused++) {
                writeFloat(out, 0f);
            }
        }
        for (int index = 0; index < names.size(); index++) {
            writeFloat(out, 0f);
            // The only one of the four that is kept: the cluster's total score.
            writeFloat(out, 0.5f + index / 100f);
            writeFloat(out, 0f);
            out.write(0);
        }
        for (final String name : names) {
            writeString(out, name + "_address");
        }
        for (final String name : names) {
            writeString(out, name);
        }
        return out.toByteArray();
    }

    static Path write(final Path file, final String manifestName, final List<String> names)
            throws IOException {
        Files.createDirectories(file.getParent());
        Files.write(file, bytes(manifestName, names));
        return file;
    }

    public static void main(final String[] args) throws Exception {
        write(Paths.get(args[0]), "fixture.bpm", List.of("rs1", "rs2", "rs3", "rs4"));
        System.out.println("wrote " + args[0]);
    }
}
