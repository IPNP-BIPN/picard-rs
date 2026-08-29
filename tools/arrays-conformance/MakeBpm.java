/*
 * A bead pool manifest, written byte by byte.
 *
 * Picard's genotyping-array tools read Illumina's binary manifests, and the reference's own test
 * data is not in the pinned clone, so the fixture is built here from the format its parser
 * documents by reading it:
 *
 *   - the file opens with `BPM`, a version byte of one, and an int version between three and five;
 *   - then two strings, the manifest's name and its control configuration;
 *   - then the number of loci, an index block of four bytes each which the parser SKIPS, the
 *     names in order, and one normalization id per locus;
 *   - then a locus entry apiece, each opening with its own version between six and eight.
 *
 * Every integer is little-endian and every string is a varint length followed by its bytes, which
 * is what `InfiniumDataFile.parseString` reads.
 *
 * Two rules the parser enforces and this writer therefore obeys: a locus entry's index has to
 * match the position of its name, and `assayType` and `addressB` have to agree, an assay type of
 * zero meaning no B address and anything else meaning one.
 *
 * Usage: MakeBpm <file>
 */

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class MakeBpm {

    /** One locus, reduced to what a manifest carries about it. */
    record Locus(String name, String snp, String chrom, int position, int addressA, int addressB,
                 int assayType, int normalizationId, String ilmnStrand, String sourceStrand,
                 String refStrand) {}

    /** The four loci the fixtures use: two assay types, two chromosomes, two normalization ids. */
    static final List<Locus> LOCI = List.of(
            new Locus("rs1", "[A/G]", "1", 1000, 11, 0, 0, 1, "TOP", "TOP", "+"),
            new Locus("rs2", "[T/C]", "1", 2000, 12, 13, 1, 1, "BOT", "BOT", "-"),
            new Locus("rs3", "[A/C]", "2", 3000, 14, 0, 0, 2, "TOP", "PLUS", "+"),
            new Locus("rs4", "[A/T]", "2", 4000, 15, 16, 2, 2, "PLUS", "TOP", "+"));

    static void writeString(final ByteArrayOutputStream out, final String text) {
        final byte[] bytes = text.getBytes(StandardCharsets.UTF_8);
        // The length is a varint: seven bits a byte, the high bit saying another follows.
        int length = bytes.length;
        while (length >= 0x80) {
            out.write((length & 0x7F) | 0x80);
            length >>= 7;
        }
        out.write(length);
        out.write(bytes, 0, bytes.length);
    }

    static void writeInt(final ByteArrayOutputStream out, final int value) {
        final ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(value);
        out.write(buffer.array(), 0, 4);
    }

    static void writeFloat(final ByteArrayOutputStream out, final float value) {
        final ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putFloat(value);
        out.write(buffer.array(), 0, 4);
    }

    /** One locus entry, at version eight, which is the one that carries a reference strand. */
    static void writeLocus(final ByteArrayOutputStream out, final Locus locus, final int index) {
        writeInt(out, 8);
        writeString(out, locus.name() + "_ilmn");
        writeString(out, locus.name());
        writeString(out, "");
        writeString(out, "");
        writeString(out, "");
        // The index the parser reads is one-based.
        writeInt(out, index + 1);
        writeString(out, "");
        writeString(out, locus.ilmnStrand());
        writeString(out, locus.snp());
        writeString(out, locus.chrom());
        writeString(out, "diploid");
        writeString(out, "Homo sapiens");
        writeString(out, String.valueOf(locus.position()));
        writeString(out, "ACGT");
        writeString(out, "TOP");
        writeInt(out, locus.addressA());
        writeInt(out, locus.addressB());
        writeString(out, "ACGTACGT");
        writeString(out, locus.assayType() == 0 ? "" : "ACGTACGA");
        writeString(out, "37");
        writeString(out, "source");
        writeString(out, "1");
        writeString(out, locus.sourceStrand());
        writeString(out, "ACGTACGTACGT");
        out.write(0);
        out.write(3);
        out.write(0);
        out.write(locus.assayType());
        writeFloat(out, 0.25f);
        writeFloat(out, 0.25f);
        writeFloat(out, 0.25f);
        writeFloat(out, 0.25f);
        writeString(out, locus.refStrand());
    }

    /** The whole file. */
    static byte[] bytes(final String manifestName, final List<Locus> loci) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        out.write('B');
        out.write('P');
        out.write('M');
        out.write(1);
        writeInt(out, 4);
        writeString(out, manifestName);
        writeString(out, "control,config");
        writeInt(out, loci.size());
        // The index block, which the parser skips: four bytes a locus.
        for (int index = 0; index < loci.size(); index++) {
            writeInt(out, index);
        }
        for (final Locus locus : loci) {
            writeString(out, locus.name());
        }
        for (final Locus locus : loci) {
            out.write(locus.normalizationId());
        }
        for (int index = 0; index < loci.size(); index++) {
            writeLocus(out, loci.get(index), index);
        }
        return out.toByteArray();
    }

    static Path write(final Path file, final String manifestName, final List<Locus> loci)
            throws IOException {
        Files.createDirectories(file.getParent());
        Files.write(file, bytes(manifestName, loci));
        return file;
    }

    public static void main(final String[] args) throws Exception {
        write(Paths.get(args[0]), "fixture.bpm", LOCI);
        System.out.println("wrote " + args[0]);
    }
}
