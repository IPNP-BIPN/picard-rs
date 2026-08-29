/*
 * A genotype call file, written byte by byte.
 *
 * The `.gtc` is what Illumina's caller writes for one sample on one chip. Unlike the manifest and
 * the cluster file it is not a stream to be read in order: it opens with a TABLE OF CONTENTS, one
 * entry per kind of data, each an id and an ABSOLUTE offset into the file, and the reader seeks.
 * So this writer lays the payloads out first, records where each landed, and writes the table over
 * the top.
 *
 * The ids are the reader's own: 1 is the number of SNPs, 10 the sample's name, 100 and 101 the
 * cluster file and the manifest it was called against, 400 the normalization transformations,
 * 1000 and 1001 the raw intensities, 1002 the genotypes, 1003 the base calls, 1004 the scores and
 * 1006 the call rate.
 *
 * Usage: MakeGtc <file>
 */

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class MakeGtc {

    /** What one sample's file says, reduced to what the reader looks at. */
    record Sample(String name, List<Integer> genotypes, List<Integer> rawX, List<Integer> rawY,
                  List<Float> scores, float callRate) {}

    /** The four loci of the shared manifest, called AA, AB, BB and no-call. */
    static Sample fixture(final String name) {
        return new Sample(name, List.of(1, 2, 3, 0), List.of(1000, 2000, 3000, 4000),
                List.of(1100, 2100, 3100, 4100), List.of(0.7f, 0.8f, 0.9f, 0.0f), 0.75f);
    }

    static void writeInt(final ByteArrayOutputStream out, final int value) {
        final ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(value);
        out.write(buffer.array(), 0, 4);
    }

    static void writeShort(final ByteArrayOutputStream out, final int value) {
        final ByteBuffer buffer = ByteBuffer.allocate(2).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putShort((short) value);
        out.write(buffer.array(), 0, 2);
    }

    static void writeFloat(final ByteArrayOutputStream out, final float value) {
        final ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putFloat(value);
        out.write(buffer.array(), 0, 4);
    }

    static byte[] string(final String text) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        MakeBpm.writeString(out, text);
        return out.toByteArray();
    }

    static byte[] unsignedShorts(final List<Integer> values) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        writeInt(out, values.size());
        for (final int value : values) {
            writeShort(out, value);
        }
        return out.toByteArray();
    }

    static byte[] floats(final List<Float> values) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        writeInt(out, values.size());
        for (final float value : values) {
            writeFloat(out, value);
        }
        return out.toByteArray();
    }

    /** The genotypes, one byte apiece: 0 is a no-call and 1, 2 and 3 are AA, AB and BB. */
    static byte[] genotypes(final List<Integer> calls) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        writeInt(out, calls.size());
        for (final int call : calls) {
            out.write(call);
        }
        return out.toByteArray();
    }

    /** The base calls, two bytes apiece; a zero byte is read back as a dash. */
    static byte[] baseCalls(final List<Integer> calls) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        writeInt(out, calls.size());
        for (final int call : calls) {
            switch (call) {
                case 1 -> { out.write('A'); out.write('A'); }
                case 2 -> { out.write('A'); out.write('B'); }
                case 3 -> { out.write('B'); out.write('B'); }
                default -> { out.write(0); out.write(0); }
            }
        }
        return out.toByteArray();
    }

    /** One normalization transformation per normalization id the manifest declares. */
    static byte[] transformations(final int count) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        writeInt(out, count);
        for (int index = 0; index < count; index++) {
            writeInt(out, 1);
            writeFloat(out, 10f);
            writeFloat(out, 20f);
            writeFloat(out, 1f);
            writeFloat(out, 1f);
            writeFloat(out, 0f);
            writeFloat(out, 0f);
            for (int reserved = 0; reserved < 6; reserved++) {
                writeFloat(out, 0f);
            }
        }
        return out.toByteArray();
    }

    /** Three bare unsigned shorts, with no length in front of them. */
    static byte[] shorts(final List<Integer> values) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (final int value : values) {
            writeShort(out, value);
        }
        return out.toByteArray();
    }

    /** A B allele frequency per locus, which a homozygous call puts at nought or one. */
    static List<Float> bAlleleFreqs(final Sample sample) {
        final List<Float> values = new java.util.ArrayList<>();
        for (final int call : sample.genotypes()) {
            values.add(call == 1 ? 0.0f : call == 2 ? 0.5f : call == 3 ? 1.0f : Float.NaN);
        }
        return values;
    }

    /** A log R ratio per locus, which says how much signal there was. */
    static List<Float> logRRatios(final Sample sample) {
        final List<Float> values = new java.util.ArrayList<>();
        for (int index = 0; index < sample.genotypes().size(); index++) {
            values.add(0.1f * index);
        }
        return values;
    }

    static byte[] callRate(final float rate) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        writeFloat(out, rate);
        return out.toByteArray();
    }

    /**
     * The whole file: a header, a table of contents, and the payloads it points at.
     *
     * The table's entries are six bytes each, so the payloads start after the header and the
     * table, and every offset is absolute.
     */
    static byte[] bytes(final Sample sample, final String clusterFile, final String manifest,
                        final int normalizations) {
        final Map<Integer, byte[]> payloads = new LinkedHashMap<>();
        payloads.put(10, string(sample.name()));
        payloads.put(100, string(clusterFile));
        payloads.put(101, string(manifest));
        payloads.put(400, transformations(normalizations));
        // The control intensities: `GtcToVcf` reads their length unconditionally, so a file
        // without them makes it throw rather than report.
        payloads.put(500, unsignedShorts(List.of(10, 20, 30)));
        payloads.put(501, unsignedShorts(List.of(11, 21, 31)));
        payloads.put(1000, unsignedShorts(sample.rawX()));
        payloads.put(1001, unsignedShorts(sample.rawY()));
        payloads.put(1002, genotypes(sample.genotypes()));
        payloads.put(1003, baseCalls(sample.genotypes()));
        payloads.put(1004, floats(sample.scores()));
        payloads.put(1006, callRate(sample.callRate()));
        // The intensity percentiles are three unsigned shorts apiece, and the comparison reads
        // them unconditionally: a file without them makes the tool throw rather than report.
        // The per-locus arrays a VCF's FORMAT fields are built from.
        payloads.put(1012, floats(bAlleleFreqs(sample)));
        payloads.put(1013, floats(logRRatios(sample)));
        payloads.put(1014, shorts(List.of(100, 500, 900)));
        payloads.put(1015, shorts(List.of(110, 510, 910)));

        // The number of SNPs is not a payload at all: its OFFSET is the value, which is what the
        // reader means by `numberOfSnps = toc.getOffset()`.
        final List<Integer> ids = new ArrayList<>(List.of(1));
        ids.addAll(payloads.keySet());

        final int headerSize = 3 + 1 + 4;
        int offset = headerSize + ids.size() * 6;
        final Map<Integer, Integer> offsets = new LinkedHashMap<>();
        offsets.put(1, sample.genotypes().size());
        for (final Map.Entry<Integer, byte[]> entry : payloads.entrySet()) {
            offsets.put(entry.getKey(), offset);
            offset += entry.getValue().length;
        }

        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        out.write('g');
        out.write('t');
        out.write('c');
        out.write(5);
        writeInt(out, ids.size());
        for (final int id : ids) {
            writeShort(out, id);
            writeInt(out, offsets.get(id));
        }
        for (final byte[] payload : payloads.values()) {
            out.write(payload, 0, payload.length);
        }
        return out.toByteArray();
    }

    static Path write(final Path file, final Sample sample, final String clusterFile,
                      final String manifest, final int normalizations) throws IOException {
        Files.createDirectories(file.getParent());
        Files.write(file, bytes(sample, clusterFile, manifest, normalizations));
        return file;
    }

    public static void main(final String[] args) throws Exception {
        write(Paths.get(args[0]), fixture("sample1"), "fixture.egt", "fixture.bpm", 3);
        System.out.println("wrote " + args[0]);
    }
}
