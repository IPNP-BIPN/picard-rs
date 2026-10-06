import java.io.ByteArrayOutputStream;
import java.io.File;
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

/**
 * The genotyping-array corpus: a GRCh37-style reference, a bead pool manifest, a cluster file, an
 * extended manifest, genotype call files, the VCFs the reference's own GtcToVcf makes of them, and
 * the text the two external tools whose output Picard parses would print.
 *
 * The binary writers below are the ones tools/arrays-conformance uses (MakeBpm, MakeEgt, MakeGtc,
 * MakeExtendedManifest), renamed so the two copies cannot meet on one classpath: the corpus is
 * compiled from this directory alone.
 */
public class MakeArraysFixtures {

    static final int LENGTH = 5000;

    static String base(final int position) {
        return String.valueOf("ACGT".charAt((position - 1) % 4));
    }

    static String fasta() {
        final StringBuilder out = new StringBuilder();
        for (final String contig : new String[]{"1", "2"}) {
            out.append('>').append(contig).append('\n');
            final StringBuilder bases = new StringBuilder();
            for (int index = 0; index < LENGTH; index++) {
                bases.append("ACGT".charAt(index % 4));
            }
            for (int index = 0; index < bases.length(); index += 60) {
                out.append(bases, index, Math.min(index + 60, bases.length())).append('\n');
            }
        }
        return out.toString();
    }

    /** One run of the reference's GtcToVcf, as the array VCFs downstream tools read. */
    static void gtcToVcf(final Path gtc, final Path out, final Path reference, final Path manifest,
                         final Path egt, final Path bpm, final String sample, final String gender)
            throws Exception {
        final int code = new picard.arrays.GtcToVcf().instanceMain(new String[]{
                "INPUT=" + gtc, "OUTPUT=" + out, "R=" + reference,
                "EXTENDED_ILLUMINA_MANIFEST=" + manifest, "CLUSTER_FILE=" + egt,
                "ILLUMINA_BEAD_POOL_MANIFEST_FILE=" + bpm, "SAMPLE_ALIAS=" + sample,
                "ANALYSIS_VERSION_NUMBER=1", "PIPELINE_VERSION=1.0", "EXPECTED_GENDER=" + gender});
        if (code != 0) {
            throw new IllegalStateException("GtcToVcf refused the fixture " + gtc);
        }
    }

    static void write(final File into) throws Exception {
        final Path dir = into.toPath();
        final Path reference = dir.resolve("arr_ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
                "R=" + reference, "O=" + dir.resolve("arr_ref.dict"),
                "GENOME_ASSEMBLY=GRCh37", "SPECIES=Homo sapiens"});
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(reference, true);

        // The files are named arr.*, but the names the files carry inside are the ones the extended
        // manifest declares (fixture.bpm), which GtcToVcf checks the call files against.
        // Four loci on two contigs, two assay types, called AA, AB, BB and no-call in sample1.
        final List<ArraysManifest.Row> rows = List.of(
                new ArraysManifest.Row("rs1", "[A/G]", "1", 1001, 11, 0, "+", base(1001), "A", "G", "PASS"),
                new ArraysManifest.Row("rs2", "[T/C]", "1", 2001, 12, 13, "+", base(2001), "T", "C", "PASS"),
                new ArraysManifest.Row("rs3", "[A/C]", "2", 3001, 14, 0, "+", base(3001), "A", "C", "PASS"),
                new ArraysManifest.Row("rs4", "[A/T]", "2", 4001, 15, 16, "+", base(4001), "A", "T", "PASS"));
        final List<ArraysBpm.Locus> loci = new ArrayList<>();
        final List<String> names = new ArrayList<>();
        final List<ArraysManifest.Row> plain = new ArrayList<>();
        for (final ArraysManifest.Row row : rows) {
            loci.add(new ArraysBpm.Locus(row.name(), row.snp(), row.chrom(), row.position(),
                    row.addressA(), row.addressB(), row.addressB() == 0 ? 0 : 1, 1,
                    "TOP", "TOP", row.refStrand()));
            names.add(row.name());
            final String[] alleles = row.snp().replace("[", "").replace("]", "").split("/");
            plain.add(new ArraysManifest.Row(row.name(), row.snp(), row.chrom(), row.position(),
                    row.addressA(), row.addressB(), row.refStrand(), "", alleles[0], alleles[1], ""));
        }
        final Path bpm = ArraysBpm.write(dir.resolve("arr.bpm"), "fixture.bpm", loci);
        final Path egt = ArraysEgt.write(dir.resolve("arr.egt"), "fixture.bpm", names);
        final Path manifest = ArraysManifest.write(dir.resolve("arr_extended.csv"), rows);
        ArraysManifest.writePlain(dir.resolve("arr_manifest.csv"), plain);
        final java.util.Set<Integer> unique = new java.util.TreeSet<>();
        for (final ArraysBpm.Locus locus : loci) {
            unique.add(locus.normalizationId() + 100 * locus.assayType());
        }

        // dbSNP over two of the four loci, indexed.
        final StringBuilder sites = new StringBuilder("##fileformat=VCFv4.2\n");
        sites.append("##contig=<ID=1,length=").append(LENGTH).append(">\n");
        sites.append("##contig=<ID=2,length=").append(LENGTH).append(">\n");
        sites.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (final int[] site : new int[][]{{1, 1001}, {2, 3001}}) {
            sites.append(site[0]).append('\t').append(site[1]).append("\trs").append(site[1])
                    .append('\t').append(base(site[1])).append("\tG\t100\tPASS\t.\n");
        }
        final Path dbsnp = dir.resolve("arr_dbsnp.vcf");
        Files.writeString(dbsnp, sites.toString(), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.writeIndex(
                htsjdk.tribble.index.IndexFactory.createLinearIndex(
                        dbsnp.toFile(), new htsjdk.variant.vcf.VCFCodec()),
                new File(dbsnp + ".idx"));

        // Three call files: sample1 as the conformance fixture has it, sample1 again with one
        // genotype and one intensity changed, and sample2 with other calls throughout.
        final ArraysGtc.Sample s1 = ArraysGtc.fixture("sample1");
        final ArraysGtc.Sample s1b = new ArraysGtc.Sample("sample1", List.of(1, 3, 3, 0),
                List.of(1000, 2000, 3500, 4000), List.of(1100, 2100, 3100, 4100),
                List.of(0.7f, 0.8f, 0.9f, 0.0f), 0.75f);
        final ArraysGtc.Sample s2 = new ArraysGtc.Sample("sample2", List.of(3, 2, 1, 1),
                List.of(900, 2200, 2800, 4300), List.of(1200, 2000, 3300, 3900),
                List.of(0.6f, 0.85f, 0.95f, 0.5f), 1.0f);
        final Path g1 = ArraysGtc.write(dir.resolve("arr_s1.gtc"), s1, "fixture.egt", "fixture.bpm", unique.size());
        ArraysGtc.write(dir.resolve("arr_s1b.gtc"), s1b, "fixture.egt", "fixture.bpm", unique.size());
        final Path g2 = ArraysGtc.write(dir.resolve("arr_s2.gtc"), s2, "fixture.egt", "fixture.bpm", unique.size());

        gtcToVcf(g1, dir.resolve("arr_s1.vcf"), reference, manifest, egt, bpm, "sample1", "Female");
        gtcToVcf(g2, dir.resolve("arr_s2.vcf"), reference, manifest, egt, bpm, "sample2", "Male");

        // What bafRegress and VerifyIDIntensity print, in the shapes Picard's parsers read.
        Files.writeString(dir.resolve("arr_bafregress.txt"),
                "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom\n"
                        + "sample1\t-0.000538639821316378\t0.00011230471350585\t-4.79623521134149"
                        + "\t1.61690798727399e-06\t0.995028901600268\t1510547\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("arr_bafregress_two.txt"),
                "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom\n"
                        + "sample1\t0.0125\t0.0031\t4.03\t5.6e-05\t0.981\t1200\n"
                        + "sample2\t-0.002\t0.0009\t-2.2\t0.0278\t0.9995\t998\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("arr_verifyidintensity_four.txt"),
                "ID\t%Mix\t\tLLK\t\tLLK0\t\n"
                        + "-----------------------------------------------------------------\n"
                        + "0\t0.214766\t157575\t177169\n"
                        + "1\t0.214767\t157576\t177170\n"
                        + "2\t0.0994769\t90260.4\t91166.7\n"
                        + "3\t0.234567\t-4703.97\t-5204.97\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("arr_verifyidintensity.txt"),
                "ID\t%Mix\t\tLLK\t\tLLK0\t\n"
                        + "-----------------------------------------------------------------\n"
                        + "0\t0.214766\t157575\t177169\n"
                        + "1\t0.0994769\t90260.4\t91166.7\n",
                StandardCharsets.UTF_8);
    }
}

final class ArraysBpm {

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
}

final class ArraysEgt {

    static void writeString(final ByteArrayOutputStream out, final String text) {
        ArraysBpm.writeString(out, text);
    }

    static void writeInt(final ByteArrayOutputStream out, final int value) {
        ArraysBpm.writeInt(out, value);
    }

    static void writeFloat(final ByteArrayOutputStream out, final float value) {
        ArraysBpm.writeFloat(out, value);
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
}

final class ArraysGtc {

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
        ArraysBpm.writeString(out, text);
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
}

final class ArraysManifest {

    /** The columns, in the order the file writes them. */
    static final List<String> COLUMNS = List.of(
            "IlmnID", "Name", "IlmnStrand", "SNP", "AddressA_ID", "AlleleA_ProbeSeq",
            "AddressB_ID", "AlleleB_ProbeSeq", "GenomeBuild", "Chr", "MapInfo", "Ploidy",
            "Species", "Source", "SourceVersion", "SourceStrand", "SourceSeq", "TopGenomicSeq",
            "BeadSetID", "Exp_Clusters", "RefStrand", "Intensity_Only",
            "build37Chr", "build37Pos", "build37RefAllele", "build37AlleleA", "build37AlleleB",
            "build37Rsid", "build37Flag");

    /** One locus's row, with the build-37 columns the extension adds. */
    record Row(String name, String snp, String chrom, int position, int addressA, int addressB,
               String refStrand, String refAllele, String alleleA, String alleleB, String flag) {}

    static String row(final Row row) {
        final List<String> values = new ArrayList<>(List.of(
                row.name() + "_ilmn", row.name(), "TOP", row.snp(),
                String.valueOf(row.addressA()), "ACGTACGT",
                row.addressB() == 0 ? "" : String.valueOf(row.addressB()),
                row.addressB() == 0 ? "" : "ACGTACGA",
                "37", row.chrom(), String.valueOf(row.position()), "diploid", "Homo sapiens",
                "source", "1", "TOP", "ACGTACGTACGT", "ACGT", "1", "3", row.refStrand(), "0",
                row.chrom(), String.valueOf(row.position()), row.refAllele(), row.alleleA(),
                row.alleleB(), row.name(), row.flag()));
        return String.join(",", values);
    }

    /** The plain manifest: the same rows without the seven columns the extension adds. */
    static String plainText(final List<Row> rows) {
        final int columns = COLUMNS.size() - 7;
        final StringBuilder out = new StringBuilder();
        out.append("Illumina, Inc.\n");
        out.append("[Heading]\n");
        out.append("Descriptor File Name,fixture.bpm\n");
        out.append("Assay Format,Infinium HTS\n");
        out.append("Date Manufactured,1/1/2020\n");
        out.append("Loci Count ,").append(rows.size()).append('\n');
        out.append("[Assay]\n");
        out.append(String.join(",", COLUMNS.subList(0, columns))).append('\n');
        for (final Row row : rows) {
            final String[] values = row(row).split(",", -1);
            out.append(String.join(",", java.util.Arrays.asList(values).subList(0, columns)))
                    .append('\n');
        }
        out.append("[Controls]\n");
        return out.toString();
    }

    static Path writePlain(final Path file, final List<Row> rows) throws IOException {
        Files.createDirectories(file.getParent());
        Files.writeString(file, plainText(rows), StandardCharsets.UTF_8);
        return file;
    }

    static String text(final List<Row> rows) {
        final StringBuilder out = new StringBuilder();
        out.append("Illumina, Inc.\n");
        out.append("[Heading]\n");
        out.append("Descriptor File Name,fixture.bpm\n");
        out.append("Assay Format,Infinium HTS\n");
        out.append("Date Manufactured,1/1/2020\n");
        out.append("Loci Count ,").append(rows.size()).append('\n');
        out.append("[Assay]\n");
        out.append(String.join(",", COLUMNS)).append('\n');
        for (final Row row : rows) {
            out.append(row(row)).append('\n');
        }
        out.append("[Controls]\n");
        return out.toString();
    }

    static Path write(final Path file, final List<Row> rows) throws IOException {
        Files.createDirectories(file.getParent());
        Files.writeString(file, text(rows), StandardCharsets.UTF_8);
        return file;
    }
}
