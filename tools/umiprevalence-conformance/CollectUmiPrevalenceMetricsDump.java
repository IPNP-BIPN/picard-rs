/*
 * CollectUmiPrevalenceMetrics' histogram, taken from the reference.
 *
 * The tool groups the reads into duplicate sets and counts how many distinct UMIs each set holds.
 * What is measured is which reads reach a set at all, which is five filters and an optional
 * sixth, and how the UMIs of a set are counted.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE HISTOGRAM IS SETS BY UMI COUNT, so three sets of one UMI each is `1 -> 3` and not
 *     `3 -> 1`;
 *   - THE UMIS OF A SET ARE A SET, so two reads carrying the same tag are one UMI;
 *   - A READ WITH NO UMI TAG IS FILTERED OUT before any set is formed;
 *   - THE BARCODE-QUALITY FILTER IS INVERTED: it keeps the reads whose barcode has a base UNDER
 *     the floor and drops the ones whose barcode is entirely above it, which is the opposite of
 *     what its name and its argument's documentation say. A file of well-formed barcodes
 *     therefore reports NOTHING, and lowering the floor makes the tool report less rather than
 *     more. Every read of this fixture carries one bad barcode base so the other nine behaviours
 *     can be seen past it;
 *   - A READ WITH NO BARCODE-QUALITY TAG IS KEPT by that same filter, the absent tag being the
 *     only way past it;
 *   - AN UNALIGNED READ IS FILTERED, and so is one under --MINIMUM_MQ;
 *   - A SECONDARY OR SUPPLEMENTARY READ IS FILTERED;
 *   - --FILTER_UNPAIRED_READS ADDS A SIXTH FILTER and is on by default, so an unpaired read
 *     reaches no set unless it is turned off;
 *   - THE TAGS ARE NAMED BY ARGUMENT, so --BARCODE_TAG and --BARCODE_BQ move which tags are read;
 *   - AND A FILE WHOSE EVERY READ IS FILTERED WRITES AN EMPTY HISTOGRAM rather than none.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     histogram\t<case>\t<the histogram section, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectUmiPrevalenceMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One read: where it sits, its flags, its mapping quality and its two barcode tags. */
    record Read(String name, int start, int flags, int mappingQuality, String barcode,
                String barcodeQuality) {}

    /**
     * A read that gets past every filter.
     *
     * Its barcode quality carries one base under the floor ON PURPOSE: the barcode-quality filter
     * is inverted, and a barcode entirely above the floor is the one it drops. A fixture of
     * well-formed barcodes produces an empty histogram in every case, which is what the first
     * version of this dump found.
     */
    static Read paired(final String name, final int start, final String barcode) {
        return new Read(name, start, 0x1 | 0x2 | 0x40, 60, barcode, "II#III");
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                if ((spec.flags() & 0x4) != 0) {
                    record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                    record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    record.setMappingQuality(0);
                    record.setCigarString("*");
                } else {
                    record.setReferenceName("chr1");
                    record.setAlignmentStart(spec.start());
                    record.setMappingQuality(spec.mappingQuality());
                    record.setCigarString("10M");
                }
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                if ((spec.flags() & 0x1) != 0) {
                    record.setMateReferenceName("chr1");
                    record.setMateAlignmentStart(spec.start() + 100);
                    // The duplicate-set iterator wants the mate's CIGAR, and refuses by name
                    // without it.
                    record.setAttribute("MC", "10M");
                }
                if (spec.barcode() != null) {
                    record.setAttribute("RX", spec.barcode());
                }
                if (spec.barcodeQuality() != null) {
                    record.setAttribute("BQ", spec.barcodeQuality());
                }
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /** The histogram section of a metrics file. */
    static String histogram(final String text) {
        final List<String> kept = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
                continue;
            }
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            if (inHistogram) {
                kept.add(line);
            }
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("umiprevalence");
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path metrics = dir.resolve("out.txt");
        final List<String> argv = new ArrayList<>(List.of("I=" + in, "O=" + metrics));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.CollectUmiPrevalenceMetrics()
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
        emit("histogram", name, histogram(Files.readString(metrics, StandardCharsets.UTF_8)));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Three sets at three positions, each of two reads sharing one UMI.
        final List<Read> oneUmiEach = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            oneUmiEach.add(paired("a" + i, 100 + i * 1000, "AAAAAA"));
            oneUmiEach.add(paired("b" + i, 100 + i * 1000, "AAAAAA"));
        }
        run("three-sets-one-umi", oneUmiEach);

        // One set of three reads carrying three different UMIs.
        run("one-set-three-umis", List.of(
                paired("a", 100, "AAAAAA"),
                paired("b", 100, "CCCCCC"),
                paired("c", 100, "GGGGGG")));

        // The same set with two reads sharing a UMI, so the count is two and not three.
        run("one-set-two-umis", List.of(
                paired("a", 100, "AAAAAA"),
                paired("b", 100, "AAAAAA"),
                paired("c", 100, "GGGGGG")));

        // A read with no UMI tag at all.
        run("no-umi-tag", List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 100, 0x1 | 0x2 | 0x40, 60, null, "IIIIII")));

        // The barcode-quality filter, whose predicate is inverted: a barcode entirely above the
        // floor is DROPPED and one with a base under it is KEPT. Two reads of one set, one of
        // each, so the set that survives holds one UMI and not two.
        run("barcode-quality-good", List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 100, 0x1 | 0x2 | 0x40, 60, "CCCCCC", "IIIIII")));
        // Every read well-formed, which leaves nothing at all.
        run("barcode-quality-all-good", List.of(
                new Read("a", 100, 0x1 | 0x2 | 0x40, 60, "AAAAAA", "IIIIII"),
                new Read("b", 100, 0x1 | 0x2 | 0x40, 60, "CCCCCC", "IIIIII")));
        // And a read with no quality tag at all, which that filter keeps.
        run("no-barcode-quality", List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 100, 0x1 | 0x2 | 0x40, 60, "CCCCCC", null)));

        // The reads the other four filters take.
        run("unaligned", List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 0, 0x1 | 0x4 | 0x40, 0, "CCCCCC", "II#III")));
        run("low-mapping-quality", List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 100, 0x1 | 0x2 | 0x40, 5, "CCCCCC", "II#III")));
        run("secondary", List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 100, 0x1 | 0x2 | 0x40 | 0x100, 60, "CCCCCC", "II#III")));

        // An unpaired read, which the default filters and --FILTER_UNPAIRED_READS=false keeps.
        final List<Read> unpaired = List.of(
                paired("a", 100, "AAAAAA"),
                new Read("b", 100, 0, 60, "CCCCCC", "II#III"));
        run("unpaired-filtered", unpaired);
        run("unpaired-kept", unpaired, "FILTER_UNPAIRED_READS=false");

        // The tags named by argument.
        run("other-tags", List.of(paired("a", 100, "AAAAAA")),
                "BARCODE_TAG=RX", "BARCODE_BQ=BQ");
        // A floor low enough that no barcode base is under it, which the inversion turns into a
        // filter that drops everything.
        run("low-barcode-floor", List.of(paired("a", 100, "AAAAAA")), "BQ=0");

        // Every read filtered.
        run("everything-filtered", List.of(
                new Read("a", 100, 0x1 | 0x2 | 0x40, 60, null, "II#III")));
        // No reads at all.
        run("empty", List.of());

        System.out.print(buf);
    }
}
