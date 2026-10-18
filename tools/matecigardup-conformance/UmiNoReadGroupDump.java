/*
 * What `UmiAwareMarkDuplicatesWithMateCigar` does with a read that has no read group.
 *
 * `UmiAwareDuplicateSetIterator.process` reads `getReadGroup().getLibrary()` off a duplicate set's
 * records, and `getReadGroup()` is null for a record with no RG tag or one naming a group the
 * header does not declare, so such a read throws a NullPointerException, single-end or paired.
 * Unlike MarkDuplicates the output has been opened by then, so the file exists; what it holds at
 * the throw is part of the answer and is recorded byte for byte.
 *
 * Output:
 *
 *     marked\t<case>\t<the output as sam, without its header, escaped>
 *     error\t<case>\t<exception class>:<message>
 *     partial\t<case>\t<the sam output's text at the throw, escaped>, a bam output's size, or "absent"
 *
 * Usage: UmiNoReadGroupDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class UmiNoReadGroupDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static SAMFileHeader header(final SAMFileHeader.SortOrder order, final boolean withGroup) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 2000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        if (withGroup) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
            group.setSample("sample1");
            group.setLibrary("lib1");
            header.addReadGroup(group);
        }
        return header;
    }

    /** name, flags, start, mate start (0 when the mate is unmapped), RG or null. */
    record Read(String name, int flags, int start, int mateStart, String group) {}

    static List<Read> pairs(final String firstGroup, final String secondGroup) {
        return List.of(
                new Read("HWI:1:FC:1:1:1:1", 99, 100, 300, firstGroup),
                new Read("HWI:1:FC:1:1:1:1", 147, 300, 100, firstGroup),
                new Read("HWI:1:FC:1:1:9000:9000", 99, 100, 300, secondGroup),
                new Read("HWI:1:FC:1:1:9000:9000", 147, 300, 100, secondGroup));
    }

    static void writeBam(final Path bam, final SAMFileHeader header, final List<Read> reads) {
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                record.setReadString("A".repeat(20));
                record.setBaseQualityString("I".repeat(20));
                if (spec.start() > 0) {
                    record.setReferenceName("chr1");
                    record.setAlignmentStart(spec.start());
                }
                if ((spec.flags() & 0x4) == 0) {
                    record.setCigarString("20M");
                    record.setMappingQuality(60);
                }
                if ((spec.flags() & 0x1) != 0) {
                    record.setMateReferenceName("chr1");
                    record.setMateAlignmentStart(spec.mateStart() == 0 ? spec.start() : spec.mateStart());
                    if ((spec.flags() & 0x8) == 0) {
                        record.setAttribute("MC", "20M");
                    }
                }
                record.setAttribute("RX", "ACGT");
                if (spec.group() != null) {
                    record.setAttribute("RG", spec.group());
                }
                writer.addAlignment(record);
            }
        }
    }

    static void run(final String name, final SAMFileHeader header, final List<Read> reads)
            throws Exception {
        run(name, header, reads, "out.sam");
    }

    static void run(final String name, final SAMFileHeader header, final List<Read> reads,
                    final String output) throws Exception {
        final Path dir = Files.createTempDirectory("uminoreadgroup");
        final Path in = dir.resolve("in.bam");
        writeBam(in, header, reads);
        final Path out = dir.resolve(output);
        final String[] argv = {"I=" + in, "O=" + out, "M=" + dir.resolve("metrics.txt"),
                "UMI_METRICS=" + dir.resolve("umi.txt")};
        try {
            final int code = new picard.sam.markduplicates.UmiAwareMarkDuplicatesWithMateCigar()
                    .instanceMain(argv);
            if (code != 0) {
                emit("error", name, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            // The writer is never closed after the throw, so this is what had reached the disk.
            emit("partial", name, !Files.exists(out) ? "absent"
                    : output.endsWith(".sam")
                    ? Files.readString(out, StandardCharsets.UTF_8).replace(dir.toString(), "<dir>")
                    : Files.size(out) + " bytes");
            return;
        }
        final StringBuilder marked = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(out.toFile())) {
            for (final SAMRecord record : reader) {
                final String line = record.getSAMString();
                marked.append(line.endsWith("\n") ? line : line + "\n");
            }
        }
        emit("marked", name, marked.toString());
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        final SAMFileHeader.SortOrder coordinate = SAMFileHeader.SortOrder.coordinate;

        run("pairs-with-rg", header(coordinate, true), pairs("rg1", "rg1"));
        run("pairs-without-rg", header(coordinate, true), pairs(null, null));
        run("pairs-without-rg-and-no-rg-lines", header(coordinate, false), pairs(null, null));
        run("one-pair-with-rg-one-without", header(coordinate, true), pairs("rg1", null));
        run("one-pair-without-rg-one-with", header(coordinate, true), pairs(null, "rg1"));
        run("undeclared-rg", header(coordinate, true), pairs("rg2", "rg2"));
        run("singles-without-rg", header(coordinate, true), List.of(
                new Read("HWI:1:FC:1:1:1:1", 0, 100, 0, null),
                new Read("HWI:1:FC:1:1:9000:9000", 0, 100, 0, null)));
        // Reads with a group before the first one without: are they on disk at the throw?
        run("rg-reads-then-one-without", header(coordinate, true), List.of(
                new Read("HWI:1:FC:1:1:1:1", 0, 100, 0, "rg1"),
                new Read("HWI:1:FC:1:1:2:2", 0, 200, 0, "rg1"),
                new Read("HWI:1:FC:1:1:3:3", 0, 300, 0, "rg1"),
                new Read("HWI:1:FC:1:1:4:4", 0, 900, 0, null)));
        run("mate-unmapped-without-rg", header(coordinate, true), List.of(
                new Read("HWI:1:FC:1:1:1:1", 0x1 | 0x8 | 0x40, 100, 0, null),
                new Read("HWI:1:FC:1:1:1:1", 0x1 | 0x4 | 0x80, 100, 0, null)));
        run("unplaced-unmapped-without-rg", header(coordinate, true), List.of(
                new Read("HWI:1:FC:1:1:1:1", 0, 100, 0, "rg1"),
                new Read("HWI:1:FC:1:1:2:2", 0x4, 0, 0, null)));
        // A BAM output is block-compressed: is the header's block on disk at the throw?
        run("pairs-without-rg-to-bam", header(coordinate, true), pairs(null, null), "out.bam");
        run("queryname-pairs-without-rg", header(SAMFileHeader.SortOrder.queryname, true),
                pairs(null, null));

        System.out.print(buf);
    }
}
