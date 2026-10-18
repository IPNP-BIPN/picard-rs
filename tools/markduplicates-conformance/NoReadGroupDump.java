/*
 * What the duplicate markers do with a read that has no read group.
 *
 * MarkDuplicates' first pass keys every primary mapped read whose mate is mapped on
 * `rec.getReadGroup().getReadGroupId()`, and `getReadGroup()` is null for a record with no RG tag,
 * so a mapped pair without one throws a NullPointerException before anything is written. That
 * holds whether or not the header declares any @RG line, and in a queryname-sorted file too. A
 * single-end read, or one whose mate is unmapped, never reaches the key and is marked as usual.
 * The two mate-cigar markers never build that key, so they mark the same files without complaint.
 *
 * Output:
 *
 *     marked\t<tool>/<case>\t<the output as sam, without its header, escaped>
 *     error\t<tool>/<case>\t<exception class>:<message>
 *
 * Usage: NoReadGroupDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import picard.cmdline.CommandLineProgram;
import picard.sam.markduplicates.MarkDuplicates;
import picard.sam.markduplicates.MarkDuplicatesWithMateCigar;
import picard.sam.markduplicates.SimpleMarkDuplicatesWithMateCigar;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.function.Supplier;

public class NoReadGroupDump {

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
                record.setReferenceName("chr1");
                record.setReadString("A".repeat(20));
                record.setBaseQualityString("I".repeat(20));
                // An unmapped mate sits at its partner's position, as aligners place it.
                record.setAlignmentStart(spec.start());
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
                if (spec.group() != null) {
                    record.setAttribute("RG", spec.group());
                }
                writer.addAlignment(record);
            }
        }
    }

    static void run(final String tool, final Supplier<CommandLineProgram> program,
                    final String name, final SAMFileHeader header, final List<Read> reads)
            throws Exception {
        final String key = tool + "/" + name;
        final Path dir = Files.createTempDirectory("noreadgroup");
        final Path in = dir.resolve("in.bam");
        writeBam(in, header, reads);
        final Path out = dir.resolve("out.sam");
        final String[] argv = {"I=" + in, "O=" + out, "M=" + dir.resolve("metrics.txt")};
        try {
            final int code = program.get().instanceMain(argv);
            if (code != 0) {
                emit("error", key, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", key, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        }
        final StringBuilder marked = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(out.toFile())) {
            for (final SAMRecord record : reader) {
                final String line = record.getSAMString();
                marked.append(line.endsWith("\n") ? line : line + "\n");
            }
        }
        emit("marked", key, marked.toString());
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        final SAMFileHeader.SortOrder coordinate = SAMFileHeader.SortOrder.coordinate;

        final List<Object[]> tools = List.of(
                new Object[]{"MarkDuplicates", (Supplier<CommandLineProgram>) MarkDuplicates::new},
                new Object[]{"MarkDuplicatesWithMateCigar",
                        (Supplier<CommandLineProgram>) MarkDuplicatesWithMateCigar::new},
                new Object[]{"SimpleMarkDuplicatesWithMateCigar",
                        (Supplier<CommandLineProgram>) SimpleMarkDuplicatesWithMateCigar::new});
        for (final Object[] one : tools) {
            final String tool = (String) one[0];
            @SuppressWarnings("unchecked")
            final Supplier<CommandLineProgram> program = (Supplier<CommandLineProgram>) one[1];
            run(tool, program, "pairs-without-rg", header(coordinate, true), pairs(null, null));
            run(tool, program, "pairs-without-rg-and-no-rg-lines", header(coordinate, false),
                    pairs(null, null));
            run(tool, program, "one-pair-with-rg-one-without", header(coordinate, true),
                    pairs("rg1", null));
            run(tool, program, "singles-without-rg", header(coordinate, true), List.of(
                    new Read("HWI:1:FC:1:1:1:1", 0, 100, 0, null),
                    new Read("HWI:1:FC:1:1:9000:9000", 0, 100, 0, null)));
            run(tool, program, "mate-unmapped-without-rg", header(coordinate, true), List.of(
                    new Read("HWI:1:FC:1:1:1:1", 0x1 | 0x8 | 0x40, 100, 0, null),
                    new Read("HWI:1:FC:1:1:1:1", 0x1 | 0x4 | 0x80, 100, 0, null)));
        }
        // Only MarkDuplicates reads a queryname-sorted file; the other two refuse it by name.
        run("MarkDuplicates", MarkDuplicates::new, "queryname-pairs-without-rg",
                header(SAMFileHeader.SortOrder.queryname, true), pairs(null, null));

        System.out.print(buf);
    }
}
