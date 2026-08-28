/*
 * PositionBasedDownsampleSam's selection, taken from the reference.
 *
 * The tool keeps the reads whose flowcell position falls outside a circle, so that neighbouring
 * reads are kept or dropped together. There is no randomness in it at all: which reads survive is
 * decided by the tile's own extent and the read's coordinates in it.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE SELECTION IS PURELY POSITIONAL AND CARRIES NO SEED, so the same file downsampled twice
 *     keeps the same reads;
 *   - A FRACTION OVER A HALF INVERTS THE MASK rather than growing the circle: the circle is built
 *     for `1 - FRACTION` and the reads INSIDE it are kept;
 *   - THE TILE'S EXTENT IS TAKEN FROM THE READS THEMSELVES, so the same read is kept or dropped
 *     depending on what else is in its tile;
 *   - THE EXTENT STARTS AT ZERO AND NOT AT THE FIRST READ, the accumulator being built on zeroes
 *     rather than on the first read's coordinates, so the mask is never centred on the reads;
 *   - THE EXTENT IS THEN WIDENED BY ITS OWN SPAN OVER THE READ COUNT, which moves the boundary
 *     for a tile of few reads and barely at all for a tile of many;
 *   - EACH TILE IS MASKED SEPARATELY, so two tiles of different extents keep different reads at
 *     the same coordinates;
 *   - --REMOVE_DUPLICATE_INFORMATION CLEARS THE DUPLICATE FLAG on the reads it keeps, and turning
 *     it off leaves the flag alone;
 *   - --STOP_AFTER LIMITS BOTH PASSES, so it changes the tile extents as well as the read count;
 *   - RUNNING THE TOOL TWICE IS REFUSED by a message naming the previous run, unless
 *     --ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS says otherwise;
 *   - A FRACTION OUTSIDE [0,1] IS REFUSED by the argument validation, which is an exit code of
 *     one rather than an exception;
 *   - AND A READ NAME THE REGEX CANNOT PARSE IS NOT REFUSED: it is given a tile of -1 and masked
 *     with every other unparseable read.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     kept\t<case>\t<the kept read names, comma-separated>
 *     flags\t<case>\t<the kept reads' flags, comma-separated>
 *     pg\t<case>\t<the @PG lines of the output, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PositionBasedDownsampleSamDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One read, named the Illumina way so the parser finds a tile and two coordinates in it. */
    record Read(String name, int flags) {}

    static Read read(final int tile, final int x, final int y, final int flags) {
        return new Read("RUN:1:" + tile + ":" + x + ":" + y, flags);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, true, bam.toFile())) {
            int position = 100;
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                record.setReferenceName("chr1");
                record.setAlignmentStart(position);
                record.setMappingQuality(60);
                record.setCigarString("10M");
                record.setReadBases("ACGTACGTAC".getBytes());
                record.setBaseQualities("IIIIIIIIII".getBytes());
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
                position += 10;
            }
        }
    }

    static void run(final String name, final List<Read> reads, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("positiondownsample");
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());
        runOn(dir, name, in, extra);
    }

    /** Runs the tool over an input already on disk, so a second run can be measured too. */
    static void runOn(final Path dir, final String name, final Path in, final String... extra)
            throws Exception {
        final Path out = dir.resolve("out-" + name + ".bam");
        final List<String> argv = new ArrayList<>(List.of("I=" + in, "O=" + out));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.sam.PositionBasedDownsampleSam()
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
        final List<String> names = new ArrayList<>();
        final List<String> flags = new ArrayList<>();
        final StringBuilder programs = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(out.toFile())) {
            for (final SAMProgramRecord program : reader.getFileHeader().getProgramRecords()) {
                programs.append(program.getId()).append('=').append(program.getProgramName())
                        .append('\n');
            }
            for (final SAMRecord record : reader) {
                names.add(record.getReadName());
                flags.add(Integer.toString(record.getFlags()));
            }
        }
        emit("kept", name, String.join(",", names));
        emit("flags", name, String.join(",", flags));
        emit("pg", name, programs.toString());
    }

    /** A grid of reads over one tile, which is enough to see the mask's shape. */
    static List<Read> grid(final int tile, final int flags) {
        final List<Read> reads = new ArrayList<>();
        for (int x = 1000; x <= 5000; x += 1000) {
            for (int y = 1000; y <= 5000; y += 1000) {
                reads.add(read(tile, x, y, flags));
            }
        }
        return reads;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final List<Read> oneTile = grid(1101, 0);

        run("fraction-one-tenth", oneTile, "F=0.1");
        run("fraction-one-half", oneTile, "F=0.5");
        // Over a half the mask inverts: the circle is built for 1-F and its inside is kept.
        run("fraction-nine-tenths", oneTile, "F=0.9");
        run("fraction-one", oneTile, "F=1.0");
        run("fraction-zero", oneTile, "F=0.0");

        // The same reads twice over, which must keep exactly the same ones.
        run("repeatable", oneTile, "F=0.3");
        run("repeatable-again", oneTile, "F=0.3");

        // Two tiles of different extents, at the same coordinates.
        final List<Read> twoTiles = new ArrayList<>(grid(1101, 0));
        twoTiles.addAll(grid(1102, 0));
        twoTiles.add(read(1102, 9000, 9000, 0));
        run("two-tiles", twoTiles, "F=0.3");

        // Duplicate-flagged reads, with the flag cleared and left alone.
        final List<Read> duplicates = grid(1101, 1024);
        run("duplicates-cleared", duplicates, "F=0.3");
        run("duplicates-kept", duplicates, "F=0.3", "REMOVE_DUPLICATE_INFORMATION=false");

        // STOP_AFTER, which limits the first pass as well and so moves the tile's extent.
        run("stop-after-ten", oneTile, "F=0.3", "STOP_AFTER=10");

        // A read name the regex cannot parse, beside ones it can.
        final List<Read> unparseable = new ArrayList<>(grid(1101, 0));
        unparseable.add(new Read("no-colons-at-all", 0));
        run("unparseable-name", unparseable, "F=0.3");

        // A fraction outside the range.
        run("fraction-too-big", oneTile, "F=1.5");
        run("fraction-negative", oneTile, "F=-0.1");

        // A second run over an already-downsampled file, refused and then allowed.
        final Path dir = Files.createTempDirectory("positiondownsample-twice");
        final Path in = dir.resolve("in.bam");
        writeBam(in, oneTile);
        final Path once = dir.resolve("out-once.bam");
        new picard.sam.PositionBasedDownsampleSam()
                .instanceMain(new String[]{"I=" + in, "O=" + once, "F=0.5"});
        runOn(dir, "downsampled-twice", once, "F=0.5");
        runOn(dir, "downsampled-twice-allowed", once, "F=0.5",
                "ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS=true");

        System.out.print(buf);
    }
}
