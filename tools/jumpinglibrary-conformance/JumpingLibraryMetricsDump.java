/*
 * CollectJumpingLibraryMetrics' metrics file, taken from the reference.
 *
 * A jumping library is one whose pairs point outwards, so the tool counts the pairs that do, the
 * pairs that do not, and the ones it calls chimeric. What is measured is which pair lands in which
 * bucket, and what the arithmetic on top of those counts produces.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE ORIENTATION IS READ OFF THE TWO STRANDS AND THE TWO POSITIONS, not off the sign of the
 *     insert size: a read on the reverse strand with its mate forward and to its right is a jump,
 *     and the same two strands the other way round are not;
 *   - ONLY THE FIRST READ OF EACH PAIR IS COUNTED, so writing the mates as well changes nothing;
 *   - A PAIR WITH ONE END UNMAPPED IS A FRAGMENT, counted in neither bucket;
 *   - A PAIR WITH BOTH ENDS UNMAPPED AND UNPLACED ENDS THE FILE, the loop breaking rather than
 *     continuing, though a coordinate-sorted file puts those reads LAST so the break never cuts
 *     an aligned pair short;
 *   - THE CHIMERA THRESHOLD IS THE GREATER OF --CHIMERA_KB_MIN AND THE OUTWARD MODE, and the
 *     mode is taken in a first pass over the same file, so lowering the argument to one still
 *     leaves the mode as the floor;
 *   - THE THREE CHIMERA KINDS ARE COUNTED SEPARATELY AND REPORTED TOGETHER: an oversized insert,
 *     a tandem pair and a cross-chromosome pair all land in CHIMERIC_PAIRS;
 *   - THE ORDER OF THOSE TESTS DECIDES which counter a pair that is two of them at once lands in,
 *     though the reported total cannot tell them apart;
 *   - MQ IS CONSULTED ONLY WHEN PRESENT, so a pair with no MQ tag passes the floor on its own
 *     mapping quality alone;
 *   - THE INSERT HISTOGRAM IS TRIMMED BEFORE THE MEAN, AND THE TRIM KEEPS ONLY CONSECUTIVE BINS:
 *     it walks to the mode, then forward only while each bin follows the last by exactly one and
 *     holds at least the mode's count over --TAIL_LIMIT. Inserts a hundred apart are therefore
 *     cut back to the mode alone whatever the limit says, and a mean of 1900 comes out of a set
 *     whose three inserts average 2000;
 *   - THE LIBRARY SIZE IS ZERO WHEN THERE ARE NO DUPLICATES rather than being estimated;
 *   - EVERY PERCENTAGE IS ZERO WHEN ITS DENOMINATOR IS, rather than a NaN;
 *   - AND A FILE THAT IS NOT COORDINATE-SORTED IS REFUSED, by a message that misspells
 *     "coordinate".
 *
 * Output:
 *
 *     metrics\t<case>\t<the metrics table without its header comments, escaped>
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

public class JumpingLibraryMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One pair, described by everything the tool reads off its first read. */
    record Pair(String name, int refIndex, int start, boolean reverse, int mateRefIndex,
                int mateStart, boolean mateReverse, int insertSize, boolean duplicate,
                Integer mateQuality, int mappingQuality, boolean unmapped, boolean mateUnmapped) { }

    static Pair pair(final String name, final int refIndex, final int start, final boolean reverse,
                     final int mateRefIndex, final boolean mateReverse, final int insertSize) {
        // The mate always sits to the RIGHT of the read, so the orientation is decided by the two
        // strands alone: reverse-then-forward is outward (RF) and forward-then-reverse is inward
        // (FR). The inferred insert size is carried separately and may be signed either way.
        return new Pair(name, refIndex, start, reverse, mateRefIndex,
                start + Math.abs(insertSize), mateReverse, insertSize, false, null, 60, false,
                false);
    }

    static SAMRecord record(final SAMFileHeader header, final Pair pair) {
        final SAMRecord read = new SAMRecord(header);
        read.setReadName(pair.name());
        read.setReadPairedFlag(true);
        read.setFirstOfPairFlag(true);
        read.setReadBases("ACGTACGTAC".getBytes());
        final byte[] quality = new byte[10];
        Arrays.fill(quality, (byte) 30);
        read.setBaseQualities(quality);
        if (pair.unmapped()) {
            read.setReadUnmappedFlag(true);
            read.setReferenceIndex(pair.refIndex());
            if (pair.refIndex() >= 0) {
                read.setAlignmentStart(pair.start());
            }
        } else {
            read.setReferenceIndex(pair.refIndex());
            read.setAlignmentStart(pair.start());
            read.setCigarString("10M");
            read.setMappingQuality(pair.mappingQuality());
            read.setReadNegativeStrandFlag(pair.reverse());
        }
        read.setMateUnmappedFlag(pair.mateUnmapped());
        if (!pair.mateUnmapped()) {
            read.setMateReferenceIndex(pair.mateRefIndex());
            read.setMateAlignmentStart(Math.max(1, pair.mateStart()));
            read.setMateNegativeStrandFlag(pair.mateReverse());
        } else if (pair.unmapped()) {
            // Both ends unmapped and unplaced: the mate must carry no position at all, or the
            // writer refuses the record before the tool ever sees it.
            read.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            read.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        } else {
            read.setMateReferenceIndex(pair.refIndex());
            read.setMateAlignmentStart(Math.max(1, pair.start()));
        }
        read.setInferredInsertSize(pair.insertSize());
        read.setDuplicateReadFlag(pair.duplicate());
        if (pair.mateQuality() != null) {
            read.setAttribute(SAMTag.MQ.name(), pair.mateQuality());
        }
        return read;
    }

    /** The second read of a pair, which the tool skips: it is here to prove that it does. */
    static SAMRecord mate(final SAMFileHeader header, final Pair pair) {
        final SAMRecord read = record(header, pair);
        read.setFirstOfPairFlag(false);
        read.setSecondOfPairFlag(true);
        read.setReadName(pair.name());
        return read;
    }

    static void run(final String name, final List<Pair> pairs, final boolean withMates,
                    final SAMFileHeader.SortOrder order, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("jump");
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 300000000));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 200000000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        final File bam = new File(dir.toFile(), "in.bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setUseAsyncIo(false)
                .makeBAMWriter(header, false, bam)) {
            for (final Pair pair : pairs) {
                writer.addAlignment(record(header, pair));
                if (withMates) {
                    writer.addAlignment(mate(header, pair));
                }
            }
        }
        final File out = new File(dir.toFile(), "metrics.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + bam.getAbsolutePath(), "O=" + out.getAbsolutePath()));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.CollectJumpingLibraryMetrics()
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
            emit("error", name, cause.getClass().getName() + ":" + cause.getMessage());
            return;
        }
        // The header lines carry a timestamp and the command line, so only the table is kept.
        final StringBuilder table = new StringBuilder();
        for (final String line : Files.readString(out.toPath()).split("\n", -1)) {
            if (!line.startsWith("#") && !line.isEmpty()) {
                table.append(line).append('\n');
            }
        }
        emit("metrics", name, table.toString());
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Outward-facing pairs, which are the jumps: the first read is reverse and its mate is
        // forward, and the insert is negative for the leftmost read of an RF pair.
        final List<Pair> jumps = List.of(
                pair("j1", 0, 1000, true, 0, false, -2000),
                pair("j2", 0, 2000, true, 0, false, -2100),
                pair("j3", 0, 3000, true, 0, false, -1900));
        run("jumps-only", jumps, false, SAMFileHeader.SortOrder.coordinate);
        // The same pairs with their mates written too, which changes nothing.
        run("jumps-with-mates", jumps, true, SAMFileHeader.SortOrder.coordinate);

        // Inward-facing pairs, which are the non-jumps.
        final List<Pair> innies = List.of(
                pair("i1", 0, 1000, false, 0, true, 300),
                pair("i2", 0, 2000, false, 0, true, 350));
        run("innies-only", innies, false, SAMFileHeader.SortOrder.coordinate);

        // The three chimera kinds, one of each.
        final List<Pair> chimeras = new ArrayList<>(jumps);
        // Oversized: past the hundred-kilobase floor.
        chimeras.add(pair("over", 0, 10000, true, 0, false, -500000));
        // Tandem: both ends on the same strand.
        chimeras.add(pair("tandem", 0, 20000, true, 0, true, -2000));
        // Cross-chromosome.
        chimeras.add(pair("cross", 0, 30000, true, 1, false, 0));
        run("chimeras", chimeras, false, SAMFileHeader.SortOrder.coordinate);

        // A pair that is BOTH oversized and tandem, and one that is BOTH tandem and
        // cross-chromosome: the order of the tests decides which counter each lands in, though
        // the reported total cannot tell them apart.
        final List<Pair> overlapping = new ArrayList<>(jumps);
        overlapping.add(pair("over-and-tandem", 0, 40000, true, 0, true, -500000));
        overlapping.add(pair("tandem-and-cross", 0, 50000, true, 1, true, 0));
        run("overlapping-chimeras", overlapping, false, SAMFileHeader.SortOrder.coordinate);

        // Duplicates, which move the duplicate counts and the library size.
        final List<Pair> withDuplicates = new ArrayList<>();
        for (final Pair jump : jumps) {
            withDuplicates.add(jump);
        }
        withDuplicates.add(new Pair("dup", 0, 4000, true, 0, 6000, false, -2000, true, null, 60,
                false, false));
        run("duplicates", withDuplicates, false, SAMFileHeader.SortOrder.coordinate);

        // A fragment: one end mapped and the other not.
        final List<Pair> withFragment = new ArrayList<>(jumps);
        withFragment.add(new Pair("frag", 0, 5000, true, 0, 7000, false, 0, false, null, 60,
                false, true));
        run("fragment", withFragment, false, SAMFileHeader.SortOrder.coordinate);

        // A pair with BOTH ends unmapped and no reference index, which ENDS the file: the jump
        // written after it is never counted.
        final List<Pair> withTerminator = new ArrayList<>(jumps);
        withTerminator.add(new Pair("stop", -1, 0, false, -1, 0, false, 0, false, null, 0, true,
                true));
        withTerminator.add(pair("after-the-stop", 0, 6000, true, 0, false, -2000));
        run("unmapped-terminator", withTerminator, false, SAMFileHeader.SortOrder.coordinate);

        // The mapping-quality floor, with and without an MQ tag.
        final List<Pair> qualities = new ArrayList<>(jumps);
        qualities.add(new Pair("low-mq", 0, 7000, true, 0, 9000, false, -2000, false, null, 5,
                false, false));
        qualities.add(new Pair("low-mate-mq", 0, 8000, true, 0, 10000, false, -2000, false, 5, 60,
                false, false));
        run("quality-floor-off", qualities, false, SAMFileHeader.SortOrder.coordinate);
        run("quality-floor-thirty", qualities, false, SAMFileHeader.SortOrder.coordinate,
                "MINIMUM_MAPPING_QUALITY=30");

        // The chimera floor lowered, which turns the ordinary jumps into chimeras.
        run("chimera-floor-one", jumps, false, SAMFileHeader.SortOrder.coordinate,
                "CHIMERA_KB_MIN=1");

        // The tail limit, which trims the histogram before the mean is taken. On a distribution
        // of one observation per bin it can change nothing, so a second fixture gives the
        // histogram a mode to hold on to and a lone outlier past it.
        run("tail-limit-one", jumps, false, SAMFileHeader.SortOrder.coordinate, "T=1");
        final List<Pair> tailHeavy = new ArrayList<>();
        for (int i = 0; i < 5; i++) {
            tailHeavy.add(pair("mode" + i, 0, 1000 + i * 100, true, 0, false, -2000));
        }
        tailHeavy.add(pair("outlier", 0, 9000, true, 0, false, -60000));
        run("tail-heavy", tailHeavy, false, SAMFileHeader.SortOrder.coordinate);
        run("tail-heavy-trimmed", tailHeavy, false, SAMFileHeader.SortOrder.coordinate, "T=2");
        // And a third whose inserts ARE consecutive, which is the only shape the trim keeps past
        // the mode: it walks forward only while each bin follows the last by exactly one.
        final List<Pair> consecutive = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            consecutive.add(pair("mode" + i, 0, 1000 + i * 100, true, 0, false, -2000));
        }
        consecutive.add(pair("next1", 0, 4000, true, 0, false, -2001));
        consecutive.add(pair("next2", 0, 5000, true, 0, false, -2001));
        consecutive.add(pair("next3", 0, 6000, true, 0, false, -2002));
        run("consecutive-bins", consecutive, false, SAMFileHeader.SortOrder.coordinate);

        // No pairs at all, where every percentage divides by zero.
        run("no-pairs", List.of(), false, SAMFileHeader.SortOrder.coordinate);

        // A file that is not coordinate-sorted.
        run("unsorted", jumps, false, SAMFileHeader.SortOrder.queryname);

        System.out.print(buf);
    }
}
