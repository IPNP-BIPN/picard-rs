/*
 * CollectDuplicateMetrics' metrics, taken from the reference.
 *
 * The tool reads a file that has already been duplicate-marked and tallies the marks. It never
 * marks anything itself, so what is measured is which reads reach which counter and what the
 * derived fields make of them.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE FOUR COUNTERS ARE A CHAIN AND NOT A SET: an unmapped read is counted as unmapped and
 *     nothing else, a secondary one as secondary and nothing else, and only what falls past both
 *     is examined as paired or unpaired;
 *   - A PAIRED READ WHOSE MATE IS UNMAPPED IS COUNTED AS UNPAIRED, so a half-mapped pair
 *     contributes one unpaired read and one unmapped one;
 *   - READ_PAIRS_EXAMINED AND READ_PAIR_DUPLICATES ARE HALVED AT THE END, counting reads while
 *     the walk runs and pairs in the file;
 *   - THE HALVING IS AN INTEGER DIVISION, so a lone paired read reports NO pairs examined and
 *     five paired reads report the same two as four do;
 *   - A DUPLICATE THAT IS UNMAPPED OR SECONDARY COUNTS NOWHERE, the duplicate tally having its
 *     own guard rather than reusing the first;
 *   - READ_PAIR_OPTICAL_DUPLICATES IS ALWAYS ZERO, which is what the tool's own summary warns
 *     about, and ESTIMATED_LIBRARY_SIZE is computed from that zero;
 *   - PERCENT_DUPLICATION WEIGHTS THE PAIRS BY TWO on both sides of its fraction, and is zero
 *     rather than absent when nothing was examined;
 *   - THERE IS ONE ROW PER LIBRARY NAMED IN THE HEADER, whether a read ever used it or not, so a
 *     file with no reads at all still writes a row of zeros;
 *   - A READ GROUP WITH NO LB AT ALL FALLS UNDER `Unknown Library`;
 *   - THE HISTOGRAM IS WRITTEN ONLY FOR A ONE-LIBRARY FILE whose estimated size is not null, so a
 *     file of unpaired reads alone has none;
 *   - AND AN ESTIMATED SIZE OF ZERO MAKES EVERY HISTOGRAM BIN NaN, which the writer renders as
 *     `?` rather than as a number or an empty field.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
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

public class CollectDuplicateMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One read: its name, its flags and the read group it belongs to. */
    record Read(String name, int flags, String group) {}

    static Read read(final String name, final int flags) {
        return new Read(name, flags, "rg1");
    }

    static final int PAIRED = 0x1;
    static final int PROPER = 0x2;
    static final int UNMAPPED = 0x4;
    static final int MATE_UNMAPPED = 0x8;
    static final int FIRST = 0x40;
    static final int SECOND = 0x80;
    static final int SECONDARY = 0x100;
    static final int DUPLICATE = 0x400;
    static final int SUPPLEMENTARY = 0x800;

    /** A header whose read groups carry the libraries named, a null naming none. */
    static SAMFileHeader header(final List<String[]> groups) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String[] group : groups) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample("sample1");
            if (group[1] != null) {
                record.setLibrary(group[1]);
            }
            header.addReadGroup(record);
        }
        return header;
    }

    static void writeBam(final Path bam, final SAMFileHeader header, final List<Read> reads) {
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, false, bam.toFile())) {
            int position = 100;
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                if ((spec.flags() & UNMAPPED) != 0) {
                    record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                    record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    record.setMappingQuality(0);
                    record.setCigarString("*");
                } else {
                    record.setReferenceName("chr1");
                    record.setAlignmentStart(position);
                    record.setMappingQuality(60);
                    record.setCigarString("10M");
                    position += 10;
                }
                if ((spec.flags() & PAIRED) != 0) {
                    if ((spec.flags() & MATE_UNMAPPED) != 0) {
                        record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                        record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    } else if ((spec.flags() & UNMAPPED) != 0) {
                        // An unmapped read whose mate is mapped is placed at its mate.
                        record.setReferenceName("chr1");
                        record.setAlignmentStart(100);
                        record.setMateReferenceName("chr1");
                        record.setMateAlignmentStart(100);
                    } else {
                        record.setMateReferenceName("chr1");
                        record.setMateAlignmentStart(200);
                    }
                }
                record.setReadBases("ACGTACGTAC".getBytes());
                record.setBaseQualities("IIIIIIIIII".getBytes());
                record.setAttribute("RG", spec.group());
                writer.addAlignment(record);
            }
        }
    }

    /** The metrics table without its comment lines, and its histogram section apart. */
    static String[] split(final String text) {
        final List<String> table = new ArrayList<>();
        final List<String> histogram = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
                continue;
            }
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            (inHistogram ? histogram : table).add(line);
        }
        return new String[]{String.join("\n", table), String.join("\n", histogram)};
    }

    static void run(final String name, final List<String[]> groups, final List<Read> reads,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("duplicatemetrics");
        final Path in = dir.resolve("in.bam");
        final SAMFileHeader header = header(groups);
        writeBam(in, header, reads);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path metrics = dir.resolve("out.txt");
        final List<String> argv = new ArrayList<>(List.of("I=" + in, "M=" + metrics));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.sam.markduplicates.CollectDuplicateMetrics()
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
        final String[] parts = split(Files.readString(metrics, StandardCharsets.UTF_8));
        emit("metrics", name, parts[0]);
        emit("histogram", name, parts[1]);
    }

    static final List<String[]> ONE_LIBRARY = List.<String[]>of(new String[]{"rg1", "lib1"});

    /** A pair of reads of one name, both mapped, optionally marked duplicate. */
    static List<Read> pair(final String name, final boolean duplicate) {
        final int extra = duplicate ? DUPLICATE : 0;
        return List.of(
                read(name, PAIRED | PROPER | FIRST | extra),
                read(name, PAIRED | PROPER | SECOND | extra));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Two pairs, one of them marked duplicate.
        final List<Read> pairs = new ArrayList<>(pair("a", false));
        pairs.addAll(pair("b", true));
        run("two-pairs-one-duplicate", ONE_LIBRARY, pairs);

        // An odd number of paired reads, so the halving loses one.
        final List<Read> odd = new ArrayList<>(pairs);
        odd.add(read("c", PAIRED | PROPER | FIRST));
        run("odd-paired-count", ONE_LIBRARY, odd);

        // An unpaired read, marked duplicate and not.
        run("unpaired", ONE_LIBRARY, List.of(read("a", 0), read("b", DUPLICATE)));

        // A pair whose mate is unmapped: one unpaired read and one unmapped one.
        run("mate-unmapped", ONE_LIBRARY, List.of(
                read("a", PAIRED | FIRST | MATE_UNMAPPED),
                read("a", PAIRED | SECOND | UNMAPPED)));

        // An unmapped read that is marked duplicate, which counts nowhere.
        run("unmapped-duplicate", ONE_LIBRARY, List.of(
                read("a", UNMAPPED | DUPLICATE),
                read("b", 0)));

        // A secondary and a supplementary read, both marked duplicate.
        run("secondary-and-supplementary", ONE_LIBRARY, List.of(
                read("a", PAIRED | PROPER | FIRST),
                read("a", PAIRED | PROPER | FIRST | SECONDARY | DUPLICATE),
                read("a", PAIRED | PROPER | FIRST | SUPPLEMENTARY | DUPLICATE)));

        // Two libraries, so two rows and no histogram.
        final List<Read> twoLibraries = new ArrayList<>();
        twoLibraries.add(new Read("a", PAIRED | PROPER | FIRST, "rg1"));
        twoLibraries.add(new Read("a", PAIRED | PROPER | SECOND, "rg1"));
        twoLibraries.add(new Read("b", PAIRED | PROPER | FIRST | DUPLICATE, "rg2"));
        twoLibraries.add(new Read("b", PAIRED | PROPER | SECOND | DUPLICATE, "rg2"));
        run("two-libraries", List.<String[]>of(
                new String[]{"rg1", "lib1"}, new String[]{"rg2", "lib2"}), twoLibraries);

        // A read group with no library at all.
        run("no-library", List.<String[]>of(new String[]{"rg1", null}), pairs);

        // Every pair a duplicate, which is where the estimate has the least to work with.
        final List<Read> allDuplicates = new ArrayList<>(pair("a", true));
        allDuplicates.addAll(pair("b", true));
        run("all-duplicates", ONE_LIBRARY, allDuplicates);

        // Nothing examined at all: only unmapped reads.
        run("only-unmapped", ONE_LIBRARY, List.of(read("a", UNMAPPED), read("b", UNMAPPED)));

        // A file with no reads.
        run("empty", ONE_LIBRARY, List.of());

        System.out.print(buf);
    }
}
