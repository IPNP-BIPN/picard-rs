/*
 * The two mate-cigar duplicate markers, taken from the reference.
 *
 * `MarkDuplicates` decides a pair's position from the two ends it has seen, which means holding
 * every unpaired end until its mate arrives. The two tools here decide it from the MATE CIGAR tag
 * instead, so one pass over a coordinate-sorted file is enough, and they differ from
 * `MarkDuplicates` and from each other in what that buys and what it costs.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE MATE CIGAR IS THE POINT: a pair whose mate is soft-clipped is keyed on the mate's
 *     UNCLIPPED end, which the `MC` tag carries, so two pairs whose mates clip differently are
 *     duplicates where their aligned positions say they are not;
 *   - A PAIR WITH NO `MC` IS SKIPPED, not guessed at: `SKIP_PAIRS_WITH_NO_MATE_CIGAR` is true by
 *     default and the record is written unmarked, and turning it off is an error rather than a
 *     second algorithm;
 *   - `MINIMUM_DISTANCE` IS THE WINDOW, and `-1` means twice the first read's length or a hundred,
 *     whichever is smaller, so a fixture whose duplicates sit further apart than that is not
 *     marked at all;
 *   - THE TWO TOOLS BREAK TIES DIFFERENTLY from `MarkDuplicates`, which is the documented
 *     difference and the reason both outputs are recorded for every case;
 *   - `SimpleMarkDuplicatesWithMateCigar` IS A `MarkDuplicates` SUBCLASS driven by
 *     `DuplicateSetIterator`, so its metrics are counted in its own loop rather than in the
 *     writing pass;
 *   - AN UNPAIRED READ IS STILL ITS OWN SET in both;
 *   - OPTICAL DUPLICATES ARE TRACKED ON THE FIRST END ONLY in the simple one, which is stricter
 *     than the pair-level tracking `MarkDuplicates` does;
 *   - AND A COORDINATE SORT IS REQUIRED, which the simple one refuses by name.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     marked\t<case>.<tool>\t<the output as sam, without its header, escaped>
 *     metrics\t<case>.<tool>\t<the metrics table without its comments, escaped>
 *     error\t<case>.<tool>\t<exception class>:<message>
 *
 * Usage: MateCigarDuplicatesDump
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

public class MateCigarDuplicatesDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One record of a fixture, with its mate's cigar spelled out rather than assumed. */
    record Read(String name, int start, String cigar, int flags, int mateStart, String mateCigar,
                int length) {}

    static SAMFileHeader header(final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 4000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads,
                         final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = header(order);
        final List<SAMRecord> records = new ArrayList<>();
        for (final Read spec : reads) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(spec.name());
            record.setFlags(spec.flags());
            record.setReferenceName("chr1");
            record.setAlignmentStart(spec.start());
            record.setCigarString(spec.cigar());
            record.setMappingQuality(60);
            record.setReadString("A".repeat(spec.length()));
            record.setBaseQualityString("I".repeat(spec.length()));
            record.setAttribute("RG", "rg1");
            if ((spec.flags() & 0x1) != 0) {
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(spec.mateStart());
                if (spec.mateCigar() != null) {
                    record.setAttribute("MC", spec.mateCigar());
                }
            }
            records.add(record);
        }
        if (order == SAMFileHeader.SortOrder.coordinate) {
            records.sort((a, b) -> {
                final int byStart = Integer.compare(a.getAlignmentStart(), b.getAlignmentStart());
                return byStart != 0 ? byStart : a.getReadName().compareTo(b.getReadName());
            });
        }
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header(order), true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /** A pair whose two ends carry each other's cigar, which is what these tools read. */
    static List<Read> pair(final String name, final int first, final String firstCigar,
                           final int second, final String secondCigar, final int length) {
        return List.of(
                new Read(name, first, firstCigar, 0x1 | 0x2 | 0x40 | 0x20, second, secondCigar,
                        length),
                new Read(name, second, secondCigar, 0x1 | 0x2 | 0x80 | 0x10, first, firstCigar,
                        length));
    }

    static String[] split(final Path file) throws Exception {
        final List<String> table = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
            }
            if (inHistogram || line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            table.add(line);
        }
        return new String[]{String.join("\n", table)};
    }

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        run(name, reads, SAMFileHeader.SortOrder.coordinate, extra);
    }

    static void run(final String name, final List<Read> reads,
                    final SAMFileHeader.SortOrder order, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("matecigar");
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads, order);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        for (final String tool : new String[]{"withmatecigar", "simple", "markduplicates"}) {
            final Path out = dir.resolve(tool + ".bam");
            final Path metrics = dir.resolve(tool + ".metrics");
            final List<String> argv = new ArrayList<>(List.of(
                    "I=" + in, "O=" + out, "M=" + metrics));
            argv.addAll(Arrays.asList(extra));
            final picard.cmdline.CommandLineProgram program = switch (tool) {
                case "withmatecigar" ->
                        new picard.sam.markduplicates.MarkDuplicatesWithMateCigar();
                case "simple" ->
                        new picard.sam.markduplicates.SimpleMarkDuplicatesWithMateCigar();
                default -> new picard.sam.markduplicates.MarkDuplicates();
            };
            final String label = name + "." + tool;
            try {
                final int code = program.instanceMain(argv.toArray(new String[0]));
                if (code != 0) {
                    emit("error", label, "exit " + code);
                    continue;
                }
            } catch (final Exception e) {
                Throwable cause = e;
                while (cause.getCause() != null && cause.getCause() != cause) {
                    cause = cause.getCause();
                }
                emit("error", label, cause.getClass().getName() + ":"
                        + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
                continue;
            }
            final StringBuilder marked = new StringBuilder();
            try (final SamReader reader = SamReaderFactory.makeDefault().open(out.toFile())) {
                for (final SAMRecord record : reader) {
                    marked.append(record.getSAMString());
                }
            }
            emit("marked", label, marked.toString());
            emit("metrics", label, split(metrics)[0]);
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Two pairs at one position, which all three tools agree about.
        final List<Read> twoPairs = new ArrayList<>();
        twoPairs.addAll(pair("HWI:1:FC:1:1:1:1", 100, "20M", 300, "20M", 20));
        twoPairs.addAll(pair("HWI:1:FC:1:1:9000:9000", 100, "20M", 300, "20M", 20));
        run("two-pairs", twoPairs);

        // The mate of the second pair is soft-clipped, so its ALIGNED start differs and its
        // UNCLIPPED one does not. Only a tool that reads the mate cigar can see that.
        final List<Read> clippedMate = new ArrayList<>();
        clippedMate.addAll(pair("HWI:1:FC:1:1:1:1", 100, "20M", 300, "20M", 20));
        clippedMate.addAll(pair("HWI:1:FC:1:1:9000:9000", 100, "20M", 305, "5S15M", 20));
        run("a-soft-clipped-mate", clippedMate);

        // The first end is clipped instead, which moves the position the set is keyed on.
        final List<Read> clippedFirst = new ArrayList<>();
        clippedFirst.addAll(pair("HWI:1:FC:1:1:1:1", 100, "20M", 300, "20M", 20));
        clippedFirst.addAll(pair("HWI:1:FC:1:1:9000:9000", 105, "5S15M", 300, "20M", 20));
        run("a-soft-clipped-first-end", clippedFirst);

        // A pair with no mate cigar at all.
        final List<Read> noMateCigar = new ArrayList<>();
        noMateCigar.addAll(pair("HWI:1:FC:1:1:1:1", 100, "20M", 300, "20M", 20));
        noMateCigar.add(new Read("HWI:1:FC:1:1:9000:9000", 100, "20M", 0x1 | 0x2 | 0x40 | 0x20,
                300, null, 20));
        noMateCigar.add(new Read("HWI:1:FC:1:1:9000:9000", 300, "20M", 0x1 | 0x2 | 0x80 | 0x10,
                100, null, 20));
        run("no-mate-cigar", noMateCigar);
        run("no-mate-cigar-not-skipped", noMateCigar, "SKIP_PAIRS_WITH_NO_MATE_CIGAR=false");

        // The window: two duplicate pairs whose 5' ends are a thousand bases apart cannot both be
        // in the buffer at the default distance.
        final List<Read> distant = new ArrayList<>();
        distant.addAll(pair("HWI:1:FC:1:1:1:1", 100, "20M", 2000, "20M", 20));
        distant.addAll(pair("HWI:1:FC:1:1:9000:9000", 100, "20M", 2000, "20M", 20));
        run("a-distant-mate", distant);
        run("a-distant-mate-with-a-wide-window", distant, "MINIMUM_DISTANCE=3000");

        // Unpaired reads, and optical duplicates.
        run("three-singles", List.of(
                new Read("HWI:1:FC:1:1:1:1", 100, "10M", 0, 0, null, 10),
                new Read("HWI:1:FC:1:1:2:2", 100, "10M", 0, 0, null, 10),
                new Read("HWI:1:FC:1:1:9000:9000", 100, "10M", 0, 0, null, 10)));
        final List<Read> optical = new ArrayList<>();
        optical.addAll(pair("HWI:1:FC:1:1:1:1", 100, "20M", 300, "20M", 20));
        optical.addAll(pair("HWI:1:FC:1:1:5:5", 100, "20M", 300, "20M", 20));
        run("optical-duplicates", optical);

        // Removal, and a file that is not coordinate sorted.
        run("remove-duplicates", twoPairs, "REMOVE_DUPLICATES=true");
        run("a-queryname-sorted-file", twoPairs, SAMFileHeader.SortOrder.queryname);

        System.out.print(buf);
    }
}
