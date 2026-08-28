/*
 * EstimateLibraryComplexity's histogram and its estimate, taken from the reference.
 *
 * The tool groups read pairs by their first bases, calls two pairs duplicates when the rest of
 * their sequence agrees closely enough, and reads a library size off the resulting histogram. What
 * is measured is which pairs reach a group, which pairs a group calls duplicates, and what the two
 * files carry.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE GROUPING IS ON THE FIRST --MIN_IDENTICAL_BASES OF BOTH ENDS, so two pairs whose sixth
 *     base differs are still one group at the default of five and two groups at six;
 *   - --MAX_DIFF_RATE IS A RATE OVER BOTH READS' COMPARED LENGTH, floored: a pair of thirty-base
 *     ends compared over fifty bases allows one error at 0.03 and none at 0.01;
 *   - AND THE COMPARISON SKIPS THE IDENTICAL PREFIX, so a difference inside it is not an error at
 *     all: it is what put the two pairs in the same group;
 *   - --MIN_MEAN_QUALITY DROPS A PAIR BEFORE ANY GROUPING, and the mean is an INTEGER division
 *     over the read's length, so a read of qualities averaging 19.9 is dropped at twenty;
 *   - AN `N` IN THE FIRST --MIN_IDENTICAL_BASES DROPS THE PAIR whatever its qualities say;
 *   - A READ SHORTER THAN --MIN_IDENTICAL_BASES DROPS THE PAIR;
 *   - THE HISTOGRAM IS ONE ROW PER DUPLICATE-SET SIZE PER LIBRARY, so three pairs in one set is a
 *     single row at size three and not three rows;
 *   - A BIN HOLDING FEWER THAN --MIN_GROUP_COUNT GROUPS IS DROPPED FROM THE METRICS AND KEPT IN
 *     THE HISTOGRAM, and the default is TWO, so one duplicate pair on its own reports nothing
 *     examined and a histogram that says otherwise;
 *   - THE ESTIMATE IS LANDER-WATERMAN OVER THE PAIRS LESS THE OPTICAL ONES AND THE UNIQUE PAIRS,
 *     so ten pairs with three duplicates answer thirteen;
 *   - AND AN OPTICAL DUPLICATE IS SUBTRACTED BEFORE THE ESTIMATE RUNS: two duplicates ten pixels
 *     apart leave nothing to estimate from and the column comes out empty;
 *   - THE READ GROUPS' LIBRARIES ARE COUNTED APART, so two libraries are two rows and a histogram
 *     of two columns, but only when the READ NAME parses as a flowcell location: the read group is
 *     recorded inside that branch, so a pair called `a` lands in `Unknown` however carefully the
 *     header declares its library;
 *   - AND --MAX_READ_LENGTH TRUNCATES BOTH THE QUALITY MEAN AND THE COMPARISON, so a pair that
 *     differs past it is a duplicate.
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

public class EstimateLibraryComplexityDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /**
     * A read name in the Illumina shape, which is what makes the read group visible.
     *
     * `ELC` only records a pair's read group when the optical-duplicate finder can parse a
     * location out of the read's NAME. A pair called `a` therefore lands in the library `Unknown`
     * however carefully the header declares one, and two libraries come out as one row.
     */
    static String readName(final int index) {
        // Far apart on the flowcell: two duplicates within the optical distance are counted as
        // OPTICAL duplicates, and the estimate subtracts those before it runs, so a fixture whose
        // reads all sit in one corner has no estimate at all.
        return "H0164ALXX140820:2:1101:" + (1000 + index * 5000) + ":" + (2000 + index * 5000);
    }

    /** Two names ten pixels apart, which the optical finder calls one cluster. */
    static String neighbouringName(final int index) {
        return "H0164ALXX140820:2:2101:" + (1000 + index * 10) + ":" + (2000 + index * 10);
    }

    /** One pair: the two ends' bases and qualities, and which library it came from. */
    record Pair(String name, String readOne, String readTwo, String qualityOne, String qualityTwo,
                String library) {}

    static int named = 0;

    static Pair pair(final String name, final String one, final String two) {
        return new Pair(readName(named++), one, two, null, null, "lib1");
    }

    static Pair inLibrary(final String library, final String one, final String two) {
        return new Pair(readName(named++), one, two, null, null, library);
    }

    static Pair withQualities(final String one, final String two, final String qualities) {
        return new Pair(readName(named++), one, two, qualities, qualities, "lib1");
    }

    static final String BASES_ONE = "ACGTACGTACGTACGTACGTACGTACGTAC";
    static final String BASES_TWO = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAA";

    /** The same pair with one base changed at the given offset of read one. */
    static Pair differing(final String name, final int offset, final char base) {
        final StringBuilder one = new StringBuilder(BASES_ONE);
        one.setCharAt(offset, base);
        return pair(name, one.toString(), BASES_TWO);
    }

    static SAMFileHeader header(final List<String> libraries) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 10000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        for (final String library : libraries) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg-" + library);
            group.setSample("sample1");
            group.setLibrary(library);
            header.addReadGroup(group);
        }
        return header;
    }

    static void writeBam(final Path bam, final List<Pair> pairs, final List<String> libraries) {
        final SAMFileHeader head = header(libraries);
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(head, false, bam.toFile())) {
            for (final Pair spec : pairs) {
                writer.addAlignment(record(head, spec, true));
                writer.addAlignment(record(head, spec, false));
            }
        }
    }

    /** One end of a pair, unmapped: the tool reads sequence and never coordinates. */
    static SAMRecord record(final SAMFileHeader head, final Pair spec, final boolean first) {
        final SAMRecord record = new SAMRecord(head);
        record.setReadName(spec.name());
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(first);
        record.setSecondOfPairFlag(!first);
        record.setReadUnmappedFlag(true);
        record.setMateUnmappedFlag(true);
        record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        final String bases = first ? spec.readOne() : spec.readTwo();
        record.setReadString(bases);
        final String qualities = first ? spec.qualityOne() : spec.qualityTwo();
        if (qualities == null) {
            record.setBaseQualityString("I".repeat(bases.length()));
        } else {
            record.setBaseQualityString(qualities);
        }
        record.setAttribute("RG", "rg-" + spec.library());
        return record;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        final List<String> one = List.of("lib1");

        // Two identical pairs: one duplicate set of two.
        run("two-identical", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", BASES_ONE, BASES_TWO)), one);

        // Three identical pairs: ONE row at size three, not three rows.
        run("three-identical", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", BASES_ONE, BASES_TWO),
                pair("c", BASES_ONE, BASES_TWO)), one);

        // Two pairs that differ at the sixth base: one group at the default seed, two at six.
        final List<Pair> sixth = List.of(
                pair("a", BASES_ONE, BASES_TWO),
                differing("b", 5, 'T'));
        run("differ-at-the-sixth-base", sixth, one);
        run("seed-of-six", sixth, one, "MIN_IDENTICAL_BASES=6");

        // Two pairs that differ inside the seed: two groups whatever the rate says.
        run("differ-inside-the-seed", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                differing("b", 2, 'T')), one);

        // One difference over sixty compared bases: inside 0.03 and outside 0.01.
        run("one-difference", sixth, one);
        run("one-difference-strict", sixth, one, "MAX_DIFF_RATE=0.01");

        // A pair whose mean quality is under the floor, and the same pair with the floor lowered.
        final String lowQuality = "4".repeat(30);
        final List<Pair> lowQualityPairs = List.of(
                pair("a", BASES_ONE, BASES_TWO),
                withQualities(BASES_ONE, BASES_TWO, lowQuality));
        run("low-mean-quality", lowQualityPairs, one);
        run("low-mean-quality-allowed", lowQualityPairs, one, "MIN_MEAN_QUALITY=1");

        // An N inside the seed, which drops the pair whatever its qualities say.
        run("an-n-in-the-seed", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", "ACNTACGTACGTACGTACGTACGTACGTAC", BASES_TWO)), one);

        // A read shorter than the seed.
        run("shorter-than-the-seed", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", "ACG", BASES_TWO)), one, "MIN_IDENTICAL_BASES=5");

        // Two duplicate groups of two, which is what the metrics' own floor wants: a bin holding
        // ONE group is dropped from the metrics and kept in the histogram.
        final String otherOne = "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTT";
        final List<Pair> twoGroups = List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", BASES_ONE, BASES_TWO),
                pair("c", otherOne, BASES_TWO),
                pair("d", otherOne, BASES_TWO));
        run("two-duplicate-groups", twoGroups, one);
        run("one-group-counted", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", BASES_ONE, BASES_TWO)), one, "MIN_GROUP_COUNT=1");

        // A library the estimate can actually work on: six singletons and two pairs of duplicates,
        // which is a fifth duplication rather than a half.
        final List<Pair> deep = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            final StringBuilder distinct = new StringBuilder(BASES_ONE);
            distinct.setCharAt(0, "ACGT".charAt(i % 4));
            distinct.setCharAt(1, "ACGT".charAt(i / 4));
            deep.add(pair("s" + i, distinct.toString(), BASES_TWO));
        }
        deep.add(pair("d1", BASES_ONE, BASES_TWO));
        deep.add(pair("d2", BASES_ONE, BASES_TWO));
        deep.add(pair("e1", otherOne, BASES_TWO));
        deep.add(pair("e2", otherOne, BASES_TWO));
        run("an-estimate", deep, one, "MIN_GROUP_COUNT=1");

        // Two duplicates ten pixels apart, which the optical finder calls one cluster: the
        // estimate subtracts them and has nothing left to work on.
        run("optical-duplicates", List.of(
                new Pair(neighbouringName(0), BASES_ONE, BASES_TWO, null, null, "lib1"),
                new Pair(neighbouringName(1), BASES_ONE, BASES_TWO, null, null, "lib1")),
                one, "MIN_GROUP_COUNT=1");

        // Two libraries, counted apart.
        run("two-libraries", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", BASES_ONE, BASES_TWO),
                inLibrary("lib2", BASES_ONE, BASES_TWO),
                inLibrary("lib2", BASES_ONE, BASES_TWO)),
                List.of("lib1", "lib2"));

        // No duplicate at all: the estimate has nothing to work with.
        run("no-duplicates", List.of(
                pair("a", BASES_ONE, BASES_TWO),
                pair("b", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", BASES_TWO)), one);

        // A difference past the truncation, which is a duplicate once the reads are cut short.
        final List<Pair> late = List.of(
                pair("a", BASES_ONE, BASES_TWO),
                differing("b", 25, 'G'));
        run("difference-past-the-truncation", late, one, "MAX_READ_LENGTH=20",
                "MAX_DIFF_RATE=0.0");
        run("difference-inside-the-window", late, one, "MAX_DIFF_RATE=0.0");

        System.out.print(buf);
    }

    static void run(final String name, final List<Pair> pairs, final List<String> libraries,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("libcomplexity");
        final Path in = dir.resolve("in.bam");
        writeBam(in, pairs, libraries);
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
            final int code = new picard.sam.markduplicates.EstimateLibraryComplexity()
                    .instanceMain(argv.toArray(new String[0]));
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
            return;
        }
        final List<String> table = new ArrayList<>();
        final List<String> histogram = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : Files.readString(metrics, StandardCharsets.UTF_8)
                .split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
                continue;
            }
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            (inHistogram ? histogram : table).add(line);
        }
        emit("metrics", name, String.join("\n", table));
        emit("histogram", name, String.join("\n", histogram));
    }
}
