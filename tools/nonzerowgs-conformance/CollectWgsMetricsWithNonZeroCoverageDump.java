/*
 * CollectWgsMetricsWithNonZeroCoverage's two rows, taken from the reference.
 *
 * The tool is CollectWgsMetrics run twice over one traversal: once over the whole territory, and
 * once over the same counts with the depth-zero bin set to zero. What is measured is what the
 * second pass changes, which columns it leaves alone, and what the chart argument does.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE OUTPUT IS TWO ROWS AND NOT ONE, carrying a CATEGORY column the parent has no trace of,
 *     WHOLE_GENOME first and NON_ZERO_REGIONS second;
 *   - THE SECOND ROW'S TERRITORY IS THE COVERED BASES, the depth-zero bin being dropped rather
 *     than the loci being walked again, so an uncovered half halves GENOME_TERRITORY;
 *   - AND ITS MEAN_COVERAGE RISES BY EXACTLY THAT RATIO, the same bases over a smaller territory;
 *   - THE EXCLUSION PERCENTAGES ARE RECOMPUTED AND NOT COPIED, their denominator being the second
 *     row's own total;
 *   - A FULLY COVERED REFERENCE MAKES THE TWO ROWS IDENTICAL apart from the category, which is
 *     what says the second row is a recomputation and not a different traversal;
 *   - THE HISTOGRAM IS TWO COLUMNS, `count_WHOLE_GENOME` and `count_NON_ZERO_REGIONS`, in one
 *     table rather than two sections, and the second column's depth-zero cell is zero;
 *   - --INCLUDE_BQ_HISTOGRAM ADDS A THIRD COLUMN to that same table;
 *   - --CHART_OUTPUT IS REQUIRED, the run refusing to start without it, which is the one argument
 *     the parent does not have;
 *   - THE CHART IS WRITTEN EVEN WHEN NOTHING IS COVERED: the emptiness test asks the histogram
 *     whether it has bins and not whether its counts are non-zero, and the bins are created for
 *     every depth up to the cap, so the "no valid bases" warning is unreachable;
 *   - AND A FILE WHOSE READS ARE ALL EXCLUDED IS THE ZERO-TERRITORY CASE, the second row dividing
 *     by a territory of nought.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     histogram\t<case>\t<the histogram section, escaped>
 *     chart\t<case>\t<the chart file's first line, or `absent`>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectWgsMetricsWithNonZeroCoverageDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 100;
    static final int TERRITORY = CONTIG_LENGTH - 10;

    /** The reference: a repeating pattern, with ten Ns at the end so the territory is smaller. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < TERRITORY; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        bases.append("NNNNNNNNNN");
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** One read: where it sits, how long it is, and what its flags say. */
    record Read(String name, int start, int length, int flags, int mappingQuality, int mateStart) {}

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader head = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(head, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(head);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setMappingQuality(spec.mappingQuality());
                record.setCigarString(spec.length() + "M");
                final StringBuilder bases = new StringBuilder();
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    bases.append('A');
                    quals.append('I');
                }
                record.setReadString(bases.toString());
                record.setBaseQualityString(quals.toString());
                if ((spec.flags() & 0x1) != 0) {
                    record.setMateReferenceName("chr1");
                    record.setMateAlignmentStart(spec.mateStart());
                }
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /** The metrics table and the histogram section, each without its comment lines. */
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

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("nonzerowgs");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
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
        final Path chart = dir.resolve("chart.pdf");
        final List<String> extras = Arrays.asList(extra);
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + metrics, "R=" + reference));
        if (!extras.contains("NO_CHART")) {
            argv.add("CHART=" + chart);
        }
        for (final String value : extras) {
            if (!value.equals("NO_CHART")) {
                argv.add(value);
            }
        }
        try {
            final int code = new picard.analysis.CollectWgsMetricsWithNonZeroCoverage()
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
        // The chart is a PDF an R script wrote, whose bytes carry the date it was written on.
        // Only its presence and its version line are stable, which is what the golden holds.
        if (Files.exists(chart)) {
            final String first = Files.readAllLines(chart, StandardCharsets.ISO_8859_1).isEmpty()
                    ? "" : Files.readAllLines(chart, StandardCharsets.ISO_8859_1).get(0);
            emit("chart", name, first);
        } else {
            emit("chart", name, "absent");
        }
    }

    static Read read(final String name, final int start, final int length) {
        return new Read(name, start, length, 0, 60, 0);
    }

    /** A pair, whose ends the defaults count where an unpaired read is excluded. */
    static List<Read> pair(final String name, final int first, final int second, final int length) {
        return List.of(
                new Read(name, first, length, 0x1 | 0x2 | 0x40, 60, second),
                new Read(name, second, length, 0x1 | 0x2 | 0x80, 60, first));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // A pair covering forty of the ninety bases, so the two rows differ.
        run("partly-covered", pair("a", 1, 21, 20));

        // Reads covering every base of the territory, so the two rows agree.
        final List<Read> whole = new ArrayList<>();
        for (int start = 1; start <= TERRITORY; start += 10) {
            whole.addAll(pair("w" + start, start, start, Math.min(10, TERRITORY - start + 1)));
        }
        run("fully-covered", whole);

        // A single unpaired read, which the defaults exclude, and the same read counted.
        run("all-excluded", List.of(read("a", 1, 20)));
        run("unpaired-counted", List.of(read("a", 1, 20)), "COUNT_UNPAIRED=true");

        // A file with no reads at all.
        run("no-reads", List.of());

        // Ten pairs over the same twenty bases, under a cap of two.
        final List<Read> deep = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            deep.addAll(pair("d" + i, 1, 1, 20));
        }
        run("deep-uncapped", deep);
        run("deep-capped", deep, "COVERAGE_CAP=2");

        // The base-quality histogram, which becomes a third column of the same table.
        run("with-bq-histogram", pair("a", 1, 21, 20), "INCLUDE_BQ_HISTOGRAM=true");

        // Every read a duplicate, so reads exist and nothing is covered.
        final List<Read> duplicates = new ArrayList<>();
        for (final Read one : pair("a", 1, 21, 20)) {
            duplicates.add(new Read(one.name(), one.start(), one.length(),
                    one.flags() | 0x400, one.mappingQuality(), one.mateStart()));
        }
        run("all-duplicates", duplicates);

        // The chart argument left out.
        run("no-chart", pair("a", 1, 21, 20), "NO_CHART");

        System.out.print(buf);
    }
}
