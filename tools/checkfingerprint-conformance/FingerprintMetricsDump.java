/*
 * `CalculateFingerprintMetrics` and `CrosscheckReadGroupFingerprints`, taken from the reference.
 *
 * The first asks how good a fingerprint IS, with nothing to compare it against: it samples the
 * haplotype blocks and reports how far the observed genotypes are from what the panel's allele
 * frequencies would predict, which is what catches a contaminated or a mixed-up file before any
 * comparison is made. The second is `CrosscheckFingerprints` with its defaults moved, and is
 * deprecated in the reference, which is itself part of what is measured: a tool that refuses a
 * combination its parent accepts.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE METRICS ARE PER FINGERPRINT, not per pair, so one file of two read groups is two rows;
 *   - `--CALCULATE_BY` DECIDES WHAT A ROW IS, the same way `--CROSSCHECK_BY` does next door;
 *   - THE CHI-SQUARED AND THE LOD COLUMNS ARE SAMPLED, a hundred times by a constant the tool does
 *     not expose, so what a port has to reproduce is a number that came out of a random walk with
 *     a fixed seed;
 *   - THE COUNTS ARE NOT: the number of haplotypes, how many were genotyped, and how many were
 *     called homozygous or heterozygous are arithmetic over the same fingerprint;
 *   - `CrosscheckReadGroupFingerprints` ROLLS UP with `--CROSSCHECK_SAMPLES` and
 *     `--CROSSCHECK_LIBRARIES` rather than with `--CROSSCHECK_BY`;
 *   - IT REFUSES `--CROSSCHECK_BY` ALTOGETHER, which its parent's own validation would have
 *     accepted, and the refusal is a command line one;
 *   - AND `--EXPECT_ALL_READ_GROUPS_TO_MATCH` IS MUTUALLY EXCLUSIVE with the parent's
 *     `--EXPECT_ALL_GROUPS_TO_MATCH`, which is a different refusal again.
 *
 * The fixture is `CheckFingerprintDump`'s, transcribed rather than shared: the same reference, the
 * same three-site haplotype map and the same BAM writer.
 *
 * Output:
 *
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     code\t<case>\t<the exit status>
 *     refusal\t<case>\t<the lines a refused command line printed, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: FingerprintMetricsDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FingerprintMetricsDump {

    static final StringBuilder buf = new StringBuilder();
    static final int CONTIG_LENGTH = 600;

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** chr1 is `ACGT` repeating, so the base at a position is known by arithmetic. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** Three sites, the first two in one block and the third alone. */
    static String database() {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6\tSO:coordinate");
        lines.add("@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH);
        lines.add("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS");
        lines.add("chr1\t101\trs1\tA\tC\t0.4\trs1\t");
        lines.add("chr1\t105\trs2\tA\tC\t0.4\trs1\t");
        lines.add("chr1\t201\trs3\tA\tC\t0.3\trs3\t");
        return String.join("\n", lines) + "\n";
    }

    record Read(String name, int start, String bases, String qualities) {}

    static void writeBam(final Path bam, final List<String> samples, final List<Read> reads) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        int index = 0;
        for (final String sample : samples) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg" + (++index));
            group.setSample(sample);
            group.setLibrary("lib" + index);
            group.setPlatformUnit("unit" + index);
            header.addReadGroup(group);
        }
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            int position = 0;
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(0);
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setMappingQuality(60);
                record.setCigarString(spec.bases().length() + "M");
                record.setReadString(spec.bases());
                record.setBaseQualityString(spec.qualities());
                record.setAttribute("RG", "rg" + ((position++ % samples.size()) + 1));
                writer.addAlignment(record);
            }
        }
    }

    /** Reads over the map's sites, carrying the allele asked for. */
    static List<Read> reads(final String allele) {
        final List<Read> reads = new ArrayList<>();
        for (int copy = 0; copy < 30; copy++) {
            final StringBuilder first = new StringBuilder();
            for (int offset = 0; offset < 12; offset++) {
                first.append("ACGT".charAt((98 + offset) % 4));
            }
            first.setCharAt(2, allele.charAt(0));
            first.setCharAt(6, allele.charAt(0));
            reads.add(new Read("r" + allele + copy, 99, first.toString(), "I".repeat(12)));
            final StringBuilder second = new StringBuilder();
            for (int offset = 0; offset < 12; offset++) {
                second.append("ACGT".charAt((198 + offset) % 4));
            }
            second.setCharAt(2, allele.charAt(0));
            reads.add(new Read("s" + allele + copy, 199, second.toString(), "I".repeat(12)));
        }
        return reads;
    }

    static String table(final Path file, final Path dir) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line.replace(dir.toString(), "<dir>"));
        }
        if (kept.size() > 1) {
            // The rows are sorted and the header is not: the tool walks its fingerprints out of a
            // hash-ordered map, so the order moves between runs while the values do not.
            final List<String> rows = new ArrayList<>(kept.subList(1, kept.size()));
            java.util.Collections.sort(rows);
            kept.subList(1, kept.size()).clear();
            kept.addAll(rows);
        }
        // A ROLL-UP is a matrix rather than a table: the header is one column per fingerprint and
        // that order is the same hash-ordered walk, so the columns are sorted too and every row is
        // permuted the same way. A matrix is told from a table by its header, whose columns are the
        // rows' own names.
        if (kept.size() > 1) {
            final String[] header = kept.get(0).split("\t");
            final List<String> columns =
                    new ArrayList<>(List.of(header).subList(1, header.length));
            boolean matrix = !columns.isEmpty();
            for (final String row : kept.subList(1, kept.size())) {
                matrix = matrix && columns.contains(row.split("\t")[0]);
            }
            if (matrix) {
                final List<String> sorted = new ArrayList<>(columns);
                java.util.Collections.sort(sorted);
                final List<String> out = new ArrayList<>();
                out.add(header[0] + "\t" + String.join("\t", sorted));
                for (final String row : kept.subList(1, kept.size())) {
                    final String[] values = row.split("\t");
                    final StringBuilder rebuilt = new StringBuilder(values[0]);
                    for (final String column : sorted) {
                        rebuilt.append('\t').append(values[columns.indexOf(column) + 1]);
                    }
                    out.add(rebuilt.toString());
                }
                return String.join("\n", out);
            }
        }
        return String.join("\n", kept);
    }

    /** One run of either tool over one input file. */
    static void run(final String name, final boolean crosscheck, final List<String> samples,
                    final String allele, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("fingerprintmetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path map = dir.resolve("map.txt");
        Files.writeString(map, database(), StandardCharsets.UTF_8);
        final Path in = dir.resolve("in.bam");
        writeBam(in, samples, reads(allele));

        final Path out = dir.resolve("metrics.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "H=" + map, "R=" + reference, "O=" + out));
        argv.addAll(Arrays.asList(extra));

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream errors = new ByteArrayOutputStream();
        final int code;
        try {
            // A command line the tool refuses is refused with a status rather than an exception,
            // and the reason is printed. Both streams are captured because the two refusals here
            // are not printed to the same one.
            System.setOut(new PrintStream(errors, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(errors, true, StandardCharsets.UTF_8));
            code = crosscheck
                    ? new picard.fingerprint.CrosscheckReadGroupFingerprints()
                            .instanceMain(argv.toArray(new String[0]))
                    : new picard.fingerprint.CalculateFingerprintMetrics()
                            .instanceMain(argv.toArray(new String[0]));
        } catch (final Exception e) {
            System.setOut(original);
            System.setErr(originalError);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        } finally {
            System.setOut(original);
            System.setErr(originalError);
        }
        emit("code", name, String.valueOf(code));
        if (code != 0) {
            // The refusal itself.
            final List<String> said = new ArrayList<>();
            for (final String line : errors.toString(StandardCharsets.UTF_8).split("\n", -1)) {
                // The stream carries the usage and the log as well as the refusal. The log lines
                // carry a clock, the usage is measured next door, and what is left is the reason.
                // The stream carries the whole usage as well as the refusal, and the usage is
                // measured next door. What a refusal IS is the `ERROR:` line, and a command line
                // refused for a reason the tool states has exactly one.
                // Two shapes: the parser's own `ERROR:` for a mutex pair, and the TOOL's own
                // sentence for a combination it declines by hand, which carries no prefix at all.
                if (!line.startsWith("ERROR:") && !line.contains("please refrain from supplying")
                        && !line.startsWith("(Found value")
                        && !line.startsWith("Use CrosscheckFingerprints")) {
                    continue;
                }
                said.add(line.replace(dir.toString(), "<dir>"));
            }
            emit("refusal", name, String.join("\n", said));
        }
        if (Files.exists(out)) {
            emit("metrics", name, table(out, dir));
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // CalculateFingerprintMetrics: one row per fingerprint, and what a row is.
        run("one-read-group", false, List.of("sample1"), "C");
        run("two-read-groups", false, List.of("sample1", "sample1"), "C");
        run("two-read-groups-by-sample", false, List.of("sample1", "sample1"), "C",
                "CALCULATE_BY=SAMPLE");
        run("two-read-groups-by-file", false, List.of("sample1", "sample1"), "C",
                "CALCULATE_BY=FILE");
        run("two-samples", false, List.of("sample1", "sample2"), "C");
        // The major allele instead of the minor one, which is a different fingerprint entirely.
        run("the-major-allele", false, List.of("sample1"), "A");

        // CrosscheckReadGroupFingerprints: the deprecated wrapper, and what it refuses.
        run("crosscheck-read-groups", true, List.of("sample1", "sample1"), "C");
        run("crosscheck-rolled-up-to-samples", true, List.of("sample1", "sample1"), "C",
                "CROSSCHECK_SAMPLES=true");
        run("crosscheck-rolled-up-to-libraries", true, List.of("sample1", "sample1"), "C",
                "CROSSCHECK_LIBRARIES=true");
        run("crosscheck-by-is-refused", true, List.of("sample1", "sample1"), "C",
                "CROSSCHECK_BY=SAMPLE");
        run("the-two-expectations-are-exclusive", true, List.of("sample1", "sample1"), "C",
                "EXPECT_ALL_READ_GROUPS_TO_MATCH=true", "EXPECT_ALL_GROUPS_TO_MATCH=true");

        System.out.print(buf);
    }
}
