/*
 * `CrosscheckFingerprints`, taken from the reference.
 *
 * `CheckFingerprint` asks whether one file's sample is the sample a set of genotypes says it is.
 * This one asks the question of every PAIR of inputs, with nothing to compare against but each
 * other, and answers with a matrix of LOD scores and a verdict per pair.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE ANSWER IS A PAIR-WISE LOD, positive where two fingerprints agree and negative where they
 *     do not, with a verdict of `EXPECTED_MATCH`, `UNEXPECTED_MATCH`, `EXPECTED_MISMATCH` or
 *     `UNEXPECTED_MISMATCH` on top of it;
 *   - WHAT IS EXPECTED COMES FROM THE SAMPLE NAMES, not from the reads: two read groups of one
 *     sample are expected to match, and a match between two samples is UNEXPECTED however good
 *     the LOD is;
 *   - `--CROSSCHECK_BY` DECIDES WHAT A ROW IS: read groups by default, so one file with two groups
 *     is two fingerprints and one comparison, and `SAMPLE` collapses them into one;
 *   - `--SECOND_INPUT` MAKES IT TWO SETS rather than one: every left input against every right
 *     one, and no comparison inside either;
 *   - `--CROSSCHECK_MODE` CHOOSES WHICH PAIRS ARE CHECKED AT ALL, `CHECK_SAME_SAMPLE` being the
 *     default and `CHECK_ALL_OTHERS` the one that compares across samples;
 *   - `--OUTPUT_ERRORS_ONLY` DROPS THE ROWS THAT AGREED, so a run where everything matched writes
 *     a file with a header and no rows;
 *   - `--MATRIX_OUTPUT` IS A SECOND SHAPE of the same numbers, one row per fingerprint and one
 *     column per fingerprint;
 *   - `--EXPECT_ALL_GROUPS_TO_MATCH` TURNS AN UNEXPECTED MISMATCH INTO AN EXIT CODE, which is
 *     `--EXIT_CODE_WHEN_MISMATCH`;
 *   - AND A RUN WITH NOTHING TO COMPARE RETURNS `--EXIT_CODE_WHEN_NO_VALID_CHECKS`.
 *
 * The fixture is `CheckFingerprintDump`'s, transcribed rather than shared: the same reference, the
 * same three-site haplotype map and the same BAM writer, so a case here and a case there ask about
 * the same sites while each dump still compiles on its own.
 *
 * Output:
 *
 *     sam\t<case>.<n>\t<each input as sam, without its header, escaped>
 *     metrics\t<case>\t<the crosscheck table without its comments, escaped>
 *     matrix\t<case>\t<the matrix file without its comments, escaped>
 *     code\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CrosscheckFingerprintsDump
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

public class CrosscheckFingerprintsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 600;

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

    /**
     * The haplotype map: three sites, the first two in one block and the third alone.
     *
     * Position 101 is A and position 105 is A and position 201 is A in the reference, so each
     * site's major allele is its reference base and its minor allele is C.
     */
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

    static SAMFileHeader header(final List<String> samples) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        int i = 0;
        for (final String sample : samples) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg" + (++i));
            group.setSample(sample);
            group.setLibrary("lib" + i);
            group.setPlatformUnit("unit" + i);
            header.addReadGroup(group);
        }
        return header;
    }

    /**
     * A BAM whose read groups are named after the FILE as well as the sample.
     *
     * Two files whose groups share an id are one fingerprint to a crosscheck by read group, which
     * is a merge rather than a comparison: the group is the key. So the id, the platform unit and
     * the library all carry the file's index, and only the sample is shared where a case wants it
     * shared.
     */
    static void writeNamedBam(final Path bam, final int file, final List<String> samples,
                              final List<Read> reads) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        int index = 0;
        for (final String sample : samples) {
            final SAMReadGroupRecord group =
                    new SAMReadGroupRecord("f" + file + "rg" + (++index));
            group.setSample(sample);
            group.setLibrary("f" + file + "lib" + index);
            group.setPlatformUnit("f" + file + "unit" + index);
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
                record.setAttribute("RG",
                        "f" + file + "rg" + ((position++ % samples.size()) + 1));
                writer.addAlignment(record);
            }
        }
    }

    static void writeBam(final Path bam, final List<String> samples, final List<Read> reads) {
        final SAMFileHeader header = header(samples);
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            int index = 0;
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
                // Round-robin the read groups, so a two-sample file has reads of each.
                record.setAttribute("RG", "rg" + ((index++ % samples.size()) + 1));
                writer.addAlignment(record);
            }
        }
    }
    /** Reads covering the map's three sites, carrying the allele asked for at each. */
    static List<Read> reads(final String allele) {
        final List<Read> reads = new ArrayList<>();
        // Thirty copies rather than a handful: a fingerprint is a likelihood ratio, and six reads
        // leave it near zero whatever the alleles are, so a case named for a disagreement has to
        // carry enough depth to state one.
        for (int copy = 0; copy < 30; copy++) {
            // Position 101 and 105 are in one block, 201 in another. The reference is `ACGT`
            // repeating, so a read starting at 99 carries the reference base unless it is edited.
            final StringBuilder bases = new StringBuilder();
            for (int offset = 0; offset < 12; offset++) {
                bases.append("ACGT".charAt((98 + offset) % 4));
            }
            bases.setCharAt(2, allele.charAt(0));
            bases.setCharAt(6, allele.charAt(0));
            reads.add(new Read(
                    "r" + allele + copy, 99, bases.toString(), "I".repeat(12)));
            final StringBuilder second = new StringBuilder();
            for (int offset = 0; offset < 12; offset++) {
                second.append("ACGT".charAt((198 + offset) % 4));
            }
            second.setCharAt(2, allele.charAt(0));
            reads.add(new Read(
                    "s" + allele + copy, 199, second.toString(), "I".repeat(12)));
        }
        return reads;
    }

    /**
     * A metrics table, without its comments and with the run's own directory taken out.
     *
     * The `LEFT_FILE` and `RIGHT_FILE` columns carry the input's URI, and the input lives in a
     * directory named after the clock. Leaving it in would make the dump differ from itself.
     */
    static String table(final Path file, final Path dir) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line.replace(dir.toString(), "<dir>"));
        }
        if (kept.size() > 1) {
            // The ROWS are sorted and the header is not. The tool walks its fingerprints out of a
            // hash-ordered map, so two runs of the same input write the same rows in different
            // orders; the values are what is being measured, and an order that moves between runs
            // would make the dump differ from itself.
            final List<String> rows = new ArrayList<>(kept.subList(1, kept.size()));
            java.util.Collections.sort(rows);
            kept.subList(1, kept.size()).clear();
            kept.addAll(rows);
        }
        return String.join("\n", kept);
    }

    /** One input file: its samples, and the allele its reads carry. */
    record Input(List<String> samples, String allele) {}

    static void run(final String name, final List<Input> inputs, final int secondFrom,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("crosscheck");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path map = dir.resolve("map.txt");
        Files.writeString(map, database(), StandardCharsets.UTF_8);

        final List<Path> files = new ArrayList<>();
        for (int index = 0; index < inputs.size(); index++) {
            final Path bam = dir.resolve("in" + index + ".bam");
            writeNamedBam(bam, index, inputs.get(index).samples(),
                    reads(inputs.get(index).allele()));
            files.add(bam);
            final StringBuilder sam = new StringBuilder();
            try (final SamReader reader = SamReaderFactory.makeDefault().open(bam.toFile())) {
                for (final SAMRecord record : reader) {
                    sam.append(record.getSAMString());
                }
            }
            emit("sam", name + "." + index, sam.toString());
        }

        final Path out = dir.resolve("crosscheck.txt");
        final Path matrix = dir.resolve("matrix.txt");
        final List<String> argv = new ArrayList<>(List.of("H=" + map, "R=" + reference));
        for (int index = 0; index < files.size(); index++) {
            argv.add((secondFrom >= 0 && index >= secondFrom ? "SI=" : "I=") + files.get(index));
        }
        final List<String> tail = Arrays.asList(extra);
        if (tail.contains("MATRIX")) {
            argv.add("MATRIX_OUTPUT=" + matrix);
        } else {
            argv.add("O=" + out);
        }
        for (final String argument : tail) {
            if (!argument.equals("MATRIX")) {
                argv.add(argument);
            }
        }

        final int code;
        // The tool writes its own command line to stdout as a banner, and that line carries the
        // temporary directory. What is measured is the files it writes, so the banner goes to a
        // sink rather than into the dump.
        final java.io.PrintStream original = System.out;
        try {
            System.setOut(new java.io.PrintStream(new java.io.ByteArrayOutputStream(), true,
                    StandardCharsets.UTF_8));
            code = new picard.fingerprint.CrosscheckFingerprints()
                    .instanceMain(argv.toArray(new String[0]));
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            System.setOut(original);
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        } finally {
            System.setOut(original);
        }
        emit("code", name, String.valueOf(code));
        if (Files.exists(out)) {
            emit("metrics", name, table(out, dir));
        }
        if (Files.exists(matrix)) {
            emit("matrix", name, table(matrix, dir));
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Input sampleOneMinor = new Input(List.of("sample1"), "C");
        final Input sampleOneMajor = new Input(List.of("sample1"), "A");
        final Input sampleTwoMinor = new Input(List.of("sample2"), "C");

        // One sample in two files, whose reads carry the same allele: a match that is expected.
        run("one-sample-agreeing", List.of(sampleOneMinor, sampleOneMinor), -1);
        // The same sample, disagreeing: a MISMATCH that is not expected.
        run("one-sample-disagreeing", List.of(sampleOneMinor, sampleOneMajor), -1);
        // Two samples that agree, which is a match nobody expected.
        run("two-samples-agreeing", List.of(sampleOneMinor, sampleTwoMinor), -1,
                "CROSSCHECK_MODE=CHECK_ALL_OTHERS");
        // Two samples, checked the default way: nothing to compare, and a code for it.
        run("two-samples-default-mode", List.of(sampleOneMinor, sampleTwoMinor), -1);

        // What a row is: two read groups of one sample in a single file.
        final Input twoGroups = new Input(List.of("sample1", "sample1"), "C");
        run("two-read-groups", List.of(twoGroups), -1);
        run("two-read-groups-by-sample", List.of(twoGroups), -1, "CROSSCHECK_BY=SAMPLE");
        run("two-read-groups-by-file", List.of(twoGroups), -1, "CROSSCHECK_BY=FILE");

        // Two sets rather than one.
        run("a-second-input", List.of(sampleOneMinor, sampleOneMajor), 1);

        // The other shapes and the exit codes.
        run("errors-only", List.of(sampleOneMinor, sampleOneMinor), -1, "OUTPUT_ERRORS_ONLY=true");
        run("a-matrix", List.of(sampleOneMinor, sampleOneMajor), -1, "MATRIX");
        run("expect-all-to-match", List.of(sampleOneMinor, sampleOneMajor), -1,
                "EXPECT_ALL_GROUPS_TO_MATCH=true");
        run("expect-all-to-match-with-a-code", List.of(sampleOneMinor, sampleOneMajor), -1,
                "EXPECT_ALL_GROUPS_TO_MATCH=true", "EXIT_CODE_WHEN_MISMATCH=7");
        run("a-lod-threshold", List.of(sampleOneMinor, sampleOneMajor), -1, "LOD_THRESHOLD=-100");

        System.out.print(buf);
    }
}
