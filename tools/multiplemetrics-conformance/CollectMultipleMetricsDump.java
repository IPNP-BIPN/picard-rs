/*
 * CollectMultipleMetrics' single pass, taken from the reference.
 *
 * The tool is a dispatcher: it builds one instance per PROGRAM, hands every one of them the same
 * records, and lets each write its own files. What is measured is therefore the dispatch and not
 * the arithmetic, which the per-program suites already hold: which programs run, which files they
 * land on, what a missing reference costs, and how an argument reaches one program and not the
 * others.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE DEFAULT SET IS FIVE PROGRAMS and not the eight the enum declares, so a plain run writes
 *     five metrics files and five charts and nothing for GC bias, RNA-seq or the artifacts;
 *   - `PROGRAM=null` EMPTIES THE SET, which is how a run names its own programs, and an empty set
 *     is refused by name rather than running nothing;
 *   - EVERY PROGRAM IN THE DEFAULT SET WRITES A PDF, and the pdf is R's output: two runs of one
 *     fixture do not agree byte for byte, so the bytes are not a claim this golden can make;
 *   - FILE_EXTENSION IS APPENDED TO THE METRICS FILE AND NOT TO THE CHART, so `EXT=.txt` renames
 *     five files and leaves five alone;
 *   - A PROGRAM THAT NEEDS THE REFERENCE IS REFUSED WITHOUT ONE, before any record is read, and
 *     the refusal names the program;
 *   - RnaSeqMetrics NEEDS REF_FLAT the same way, which is a second and separate check;
 *   - EXTRA_ARGUMENT REACHES ONE PROGRAM, `<PROGRAM>::<ARGUMENT_AND_VALUE>`, and changes that
 *     program's file and no other;
 *   - AN EXTRA_ARGUMENT FOR A PROGRAM THAT IS NOT RUNNING IS AN ERROR, raised after the programs
 *     are built and not when it is parsed;
 *   - AND ONE THAT NAMES INPUT OR STOP_AFTER IS NOT, the four arguments the pass owns being
 *     documented as silently ignored;
 *   - STOP_AFTER TRUNCATES THE PASS FOR EVERY PROGRAM AT ONCE, which is what makes this one pass
 *     rather than five;
 *   - AND THE NUMBERS ARE THE STANDALONE TOOL'S, so the same fixture through
 *     CollectAlignmentSummaryMetrics alone gives the file this tool writes.
 *
 * Output:
 *
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     metrics\t<case>/<suffix>\t<the metrics table without its comments, escaped>
 *     chart\t<case>\t<name>=<empty|non-empty>
 *     chart-stability\t<case>\t<name>=<same|differs>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CollectMultipleMetricsDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CollectMultipleMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 400;

    /** The reference: a repeating pattern, long enough for a pair to sit inside it. */
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

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        header.addReadGroup(group);
        return header;
    }

    /**
     * Ten pairs, each a hundred bases apart, every base matching the reference.
     *
     * The insert sizes are deliberately not all equal: CollectInsertSizeMetrics refuses to write a
     * histogram it cannot summarise, and one width of one is not a distribution.
     */
    static void writeBam(final Path bam) {
        final SAMFileHeader header = header();
        final List<SAMRecord> records = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            final int first = 1 + i * 20;
            final int second = first + 100 + (i % 3) * 10;
            records.add(read(header, "p" + i, first, second, true));
            records.add(read(header, "p" + i, second, first, false));
        }
        records.sort((a, b) -> Integer.compare(a.getAlignmentStart(), b.getAlignmentStart()));
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final int mateStart, final boolean first) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setMappingQuality(60);
        record.setCigarString("20M");
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < 20; i++) {
            bases.append("ACGT".charAt((start - 1 + i) % 4));
        }
        record.setReadString(bases.toString());
        record.setBaseQualityString("IIIIIIIIIIIIIIIIIIII");
        record.setFlags(0x1 | 0x2 | (first ? 0x40 : 0x80) | (first ? 0x20 : 0x10));
        record.setMateReferenceName("chr1");
        record.setMateAlignmentStart(mateStart);
        record.setInferredInsertSize((first ? 1 : -1) * (Math.abs(mateStart - start) + 20));
        record.setAttribute("RG", "rg1");
        return record;
    }

    /** A metrics table without its comment lines, which carry the command line and the clock. */
    static String table(final Path file) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static Path fixture() throws Exception {
        final Path dir = Files.createTempDirectory("multiplemetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        writeBam(dir.resolve("in.bam"));
        return dir;
    }

    /** One run, into its own directory, so the files it lands on are the files it wrote. */
    static Path run(final String name, final boolean withReference, final String... extra)
            throws Exception {
        final Path dir = fixture();
        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + dir.resolve("in.bam"), "O=" + out.resolve("m")));
        if (withReference) {
            argv.add("R=" + dir.resolve("ref.fasta"));
        }
        argv.addAll(Arrays.asList(extra));
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        try {
            final int code;
            try {
                System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
                code = new picard.analysis.CollectMultipleMetrics()
                        .instanceMain(argv.toArray(new String[0]));
            } finally {
                System.err.flush();
                System.setErr(realErr);
            }
            if (code != 0) {
                // A validation failure is a return code and not an exception, and the reason for
                // it is on the error stream under a usage the golden has no room for.
                emit("error", name, "exit " + code + " " + refusal(errBytes.toString(StandardCharsets.UTF_8)));
                return null;
            }
        } catch (final Exception e) {
            System.setErr(realErr);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage())
                            .replace(dir.toString(), "<dir>"));
            return null;
        }
        final List<String> written = new ArrayList<>();
        for (final File file : listed(out.toFile().listFiles())) {
            written.add(file.getName());
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            final Path path = out.resolve(file);
            if (file.endsWith(".pdf")) {
                emit("chart", name, file + "=" + (Files.size(path) > 0 ? "non-empty" : "empty"));
            } else {
                emit("metrics", name + "/" + file.substring("m".length()), table(path));
            }
        }
        return out;
    }

    /** The reason under the usage: the lines a refusal writes, without the tool's own manual. */
    static String refusal(final String stderr) {
        final List<String> kept = new ArrayList<>();
        for (final String line : stderr.split("\n", -1)) {
            final String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("*")) {
                continue;
            }
            kept.add(trimmed);
        }
        // The reason is printed UNDER the usage, as bare lines with nothing marking them, so what
        // is kept is the tail rather than a line matched by its shape.
        return kept.isEmpty() ? "" : kept.get(kept.size() - 1);
    }

    static File[] listed(final File[] files) {
        return files == null ? new File[0] : files;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // The default set, which is five of the eight programs the enum declares.
        final Path plain = run("default-programs", true);

        // The same fixture again, to ask whether the charts are bytes a golden could hold.
        final Path again = run("default-programs-again", true);
        if (plain != null && again != null) {
            for (final File file : listed(plain.toFile().listFiles())) {
                if (!file.getName().endsWith(".pdf")) {
                    continue;
                }
                final byte[] first = Files.readAllBytes(file.toPath());
                final byte[] second = Files.readAllBytes(again.resolve(file.getName()));
                emit("chart-stability", "default-programs",
                        file.getName() + "=" + (Arrays.equals(first, second) ? "same" : "differs"));
            }
        }

        // A run that names its own programs, which starts by emptying the default set.
        run("one-program", true, "PROGRAM=null", "PROGRAM=QualityScoreDistribution");
        run("two-programs", true, "PROGRAM=null", "PROGRAM=QualityScoreDistribution",
                "PROGRAM=MeanQualityByCycle");
        run("no-programs", true, "PROGRAM=null");

        // The extension, which lands on the metrics files and not on the charts.
        run("file-extension", true, "EXT=.txt");

        // The programs that need more than the reads.
        run("gc-bias-without-a-reference", false, "PROGRAM=null", "PROGRAM=CollectGcBiasMetrics");
        run("gc-bias", true, "PROGRAM=null", "PROGRAM=CollectGcBiasMetrics");
        run("artifacts-without-a-reference", false, "PROGRAM=null",
                "PROGRAM=CollectSequencingArtifactMetrics");
        run("rna-seq-without-a-refflat", true, "PROGRAM=null", "PROGRAM=RnaSeqMetrics");

        // The default set without a reference at all, which none of the five needs.
        run("default-programs-without-a-reference", false);

        // One argument, to one program.
        run("extra-argument", true, "PROGRAM=null", "PROGRAM=CollectInsertSizeMetrics",
                "EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200");
        run("extra-argument-for-another-program", true, "PROGRAM=null",
                "PROGRAM=QualityScoreDistribution",
                "EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200");
        run("extra-argument-malformed", true, "EXTRA_ARGUMENT=HISTOGRAM_WIDTH=200");
        run("extra-argument-unknown-program", true, "EXTRA_ARGUMENT=NoSuchProgram::X=1");
        run("extra-argument-the-pass-owns", true, "PROGRAM=null",
                "PROGRAM=QualityScoreDistribution",
                "EXTRA_ARGUMENT=QualityScoreDistribution::STOP_AFTER=1");

        // The pass, truncated for every program at once.
        run("stop-after", true, "STOP_AFTER=4");

        // A level the two accumulating programs support and the other three do not.
        run("accumulation-level", true, "LEVEL=ALL_READS", "LEVEL=SAMPLE");

        // And the same fixture through one of the programs on its own, which is what says the
        // dispatcher changes no number.
        final Path dir = fixture();
        final Path alone = dir.resolve("alone.alignment_summary_metrics");
        new picard.analysis.CollectAlignmentSummaryMetrics().instanceMain(new String[]{
                "I=" + dir.resolve("in.bam"), "O=" + alone, "R=" + dir.resolve("ref.fasta")});
        emit("metrics", "standalone/.alignment_summary_metrics", table(alone));

        System.out.print(buf);
    }
}
