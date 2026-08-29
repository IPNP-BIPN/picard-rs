/*
 * `CollectSamErrorMetrics`, taken from the reference.
 *
 * The tool counts how often a base disagrees with the reference where it has no business
 * disagreeing, and reports that as an empirical quality. What makes it more than a mismatch
 * counter is everything it refuses to count: a base at a site the sample is known to be polymorphic
 * at, a base below a quality, a read below a mapping quality, and the second observation of a pair
 * that overlaps itself.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE ERROR RATE IS BAYESIAN, not a ratio: `PRIOR_Q` is a pseudo-count in phred space, so a
 *     file with no errors at all reports a finite quality rather than an infinite one;
 *   - THE OUTPUT IS ONE FILE PER METRIC, named `<basename>.<metric>`, so `--ERROR_METRICS` decides
 *     how many files a run writes;
 *   - A STRATIFIER SPLITS THE ROWS and nothing else: the same bases, counted per bin;
 *   - `--VCF` IS SUBTRACTIVE: a mismatch at a site the VCF carries is not an error, and removing
 *     it changes the count rather than the rows;
 *   - `--MIN_BASE_Q` AND `--MIN_MAPPING_Q` DROP OBSERVATIONS, one by base and one by read;
 *   - AN OVERLAPPING PAIR IS COUNTED ONCE, and `OVERLAPPING_ERROR` is where the second
 *     observation goes;
 *   - AN INDEL IS ITS OWN METRIC with its own denominator;
 *   - `--INTERVALS` LIMITS THE LOCI, and the rows shrink with them;
 *   - `--MAX_LOCI` STOPS EARLY, which is a different kind of subset again;
 *   - AND `--ERROR_METRICS` IS APPENDED TO rather than replaced: its default is twenty-seven
 *     entries, so naming one without clearing the list asks for it twice, and the run is REFUSED
 *     for a duplicated suffix rather than writing a file twice.
 *
 * Output:
 *
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     metrics\t<case>.<file>\t<the table without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: SamErrorMetricsDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SamErrorMetricsDump {

    static final StringBuilder buf = new StringBuilder();
    static final int CONTIG_LENGTH = 600;

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** chr1 is `ACGT` repeating, so a mismatch is made by writing any other base. */
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

    /** The reference's own bases over a window, which a read copies before it is edited. */
    static String reference(final int start, final int length) {
        final StringBuilder bases = new StringBuilder();
        for (int offset = 0; offset < length; offset++) {
            bases.append("ACGT".charAt((start - 1 + offset) % 4));
        }
        return bases.toString();
    }

    record Read(String name, int start, String bases, String qualities, int flags, int mateStart,
                String cigar, int mappingQuality) {}

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

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        final List<SAMRecord> records = new ArrayList<>();
        for (final Read spec : reads) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(spec.name());
            record.setFlags(spec.flags());
            record.setReferenceName("chr1");
            record.setAlignmentStart(spec.start());
            record.setCigarString(spec.cigar());
            record.setMappingQuality(spec.mappingQuality());
            record.setReadString(spec.bases());
            record.setBaseQualityString(spec.qualities());
            record.setAttribute("RG", "rg1");
            if ((spec.flags() & 0x1) != 0) {
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(spec.mateStart());
                record.setInferredInsertSize(
                        (spec.flags() & 0x10) != 0 ? -100 : 100);
            }
            records.add(record);
        }
        records.sort((a, b) -> {
            final int byStart = Integer.compare(a.getAlignmentStart(), b.getAlignmentStart());
            return byStart != 0 ? byStart : a.getReadName().compareTo(b.getReadName());
        });
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header(), false, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /** A VCF of the sites the sample is known to be polymorphic at. */
    static String vcf(final List<Integer> positions) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=").append(CONTIG_LENGTH).append(">\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
        for (final int position : positions) {
            final char base = "ACGT".charAt((position - 1) % 4);
            final char alternate = base == 'A' ? 'C' : 'A';
            text.append("chr1\t").append(position).append("\t.\t").append(base).append('\t')
                    .append(alternate).append("\t100\tPASS\t.\tGT\t0/1\n");
        }
        return text.toString();
    }

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

    /** One run: the reads given, the known sites given, and whatever else the case asks for. */
    static void run(final String name, final List<Read> reads, final List<Integer> known,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("samerrormetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads);
        final Path sites = dir.resolve("known.vcf");
        Files.writeString(sites, vcf(known), StandardCharsets.UTF_8);
        // The VCF is QUERIED rather than streamed, so it needs its Tribble index beside it: the
        // tool asks the file for the sites overlapping each locus.
        htsjdk.tribble.index.IndexFactory.writeIndex(
                htsjdk.tribble.index.IndexFactory.createLinearIndex(
                        sites.toFile(), new htsjdk.variant.vcf.VCFCodec()),
                new File(sites + ".idx"));

        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "R=" + reference, "V=" + sites, "O=" + out.resolve("errors")));
        // `--ERROR_METRICS` is a collection whose DEFAULT is twenty-seven entries, and Picard's
        // parser APPENDS to a collection rather than replacing it, so naming one metric without
        // clearing the list first asks for it twice and the run is refused. Every case that names
        // metrics clears the list, and `the-default-list-is-appended-to` is the case that does not.
        if (Arrays.stream(extra).anyMatch(argument -> argument.startsWith("ERROR_METRICS="))
                && !name.equals("the-default-list-is-appended-to")) {
            argv.add("ERROR_METRICS=null");
        }
        argv.addAll(Arrays.asList(extra));

        final PrintStream original = System.out;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = new picard.sam.SamErrorMetric.CollectSamErrorMetrics()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                System.setOut(original);
                emit("error", name, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            System.setOut(original);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        } finally {
            System.setOut(original);
        }

        final List<String> written = new ArrayList<>();
        for (final File file : out.toFile().listFiles()) {
            written.add(file.getName());
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            emit("metrics", name + "." + file, table(out.resolve(file)));
        }
    }

    /** Reads over one window, the given number of them carrying a mismatch at one position. */
    static List<Read> reads(final int mismatches, final int mismatchAt, final char quality) {
        final List<Read> reads = new ArrayList<>();
        for (int index = 0; index < 8; index++) {
            final StringBuilder bases = new StringBuilder(reference(101, 20));
            final StringBuilder qualities = new StringBuilder("I".repeat(20));
            if (index < mismatches) {
                final int offset = mismatchAt - 101;
                bases.setCharAt(offset, bases.charAt(offset) == 'A' ? 'C' : 'A');
                qualities.setCharAt(offset, quality);
            }
            reads.add(new Read("r" + index, 101, bases.toString(), qualities.toString(), 0, 0,
                    "20M", 60));
        }
        return reads;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final List<Read> plain = reads(2, 105, 'I');
        run("two-mismatches", plain, List.of(), "ERROR_METRICS=ERROR");
        run("no-mismatches", reads(0, 105, 'I'), List.of(), "ERROR_METRICS=ERROR");
        // The prior is a pseudo-count, so a file with no errors still reports a finite quality and
        // the number moves when the prior does.
        run("no-mismatches-with-another-prior", reads(0, 105, 'I'), List.of(),
                "ERROR_METRICS=ERROR", "PRIOR_Q=10");

        // The same bases, split by a stratifier: more rows, the same totals.
        run("stratified-by-base-quality", plain, List.of(), "ERROR_METRICS=ERROR:BASE_QUALITY");
        run("stratified-by-cycle", plain, List.of(), "ERROR_METRICS=ERROR:CYCLE");
        run("two-metrics-two-files", plain, List.of(),
                "ERROR_METRICS=ERROR", "ERROR_METRICS=ERROR:GC_CONTENT");

        // The VCF is subtractive: the same mismatch at a known site is not an error.
        run("a-known-site", plain, List.of(105), "ERROR_METRICS=ERROR");

        // The two thresholds, one per base and one per read.
        run("a-low-quality-mismatch", reads(2, 105, '#'), List.of(), "ERROR_METRICS=ERROR");
        run("a-low-quality-mismatch-with-a-lower-threshold", reads(2, 105, '#'), List.of(),
                "ERROR_METRICS=ERROR", "MIN_BASE_Q=2");
        final List<Read> poorlyMapped = new ArrayList<>();
        for (final Read read : reads(2, 105, 'I')) {
            poorlyMapped.add(new Read(read.name(), read.start(), read.bases(), read.qualities(),
                    read.flags(), read.mateStart(), read.cigar(), 5));
        }
        run("poorly-mapped-reads", poorlyMapped, List.of(), "ERROR_METRICS=ERROR");
        run("poorly-mapped-reads-with-a-lower-threshold", poorlyMapped, List.of(),
                "ERROR_METRICS=ERROR", "MIN_MAPPING_Q=1");

        // A pair that overlaps itself, which is where OVERLAPPING_ERROR comes from.
        final List<Read> overlapping = new ArrayList<>();
        final String bases = reference(101, 20);
        final StringBuilder edited = new StringBuilder(bases);
        edited.setCharAt(4, edited.charAt(4) == 'A' ? 'C' : 'A');
        overlapping.add(new Read("p1", 101, edited.toString(), "I".repeat(20),
                0x1 | 0x2 | 0x40 | 0x20, 106, "20M", 60));
        overlapping.add(new Read("p1", 106, reference(106, 20), "I".repeat(20),
                0x1 | 0x2 | 0x80 | 0x10, 101, "20M", 60));
        run("an-overlapping-pair", overlapping, List.of(),
                "ERROR_METRICS=ERROR", "ERROR_METRICS=OVERLAPPING_ERROR");

        // An indel, which has a metric and a denominator of its own.
        final List<Read> indels = new ArrayList<>();
        for (int index = 0; index < 4; index++) {
            indels.add(new Read("d" + index, 101,
                    reference(101, 10) + reference(113, 8), "I".repeat(18), 0, 0, "10M2D8M", 60));
        }
        run("a-deletion", indels, List.of(), "ERROR_METRICS=INDEL_ERROR");

        // The loci themselves: a cap.
        run("a-cap-on-the-loci", plain, List.of(), "ERROR_METRICS=ERROR", "MAX_LOCI=5");

        // And the collection itself: naming a metric the default list already carries, without
        // clearing it, is a refusal rather than a duplicate row.
        run("the-default-list-is-appended-to", plain, List.of(), "ERROR_METRICS=ERROR");

        System.out.print(buf);
    }
}
