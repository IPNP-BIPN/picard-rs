/*
 * CollectWgsMetrics' metrics, taken from the reference.
 *
 * The tool walks every base of the reference and counts, for each one, how many reads cover it
 * and why the rest do not. What is measured is which base reaches the depth histogram and which
 * of the seven exclusion counters takes the ones that do not.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - GENOME_TERRITORY IS THE NON-N BASES OF THE REFERENCE and not the covered ones, so a
 *     reference padded with Ns has a territory smaller than its length;
 *   - THE EXCLUSIONS ARE A PARTITION OF THE BASES THAT DID NOT COUNT, and PCT_EXC_TOTAL is their
 *     sum, so the seven and the depth account for every base of every read over the territory;
 *   - --MINIMUM_MAPPING_QUALITY TAKES A WHOLE READ, its bases going to PCT_EXC_MAPQ;
 *   - --MINIMUM_BASE_QUALITY TAKES SINGLE BASES, so a read may contribute some of its bases and
 *     not others;
 *   - AN `N` BASE IS EXCLUDED BY QUALITY whatever its quality says, being treated as quality
 *     zero;
 *   - A DUPLICATE READ IS EXCLUDED WHOLE, and separately from the mapping-quality one;
 *   - AN UNPAIRED READ IS EXCLUDED UNLESS --COUNT_UNPAIRED, which is what PCT_EXC_UNPAIRED
 *     counts, and a pair with one end unmapped is unpaired for this purpose;
 *   - THE OVERLAP OF A PAIR IS COUNTED ONCE, the second end's bases going to PCT_EXC_OVERLAP,
 *     so a pair whose ends overlap entirely covers its span at depth one and not two;
 *   - --COVERAGE_CAP TRUNCATES THE DEPTH AND COUNTS THE REMAINDER, so a base covered ten times
 *     under a cap of two is depth two and eight bases of PCT_EXC_CAPPED;
 *   - MEAN_COVERAGE IS OVER THE TERRITORY AND NOT OVER THE COVERED BASES, so an uncovered half
 *     halves it;
 *   - --INCLUDE_BQ_HISTOGRAM ADDS A THIRD COLUMN TO THE HISTOGRAM, `unfiltered_baseq_count`,
 *     rather than a second table, and changes no metric;
 *   - AND A FILE WITH NO READS STILL REPORTS ITS TERRITORY, every depth being zero.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     histogram\t<case>\t<the histogram section, escaped>
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

public class CollectWgsMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 200;

    /** The reference: a repeating pattern, with ten Ns at the end so the territory is smaller. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH - 10; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        bases.append("NNNNNNNNNN");
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** One read: where it sits, how long it is, and what its flags and qualities say. */
    record Read(String name, int start, int length, int flags, int mappingQuality,
                String qualities, int mateStart, String bases) {}

    static Read read(final String name, final int start, final int length) {
        return new Read(name, start, length, 0, 60, null, 0, null);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setMappingQuality(spec.mappingQuality());
                record.setCigarString(spec.length() + "M");
                final StringBuilder bases = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    bases.append(spec.bases() == null ? 'A' : spec.bases().charAt(i));
                }
                record.setReadString(bases.toString());
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    quals.append(spec.qualities() == null ? 'I' : spec.qualities().charAt(i));
                }
                record.setBaseQualityString(quals.toString());
                if ((spec.flags() & 0x1) != 0) {
                    if ((spec.flags() & 0x8) != 0) {
                        record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                        record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    } else {
                        record.setMateReferenceName("chr1");
                        record.setMateAlignmentStart(spec.mateStart());
                    }
                }
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

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
        final Path dir = Files.createTempDirectory("wgsmetrics");
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
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + metrics, "R=" + reference));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.CollectWgsMetrics()
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

    /** A pair whose two ends are given explicitly, so their overlap can be chosen. */
    static List<Read> pair(final String name, final int first, final int second, final int length) {
        return List.of(
                new Read(name, first, length, 0x1 | 0x2 | 0x40, 60, null, second, null),
                new Read(name, second, length, 0x1 | 0x2 | 0x80, 60, null, first, null));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // One unpaired read, which the defaults exclude entirely.
        run("one-unpaired-read", List.of(read("a", 1, 20)));
        run("one-unpaired-read-counted", List.of(read("a", 1, 20)), "COUNT_UNPAIRED=true");

        // A pair whose ends do not overlap, and one whose ends overlap entirely.
        run("pair-disjoint", pair("a", 1, 50, 20));
        run("pair-overlapping", pair("a", 1, 1, 20));

        // A read under the mapping-quality floor.
        run("low-mapping-quality", List.of(
                new Read("a", 1, 20, 0, 5, null, 0, null)), "COUNT_UNPAIRED=true");
        // A read half of whose bases are under the base-quality floor.
        run("low-base-quality", List.of(
                new Read("a", 1, 20, 0, 60, "IIIIIIIIII##########", 0, null)),
                "COUNT_UNPAIRED=true");
        // A read with N bases, which are excluded whatever their quality says.
        run("n-bases", List.of(
                new Read("a", 1, 20, 0, 60, null, 0, "AAAAANNNNNAAAAAAAAAA")),
                "COUNT_UNPAIRED=true");
        // A duplicate read.
        run("duplicate", List.of(
                new Read("a", 1, 20, 0x400, 60, null, 0, null)), "COUNT_UNPAIRED=true");
        // A pair with one end unmapped, which counts as unpaired.
        run("mate-unmapped", List.of(
                new Read("a", 1, 20, 0x1 | 0x8 | 0x40, 60, null, 0, null)));

        // Ten reads over the same twenty bases, under a cap of two.
        final List<Read> deep = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            deep.add(new Read("d" + i, 1, 20, 0, 60, null, 0, null));
        }
        run("deep-uncapped", deep, "COUNT_UNPAIRED=true");
        run("deep-capped", deep, "COUNT_UNPAIRED=true", "COVERAGE_CAP=2");

        // The base-quality histogram, which changes no metric.
        run("with-histogram", pair("a", 1, 50, 20), "INCLUDE_BQ_HISTOGRAM=true");

        // A file with no reads at all.
        run("empty", List.of());

        System.out.print(buf);
    }
}
